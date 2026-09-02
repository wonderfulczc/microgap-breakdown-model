from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.source.adapters import (  # noqa: E402
    AFIVO_CONTINUITY_CONSISTENT_SOURCE_DEFINITION,
    AFIVO_STAGE_E_SOURCE_DEFINITION,
    infer_afivo_cell_widths_from_config,
    load_afivo_stage_e_source_csv,
    source_record_from_petsc_axisymmetric,
)
from streamer_rf.rf.source.audit import (  # noqa: E402
    boundary_flux,
    global_continuity_audit,
    local_continuity_on_matching_grid,
    remap_audit,
)
from streamer_rf.rf.source.derivatives import central_difference_three_point  # noqa: E402
from streamer_rf.rf.source.integrals import current_moment, current_profile_z, frequency_metadata, total_charge  # noqa: E402
from streamer_rf.rf.source.remap import conservative_remap, nearest_remap  # noqa: E402
from streamer_rf.rf.source.schema import SourceMetadata, SourceRecord, SourceSeries  # noqa: E402
from streamer_rf.rf.source.synthetic import moving_gaussian_on_edges, moving_gaussian_series  # noqa: E402


def test_source_record_schema_and_3d_integrals() -> None:
    meta = SourceMetadata(
        case_id="unit",
        solver="synthetic",
        solver_version="test",
        time_s=0.0,
        coordinate_system="cartesian_3d_native_amr",
        pressure_Pa=101325.0,
        temperature_K=300.0,
        geometry_id="box",
        voltage_state="none",
        photoionization="off",
        source_definition="unit",
    )
    cols = {
        "level": np.array([0]),
        "x_center": np.array([0.0]),
        "y_center": np.array([0.0]),
        "z_center": np.array([0.0]),
        "dx": np.array([2.0]),
        "dy": np.array([3.0]),
        "dz": np.array([4.0]),
        "cell_volume": np.array([24.0]),
        "rho": np.array([5.0]),
        "Jx": np.array([1.0]),
        "Jy": np.array([2.0]),
        "Jz": np.array([3.0]),
    }
    rec = SourceRecord(meta, cols)
    assert total_charge(rec) == 120.0
    np.testing.assert_allclose(current_moment(rec), [24.0, 48.0, 72.0])
    np.testing.assert_allclose(current_profile_z(rec, np.array([0.0])), [18.0])


def test_axisymmetric_volume_integration_and_petsc_adapter() -> None:
    r = np.array([0.5, 1.5])
    z = np.array([0.25, 0.75])
    rho = np.ones((2, 2)) * 2.0
    Jr = np.zeros_like(rho)
    Jz = np.ones_like(rho) * 4.0
    meta = SourceMetadata("petsc", "PETSc-2D", "stage-c", 0.0, "axisymmetric_rz", 1.0, 1.0, "g", "v", "off", "adapter")
    rec = source_record_from_petsc_axisymmetric(r, z, 1.0, 0.5, rho, Jr, Jz, meta)
    expected_vol = np.sum(2.0 * np.pi * r[:, None] * 1.0 * 0.5 * np.ones_like(rho))
    np.testing.assert_allclose(total_charge(rec), 2.0 * expected_vol)
    np.testing.assert_allclose(current_moment(rec)[2], 4.0 * expected_vol)


def test_conservative_amr_remap_preserves_Q_and_M() -> None:
    src_edges = (np.array([0.0, 0.25, 0.5, 1.0]), np.array([0.0, 1.0]), np.array([0.0, 1.0]))
    dst_edges = (np.array([0.0, 0.5, 1.0]), np.array([0.0, 1.0]), np.array([0.0, 1.0]))
    src = moving_gaussian_on_edges(*src_edges, time_s=0.0)
    dst_geom = moving_gaussian_on_edges(*dst_edges, time_s=0.0).geometry_like()
    remapped = conservative_remap(src, dst_geom)
    audit = remap_audit(src, remapped)
    assert audit.relative_charge_error < 1e-14
    assert audit.relative_M_error < 1e-14


def test_bad_remap_artifact_is_larger_than_conservative() -> None:
    ref_edges = (
        np.linspace(-3.0, 3.0, 17),
        np.linspace(-3.0, 3.0, 9),
        np.linspace(-3.0, 3.0, 9),
    )
    shifted_edges = (
        np.array([-3.0, -2.25, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.25, 3.0]),
        np.linspace(-3.0, 3.0, 9),
        np.linspace(-3.0, 3.0, 9),
    )
    mid = moving_gaussian_on_edges(*shifted_edges, time_s=0.4)
    target = moving_gaussian_on_edges(*ref_edges, time_s=0.4).geometry_like()
    cons = remap_audit(mid, conservative_remap(mid, target))
    bad = remap_audit(mid, nearest_remap(mid, target))
    assert cons.relative_charge_error < 1e-13
    assert bad.relative_charge_error > cons.relative_charge_error + 1e-3


def test_moving_gaussian_continuity_and_central_derivative() -> None:
    times = np.array([0.0, 0.4, 0.8])
    vals = np.vstack([times**2, 2.0 * times + 1.0]).T
    deriv = central_difference_three_point(times, vals)
    np.testing.assert_allclose(deriv, [0.8, 2.0], rtol=1e-14, atol=1e-14)

    edges = (np.linspace(-3.0, 3.0, 17), np.linspace(-3.0, 3.0, 13), np.linspace(-3.0, 3.0, 13))
    series = moving_gaussian_series(edges, tuple(times))
    global_audit = global_continuity_audit(series, boundary_flux_at_center=boundary_flux(series.records[1], atol=1e-12))
    assert global_audit["continuity_global_abs_A"] < 2e-3
    local = local_continuity_on_matching_grid(series)
    assert local["continuity_local_L2"] < 0.25


def test_frequency_metadata_nonuniform_and_uniform() -> None:
    edges = (np.linspace(0.0, 1.0, 3), np.linspace(0.0, 1.0, 3), np.linspace(0.0, 1.0, 3))
    uniform = frequency_metadata(moving_gaussian_series(edges, (0.0, 1.0, 2.0)))
    assert uniform["uniform_output_interval"] is True
    assert uniform["nyquist_frequency_Hz"] == 0.5
    nonuniform = frequency_metadata(moving_gaussian_series(edges, (0.0, 1.0, 3.0)))
    assert nonuniform["uniform_output_interval"] is False


def test_afivo_stage_e_adapter_smoke_if_raw_available() -> None:
    cfg = ROOT / "solver3d/afivo_reference/stage_e/configs/e3_aligned_needle_pair_500V.cfg"
    raw = ROOT / "solver3d/afivo_reference/stage_e/results_raw/e3_aligned_needle_pair_500V_source_000000.csv"
    if not raw.exists():
        return
    rec = load_afivo_stage_e_source_csv(
        raw,
        cfg,
        case_id="E3_aligned_needle_pair_500V",
        solver_version="a50b5508775086e90dfe423455fb58d812578410",
        geometry_id="E3_aligned_needle_pair",
        voltage_state="500 V constant",
        photoionization="off",
        max_rows=1024,
    )
    assert rec.n_cells == 1024
    assert "CELL_CENTERED_DRIFT_CURRENT" in rec.metadata.source_definition
    assert rec.metadata.extra["J_transport_available"] is False
    assert total_charge(rec) == total_charge(rec)
    assert np.linalg.norm(current_moment(rec)) >= 0.0
    assert AFIVO_STAGE_E_SOURCE_DEFINITION == rec.metadata.source_definition


def test_afivo_adapter_prefers_continuity_consistent_jrf_columns(tmp_path: Path) -> None:
    cfg = tmp_path / "case.cfg"
    cfg.write_text("domain_len = 1.0 1.0 1.0\nbox_size = 1\n")
    raw = tmp_path / "source.csv"
    raw.write_text(
        "time_s,x_m,y_m,z_m,cell_volume_m3,rho_Cpm3,ne_m3,"
        "Ex_Vpm,Ey_Vpm,Ez_Vpm,Jx_Apm2,Jy_Apm2,Jz_Apm2,"
        "Jrf_x_Apm2,Jrf_y_Apm2,Jrf_z_Apm2,Eabs_Vpm,lsf_m,level\n"
        "1e-12,0.5,0.5,0.5,1.0,2.0,3.0,4.0,5.0,6.0,10.0,20.0,30.0,-1.0,-2.0,-3.0,7.0,8.0,1\n"
    )
    rec = load_afivo_stage_e_source_csv(
        raw,
        cfg,
        case_id="unit",
        solver_version="afivo-test",
        geometry_id="g",
        voltage_state="v",
        photoionization="off",
    )
    assert rec.metadata.source_definition == AFIVO_CONTINUITY_CONSISTENT_SOURCE_DEFINITION
    assert rec.metadata.extra["J_transport_available"] is True
    assert rec.metadata.extra["J_export_columns"] == "Jrf_x_Apm2,Jrf_y_Apm2,Jrf_z_Apm2"
    np.testing.assert_allclose(rec.columns["Jx"], [-1.0])
    np.testing.assert_allclose(rec.columns["Jy"], [-2.0])
    np.testing.assert_allclose(rec.columns["Jz"], [-3.0])


def test_afivo_cell_width_inference_matches_known_config() -> None:
    cfg = ROOT / "solver3d/afivo_reference/stage_e/configs/e3_aligned_needle_pair_500V.cfg"
    levels = np.array([1, 5])
    dx, dy, dz = infer_afivo_cell_widths_from_config(cfg, levels)
    np.testing.assert_allclose(dx, [1.5e-5, 9.375e-7])
    np.testing.assert_allclose(dy, [1.5e-5, 9.375e-7])
    np.testing.assert_allclose(dz, [2.5e-5, 1.5625e-6])
