from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.jefimenko.axisymmetric import rotate_axisymmetric_series  # noqa: E402
from streamer_rf.rf.jefimenko.constants import C0, K_B, K_E  # noqa: E402
from streamer_rf.rf.jefimenko.diagnostics import current_moment_radiation_approx, near_far_audit  # noqa: E402
from streamer_rf.rf.jefimenko.observer import Observer  # noqa: E402
from streamer_rf.rf.jefimenko.solver import evaluate_observer, evaluate_waveform  # noqa: E402
from streamer_rf.rf.jefimenko.synthetic import (  # noqa: E402
    dipole_far_field_E,
    dipole_radiation_series,
    make_source_record_from_arrays,
    static_gaussian_charge_series,
    steady_current_element_series,
)
from streamer_rf.rf.source.adapters import load_afivo_stage_e_source_csv, source_record_from_petsc_axisymmetric  # noqa: E402
from streamer_rf.rf.source.derivatives import central_difference_three_point  # noqa: E402
from streamer_rf.rf.source.integrals import current_moment  # noqa: E402
from streamer_rf.rf.source.schema import SourceMetadata, SourceSeries  # noqa: E402


def _retarded_time(observer: np.ndarray, source_center: np.ndarray, tr: float) -> float:
    return tr + float(np.linalg.norm(observer - source_center)) / C0


def _relative(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def test_static_charge_coulomb_field_and_B_zero() -> None:
    series = static_gaussian_charge_series(total_charge_C=1e-12, times_s=(0.0, 1e-9, 2e-9, 3e-9, 4e-9))
    obs = Observer("x_far", 0.2, 0.0, 0.0)
    sample = evaluate_observer(series, obs, _retarded_time(obs.position, np.zeros(3), 2e-9))
    expected = K_E * 1e-12 / 0.2**2
    assert sample.retarded_time_valid
    assert _relative(sample.E_total[0], expected) < 4e-3
    assert abs(sample.E_total[1]) < abs(sample.E_total[0]) * 1e-12
    assert abs(sample.E_total[2]) < abs(sample.E_total[0]) * 1e-12
    assert np.linalg.norm(sample.B_total) < 1e-18


def test_static_charge_has_no_induction_or_radiation_terms() -> None:
    series = static_gaussian_charge_series(total_charge_C=1e-12, times_s=(0.0, 1e-9, 2e-9, 3e-9, 4e-9))
    obs = Observer("x_far", 0.2, 0.0, 0.0)
    sample = evaluate_observer(series, obs, _retarded_time(obs.position, np.zeros(3), 2e-9))
    assert np.linalg.norm(sample.E_drho_induction) < np.linalg.norm(sample.E_rho_near) * 1e-12
    assert np.linalg.norm(sample.E_dJ_radiation) == 0.0
    assert np.linalg.norm(sample.B_J_near) == 0.0
    assert np.linalg.norm(sample.B_dJ_radiation) == 0.0


def test_steady_current_magnetic_field_direction_and_no_radiation() -> None:
    Mz = 2e-6
    series = steady_current_element_series(current_moment_Am=Mz, times_s=(0.0, 1e-9, 2e-9, 3e-9, 4e-9))
    obs = Observer("x_far", 0.1, 0.0, 0.0)
    sample = evaluate_observer(series, obs, _retarded_time(obs.position, np.zeros(3), 2e-9))
    expected_By = K_B * Mz / 0.1**2
    assert sample.retarded_time_valid
    assert _relative(sample.B_total[1], expected_By) < 3e-4
    assert sample.B_total[1] > 0
    assert np.linalg.norm(sample.B_dJ_radiation) < abs(expected_By) * 1e-12


def test_magnetic_cross_product_right_hand_rule() -> None:
    series = steady_current_element_series(current_moment_Am=1e-6, times_s=(0.0, 1e-9, 2e-9, 3e-9, 4e-9))
    pos_x = Observer("pos_x", 0.1, 0.0, 0.0)
    neg_x = Observer("neg_x", -0.1, 0.0, 0.0)
    b_pos = evaluate_observer(series, pos_x, _retarded_time(pos_x.position, np.zeros(3), 2e-9)).B_total
    b_neg = evaluate_observer(series, neg_x, _retarded_time(neg_x.position, np.zeros(3), 2e-9)).B_total
    assert b_pos[1] > 0.0
    assert b_neg[1] < 0.0
    np.testing.assert_allclose(abs(b_pos[1]), abs(b_neg[1]), rtol=1e-12)


def test_radiation_scaling_angular_dependence_and_E_over_B() -> None:
    omega = 2.0 * math.pi * 1.0e9
    p0 = 1e-18
    tr = 0.25e-9
    series = dipole_radiation_series(p0_Cm=p0, omega_rad_s=omega, times_s=np.linspace(-1.5e-9, 1.5e-9, 121))
    radii = np.array([0.2, 0.3, 0.5])
    e_perp = []
    for R in radii:
        obs = Observer(f"R{R}", R, 0.0, 0.0)
        sample = evaluate_observer(series, obs, _retarded_time(obs.position, np.zeros(3), tr))
        assert sample.retarded_time_valid
        rhat = obs.position / np.linalg.norm(obs.position)
        evec = sample.E_total - np.dot(sample.E_total, rhat) * rhat
        bvec = sample.B_total
        e_perp.append(np.linalg.norm(evec))
        assert _relative(np.linalg.norm(evec) / np.linalg.norm(bvec), C0) < 0.06
        assert abs(np.dot(evec, rhat)) < np.linalg.norm(evec) * 1e-10
    slope = np.polyfit(np.log(radii), np.log(e_perp), 1)[0]
    assert abs(slope + 1.0) < 0.12

    pddot = -p0 * omega**2 * math.sin(omega * tr)
    angles = np.deg2rad([30.0, 60.0, 90.0])
    amps = []
    for theta in angles:
        pos = np.array([0.3 * math.sin(theta), 0.0, 0.3 * math.cos(theta)])
        sample = evaluate_observer(series, Observer(f"theta{theta}", *pos), _retarded_time(pos, np.zeros(3), tr))
        rhat = pos / np.linalg.norm(pos)
        evec = sample.E_total - np.dot(sample.E_total, rhat) * rhat
        amps.append(np.linalg.norm(evec))
    ratios = np.asarray(amps) / amps[-1]
    np.testing.assert_allclose(ratios, np.sin(angles), rtol=0.08, atol=0.02)
    expected_90 = dipole_far_field_E(math.pi / 2.0, 0.3, pddot)
    assert _relative(amps[-1], expected_90) < 0.08


def test_static_charge_near_term_has_inverse_square_scaling() -> None:
    series = static_gaussian_charge_series(total_charge_C=1e-12, times_s=(0.0, 1e-9, 2e-9, 3e-9, 4e-9))
    radii = np.array([0.08, 0.14, 0.25])
    amps = []
    for R in radii:
        obs = Observer(f"R{R}", R, 0.0, 0.0)
        sample = evaluate_observer(series, obs, _retarded_time(obs.position, np.zeros(3), 2e-9))
        amps.append(np.linalg.norm(sample.E_rho_near))
    slope = np.polyfit(np.log(radii), np.log(amps), 1)[0]
    assert abs(slope + 2.0) < 0.03


def test_propagation_delay_from_retarded_time() -> None:
    omega = 2.0 * math.pi * 1.0e9
    series = dipole_radiation_series(p0_Cm=1e-18, omega_rad_s=omega, times_s=np.linspace(-1.5e-9, 1.5e-9, 121))
    R1, R2 = 0.2, 0.35
    tr = 0.25e-9
    t1 = _retarded_time(np.array([R1, 0.0, 0.0]), np.zeros(3), tr)
    t2 = _retarded_time(np.array([R2, 0.0, 0.0]), np.zeros(3), tr)
    s1 = evaluate_observer(series, Observer("R1", R1, 0.0, 0.0), t1)
    s2 = evaluate_observer(series, Observer("R2", R2, 0.0, 0.0), t2)
    np.testing.assert_allclose(t2 - t1, (R2 - R1) / C0, rtol=1e-15)
    assert _relative(np.linalg.norm(s1.E_total) * R1, np.linalg.norm(s2.E_total) * R2) < 0.04


def test_retarded_interpolation_converges_with_output_dt() -> None:
    omega = 2.0 * math.pi * 1.0e9
    p0 = 1e-18
    tr = 0.25e-9
    obs = Observer("x", 0.3, 0.0, 0.0)
    expected = dipole_far_field_E(math.pi / 2.0, 0.3, -p0 * omega**2 * math.sin(omega * tr))
    coarse = dipole_radiation_series(p0_Cm=p0, omega_rad_s=omega, times_s=np.linspace(-1.5e-9, 1.5e-9, 61))
    fine = dipole_radiation_series(p0_Cm=p0, omega_rad_s=omega, times_s=np.linspace(-1.5e-9, 1.5e-9, 121))
    ec = np.linalg.norm(evaluate_observer(coarse, obs, _retarded_time(obs.position, np.zeros(3), tr)).E_total)
    ef = np.linalg.norm(evaluate_observer(fine, obs, _retarded_time(obs.position, np.zeros(3), tr)).E_total)
    assert _relative(ef, expected) < _relative(ec, expected)


def test_axisymmetric_phi_quadrature_converges() -> None:
    r = np.linspace(0.2e-3, 1.0e-3, 5)
    z = np.linspace(-0.4e-3, 0.4e-3, 5)
    rr, zz = np.meshgrid(r, z, indexing="ij")
    rho = 1e-7 * np.exp(-((rr - 0.6e-3) ** 2 + zz**2) / (2 * (0.25e-3) ** 2))
    Jr = np.zeros_like(rho)
    Jz = np.zeros_like(rho)
    meta = SourceMetadata("axisym", "PETSc-2D", "test", 0.0, "axisymmetric_rz", 1.0, 1.0, "ring", "none", "off", "static axisymmetric")
    base = source_record_from_petsc_axisymmetric(r, z, r[1] - r[0], z[1] - z[0], rho, Jr, Jz, meta)
    series = SourceSeries(
        tuple(
            base.with_source_columns(base.columns["rho"], np.zeros((base.n_cells, 3)), time_s=t)
            for t in (0.0, 1e-9, 2e-9, 3e-9, 4e-9)
        )
    )
    obs = Observer("off_axis", 0.01, 0.003, 0.002)
    vals = {}
    for nphi in (8, 16, 32, 128):
        rotated = rotate_axisymmetric_series(series, nphi)
        vals[nphi] = evaluate_observer(rotated, obs, _retarded_time(obs.position, np.zeros(3), 2e-9)).E_total
    err8 = np.linalg.norm(vals[8] - vals[128])
    err16 = np.linalg.norm(vals[16] - vals[128])
    err32 = np.linalg.norm(vals[32] - vals[128])
    assert err16 < err8
    assert err32 / np.linalg.norm(vals[128]) < 1e-12


def test_current_moment_far_field_approximation_improves_with_distance() -> None:
    omega = 2.0 * math.pi * 1.0e9
    p0 = 1e-18
    tr = 0.25e-9
    series = dipole_radiation_series(p0_Cm=p0, omega_rad_s=omega, sigma_m=8e-4, half_width_m=4e-3, times_s=np.linspace(-1.5e-9, 1.5e-9, 121))
    dMdt = np.array([0.0, 0.0, -p0 * omega**2 * math.sin(omega * tr)])
    errors = []
    for R in (0.12, 0.4):
        obs = Observer(f"R{R}", R, 0.0, 0.0)
        full = evaluate_observer(series, obs, _retarded_time(obs.position, np.zeros(3), tr)).E_total
        approx = current_moment_radiation_approx(dMdt, obs.position)
        errors.append(np.linalg.norm(full - approx) / np.linalg.norm(full))
    assert errors[1] < errors[0]
    assert errors[1] < 0.08


def test_observer_inside_source_validation_and_invalid_retarded_time() -> None:
    series = static_gaussian_charge_series(times_s=(0.0, 1e-9, 2e-9, 3e-9))
    inside = Observer("inside", 0.0, 0.0, 0.0)
    try:
        evaluate_observer(series, inside, 2e-9)
    except ValueError as exc:
        assert "inside a source cell" in str(exc)
    else:
        raise AssertionError("observer inside source should fail")
    outside = Observer("outside", 0.2, 0.0, 0.0)
    sample = evaluate_observer(series, outside, 0.0)
    assert not sample.retarded_time_valid
    assert sample.retarded_time_valid_fraction < 1.0


def test_observer_waveform_output_schema_and_near_far_audit() -> None:
    series = steady_current_element_series(times_s=(0.0, 1e-9, 2e-9, 3e-9, 4e-9))
    obs = Observer("x", 0.1, 0.0, 0.0)
    samples = evaluate_waveform(series, [obs], np.array([_retarded_time(obs.position, np.zeros(3), 2e-9)]), source_manifest_id="unit")
    row = samples[0].to_row()
    for key in ("Ex_total", "Ey_rho", "Ez_dJ", "Bx_total", "By_J", "Bz_dJ", "retarded_time_valid"):
        assert key in row
    audit = near_far_audit(series, obs.position, 1e-9)
    assert audit["epsilon_geom_L_over_R"] < 0.01
    assert audit["epsilon_EM_L_over_ctau"] < 0.01


def test_real_afivo_pipeline_smoke_scientific_gate() -> None:
    manifest_path = ROOT / "rf/source/audit/rf_source_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["real_afivo_source"]["J_transport_available"] is False
    assert manifest["real_afivo_source"]["USABLE_FOR_VHF_UHF_SCIENTIFIC_SPECTRUM"] is False

    cfg = ROOT / "solver3d/afivo_reference/stage_e/configs/e3_aligned_needle_pair_500V.cfg"
    files = [
        ROOT / "solver3d/afivo_reference/stage_e/results_raw/e3_aligned_needle_pair_500V_source_000002.csv",
        ROOT / "solver3d/afivo_reference/stage_e/results_raw/e3_aligned_needle_pair_500V_source_000003.csv",
        ROOT / "solver3d/afivo_reference/stage_e/results_raw/e3_aligned_needle_pair_500V_source_000004.csv",
        ROOT / "solver3d/afivo_reference/stage_e/results_raw/e3_aligned_needle_pair_500V_source_000005.csv",
    ]
    if not all(p.exists() for p in files):
        return
    selected = []
    import csv

    with files[-1].open(newline="") as handle:
        reader = csv.DictReader(handle)
        for i, row in enumerate(reader):
            if abs(float(row["rho_Cpm3"])) + abs(float(row["Jz_Apm2"])) + 1e-30 * abs(float(row["ne_m3"])) > 0.0:
                selected.append(i)
            if len(selected) >= 512:
                break
    records = [
        load_afivo_stage_e_source_csv(
            path,
            cfg,
            case_id="E3_aligned_pipeline_smoke",
            solver_version="a50b5508775086e90dfe423455fb58d812578410",
            geometry_id="E3_aligned_needle_pair",
            voltage_state="500 V constant",
            photoionization="off",
            row_indices=np.asarray(selected, dtype=int),
        )
        for path in files
    ]
    series = SourceSeries(tuple(records))
    obs = Observer("pipeline_observer", 2e-4, 0.0, 0.0)
    xyz = np.column_stack((records[0].columns["x_center"], records[0].columns["y_center"], records[0].columns["z_center"]))
    median_delay = float(np.median(np.linalg.norm(obs.position[None, :] - xyz, axis=1)) / C0)
    sample = evaluate_observer(series, obs, 3.5e-12 + median_delay, source_manifest_id="PIPELINE_ONLY")
    assert sample.source_manifest_id == "PIPELINE_ONLY"
    assert np.all(np.isfinite(sample.E_total))
    assert np.all(np.isfinite(sample.B_total))
