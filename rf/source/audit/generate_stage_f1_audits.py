#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))

from streamer_rf.rf.source.adapters import (  # noqa: E402
    AFIVO_STAGE_E_SOURCE_DEFINITION,
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
from streamer_rf.rf.source.integrals import current_moment, frequency_metadata, total_charge  # noqa: E402
from streamer_rf.rf.source.remap import conservative_remap, nearest_remap  # noqa: E402
from streamer_rf.rf.source.schema import SourceMetadata, SourceSeries  # noqa: E402
from streamer_rf.rf.source.synthetic import moving_gaussian_series  # noqa: E402


AFIVO_SHA = "a50b5508775086e90dfe423455fb58d812578410"


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def synthetic_audit() -> dict:
    ref_edges = (
        np.linspace(-3.0, 3.0, 17),
        np.linspace(-3.0, 3.0, 13),
        np.linspace(-3.0, 3.0, 13),
    )
    shifted_edges = (
        np.array([-3.0, -2.25, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.25, 3.0]),
        np.linspace(-3.0, 3.0, 13),
        np.linspace(-3.0, 3.0, 13),
    )
    times = (0.0, 0.3, 0.8)
    ref_series = moving_gaussian_series(ref_edges, times)
    shifted_mid = moving_gaussian_series(shifted_edges, (times[1],)).records[0]
    remapped_mid = conservative_remap(shifted_mid, ref_series.records[1].geometry_like())
    bad_mid = nearest_remap(shifted_mid, ref_series.records[1].geometry_like())

    cons = remap_audit(shifted_mid, remapped_mid).to_dict()
    bad = remap_audit(shifted_mid, bad_mid).to_dict()
    fixed_series = SourceSeries((ref_series.records[0], remapped_mid, ref_series.records[2]))
    continuity = global_continuity_audit(
        fixed_series, boundary_flux_at_center=boundary_flux(remapped_mid, atol=1e-12)
    )
    local = local_continuity_on_matching_grid(ref_series)

    q_cons = [total_charge(ref_series.records[0]), total_charge(remapped_mid), total_charge(ref_series.records[2])]
    q_bad = [total_charge(ref_series.records[0]), total_charge(bad_mid), total_charge(ref_series.records[2])]
    m_cons = np.vstack([current_moment(ref_series.records[0]), current_moment(remapped_mid), current_moment(ref_series.records[2])])
    m_bad = np.vstack([current_moment(ref_series.records[0]), current_moment(bad_mid), current_moment(ref_series.records[2])])
    dQdt_cons = float(central_difference_three_point(np.asarray(times), np.asarray(q_cons)))
    dQdt_bad = float(central_difference_three_point(np.asarray(times), np.asarray(q_bad)))
    dMdt_cons = central_difference_three_point(np.asarray(times), m_cons)
    dMdt_bad = central_difference_three_point(np.asarray(times), m_bad)

    return {
        "test": "moving Gaussian across AMR refinement boundary",
        "conservative_remap": cons,
        "bad_nearest_remap": bad,
        "dQdt_artifact_conservative": abs(dQdt_cons),
        "dQdt_artifact_bad_nearest": abs(dQdt_bad),
        "dMdt_artifact_conservative_norm": float(np.linalg.norm(dMdt_cons)),
        "dMdt_artifact_bad_nearest_norm": float(np.linalg.norm(dMdt_bad)),
        "global_continuity": continuity,
        "local_continuity": local,
        "frequency_metadata": frequency_metadata(ref_series),
    }


def petsc_adapter_audit() -> dict:
    r = np.array([0.5e-6, 1.5e-6, 2.5e-6])
    z = np.array([0.5e-6, 1.5e-6])
    rho = np.full((r.size, z.size), 2.0)
    Jr = np.zeros_like(rho)
    Jz = np.full_like(rho, 3.0)
    meta = SourceMetadata(
        case_id="petsc_axisymmetric_smoke",
        solver="PETSc-2D",
        solver_version="stage-c-frozen",
        time_s=0.0,
        coordinate_system="axisymmetric_rz",
        pressure_Pa=101325.0,
        temperature_K=300.0,
        geometry_id="axisymmetric_adapter_smoke",
        voltage_state="synthetic",
        photoionization="off",
        source_definition="axisymmetric rho/Jr/Jz post-processing adapter",
    )
    record = source_record_from_petsc_axisymmetric(r, z, 1.0e-6, 1.0e-6, rho, Jr, Jz, meta)
    expected_volume = float(np.sum(2.0 * np.pi * r[:, None] * 1.0e-6 * 1.0e-6 * np.ones_like(rho)))
    return {
        "Q_C": total_charge(record),
        "Q_expected_C": 2.0 * expected_volume,
        "Mz_Am": float(current_moment(record)[2]),
        "Mz_expected_Am": 3.0 * expected_volume,
        "n_cells": record.n_cells,
    }


def afivo_real_source_audit() -> dict:
    cfg = ROOT / "solver3d/afivo_reference/stage_e/configs/e3_aligned_needle_pair_500V.cfg"
    src0 = ROOT / "solver3d/afivo_reference/stage_e/results_raw/e3_aligned_needle_pair_500V_source_000000.csv"
    src1 = ROOT / "solver3d/afivo_reference/stage_e/results_raw/e3_aligned_needle_pair_500V_source_000001.csv"
    if not (src0.exists() and src1.exists()):
        return {
            "available": False,
            "reason": "Stage E raw source snapshots not present; rerun frozen E3 source export to audit real data.",
        }
    records = [
        load_afivo_stage_e_source_csv(
            p,
            cfg,
            case_id="E3_aligned_needle_pair_500V",
            solver_version=AFIVO_SHA,
            geometry_id="E3_aligned_needle_pair",
            voltage_state="500 V constant",
            photoionization="off",
        )
        for p in (src0, src1)
    ]
    freq = frequency_metadata(SourceSeries(tuple(records)))
    duration = freq["record_duration_s"]
    return {
        "available": True,
        "source_files": [str(src0), str(src1)],
        "snapshot_count": len(records),
        "actual_output_times_s": freq["actual_output_times_s"],
        "Q_C": [total_charge(r) for r in records],
        "M_Am": [[float(x) for x in current_moment(r)] for r in records],
        "boundary_flux_A_final": boundary_flux(records[-1], atol=1e-15),
        "same_partition_identity_audit": remap_audit(records[-1], records[-1]).to_dict(),
        "source_definition": AFIVO_STAGE_E_SOURCE_DEFINITION,
        "J_transport_available": False,
        "USABLE_FOR_RF_PIPELINE_TEST": True,
        "USABLE_FOR_VHF_UHF_SCIENTIFIC_SPECTRUM": False,
        "output_dt_s": freq["output_dt_s"],
        "duration_s": duration,
        "nyquist_frequency_Hz": freq["nyquist_frequency_Hz"],
        "raw_frequency_resolution_Hz": freq["raw_frequency_resolution_Hz"],
        "note": "Existing committed Stage E source cadence contains two full raw source snapshots, so real-source time derivatives are not claimed here. Synthetic three-snapshot tests validate derivative and remap logic.",
    }


def main() -> None:
    out = ROOT / "rf/source/audit"
    syn = synthetic_audit()
    petsc = petsc_adapter_audit()
    afivo = afivo_real_source_audit()
    manifest = {
        "stage": "F1",
        "CURRENT_SOURCE_PROVENANCE": "AFIVO_STAGE_E_CELL_CENTERED_DRIFT_CURRENT",
        "final_RF_J_definition": "-e * Gamma_e from the finite-volume electron transport flux",
        "NATIVE_AMR_SOURCE": "reference",
        "synthetic": syn,
        "petsc_adapter": petsc,
        "real_afivo_source": afivo,
        "usable_for_rf": bool(afivo.get("available", False)),
        "usable_for_vhf_uhf_scientific_spectrum": False,
    }
    write_json(out / "stage_f1_synthetic_audit.json", syn)
    write_json(out / "stage_f1_petsc_adapter_manifest.json", petsc)
    write_json(out / "stage_f1_real_afivo_manifest.json", afivo)
    write_json(out / "rf_source_manifest.json", manifest)


if __name__ == "__main__":
    main()
