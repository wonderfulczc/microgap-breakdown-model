#!/usr/bin/env python3
"""Run the compact H4 2.5 mm short-dipole receiver fixtures."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.ports import CurvePort, UI_data

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
from streamer_rf.fullwave.foundation import mesh_report  # noqa: E402


F_START_HZ = 2.8e9
F_STOP_HZ = 10.5e9
F_CENTER_HZ = 0.5 * (F_START_HZ + F_STOP_HZ)
FREQUENCY_HZ = np.linspace(F_START_HZ, F_STOP_HZ, 771)


def build_case(output: Path, resolution: str, excitation: str) -> None:
    started = time.perf_counter()
    if output.exists():
        raise ValueError("EXISTING_RAW_CASE_REQUIRES_EXPLICIT_REUSE")
    output.mkdir(parents=True)
    if resolution not in {"baseline", "fine"}:
        raise ValueError("INVALID_H4_RESOLUTION")
    if excitation not in {"port", "plane_parallel", "plane_orthogonal"}:
        raise ValueError("INVALID_H4_EXCITATION")

    local_step_mm = 0.25 if resolution == "baseline" else 0.1875
    max_step_mm = 1.25 if resolution == "baseline" else 0.9375
    inner_half_mm = 18.0
    plane_wave_half_mm = 12.0
    arm_mm = 1.25

    fdtd = openEMS(NrTS=40000, EndCriteria=1e-5)
    fdtd.SetGaussExcite(F_CENTER_HZ, 0.5 * (F_STOP_HZ - F_START_HZ))
    fdtd.SetBoundaryCond(["PML_8"] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    transverse = [-inner_half_mm, -local_step_mm, 0.0, local_step_mm, inner_half_mm]
    local_axial = (
        np.linspace(-arm_mm, arm_mm, 11)
        if resolution == "baseline"
        else np.linspace(-arm_mm, arm_mm, 15)
    )
    axial = np.r_[-inner_half_mm, local_axial, inner_half_mm]
    grid.AddLine("x", transverse)
    grid.AddLine("y", transverse)
    grid.AddLine("z", axial)
    grid.SmoothMeshLines("all", max_step_mm, 1.4)
    for axis in "xyz":
        lines = np.asarray(grid.GetLines(axis), dtype=float)
        left_step = lines[1] - lines[0]
        right_step = lines[-1] - lines[-2]
        grid.AddLine(
            axis,
            np.r_[lines[0] - np.arange(8, 0, -1) * left_step, lines[-1] + np.arange(1, 9) * right_step],
        )
    axes_mm = [np.asarray(grid.GetLines(axis), dtype=float) for axis in "xyz"]
    mesh = mesh_report([axis * 1e-3 for axis in axes_mm], F_STOP_HZ)

    receiver = CurvePort(
        csx,
        1,
        R=50.0,
        start=[0.0, 0.0, -arm_mm],
        stop=[0.0, 0.0, arm_mm],
        excite=1 if excitation == "port" else 0,
    )
    if excitation.startswith("plane_"):
        e_dir = [0.0, 0.0, 1.0] if excitation == "plane_parallel" else [0.0, 1.0, 0.0]
        plane = csx.AddExcitation("plane_wave", exc_type=10, exc_val=e_dir)
        plane.SetPropagationDir([1.0, 0.0, 0.0])
        plane.SetFrequency(F_CENTER_HZ)
        plane.AddBox([-plane_wave_half_mm] * 3, [plane_wave_half_mm] * 3)

    fdtd.Run(str(output), cleanup=False, numThreads=4)
    receiver.CalcPort(str(output), FREQUENCY_HZ, ref_impedance=50.0)
    result_columns = [FREQUENCY_HZ, receiver.uf_tot.real, receiver.uf_tot.imag, receiver.if_tot.real, receiver.if_tot.imag]
    header = "frequency_Hz,Vport_real,Vport_imag,Iport_real,Iport_imag"
    result = {
        "case": f"H4_{resolution.upper()}_{excitation.upper()}",
        "resolution": resolution,
        "excitation": excitation,
        "frequency_start_Hz": F_START_HZ,
        "frequency_stop_Hz": F_STOP_HZ,
        "frequency_points": int(FREQUENCY_HZ.size),
        "receiver_id": "H4_CANONICAL_SHORT_DIPOLE_FIELD_SENSOR",
        "receiver_length_m": 0.0025,
        "receiver_axis": [0.0, 0.0, 1.0],
        "receiver_load_ohm": 50.0,
        "wire_status": "OPENEMS_CURVEPORT_PHYSICAL_WIRE_RADIUS_UNDEFINED",
        "propagation_direction": [1.0, 0.0, 0.0] if excitation.startswith("plane_") else None,
        "incident_E_V_m": [0.0, 0.0, 1.0] if excitation == "plane_parallel" else ([0.0, 1.0, 0.0] if excitation == "plane_orthogonal" else None),
        **mesh,
    }
    if excitation == "port":
        s11 = receiver.uf_ref / receiver.uf_inc
        zin = receiver.uf_tot / receiver.if_tot
        result_columns.extend([s11.real, s11.imag, zin.real, zin.imag])
        header += ",S11_real,S11_imag,Zin_real_ohm,Zin_imag_ohm"
    else:
        incident = UI_data("et", str(output), freq=FREQUENCY_HZ)
        incident_spectrum = np.asarray(incident.ui_f_val[0], dtype=complex)
        if np.any(np.abs(incident_spectrum) <= np.finfo(float).tiny):
            raise RuntimeError("ZERO_CALIBRATED_INCIDENT_FIELD")
        transfer = receiver.uf_tot / incident_spectrum
        result_columns.extend([incident_spectrum.real, incident_spectrum.imag, transfer.real, transfer.imag])
        header += ",Einc_spectrum_real,Einc_spectrum_imag,H_rx_E_real_m,H_rx_E_imag_m"
    np.savetxt(output / "response.csv", np.column_stack(result_columns), delimiter=",", header=header, comments="")
    result.update(
        {
            "runtime_s": time.perf_counter() - started,
            "peak_RSS_KiB": max(
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            ),
            "raw_bytes": sum(path.stat().st_size for path in output.rglob("*") if path.is_file()),
            "response_finite": bool(np.all(np.isfinite(np.column_stack(result_columns)))),
        }
    )
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resolution", choices=("baseline", "fine"), required=True)
    parser.add_argument(
        "--excitation", choices=("port", "plane_parallel", "plane_orthogonal"), required=True
    )
    args = parser.parse_args()
    build_case(args.output, args.resolution, args.excitation)
