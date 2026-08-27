from pathlib import Path
import os
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
APP = ROOT / "build/bin/stage_c3_current"


def read_summary(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k] = v
    return out


def run_current(out: Path, ranks: int):
    if not APP.exists():
        pytest.skip("stage_c3_current is not built")
    cmd = [
        str(APP),
        str(out),
        "--case",
        f"stage-c3-mpi-{ranks}",
        "--mode",
        "case-b",
        "--voltage",
        "500",
        "--sp3",
        "0",
        "--steps",
        "20",
        "--dt-scale",
        "0.001",
        "--dt-cap",
        "1e-16",
    ]
    if ranks > 1:
        mpirun = shutil.which("mpirun") or shutil.which("mpiexec")
        if mpirun is None:
            pytest.skip("mpirun/mpiexec is not available")
        cmd = [mpirun, "-np", str(ranks), *cmd]
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    result = subprocess.run(cmd, cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if ranks > 1 and result.returncode != 0 and "no available" in result.stderr and "interfaces were found" in result.stderr:
        pytest.skip("mpirun cannot initialize an MPI interface in this sandbox")
    assert result.returncode == 0, result.stdout + result.stderr
    return pd.read_csv(out / "current_diagnostics.csv"), read_summary(out / "summary.txt")


def test_stage_c3_current_mpi_consistency(tmp_path):
    one, s1 = run_current(tmp_path / "rank1", 1)
    two, s2 = run_current(tmp_path / "rank2", 2)
    tolerances = {
        "I_cond_HV_A": (1e-10, 1e-18),
        "I_disp_HV_A": (1e-10, 2e-7),
        "I_total_HV_A": (1e-10, 2e-7),
        "Q_HV_C": (1e-10, 2e-23),
        "Gb_S": (1e-10, 1e-20),
        "Rb_ohm": (1e-10, 10.0),
    }
    assert len(one) == len(two) == 20
    for col, (rtol, atol) in tolerances.items():
        assert np.allclose(one[col].to_numpy(), two[col].to_numpy(), rtol=rtol, atol=atol), col
    summary_tolerances = {
        "I_cond_HV_max_A": (1e-10, 1e-18),
        "I_disp_HV_max_A": (1e-10, 2e-7),
        "Gb_max_S": (1e-10, 1e-20),
    }
    for key, (rtol, atol) in summary_tolerances.items():
        assert np.isclose(float(s1[key]), float(s2[key]), rtol=rtol, atol=atol), key
