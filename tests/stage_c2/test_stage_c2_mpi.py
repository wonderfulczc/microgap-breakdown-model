from pathlib import Path
import os
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
APP = ROOT / "build/bin/stage_c2_dynamic"


def read_summary(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k] = v
    return out


def run_case(out: Path, ranks: int) -> tuple[pd.DataFrame, dict[str, str]]:
    if not APP.exists():
        pytest.skip("stage_c2_dynamic is not built")
    env = os.environ.copy()
    env.setdefault("OMPI_MCA_btl", "self,vader,tcp")
    cmd = [
        str(APP),
        str(out),
        "--case",
        f"stage-c2-mpi-{ranks}",
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
        "--n0",
        "1e16",
        "--sigma",
        "3e-6",
        "--seed-offset",
        "-10e-6",
    ]
    if ranks > 1:
        mpirun = shutil.which("mpirun") or shutil.which("mpiexec")
        if mpirun is None:
            pytest.skip("mpirun/mpiexec is not available")
        cmd = [mpirun, "-np", str(ranks), *cmd]
    result = subprocess.run(cmd, cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if ranks > 1 and result.returncode != 0 and "no available" in result.stderr and "interfaces were found" in result.stderr:
        pytest.skip("mpirun cannot initialize an MPI interface in this sandbox")
    assert result.returncode == 0, result.stdout + result.stderr
    return pd.read_csv(out / "diagnostics.csv"), read_summary(out / "summary.txt")


def test_stage_c2_mpi_timestep_consistency(tmp_path):
    one, s1 = run_case(tmp_path / "rank1", 1)
    two, s2 = run_case(tmp_path / "rank2", 2)

    assert len(one) == len(two) == 20
    assert one["dt_controller"].tolist() == two["dt_controller"].tolist()
    assert np.allclose(one["dt"].to_numpy(), two["dt"].to_numpy(), rtol=1e-12, atol=0.0)

    for key in ["final_physical_time", "Emax", "ne_max", "total_electrons"]:
        assert np.isclose(float(s1[key]), float(s2[key]), rtol=1e-10, atol=1e-24), key
