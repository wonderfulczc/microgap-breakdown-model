"""Non-destructive environment diagnostics for core and optional backends."""
from __future__ import annotations

import importlib.util
from importlib import metadata
import os
import shutil
import subprocess
import sys
from pathlib import Path


def _command(name, args=("--version",)):
    path = shutil.which(name)
    if not path:
        return {"status": "REQUIRED_MISSING", "path": None, "detail": "not found"}
    try:
        text = subprocess.run([path, *args], capture_output=True, text=True, timeout=5).stdout.strip()
    except (OSError, subprocess.SubprocessError) as exc:
        text = str(exc)
    return {"status": "AVAILABLE", "path": path, "detail": text.splitlines()[0] if text else "available"}


def doctor_report(repo_root: Path) -> dict:
    dependencies = {}
    for name in ("numpy", "scipy", "pandas", "matplotlib", "yaml", "h5py"):
        distribution = "PyYAML" if name == "yaml" else name
        if importlib.util.find_spec(name) is None:
            dependencies[name] = {"status": "REQUIRED_MISSING", "version": None}
        else:
            try:
                version = metadata.version(distribution)
            except metadata.PackageNotFoundError:
                version = "unknown"
            dependencies[name] = {"status": "AVAILABLE", "version": version}
    commands = {"cmake": _command("cmake"), "cxx": _command("c++"), "mpi": _command("mpirun")}
    petsc = _command("pkg-config", ("--modversion", "PETSc"))
    commands["PETSc"] = petsc
    optional = {}
    for label, variable in (("Afivo", "AFIVO_STREAMER_ROOT"), ("openEMS", "OPENEMS_ROOT")):
        value = os.environ.get(variable)
        available = bool(value and Path(value).exists())
        optional[label] = {
            "status": "AVAILABLE" if available else "OPTIONAL_MISSING",
            "environment_variable": variable,
            "path": value,
        }
    optional["COMSOL"] = {"status": "OPTIONAL_MISSING", "detail": "proprietary external interface"}
    environment = {}
    for variable in ("PETSC_DIR", "PETSC_ARCH", "AFIVO_STREAMER_ROOT", "OPENEMS_ROOT", "OPENEMS_PYTHON"):
        value = os.environ.get(variable)
        environment[variable] = {"status": "AVAILABLE" if value else "OPTIONAL_MISSING", "value": value}
    try:
        if not (repo_root / ".git").exists():
            raise FileNotFoundError("packaged resources have no Git worktree")
        commit = subprocess.check_output(["git", "-C", str(repo_root), "rev-parse", "HEAD"], text=True).strip()
        branch = subprocess.check_output(["git", "-C", str(repo_root), "branch", "--show-current"], text=True).strip()
    except (OSError, subprocess.SubprocessError):
        import json
        metadata_paths = (
            repo_root / "packaging/release_build_metadata.json",
            repo_root / "python/streamer_rf/resources/packaging/release_build_metadata.json",
        )
        metadata_path = next((path for path in metadata_paths if path.is_file()), None)
        if metadata_path is None:
            commit, branch = "UNAVAILABLE", "UNAVAILABLE"
        else:
            record = json.loads(metadata_path.read_text())
            commit, branch = f"{record['base_checkpoint']}+{record['source_state']}", "SOURCE_ARCHIVE_OR_INSTALLED_WHEEL"
    required_missing = [name for name, item in {**dependencies, **commands}.items() if item["status"] == "REQUIRED_MISSING"]
    return {
        "python": {"status": "AVAILABLE", "version": sys.version.split()[0], "executable": sys.executable},
        "core_dependencies": dependencies,
        "toolchain": commands,
        "optional_backends": optional,
        "environment_variables": environment,
        "repository": {"path": str(repo_root), "git_available": (repo_root / ".git").exists(), "commit": commit, "branch": branch},
        "required_missing": required_missing,
        "overall": "AVAILABLE" if not required_missing else "REQUIRED_MISSING",
    }
