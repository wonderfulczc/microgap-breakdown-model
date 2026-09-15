"""Thin orchestration over frozen project modules and contracts."""
from __future__ import annotations

import hashlib
import json
import shlex
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

from streamer_rf.fullwave.transient import H3_BAND_HZ
from streamer_rf.thermal.port import PortTransientContract
from streamer_rf.validation.wp_i_b import validate_vna_input, verify_synthetic_bundle


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repository_root(start: Path | None = None) -> Path:
    current = (start or Path.cwd()).resolve()
    for candidate in (current, *current.parents):
        if (candidate / "packaging/final_project_contract.json").is_file():
            return candidate
    package_root = Path(__file__).resolve().parents[3]
    if (package_root / "packaging/final_project_contract.json").is_file():
        return package_root
    raise ValueError("REPOSITORY_ROOT_NOT_FOUND")


def _git_commit(root: Path) -> str:
    return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()


def _output_dir(root: Path, cfg: dict) -> Path:
    configured = Path(cfg["output"]["directory"])
    path = configured if configured.is_absolute() else root / configured
    path = path.resolve()
    if path.exists():
        raise ValueError(f"CASE_OUTPUT_ALREADY_EXISTS:{path}")
    return path


def _prepare(root: Path, cfg: dict, command: list[str]):
    output = _output_dir(root, cfg)
    for name in ("logs", "data", "figures", "report"):
        (output / name).mkdir(parents=True, exist_ok=True)
    input_path = Path(cfg["_config_path"])
    input_bytes = input_path.read_bytes() if input_path.is_file() else yaml.safe_dump(cfg).encode()
    (output / "config_input.yaml").write_bytes(input_bytes)
    resolved = {key: value for key, value in cfg.items() if not key.startswith("_")}
    (output / "config_resolved.yaml").write_text(yaml.safe_dump(resolved, sort_keys=False), encoding="utf-8")
    commit = _git_commit(root)
    external = json.loads((root / "packaging/external_backends.json").read_text())
    versions = {}
    for logical_name, backend_name in cfg.get("backend", {}).items():
        if backend_name == "openems":
            versions[logical_name] = external["openEMS"]["commits"]
        elif backend_name == "afivo":
            versions[logical_name] = external["Afivo"]["commit"]
        elif backend_name == "petsc":
            versions[logical_name] = "PETSc_3.24.4_REFERENCE_ENVIRONMENT"
        else:
            versions[logical_name] = commit
    manifest = {
        "case_id": cfg["case_id"],
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": commit,
        "command": shlex.join(command),
        "input_config_hash": hashlib.sha256(input_bytes).hexdigest(),
        "resolved_config_hash": sha256(output / "config_resolved.yaml"),
        "backend": cfg.get("backend", {}),
        "backend_version": versions,
        "model_version": cfg["model_version"],
        "contract_version": cfg["contract_version"],
        "schema_version": cfg["schema_version"],
        "input_hashes": {},
        "output_hashes": {},
        "scientific_status": cfg["scientific_status"],
        "runtime_s": None,
        "exit_status": "RUNNING",
    }
    return output, manifest


def _architecture_smoke(root: Path) -> dict:
    final = json.loads((root / "packaging/final_project_contract.json").read_text())
    # Importing these frozen modules verifies the principal Python handoffs.
    assert H3_BAND_HZ == (200e6, 500e6)
    assert PortTransientContract is not None
    return {
        "stage_j": final["stages"]["J"],
        "h3_band_Hz": list(H3_BAND_HZ),
        "physics_simulation_executed": False,
        "status": "PASS",
    }


def execute_run(root: Path, cfg: dict, command: list[str]) -> Path:
    output, manifest = _prepare(root, cfg, command)
    started = time.perf_counter()
    try:
        kind = cfg["workflow"]["type"]
        if kind == "validation":
            raise ValueError("USE_VALIDATE_COMMAND_FOR_VALIDATION_WORKFLOW")
        result = _architecture_smoke(root)
        if kind == "system_rf":
            transfer = root / cfg.get("inputs", {}).get("h3_contract", "fullwave/h3/h3_result_contract.json")
            if not transfer.is_file():
                raise ValueError("H3_CONTRACT_NOT_FOUND")
            manifest["input_hashes"][str(transfer.relative_to(root))] = sha256(transfer)
            result["system_band_Hz"] = [200e6, 500e6]
            result["target_Hz"] = 350e6
            result["H3_contract"] = str(transfer.relative_to(root))
        (output / "data/result.json").write_text(json.dumps(result, indent=2) + "\n")
        manifest["output_hashes"]["data/result.json"] = sha256(output / "data/result.json")
        manifest["exit_status"] = "PASS"
    except Exception:
        manifest["exit_status"] = "FAILED"
        raise
    finally:
        manifest["runtime_s"] = time.perf_counter() - started
        (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        (output / "provenance.json").write_text(json.dumps({"git_commit": manifest["git_commit"], "input_hashes": manifest["input_hashes"], "output_hashes": manifest["output_hashes"]}, indent=2) + "\n")
        (output / "status.json").write_text(json.dumps({"exit_status": manifest["exit_status"], "scientific_status": manifest["scientific_status"]}, indent=2) + "\n")
    return output


def execute_validation(root: Path, cfg: dict, command: list[str]) -> Path:
    if cfg["workflow"]["type"] != "validation":
        raise ValueError("VALIDATE_REQUIRES_VALIDATION_WORKFLOW")
    output, manifest = _prepare(root, cfg, command)
    started = time.perf_counter()
    try:
        relative = cfg.get("inputs", {}).get("dataset")
        if not relative:
            raise ValueError("VALIDATION_DATASET_REQUIRED")
        dataset = (root / relative).resolve() if not Path(relative).is_absolute() else Path(relative).resolve()
        bundle = verify_synthetic_bundle(dataset)
        for name in bundle["metadata"]["vna"]["files"]:
            validate_vna_input(dataset / name)
        manifest["input_hashes"]["sha256_manifest.json"] = sha256(dataset / "sha256_manifest.json")
        result = {
            "validation_mode": "SYNTHETIC_DRY_RUN",
            "verified_file_count": bundle["verified_file_count"],
            "scientific_validation_allowed": False,
            "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
            "status": "PASS",
        }
        (output / "data/validation_result.json").write_text(json.dumps(result, indent=2) + "\n")
        manifest["output_hashes"]["data/validation_result.json"] = sha256(output / "data/validation_result.json")
        manifest["exit_status"] = "PASS"
    except Exception:
        manifest["exit_status"] = "FAILED"
        raise
    finally:
        manifest["runtime_s"] = time.perf_counter() - started
        (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        (output / "provenance.json").write_text(json.dumps({"git_commit": manifest["git_commit"], "input_hashes": manifest["input_hashes"], "output_hashes": manifest["output_hashes"]}, indent=2) + "\n")
        (output / "status.json").write_text(json.dumps({"exit_status": manifest["exit_status"], "scientific_status": manifest["scientific_status"], "evidence": "SYNTHETIC_DRY_RUN"}, indent=2) + "\n")
    return output
