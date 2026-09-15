"""Validation and resolution of the compact RP-1 case configuration."""
from __future__ import annotations

import re
import os
from copy import deepcopy
from pathlib import Path

import yaml

from . import CONTRACT_VERSION, MODEL_VERSION, SCHEMA_VERSION


WORKFLOWS = {"smoke", "system_rf", "validation"}
BACKENDS = {"core", "petsc", "openems", "afivo", "stage_i"}
EXTERNAL_BACKENDS = {"openems": "OPENEMS_ROOT", "afivo": "AFIVO_STREAMER_ROOT"}
CASE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$")


def load_case_config(path: str | Path) -> dict:
    path = Path(path).resolve()
    if not path.is_file():
        raise ValueError(f"CONFIG_NOT_FOUND:{path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("CONFIG_MUST_BE_A_MAPPING")
    return resolve_case_config(data, config_path=path)


def resolve_case_config(data: dict, *, config_path: Path | None = None) -> dict:
    cfg = deepcopy(data)
    if not CASE_ID.fullmatch(str(cfg.get("case_id", ""))):
        raise ValueError("INVALID_OR_MISSING_CASE_ID")
    workflow = cfg.get("workflow")
    if not isinstance(workflow, dict) or workflow.get("type") not in WORKFLOWS:
        raise ValueError("UNKNOWN_OR_MISSING_WORKFLOW")
    backend = cfg.get("backend", {})
    if not isinstance(backend, dict):
        raise ValueError("BACKEND_MUST_BE_A_MAPPING")
    unknown = set(backend.values()) - BACKENDS
    if unknown:
        raise ValueError(f"UNKNOWN_BACKEND:{sorted(unknown)[0]}")
    execution_mode = cfg.get("backend_options", {}).get("execution_mode", "FROZEN_RESULT_REUSE")
    if execution_mode not in {"FROZEN_RESULT_REUSE", "LIVE"}:
        raise ValueError("UNKNOWN_BACKEND_EXECUTION_MODE")
    if execution_mode == "LIVE":
        for name in set(backend.values()) & EXTERNAL_BACKENDS.keys():
            variable = EXTERNAL_BACKENDS[name]
            if not os.environ.get(variable):
                raise ValueError(f"BACKEND_ENVIRONMENT_REQUIRED:{variable}")
    if workflow["type"] == "validation" and "stage_i" not in backend.values():
        raise ValueError("VALIDATION_WORKFLOW_REQUIRES_STAGE_I_BACKEND")

    frequency = cfg.get("frequency")
    if frequency is not None:
        if not isinstance(frequency, dict):
            raise ValueError("FREQUENCY_MUST_BE_A_MAPPING")
        band = frequency.get("analysis_band_MHz")
        target = frequency.get("target_MHz")
        if not isinstance(band, list) or len(band) != 2:
            raise ValueError("INVALID_ANALYSIS_BAND")
        low, high = map(float, band)
        if not 0 < low < high or target is None or not low <= float(target) <= high:
            raise ValueError("INVALID_FREQUENCY_RANGE")
        if workflow["type"] == "system_rf" and (low, high) != (200.0, 500.0):
            raise ValueError("SYSTEM_RF_FORMAL_BAND_MUST_BE_200_TO_500_MHZ")

    claims = cfg.get("scientific_claims", {})
    if claims.get("native_rf_trusted_at_target") is True:
        raise ValueError("NATIVE_RF_350MHZ_TRUST_NOT_RESOLVED")
    forbidden = {"VALIDATED", "VALIDATED_WITH_VNA", "SYSTEM_350MHZ_VALIDATED", "MEASURED"}
    if workflow["type"] == "validation" and claims.get("requested_status") in forbidden:
        raise ValueError("SYNTHETIC_OR_UNVERIFIED_INPUT_CANNOT_VALIDATE_SCIENCE")

    cfg.setdefault("schema_version", SCHEMA_VERSION)
    cfg.setdefault("contract_version", CONTRACT_VERSION)
    cfg.setdefault("model_version", MODEL_VERSION)
    if cfg["schema_version"] != SCHEMA_VERSION:
        raise ValueError("UNSUPPORTED_SCHEMA_VERSION")
    cfg.setdefault("output", {})
    cfg["output"].setdefault("directory", f"results/{cfg['case_id']}")
    cfg["scientific_status"] = {
        "STAGE_I_SCIENTIFIC_VALIDATION": "PENDING_REAL_EXPERIMENT",
        "PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE": False,
        "NATIVE_RF_350MHZ": "NOT_RESOLVED",
    }
    cfg["_config_path"] = str(config_path) if config_path else "IN_MEMORY"
    return cfg
