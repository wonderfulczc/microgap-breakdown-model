import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from streamer_rf.literature import impact_analysis, load_registry
from streamer_rf.release.config import resolve_case_config
from streamer_rf.release.reporting import generate_report
from streamer_rf.release.runner import execute_run, execute_validation


ROOT = Path(__file__).resolve().parents[2]


def base_config(tmp_path, case_id="test_case", workflow="smoke"):
    return {
        "case_id": case_id,
        "workflow": {"type": workflow},
        "backend": {"core": "core" if workflow != "validation" else "stage_i"},
        "output": {"directory": str(tmp_path / case_id)},
        "scientific_claims": {"native_rf_trusted_at_target": False},
    }


def invoke(*args, cwd=ROOT):
    env = {**os.environ, "PYTHONPATH": str(ROOT / "python")}
    return subprocess.run([sys.executable, "-m", "streamer_rf.cli", *args], cwd=cwd, env=env, text=True, capture_output=True)


def test_cli_help():
    result = invoke("--help")
    assert result.returncode == 0
    assert all(name in result.stdout for name in ("doctor", "run", "validate", "report"))


def test_report_rejects_path_traversal():
    result = invoke("report", "../outside")
    assert result.returncode == 2
    assert "INVALID_CASE_ID" in result.stderr


def test_doctor_reports_optional_backends_without_crashing():
    result = invoke("doctor")
    assert result.returncode in (0, 2)
    report = json.loads(result.stdout)
    assert "optional_backends" in report
    assert report["optional_backends"]["Afivo"]["status"] in {"AVAILABLE", "OPTIONAL_MISSING"}
    assert len(report["repository"]["commit"]) == 40
    assert "PETSC_DIR" in report["environment_variables"]


def test_config_validation_and_unknown_backend(tmp_path):
    resolved = resolve_case_config(base_config(tmp_path))
    assert resolved["schema_version"] == "1.0"
    assert resolved["model_version"] == "frozen-v2.0"
    invalid = base_config(tmp_path)
    invalid["backend"] = {"core": "unknown"}
    with pytest.raises(ValueError, match="UNKNOWN_BACKEND"):
        resolve_case_config(invalid)


def test_invalid_frequency_and_native_trust_are_rejected(tmp_path):
    config = base_config(tmp_path, workflow="system_rf")
    config["frequency"] = {"analysis_band_MHz": [700, 1300], "target_MHz": 1000}
    with pytest.raises(ValueError, match="200_TO_500"):
        resolve_case_config(config)
    config["frequency"] = {"analysis_band_MHz": [200, 500], "target_MHz": 350}
    config["scientific_claims"]["native_rf_trusted_at_target"] = True
    with pytest.raises(ValueError, match="NOT_RESOLVED"):
        resolve_case_config(config)


def test_live_external_backend_requires_explicit_environment(tmp_path, monkeypatch):
    config = base_config(tmp_path)
    config["backend"] = {"fullwave": "openems"}
    config["backend_options"] = {"execution_mode": "LIVE"}
    monkeypatch.delenv("OPENEMS_ROOT", raising=False)
    with pytest.raises(ValueError, match="OPENEMS_ROOT"):
        resolve_case_config(config)


def test_run_manifest_hashes_versions_and_science(tmp_path):
    source = tmp_path / "case.yaml"
    source.write_text(yaml.safe_dump(base_config(tmp_path)))
    config = resolve_case_config(yaml.safe_load(source.read_text()), config_path=source)
    output = execute_run(ROOT, config, ["microgap-rf", "run", str(source)])
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["exit_status"] == "PASS"
    assert len(manifest["input_config_hash"]) == 64
    assert len(manifest["resolved_config_hash"]) == 64
    assert len(manifest["backend_version"]["core"]) == 40
    assert (manifest["model_version"], manifest["contract_version"], manifest["schema_version"]) == ("frozen-v2.0", "rp1.1", "1.0")
    assert manifest["scientific_status"]["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"


def test_synthetic_validation_protection_and_report(tmp_path):
    source = tmp_path / "validation.yaml"
    config = base_config(tmp_path, case_id="synthetic", workflow="validation")
    config["inputs"] = {"dataset": "validation/stage_i/system_350mhz/synthetic_input"}
    config["scientific_claims"]["requested_status"] = "SYNTHETIC_DRY_RUN"
    source.write_text(yaml.safe_dump(config))
    resolved = resolve_case_config(yaml.safe_load(source.read_text()), config_path=source)
    output = execute_validation(ROOT, resolved, ["microgap-rf", "validate", str(source)])
    result = json.loads((output / "data/validation_result.json").read_text())
    assert result["scientific_validation_allowed"] is False
    assert result["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    report = generate_report(output)
    assert report.is_file() and "PENDING_REAL_EXPERIMENT" in report.read_text()


def test_forbidden_synthetic_validated_status(tmp_path):
    config = base_config(tmp_path, workflow="validation")
    config["scientific_claims"]["requested_status"] = "VALIDATED"
    with pytest.raises(ValueError, match="CANNOT_VALIDATE"):
        resolve_case_config(config)


def test_literature_registry_and_impact_never_change_code():
    registry = load_registry(ROOT / "literature/literature_registry.yaml")
    result = impact_analysis(registry, "DANGOLA_2008_EQUILIBRIUM_AIR")
    assert result["recommendation"] == "NO_CHANGE"
    assert result["automatic_code_change"] is False


def test_evidence_and_change_request_schemas_are_protected_templates():
    evidence = yaml.safe_load((ROOT / "literature/evidence_cards/EVIDENCE_CARD_TEMPLATE.zh-CN.yaml").read_text())
    change = yaml.safe_load((ROOT / "literature/change_requests/CR-TEMPLATE.yaml").read_text())
    assert evidence["modification_recommendation"] == "PENDING_EXPERT_REVIEW"
    assert change["automatic_code_change_permitted"] is False
    assert change["decision_status"] == "DRAFT"


def test_release_status_preserves_all_scientific_boundaries():
    status = json.loads((ROOT / "packaging/release_preparation_status.json").read_text())
    assert status["UNIFIED_CLI_READY"] is True
    assert status["CONFIG_DRIVEN_WORKFLOW_READY"] is True
    assert status["STAGE_I_SCIENTIFIC_VALIDATION"] == "PENDING_REAL_EXPERIMENT"
    assert status["PUBLIC_SCIENTIFIC_VALIDATION_COMPLETE"] is False
    assert status["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"
    assert status["release_created"] is False


def test_critical_project_topics_have_authoritative_chinese_docs():
    required = {
        "软件架构.md", "CLI架构说明.md", "科学状态.md", "数据政策.md",
        "实验数据重入.md", "文献维护流程.md", "科研模型变更流程.md", "Release准备说明.md",
    }
    assert required <= {path.name for path in (ROOT / "docs/zh").glob("*.md")}
    assert (ROOT / "README.zh-CN.md").is_file()


def test_release_preparation_manifest_has_no_physics_or_release_side_effects():
    manifest = json.loads((ROOT / "packaging/release_preparation_manifest.json").read_text())
    assert manifest["base_stage_j_commit"] == "825c65c"
    assert manifest["scientific_model_changed"] is False
    assert manifest["frozen_results_reformatted"] is False
    assert manifest["release_created"] is False
    assert manifest["scientific_status"]["NATIVE_RF_350MHZ"] == "NOT_RESOLVED"


def test_pyproject_includes_runtime_property_data_and_console_entry():
    text = (ROOT / "pyproject.toml").read_text()
    assert 'microgap-rf = "streamer_rf.cli:main"' in text
    assert 'streamer_rf = ["thermal/*.json"]' in text
