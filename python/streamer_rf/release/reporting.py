"""Lightweight Markdown reporting for unified run directories."""
from __future__ import annotations

import json
from pathlib import Path


def generate_report(case_dir: Path) -> Path:
    case_dir = case_dir.resolve()
    manifest_path = case_dir / "run_manifest.json"
    status_path = case_dir / "status.json"
    if not manifest_path.is_file() or not status_path.is_file():
        raise ValueError("CASE_MANIFEST_OR_STATUS_NOT_FOUND")
    manifest = json.loads(manifest_path.read_text())
    status = json.loads(status_path.read_text())
    lines = [
        f"# Case {manifest['case_id']}", "",
        "## 运行信息", "",
        f"- 后端：`{json.dumps(manifest['backend'], ensure_ascii=False)}`",
        f"- 退出状态：`{manifest['exit_status']}`",
        f"- 运行时间：`{manifest['runtime_s']:.6f} s`", "",
        "## 科学状态", "",
    ]
    lines.extend(f"- `{key}`：`{value}`" for key, value in status["scientific_status"].items())
    lines += ["", "## 来源追踪", "", f"- Git：`{manifest['git_commit']}`", f"- 输入配置哈希：`{manifest['input_config_hash']}`", f"- 解析配置哈希：`{manifest['resolved_config_hash']}`", "", "## 已知限制", "", "- 本报告由统一 CLI 生成，不改变冻结科学结论。", "- 合成验证结果不能作为实验验证。", ""]
    target = case_dir / "report/report.md"
    target.parent.mkdir(exist_ok=True)
    target.write_text("\n".join(lines), encoding="utf-8")
    return target

