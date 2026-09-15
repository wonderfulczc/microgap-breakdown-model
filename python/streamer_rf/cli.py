"""Unified command-line entry point for the frozen v2.0 toolkit."""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from .literature import impact_analysis, load_registry
from .release.config import CASE_ID, load_case_config
from .release.audit import release_audit
from .release.doctor import doctor_report
from .release.reporting import generate_report
from .release.runner import execute_run, execute_validation, repository_root, user_execution_root


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="microgap-rf", description="多物理微间隙放电与射频科研工具包")
    sub = root.add_subparsers(dest="command", required=True)
    sub.add_parser("doctor", help="检查核心工具链和可选外部后端")
    run = sub.add_parser("run", help="运行配置驱动的轻量工作流")
    run.add_argument("config", type=Path)
    validate = sub.add_parser("validate", help="运行 Stage-I 验证入口")
    validate.add_argument("config", type=Path)
    report = sub.add_parser("report", help="为现有 case 生成轻量报告")
    report.add_argument("case_id")
    report.add_argument("--results-root", type=Path)
    literature = sub.add_parser("literature", help="文献证据维护工具")
    literature_sub = literature.add_subparsers(dest="literature_command", required=True)
    impact = literature_sub.add_parser("impact", help="分析登记论文的潜在影响")
    impact.add_argument("paper_id")
    impact.add_argument("--registry", type=Path)
    release = sub.add_parser("release", help="只读候选发布审计")
    release_sub = release.add_subparsers(dest="release_command", required=True)
    release_sub.add_parser("audit", help="检查发布门，不创建发布")
    return root


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    try:
        repo = repository_root()
        command = ["microgap-rf", *(argv if argv is not None else sys.argv[1:])]
        if args.command == "doctor":
            report = doctor_report(repo)
            print(json.dumps(report, indent=2, ensure_ascii=False))
            return 0 if report["overall"] == "AVAILABLE" else 2
        if args.command == "run":
            output = execute_run(repo, load_case_config(args.config), command)
            print(output)
            return 0
        if args.command == "validate":
            output = execute_validation(repo, load_case_config(args.config), command)
            print(output)
            return 0
        if args.command == "report":
            if not CASE_ID.fullmatch(args.case_id):
                raise ValueError("INVALID_CASE_ID")
            base = (args.results_root or user_execution_root(repo) / "results").resolve()
            print(generate_report(base / args.case_id))
            return 0
        if args.command == "release":
            print(json.dumps(release_audit(repo), indent=2, ensure_ascii=False))
            return 0
        registry_path = args.registry or repo / "literature/literature_registry.yaml"
        print(json.dumps(impact_analysis(load_registry(registry_path), args.paper_id), indent=2, ensure_ascii=False))
        return 0
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
