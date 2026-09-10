"""Run or hash-check the three bounded H2 canonical configurations."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


CASES = {
    "baseline": (5.0, 1.0),
    "fine": (3.75, 1.0),
    "expanded_domain": (5.0, 4.0 / 3.0),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--case", choices=[*CASES, "all"], default="all")
    args = parser.parse_args()
    args.raw_root.mkdir(parents=True, exist_ok=True)
    runner = Path(__file__).with_name("run_h2_reference.py")
    selected = CASES if args.case == "all" else [args.case]
    for case in selected:
        result = args.raw_root / case / "result.json"
        if result.exists():
            stored = json.loads(result.read_text())
            if stored["script_sha256"] != hashlib.sha256(runner.read_bytes()).hexdigest():
                raise ValueError("EXISTING_RESULT_CODE_HASH_DIFFERS")
            print("Reusing", case, flush=True)
            continue
        spacing, domain_scale = CASES[case]
        with (args.raw_root / f"{case}.log").open("x") as log:
            subprocess.run(
                [
                    sys.executable, str(runner), "--output", str(args.raw_root / case),
                    "--case", case, "--spacing-mm", str(spacing),
                    "--domain-scale", str(domain_scale),
                ],
                stdout=log, stderr=subprocess.STDOUT, check=True, timeout=900,
            )
        print(case, result, flush=True)


if __name__ == "__main__":
    main()
