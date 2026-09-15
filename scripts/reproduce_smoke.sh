#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT/build-smoke"}
PYTHON_BIN=${PYTHON_BIN:-"$ROOT/.venv/bin/python"}

if [[ ! -x "$PYTHON_BIN" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN=$(command -v python3)
  else
    printf 'ERROR: no Python interpreter found; set PYTHON_BIN\n' >&2
    exit 2
  fi
fi

export PYTHONPATH="$ROOT/python${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR=${MPLCONFIGDIR:-"$BUILD_DIR/matplotlib"}
mkdir -p "$MPLCONFIGDIR"

"$PYTHON_BIN" -c "import numpy, scipy, pandas, h5py, yaml, streamer_rf"
cmake -S "$ROOT" -B "$BUILD_DIR" -G "${CMAKE_GENERATOR:-Ninja}"
cmake --build "$BUILD_DIR" --parallel
ctest --test-dir "$BUILD_DIR" --output-on-failure
"$PYTHON_BIN" -m pytest -q \
  "$ROOT/tests/rf/test_stage_f2_jefimenko.py" \
  "$ROOT/tests/thermal/test_stage_g3_port.py" \
  "$ROOT/tests/fullwave/test_h3_transient.py" \
  "$ROOT/tests/validation/test_stage_i_wp_e.py" \
  "$ROOT/tests/release/test_release_preparation_cli.py"

printf 'Stage-J/RP-1 smoke reproduction passed.\n'
