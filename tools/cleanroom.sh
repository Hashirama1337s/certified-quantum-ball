#!/bin/sh
# cleanroom.sh -- replay this repository from a fresh copy, using only the documented commands.
# Usage: sh tools/cleanroom.sh [quick|full|bg2] [PYTHON]     (default: quick, python)
# 1. copies the repository (without work folders or caches) into a new empty folder next to it;
# 2. runs tools/check_pins.py there, then tools/reproduce.py --tier TIER with a fresh work folder;
# 3. writes CLEANROOM-<tier>.log in that folder (the runner's JSON report is work/_out/reproduce_report.json).
# It never edits the repository it was started from.
set -eu
TIER="${1:-quick}"
PY="${2:-python}"
SRC="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${SRC}-cleanroom-${TIER}-$(date +%Y%m%d-%H%M%S)"
if [ -e "$DEST" ]; then echo "cleanroom: $DEST exists" >&2; exit 1; fi
mkdir "$DEST"
( cd "$SRC" && tar --exclude=./work --exclude=./work_full --exclude='*/__pycache__' --exclude='*.pyc' -cf - . ) | ( cd "$DEST" && tar -xf - )
cd "$DEST"
LOG="CLEANROOM-${TIER}.log"
{
  echo "cleanroom: tier $TIER, fresh copy $DEST, started $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  "$PY" --version
  "$PY" -c "import flint, sympy, mpmath, numpy; print('python-flint', flint.__version__, 'sympy', sympy.__version__, 'mpmath', mpmath.__version__, 'numpy', numpy.__version__)"
  "$PY" tools/check_pins.py | tail -n 3
  "$PY" tools/reproduce.py --tier "$TIER" --work "$DEST/work" --python "$PY"
  echo "cleanroom: finished $(date -u +%Y-%m-%dT%H:%M:%SZ); report work/_out/reproduce_report.json"
} 2>&1 | tee "$LOG"
