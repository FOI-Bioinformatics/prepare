#!/usr/bin/env bash
# Fail when print( or sys.stdout.write appears anywhere under prepare/ except prepare/cli/base.py.
#
# stdout carries a command's result and stderr carries logging. Commands write their
# result through BaseCommand.emit / emit_raw / emit_json (or emit_error_json); library modules log.
# Matches print( as a call, not as part of a longer name (pprint(, _print_table( etc.).
# Usage: scripts/check_no_print.sh [ROOT]   (ROOT defaults to the repository root)
set -euo pipefail

cd "${1:-$(dirname "$0")/..}"

hits=$(grep -rnE --include='*.py' '(^|[^.a-zA-Z0-9_])print\(|sys\.stdout\.write' prepare | grep -v '^prepare/cli/base\.py:' || true)

if [ -n "$hits" ]; then
    echo "$hits"
    echo "ERROR: print( or sys.stdout.write outside prepare/cli/base.py; use self.emit / emit_raw / emit_json in commands and logging elsewhere"
    exit 1
fi
echo "No print( or sys.stdout.write calls outside prepare/cli/base.py"
