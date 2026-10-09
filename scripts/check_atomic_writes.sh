#!/usr/bin/env bash
# Fail on a direct, non-atomic file write under prepare/.
#
# Every output file and every shared state file is written through the helpers in
# prepare/data/file_io.py (write_text_atomic, write_bytes_atomic, write_csv, open_atomic,
# atomic_path): the content goes to a unique temporary name next to the target and is moved into
# place with one rename, so a reader, a second process or an interrupted run never sees a partial
# file. This gate flags the direct forms that bypass them:
#   .write_text(   .write_bytes(   .to_csv(<path> ...)   copy2(   copyfile(   write_html(   savefig(
#   open(<path>, "w..."/"x...")   <path>.open("w..."/"x...")   open(<path>, <mode variable>)
# A to_csv call without a path (to_csv(), to_csv(sep=...)) returns a string and is allowed, as are
# read and append modes. Comment lines are ignored. A module that must write directly (the helpers
# themselves, an append-only log) is listed in scripts/atomic_writes_allowlist.txt, one path per
# line followed by the reason; an entry whose file no longer exists, or that has no reason, fails
# the gate so the list only shrinks with the code.
# Limitations: plain grep, one line at a time, so a call split across lines is not seen, and an
# open( whose path argument itself holds a parenthesis (open(Path(x), "w")) is not seen either.
# Usage: scripts/check_atomic_writes.sh [ROOT]   (ROOT defaults to the repository root)
set -euo pipefail

cd "${1:-$(dirname "$0")/..}"

allowlist_file="scripts/atomic_writes_allowlist.txt"
allowed=""
problems=""
if [ -f "$allowlist_file" ]; then
    while IFS= read -r line || [ -n "$line" ]; do
        case "$line" in "" | "#"*) continue ;; esac
        path="${line%%[[:space:]]*}"
        reason="${line#"$path"}"
        reason="${reason#"${reason%%[![:space:]]*}"}"
        if [ ! -f "$path" ]; then
            problems="${problems}${allowlist_file}: ${path} does not exist; remove the entry
"
        elif [ -z "$reason" ]; then
            problems="${problems}${allowlist_file}: ${path} has no reason
"
        fi
        allowed="${allowed}${path}
"
    done < "$allowlist_file"
fi

quote="[\"']"
mode_kw='(mode[[:space:]]*=[[:space:]]*)?'
direct_write='\.write_text\(|\.write_bytes\(|\.to_csv\(|copy2\(|copyfile\(|write_html\(|savefig\('
direct_write="${direct_write}|open\\([^)]*,[[:space:]]*${mode_kw}${quote}[wx]"
direct_write="${direct_write}|\\.open\\([[:space:]]*${mode_kw}${quote}[wx]"
# open(<path>, mode) with the mode in a variable; os.open, gzip.open and the like are not matched.
direct_write="${direct_write}|(^|[^.A-Za-z_])open\\([^),]*,[[:space:]]*${mode_kw}[A-Za-z_][A-Za-z_0-9]*[[:space:]]*[,)]"
# A to_csv call whose first argument is a keyword other than path_or_buf, or that has none, writes
# no file: it returns the CSV text.
to_csv_no_path='\.to_csv\([[:space:]]*(\)|[A-Za-z_]+[[:space:]]*=)'
to_csv_path_kw='\.to_csv\([[:space:]]*path_or_buf[[:space:]]*='

checked=0
hits=""
while IFS= read -r file; do
    checked=$((checked + 1))
    if printf '%s' "$allowed" | LC_ALL=C grep -qxF "$file"; then
        continue
    fi
    match=$(LC_ALL=C grep -nE "$direct_write" "$file" | LC_ALL=C grep -vE '^[0-9]+:[[:space:]]*#' || true)
    if [ -n "$match" ]; then
        match=$(printf '%s\n' "$match" | while IFS= read -r hit; do
            if printf '%s' "$hit" | LC_ALL=C grep -qE "$to_csv_no_path" \
                && ! printf '%s' "$hit" | LC_ALL=C grep -qE "$to_csv_path_kw" \
                && ! printf '%s' "$hit" | LC_ALL=C grep -qE '\.write_text\(|\.write_bytes\(|open\('; then
                continue
            fi
            printf '%s\n' "$hit"
        done)
    fi
    if [ -n "$match" ]; then
        hits="${hits}${file}:
${match}
"
    fi
done < <(find prepare -name '*.py' -type f | LC_ALL=C sort)

if [ "$checked" -eq 0 ]; then
    echo "ERROR: no Python files found under $(pwd)/prepare; nothing was checked"
    exit 1
fi

status=0
if [ -n "$problems" ]; then
    printf '%s' "$problems"
    status=1
fi
if [ -n "$hits" ]; then
    printf '%s' "$hits"
    echo "ERROR: direct file write(s) found; use write_text_atomic, write_csv, open_atomic or atomic_path"
    echo "       from prepare/data/file_io.py, or list the module in ${allowlist_file} with a reason"
    status=1
fi
if [ "$status" -eq 0 ]; then
    echo "No direct file writes outside ${allowlist_file}"
fi
exit "$status"
