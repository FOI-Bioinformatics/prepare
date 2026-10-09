#!/bin/bash
# local_test.sh: end-to-end CLI walkthrough on the bundled Branchwater sample.
# Runs in build/local_test so tracked files under test_data/ are never rewritten.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
WORK="$ROOT/build/local_test"
rm -rf "$WORK"
mkdir -p "$WORK/branchwater" "$WORK/genomes"
cp "$ROOT/test_data/branchwater/salmonella_subset.csv" "$WORK/branchwater/"
if [ -f "$ROOT/test_data/GCF_000008985.1.fasta" ]; then
    cp "$ROOT/test_data/GCF_000008985.1.fasta" "$WORK/genomes/GCF_000008985.1.fna"
fi
cd "$WORK"

check() { if [ -e "$1" ]; then echo "ok   $1"; else echo "FAIL $1 missing"; exit 1; fi; }

echo "use_branchwater"
python -m prepare.cli.main use_branchwater --branchwater-folder branchwater --matches-folder matches
check matches/salmonella_subset.csv

echo "parse_containment"
python -m prepare.cli.main parse_containment --matches-folder matches \
    --parsed-containment-file parsed_containment.txt \
    --summary-containment-file summary_containment.txt --step-size 0.1
check parsed_containment.txt
check summary_containment.txt

echo "extract_branchwater_metadata"
python -m prepare.cli.main extract_branchwater_metadata --branchwater-folder branchwater --metadata-folder metadata
check metadata/branchwater_metadata.txt

echo "count_metadata (metadata table resolved from metadata/branchwater_metadata.txt)"
python -m prepare.cli.main count_metadata --metadata-column Sample_Scientific_Name \
    --threshold 0.9 --output-file metadata_counts.txt
check metadata_counts.txt

echo "plot_metadata_counts"
python -m prepare.cli.main plot_metadata_counts --file-path metadata_counts.txt --plot-type bar --save-format png
check metadata_counts_bar.png

echo "plot_containment"
python -m prepare.cli.main plot_containment --file-path parsed_containment.txt \
    --column max_containment --plot-type rank --save-format png
check parsed_containment_rank_max_containment.png

echo "select_datasets"
python -m prepare.cli.main select_datasets --threshold 0.95 --output accessions.txt
check accessions.txt
test "$(wc -l < accessions.txt)" -gt 0

# use_branchwater/parse_containment/select_datasets already recorded screening and selection in
# prepare_registry.json as they ran. Rename it aside to simulate a project that was created
# before the registry existed, so status --init can rebuild one from disk (matches folder, parsed
# containment table, accessions file); the renamed copy is kept so the two can be compared below.
echo "status --init"
mv prepare_registry.json prepare_registry.recorded.json
python -m prepare.cli.main status --init --parsed-containment parsed_containment.txt --accessions-file accessions.txt
check prepare_registry.json

echo "blacklist"
python -m prepare.cli.main blacklist --add SRR31320538 --reason "walkthrough example"

echo "status --stage excluded"
python -m prepare.cli.main status --stage excluded

echo "status"
python -m prepare.cli.main status --parsed-containment parsed_containment.txt --list-missing

if [ -f genomes/GCF_000008985.1.fna ]; then
    echo "extract_target_reads --dry-run"
    python -m prepare.cli.main extract_target_reads \
        --parsed-containment parsed_containment.txt --genome-id salmonella_subset \
        --genome-fasta genomes/GCF_000008985.1.fna --threshold 0.95 --dry-run
fi

echo "status --json"
python -m prepare.cli.main status --json > status.json
python -c "
import json
with open('status.json') as f:
    data = json.load(f)
selected = data['stages']['selected']['count']
excluded = data['stages']['excluded']['count']
assert selected > 0, f'expected stages.selected.count > 0, got {selected}'
assert excluded == 1, f'expected stages.excluded.count == 1, got {excluded}'
with open('prepare_registry.recorded.json') as f:
    recorded = json.load(f)
with open('prepare_registry.json') as f:
    rebuilt = json.load(f)
recorded_count = len(recorded['datasets'])
rebuilt_count = len(rebuilt['datasets'])
assert recorded_count == rebuilt_count, (
    f'expected the rebuilt registry to agree with the recorded one, '
    f'got recorded={recorded_count} rebuilt={rebuilt_count}'
)
print(f'status.json ok: selected={selected} excluded={excluded}')
print(f'registry rebuild ok: recorded datasets={recorded_count} rebuilt datasets={rebuilt_count}')
" || { echo "FAIL status.json stage counts"; exit 1; }

# Shared data store walkthrough. Runs with HOME redirected under $WORK and PREPARE_DATA unset, so
# store discovery never falls back to the real user config or a real store; --data-root is explicit
# on every call instead.
STORE_HOME="$WORK/store_home"
mkdir -p "$STORE_HOME"
STORE_ROOT="$WORK/store"
N_SELECTED="$(wc -l < accessions.txt | tr -d ' ')"

echo "store_init"
(
    unset PREPARE_DATA
    export HOME="$STORE_HOME"
    python -m prepare.cli.main store_init --data-root "$STORE_ROOT"
)

echo "download_sra --dry-run --data-root (store)"
(
    unset PREPARE_DATA
    export HOME="$STORE_HOME"
    python -m prepare.cli.main download_sra --accessions-file accessions.txt --dry-run \
        --data-root "$STORE_ROOT" 2>&1 \
        | tee store_download_dry_run.log
)
grep -q "would download $N_SELECTED of $N_SELECTED datasets" store_download_dry_run.log \
    || { echo "FAIL store dry-run did not report would download $N_SELECTED of $N_SELECTED datasets"; exit 1; }

echo "store_status --json"
(
    unset PREPARE_DATA
    export HOME="$STORE_HOME"
    python -m prepare.cli.main store_status --json --data-root "$STORE_ROOT" | tee store_status.json
)
grep -q '"projects": 1' store_status.json || { echo "FAIL store_status --json: expected \"projects\": 1"; exit 1; }

echo "select_datasets --top-n 5"
python -m prepare.cli.main select_datasets --threshold 0.95 --top-n 5 --output accessions_top5.txt
check accessions_top5.txt
test "$(wc -l < accessions_top5.txt)" -eq 5 || { echo "FAIL accessions_top5.txt: expected 5 lines"; exit 1; }

echo "All steps passed (outputs in $WORK)"
