#!/usr/bin/env bash
# Analysis runner used by generated run_tcr_sample.sh files.
# Activate the scNanoGPS conda environment before running this script.
set -euo pipefail

TCR_REPO="${TCR_REPO:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)}"
FASTQ_DIR="${FASTQ_DIR:-/data/nanopore/2026/NP359_TCR_Enrich_CGF021-022_PBMC_OCM_Prom114_WANG_07022026/20260702_1712_2A_PBM45682_e171202b/fastq_pass}"
WORK_DIR="${WORK_DIR:-/data/Extra_data_USERs/Tim/HR_transplant/CGF/tcr_enrichment/NP359}"
WHITELIST="${WHITELIST:-$WORK_DIR/reference/whitelist_ob.txt.gz}"
SCANNER_DIR="${SCANNER_DIR:-$WORK_DIR/scanner}"
EXPORT_DIR="${EXPORT_DIR:-$WORK_DIR/export}"
NCORES="${NCORES:-8}"
BARCODE_SUBSTITUTIONS="${BARCODE_SUBSTITUTIONS:-0}"
PYTHON_BIN="${PYTHON_BIN:-python}"

usage() {
    cat <<'HELP'
Usage: bash run_tcr_sample.sh [all|scan|export]
  all      Run Scanner, then export using the barcode whitelist (default).
  scan     Run only Scanner, writing into a new scanner directory.
  export   Export existing Scanner output using the barcode whitelist.

Activate the scNanoGPS conda environment first. Edit paths and settings at the
top of the generated run_tcr_sample.sh. Supply either an OB or 10x
whitelist with WHITELIST. It has one 16-base barcode per line and no header;
gzip is supported. The same commands work for either list. There is no OB mode
or sample assignment, and no R/OB.rda/sample-map input is consumed.
All analysis outputs go under WORK_DIR by default. Raw FASTQs are read-only.
This wrapper does not split concatenated reads or run UMI consensus/MiXCR.
HELP
}

die() { printf 'Error: %s\n' "$*" >&2; exit 1; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

require_whitelist() {
    [[ -f "$WHITELIST" && -r "$WHITELIST" && -s "$WHITELIST" ]] ||
        die "Cannot read whitelist: $WHITELIST. Copy the whitelist here or set WHITELIST=/path/to/barcodes.txt.gz."
    printf 'Barcode whitelist: %s\n' "$WHITELIST"
}

run_scanner() {
    [[ -d "$FASTQ_DIR" ]] || die "FASTQ directory not found: $FASTQ_DIR"
    [[ ! -e "$SCANNER_DIR" ]] || die "Scanner directory already exists: $SCANNER_DIR. Use 'export' for completed Scanner output or choose a new SCANNER_DIR."
    [[ -r "$TCR_REPO/scanner.py" ]] || die "Cannot find scanner.py in $TCR_REPO"
    [[ "$NCORES" =~ ^[1-9][0-9]*$ ]] || die "NCORES must be a positive integer"
    require_command "$PYTHON_BIN"
    "$PYTHON_BIN" -c 'import Bio, numpy'
    mkdir -p "$SCANNER_DIR" "$WORK_DIR/logs"
    printf 'Input FASTQs: %s\nScanner output: %s\n' "$FASTQ_DIR" "$SCANNER_DIR"
    printf '%s\n' 'Note: concatenated reads are not split; the exporter quarantines detected internal Read 1 adapters.'
    (
        cd "$WORK_DIR"
        "$PYTHON_BIN" "$TCR_REPO/scanner.py" \
            -i "$FASTQ_DIR" \
            -d "$SCANNER_DIR" \
            --log "$WORK_DIR/logs/scanner.log.txt" \
            --lCB 16 --lUMI 12 \
            -t "$NCORES"
    ) 2>&1 | tee -a "$WORK_DIR/logs/scanner.console.log"
}

export_fastqs() {
    require_whitelist
    [[ ! -e "$EXPORT_DIR" ]] || die "Export directory already exists: $EXPORT_DIR. Choose a new EXPORT_DIR."
    [[ -r "$SCANNER_DIR/processed.fastq.gz" && -r "$SCANNER_DIR/barcode_list.tsv.gz" ]] ||
        die "Paired Scanner outputs not found in $SCANNER_DIR"
    [[ "$BARCODE_SUBSTITUTIONS" == 0 || "$BARCODE_SUBSTITUTIONS" == 1 ]] || die "BARCODE_SUBSTITUTIONS must be 0 or 1"
    [[ -r "$TCR_REPO/tcr_export.py" ]] || die "Cannot find tcr_export.py in $TCR_REPO"
    require_command "$PYTHON_BIN"
    mkdir -p "$WORK_DIR/logs"
    "$PYTHON_BIN" "$TCR_REPO/tcr_export.py" \
        --scanner-dir "$SCANNER_DIR" \
        --whitelist "$WHITELIST" \
        --max-barcode-substitutions "$BARCODE_SUBSTITUTIONS" \
        --output "$EXPORT_DIR" \
        2>&1 | tee -a "$WORK_DIR/logs/export.console.log"
}

[[ $# -le 1 ]] || { usage >&2; exit 2; }
case "${1:-all}" in
    -h|--help) usage ;;
    scan) run_scanner ;;
    export) export_fastqs ;;
    all)
        [[ ! -e "$SCANNER_DIR" ]] || die "Scanner directory already exists. Use 'export' or choose a new SCANNER_DIR."
        [[ ! -e "$EXPORT_DIR" ]] || die "Export directory already exists. Choose a new EXPORT_DIR."
        [[ -d "$FASTQ_DIR" ]] || die "FASTQ directory not found: $FASTQ_DIR"
        require_whitelist
        run_scanner
        export_fastqs
        ;;
    *) usage >&2; exit 2 ;;
esac
