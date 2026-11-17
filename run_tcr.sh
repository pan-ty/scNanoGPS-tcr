#!/usr/bin/env bash
# Generate an editable sample run script in the caller's working directory.
set -euo pipefail

usage() {
    cat <<'HELP'
Usage: bash /path/to/scNanoGPS-tcr/run_tcr.sh [--output run_file.sh]

Create ./run_tcr_sample.sh with editable paths at the top. The generator does
not start analysis or require an activated conda environment. Edit FASTQ_DIR,
WORK_DIR, and WHITELIST in the generated file, activate the scNanoGPS conda
environment, then run: bash run_tcr_sample.sh [all|scan|export]

Sample paths default to empty. FASTQ_DIR, WORK_DIR, WHITELIST, GPS_REPO, NCORES,
BARCODE_SUBSTITUTIONS, and PYTHON_BIN may also be set when generating the file.
Existing files are never overwritten; use --output to choose another name.
HELP
}

OUTPUT="$PWD/run_tcr_sample.sh"
case "$#" in
    0) ;;
    1)
        if [[ "$1" == -h || "$1" == --help ]]; then usage; exit 0; fi
        usage >&2; exit 2
        ;;
    2)
        [[ "$1" == --output && -n "$2" ]] || { usage >&2; exit 2; }
        OUTPUT="$2"
        ;;
    *) usage >&2; exit 2 ;;
esac

[[ ! -e "$OUTPUT" && ! -L "$OUTPUT" ]] || {
    printf 'Error: %s already exists; choose another --output filename.\n' "$OUTPUT" >&2
    exit 1
}

GPS_REPO="${GPS_REPO:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)}"
FASTQ_DIR="${FASTQ_DIR:-}"
WORK_DIR="${WORK_DIR:-}"
WHITELIST="${WHITELIST:-}"
NCORES="${NCORES:-8}"
BARCODE_SUBSTITUTIONS="${BARCODE_SUBSTITUTIONS:-0}"
PYTHON_BIN="${PYTHON_BIN:-python}"

# %q preserves literal paths, including spaces and shell metacharacters.
(
    set -o noclobber
    {
        printf '%s\n' '#!/usr/bin/env bash' '# Edit these paths for this sequencing library.'
        printf 'FASTQ_DIR=%q\n' "$FASTQ_DIR"
        printf 'WORK_DIR=%q\n' "$WORK_DIR"
        printf 'WHITELIST=%q\n' "$WHITELIST"
        printf '\n%s\n' '# Fork location and run settings; use the existing scNanoGPS conda environment.'
        printf 'GPS_REPO=%q\nNCORES=%q\nBARCODE_SUBSTITUTIONS=%q\nPYTHON_BIN=%q\n' \
            "$GPS_REPO" "$NCORES" "$BARCODE_SUBSTITUTIONS" "$PYTHON_BIN"
        cat <<'RUN'

set -euo pipefail
export FASTQ_DIR WORK_DIR WHITELIST GPS_REPO NCORES BARCODE_SUBSTITUTIONS PYTHON_BIN
exec bash "$GPS_REPO/run_tcr_pipeline.sh" "$@"
RUN
    } > "$OUTPUT"
)
chmod u+x "$OUTPUT"
printf 'Created: %s\n' "$OUTPUT"
printf '%s\n' 'Edit FASTQ_DIR, WORK_DIR, and WHITELIST at the top, then activate scNanoGPS.'
printf 'Run: bash %q\n' "$OUTPUT"
