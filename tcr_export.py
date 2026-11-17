#!/usr/bin/env python3
"""Export Scanner-processed TCR reads by whitelist barcode, without cell calling.

Run scanner.py first. This stage never invokes Assigner, Curator, or Reporter.
Accept a combined OB whitelist or a non-OCM chemistry-matched whitelist.
Sample assignment is performed outside this pipeline.
"""

import argparse
from collections import Counter, OrderedDict
import csv
from functools import lru_cache
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile


DNA = "ACGT"
READ1 = "CTACACGACGCTCTTCCGATCT"
READ1_RC = READ1.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def open_text(path, mode="rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="ascii", newline="")
    return path.open(mode, encoding="ascii", newline="")


def encode_barcode(barcode):
    code = 0
    for base in barcode:
        code = (code << 2) | DNA.index(base)
    return code


def decode_barcode(code):
    return "".join(DNA[(code >> shift) & 3] for shift in range(30, -1, -2))


def load_reference(path):
    """Use integer keys to limit memory for multimillion-barcode references."""
    reference = set()
    with open_text(path) as stream:
        for line in stream:
            barcode = line.strip()
            if not barcode:
                continue
            if not re.fullmatch(r"[ACGT]{16}", barcode):
                raise ValueError("Reference requires full 16-base A/C/G/T barcodes")
            code = encode_barcode(barcode)
            if code in reference:
                raise ValueError("Repeated reference barcode: " + barcode)
            reference.add(code)
    if not reference:
        raise ValueError("Empty barcode reference")
    return reference


def make_matcher(reference, substitutions):
    @lru_cache(maxsize=100000)
    def match(barcode):
        if not re.fullmatch(r"[ACGT]{16}", barcode):
            return None, "invalid_barcode"
        code = encode_barcode(barcode)
        if code in reference:
            return code, "exact"
        if not substitutions:
            return None, "no_match"
        candidate = None
        for shift in range(0, 32, 2):
            for change in (1, 2, 3):
                neighbor = code ^ (change << shift)
                if neighbor in reference:
                    if candidate is not None:
                        return None, "ambiguous"
                    candidate = neighbor
        return (candidate, "corrected") if candidate is not None else (None, "no_match")
    return match


def fastq_records(path):
    with open_text(path) as stream:
        while True:
            header = stream.readline()
            if not header:
                break
            sequence, separator, quality = (stream.readline() for _ in range(3))
            if not sequence or not separator or not quality:
                raise ValueError("Incomplete FASTQ record in " + str(path))
            header, sequence, separator, quality = (
                line.rstrip("\r\n") for line in (header, sequence, separator, quality))
            if not header.startswith("@") or not header[1:].strip() or not separator.startswith("+"):
                raise ValueError("Malformed FASTQ record in " + str(path))
            if not sequence or not re.fullmatch(r"[ACGTNacgtn]+", sequence):
                raise ValueError("Invalid FASTQ sequence in " + str(path))
            if len(sequence) != len(quality) or any(not 33 <= ord(c) <= 126 for c in quality):
                raise ValueError("Invalid FASTQ qualities in " + str(path))
            yield header, sequence, quality


class FastqWriters:
    """Bound open handles; reopened gzip streams form standard gzip members."""

    def __init__(self, root, limit):
        self.root, self.limit = root, limit
        self.handles = OrderedDict()

    def write(self, relative, header, sequence, quality):
        handle = self.handles.pop(relative, None)
        if handle is None:
            if len(self.handles) >= self.limit:
                self.handles.popitem(last=False)[1].close()
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            handle = gzip.open(path, "at", encoding="ascii", compresslevel=3)
        self.handles[relative] = handle
        handle.write("{}\n{}\n+\n{}\n".format(header, sequence, quality))

    def close(self):
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def export(args):
    output = Path(args.output).expanduser().resolve()
    if output.exists():
        raise ValueError("Output must be a new directory: " + str(output))
    directories = [Path(p).expanduser().resolve() for p in args.scanner_dir]
    if len(set(directories)) != len(directories):
        raise ValueError("Repeated Scanner directory would duplicate reads")
    if any(output == d or d in output.parents for d in directories):
        raise ValueError("Keep exports outside the Scanner input directories")
    inputs = []
    for directory in directories:
        pair = (directory / args.fastq_name, directory / args.barcode_name)
        if not all(p.is_file() for p in pair):
            raise ValueError("Missing paired Scanner outputs in " + str(directory))
        inputs.append(pair)
    reference_path = Path(args.whitelist).expanduser().resolve()
    reference = load_reference(reference_path)
    match = make_matcher(reference, args.max_barcode_substitutions)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="." + output.name + ".", dir=output.parent))
    writers = FastqWriters(staging, args.max_open_files)
    counts, barcode_counts = Counter(), Counter()
    fields = ["source", "read_id", "raw_barcode", "barcode", "umi", "status"]
    try:
        with gzip.open(staging / "assignments.tsv.gz", "wt", encoding="ascii") as audit:
            audit_writer = csv.writer(audit, delimiter="\t", lineterminator="\n")
            audit_writer.writerow(fields)
            for fastq, table in inputs:
                with open_text(table) as stream:
                    rows = csv.DictReader(stream, delimiter="\t")
                    if not {"rid", "BC", "UMI"}.issubset(rows.fieldnames or []):
                        raise ValueError("Unrecognized Scanner barcode table: " + str(table))
                    for record, row in zip_longest(fastq_records(fastq), rows):
                        if record is None or row is None:
                            raise ValueError("Scanner FASTQ/table record counts differ: " + str(fastq))
                        if None in row or any(value is None for value in row.values()):
                            raise ValueError("Malformed Scanner barcode table row: " + str(table))
                        header, sequence, quality = record
                        identifier = header[1:].split()[0]
                        if row["rid"] != identifier:
                            raise ValueError("Scanner FASTQ/table read IDs differ: " + identifier)
                        barcode, umi = row["BC"], row["UMI"]
                        corrected, status = match(barcode)
                        corrected_barcode = decode_barcode(corrected) if corrected is not None else ""
                        if not re.fullmatch(r"[ACGTN]{12}", umi):
                            status = "invalid_umi"
                        elif any(motif in sequence.upper() for motif in (READ1, READ1_RC)):
                            # Scanner already removed the assigned R1 end. Another R1
                            # indicates possible extra molecules, so do not assign it.
                            status = "internal_adapter"
                        counts["input_reads"] += 1
                        counts[status] += 1
                        audit_writer.writerow([str(fastq.parent), identifier, barcode,
                                               corrected_barcode, umi, status])
                        if status in {"exact", "corrected"}:
                            tagged = "@{}_{} CB:Z:{} CR:Z:{} RX:Z:{}".format(
                                identifier, umi, corrected_barcode, barcode, umi)
                            relative = Path("fastq") / (corrected_barcode + ".fastq.gz")
                            writers.write(relative, tagged, sequence, quality)
                            barcode_counts[corrected_barcode] += 1
                            counts["assigned_reads"] += 1
                        else:
                            writers.write(Path("unassigned.fastq.gz"), header, sequence, quality)
        writers.close()
        with (staging / "barcodes.tsv").open("w", encoding="ascii", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(["barcode", "reads", "fastq"])
            for barcode, count in sorted(barcode_counts.items()):
                writer.writerow([barcode, count, str(output / "fastq" / (barcode + ".fastq.gz"))])
        summary = {
            "sample_assignment": False,
            "cell_calling": False,
            "minimum_reads_per_barcode": None,
            "minimum_umis_per_barcode": None,
            "max_barcode_substitutions": args.max_barcode_substitutions,
            "reference": str(reference_path),
            "reference_sha256": sha256(reference_path),
            "reference_barcodes": len(reference),
            "scanner_directories": [str(p) for p in directories],
            "barcode_files": len(barcode_counts),
            "counts": dict(counts),
        }
        (staging / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        if output.exists():
            raise ValueError("Output appeared during processing; refusing to replace it")
        staging.rename(output)
        return summary
    finally:
        writers.close()
        if staging.exists():
            shutil.rmtree(staging)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scanner-dir", action="append", required=True,
                   help="Scanner directory; repeat only for chunks of the same library")
    p.add_argument("--fastq-name", default="processed.fastq.gz")
    p.add_argument("--barcode-name", default="barcode_list.tsv.gz")
    p.add_argument("--whitelist", required=True,
                   help="Combined OB or non-OCM whitelist: one 16-base barcode per line, no header; .gz supported")
    p.add_argument("--output", required=True, help="New directory outside Scanner inputs")
    p.add_argument("--max-barcode-substitutions", type=int, choices=(0, 1), default=0,
                   help="Default 0: exact matches. 1: rescue a unique one-substitution match; no indel correction")
    p.add_argument("--max-open-files", type=int, default=64)
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.max_open_files < 1:
        raise ValueError("--max-open-files must be positive")
    summary = export(args)
    print(json.dumps(summary["counts"], indent=2))
    print("Exported {} barcode FASTQs to {}".format(summary["barcode_files"], args.output))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, EOFError) as error:
        sys.exit("Error: " + str(error))
