import csv
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


SPEC = importlib.util.spec_from_file_location("tcr_export", Path(__file__).resolve().parents[1] / "tcr_export.py")
exporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(exporter)


class TcrExportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.reference = self.root / "reference.tsv"
        self.reference.write_text("barcode\tocm_group\n" + "A" * 16 + "\tOB1\n" + "C" * 16 + "\tOB2\n")
        self.samples = self.root / "samples.tsv"
        self.samples.write_text("ocm_group\tsample_id\nOB1\tS1\nOB2\tS2\n")

    def scanner(self, reads, name="scanner"):
        directory = self.root / name
        directory.mkdir()
        with gzip.open(directory / "processed.fastq.gz", "wt") as fq, \
                gzip.open(directory / "barcode_list.tsv.gz", "wt") as table:
            table.write("rid\torientation\tBC_start\tBC\tUMI\tSeq_end\tmean_BC_quality\n")
            for identifier, barcode, umi, sequence in reads:
                fq.write("@{} original_comment\n{}\n+\n{}\n".format(identifier, sequence, "I" * len(sequence)))
                table.write("{}\tH\t22\t{}\t{}\t{}\t30\n".format(identifier, barcode, umi, len(sequence)))
        return directory

    def options(self, directories, **kwargs):
        argv = ["--output", str(self.root / "export"), "--ocm-reference", str(self.reference),
                "--sample-map", str(self.samples)]
        for directory in directories:
            argv.extend(["--scanner-dir", str(directory)])
        args = exporter.parser().parse_args(argv)
        for key, value in kwargs.items():
            setattr(args, key, value)
        return args

    def test_single_read_barcodes_retained_and_umis_preserved(self):
        directory = self.scanner([("r1", "A" * 16, "G" * 12, "ACGT" * 20),
                                  ("r2", "C" * 16, "T" * 12, "TGCA" * 20)])
        original = (directory / "processed.fastq.gz").read_bytes()
        result = exporter.export(self.options([directory]))
        self.assertFalse(result["cell_calling"])
        self.assertEqual(result["barcode_files"], 2)
        self.assertEqual(result["counts"]["assigned_reads"], 2)
        output = self.root / "export/fastq/S1" / ("A" * 16 + ".fastq.gz")
        with gzip.open(output, "rt") as stream:
            records = stream.read().splitlines()
        self.assertTrue(records[0].startswith("@r1_" + "G" * 12))
        self.assertIn("RX:Z:" + "G" * 12, records[0])
        self.assertEqual(records[1], "ACGT" * 20)
        self.assertEqual(records[3], "I" * 80)
        self.assertEqual((directory / "processed.fastq.gz").read_bytes(), original)

    def test_substitution_requires_unique_full_barcode(self):
        code, status = exporter.make_matcher(exporter.load_reference(self.reference, True), 1)("T" + "A" * 15)
        self.assertEqual(exporter.decode_barcode(code), "A" * 16)
        self.assertEqual(status, "corrected")
        self.reference.write_text(self.reference.read_text() + "G" + "A" * 15 + "\tOB2\n")
        match = exporter.make_matcher(exporter.load_reference(self.reference, True), 1)
        self.assertEqual(match("T" + "A" * 15), (None, "ambiguous"))
        self.assertEqual(match("A" * 16), (exporter.encode_barcode("A" * 16), "exact"))

    def test_exact_default_does_not_rescue_errors(self):
        match = exporter.make_matcher(exporter.load_reference(self.reference, True), 0)
        self.assertEqual(match("T" + "A" * 15), (None, "no_match"))
        self.assertEqual(match("N" + "A" * 15), (None, "invalid_barcode"))

    def test_multiple_chunks_and_lru_reopening_preserve_all_reads(self):
        reads = [("r1", "A" * 16, "G" * 12, "ACGT" * 20),
                 ("r2", "C" * 16, "T" * 12, "ACGT" * 20),
                 ("r3", "A" * 16, "C" * 12, "ACGT" * 20)]
        first = self.scanner(reads)
        second = self.scanner([("r4", "A" * 16, "G" * 12, "TGCA" * 20)], "scanner2")
        result = exporter.export(self.options([first, second], max_open_files=1))
        self.assertEqual(result["counts"]["assigned_reads"], 4)
        output = self.root / "export/fastq/S1" / ("A" * 16 + ".fastq.gz")
        self.assertEqual(len(list(exporter.fastq_records(output))), 3)

    def test_rejected_reads_and_internal_adapters_remain_auditable(self):
        directory = self.scanner([("unknown", "G" * 16, "A" * 12, "ACGT" * 20),
                                  ("joined", "A" * 16, "A" * 12, "ACGT" * 20 + exporter.READ1 + "TGCA" * 20),
                                  ("umi", "A" * 16, "A" * 10, "TGCA" * 20)])
        result = exporter.export(self.options([directory]))
        self.assertEqual(result["barcode_files"], 0)
        self.assertEqual(result["counts"]["no_match"], 1)
        self.assertEqual(result["counts"]["internal_adapter"], 1)
        self.assertEqual(result["counts"]["invalid_umi"], 1)
        self.assertEqual(len(list(exporter.fastq_records(self.root / "export/unassigned.fastq.gz"))), 3)

    def test_non_ocm_uses_fixed_sample_without_ob_assignment(self):
        directory = self.scanner([("r1", "A" * 16, "C" * 12, "TGCA" * 20)])
        whitelist = self.root / "whitelist.txt"
        whitelist.write_text("A" * 16 + "\n")
        args = self.options([directory], ocm_reference=None, sample_map=None,
                            whitelist=str(whitelist), sample_id="SingleSample")
        result = exporter.export(args)
        self.assertEqual(result["mode"], "non_ocm")
        self.assertTrue((self.root / "export/fastq/SingleSample" / ("A" * 16 + ".fastq.gz")).exists())

    def test_mismatched_ids_fail_without_publishing(self):
        directory = self.scanner([("r1", "A" * 16, "A" * 12, "ACGT" * 20)])
        with gzip.open(directory / "processed.fastq.gz", "wt") as stream:
            stream.write("@different\nACGT\n+\nIIII\n")
        with self.assertRaisesRegex(ValueError, "read IDs differ"):
            exporter.export(self.options([directory]))
        self.assertFalse((self.root / "export").exists())
        self.assertFalse(list(self.root.glob(".export.*")))

    def test_missing_records_and_malformed_fastq_fail(self):
        directory = self.scanner([("r1", "A" * 16, "A" * 12, "ACGT" * 20)])
        with gzip.open(directory / "processed.fastq.gz", "wt") as stream:
            stream.write("")
        with self.assertRaisesRegex(ValueError, "record counts differ"):
            exporter.export(self.options([directory]))
        with gzip.open(directory / "processed.fastq.gz", "wt") as stream:
            stream.write("@r1\nACGT\n+\nIII\n")
        with self.assertRaisesRegex(ValueError, "qualities"):
            exporter.export(self.options([directory]))

    def test_existing_output_and_duplicated_inputs_rejected(self):
        directory = self.scanner([])
        with self.assertRaisesRegex(ValueError, "Repeated Scanner"):
            exporter.export(self.options([directory, directory]))
        output = self.root / "export"
        output.mkdir()
        (output / "keep").write_text("unchanged")
        with self.assertRaisesRegex(ValueError, "new directory"):
            exporter.export(self.options([directory]))
        self.assertEqual((output / "keep").read_text(), "unchanged")

    def test_whitelist_conflicts_rejected(self):
        self.reference.write_text(self.reference.read_text() + "A" * 16 + "\tOB2\n")
        with self.assertRaisesRegex(ValueError, "conflicting"):
            exporter.load_reference(self.reference, True)

    def test_unmapped_ob_group_not_assigned_to_another_sample(self):
        self.samples.write_text("ocm_group\tsample_id\nOB1\tS1\n")
        directory = self.scanner([("r1", "C" * 16, "G" * 12, "ACGT" * 20)])
        result = exporter.export(self.options([directory]))
        self.assertEqual(result["counts"]["unmapped_ocm_group"], 1)
        self.assertEqual(result["barcode_files"], 0)


if __name__ == "__main__":
    unittest.main()
