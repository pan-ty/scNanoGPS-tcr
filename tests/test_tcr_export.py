import gzip
import importlib.util
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
        self.reference = self.root / "barcodes.txt"
        self.reference.write_text("A" * 16 + "\n" + "C" * 16 + "\n")

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
        argv = ["--output", str(self.root / "export"), "--whitelist", str(self.reference)]
        for directory in directories:
            argv.extend(["--scanner-dir", str(directory)])
        args = exporter.parser().parse_args(argv)
        for key, value in kwargs.items():
            setattr(args, key, value)
        return args

    def test_single_read_barcodes_retained_qualities_and_umis_preserved(self):
        directory = self.scanner([("r1", "A" * 16, "G" * 12, "ACGT" * 20),
                                  ("r2", "C" * 16, "T" * 12, "TGCA" * 20)])
        original = (directory / "processed.fastq.gz").read_bytes()
        result = exporter.export(self.options([directory]))
        self.assertFalse(result["cell_calling"])
        self.assertFalse(result["sample_assignment"])
        self.assertEqual(result["barcode_files"], 2)
        output = self.root / "export/fastq" / ("A" * 16 + ".fastq.gz")
        with gzip.open(output, "rt") as stream:
            records = stream.read().splitlines()
        self.assertTrue(records[0].startswith("@r1_" + "G" * 12))
        self.assertIn("RX:Z:" + "G" * 12, records[0])
        self.assertNotIn("OG:Z:", records[0])
        self.assertEqual(records[1], "ACGT" * 20)
        self.assertEqual(records[3], "I" * 80)
        self.assertEqual((directory / "processed.fastq.gz").read_bytes(), original)

    def test_substitution_requires_unique_full_barcode(self):
        code, status = exporter.make_matcher(exporter.load_reference(self.reference), 1)("T" + "A" * 15)
        self.assertEqual(exporter.decode_barcode(code), "A" * 16)
        self.assertEqual(status, "corrected")
        self.reference.write_text(self.reference.read_text() + "G" + "A" * 15 + "\n")
        match = exporter.make_matcher(exporter.load_reference(self.reference), 1)
        self.assertEqual(match("T" + "A" * 15), (None, "ambiguous"))
        self.assertEqual(match("A" * 16), (exporter.encode_barcode("A" * 16), "exact"))

    def test_exact_default_does_not_rescue_errors(self):
        match = exporter.make_matcher(exporter.load_reference(self.reference), 0)
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
        output = self.root / "export/fastq" / ("A" * 16 + ".fastq.gz")
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

    def test_gzipped_whitelist_uses_same_sample_agnostic_interface(self):
        directory = self.scanner([("r1", "A" * 16, "C" * 12, "TGCA" * 20)])
        whitelist = self.root / "whitelist.txt.gz"
        with gzip.open(whitelist, "wt") as stream:
            stream.write("A" * 16 + "\n")
        result = exporter.export(self.options([directory], whitelist=str(whitelist)))
        self.assertEqual(result["barcode_files"], 1)
        self.assertFalse(result["sample_assignment"])
        self.assertNotIn("sample_map", result)
        self.assertTrue((self.root / "export/fastq" / ("A" * 16 + ".fastq.gz")).exists())

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

    def test_repeated_reference_barcode_rejected(self):
        self.reference.write_text(self.reference.read_text() + "A" * 16 + "\n")
        with self.assertRaisesRegex(ValueError, "Repeated reference"):
            exporter.load_reference(self.reference)

    def test_corrected_barcode_preserves_raw_barcode_in_header(self):
        raw = "T" + "A" * 15
        directory = self.scanner([("r1", raw, "G" * 12, "ACGT" * 20)])
        result = exporter.export(self.options([directory], max_barcode_substitutions=1))
        self.assertEqual(result["counts"]["corrected"], 1)
        output = self.root / "export/fastq" / ("A" * 16 + ".fastq.gz")
        with gzip.open(output, "rt") as stream:
            header = stream.readline()
        self.assertIn("CB:Z:" + "A" * 16, header)
        self.assertIn("CR:Z:" + raw, header)


if __name__ == "__main__":
    unittest.main()
