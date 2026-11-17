import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest


GENERATOR = Path(__file__).resolve().parents[1] / "run_tcr.sh"


class TcrWrapperTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="tcr wrapper ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.environment = dict(os.environ)
        for key in ("FASTQ_DIR", "WORK_DIR", "WHITELIST", "GPS_REPO", "NCORES",
                    "BARCODE_SUBSTITUTIONS", "PYTHON_BIN"):
            self.environment.pop(key, None)

    def generate(self, *args, **settings):
        return subprocess.run(["bash", str(GENERATOR), *args], cwd=self.root,
                              env=dict(self.environment, **settings), text=True, capture_output=True)

    def stub_runner(self):
        repo = self.root / "tools with spaces"
        repo.mkdir()
        capture = self.root / "capture.json"
        code = ("import json,os,sys; from pathlib import Path; "
                "Path(os.environ['CAPTURE']).write_text(json.dumps({"
                "'settings':{k:os.environ[k] for k in ['FASTQ_DIR','WORK_DIR','WHITELIST','NCORES']},"
                "'args':sys.argv[1:]}))")
        (repo / "run_tcr_pipeline.sh").write_text(
            "#!/usr/bin/env bash\nexec " + shlex.quote(sys.executable) + " -c " + shlex.quote(code) + ' "$@"\n')
        return repo, capture

    def test_default_creates_editable_file_in_calling_directory_only(self):
        result = self.generate()
        self.assertEqual(result.returncode, 0, result.stderr)
        generated = self.root / "run_tcr_sample.sh"
        self.assertTrue(generated.is_file())
        lines = generated.read_text().splitlines()
        self.assertTrue(lines[2].startswith("FASTQ_DIR="))
        self.assertTrue(lines[3].startswith("WORK_DIR="))
        self.assertTrue(lines[4].startswith("WHITELIST="))
        self.assertEqual(lines[2:5], ["FASTQ_DIR=''", "WORK_DIR=''", "WHITELIST=''"])
        self.assertIn("GPS_REPO=", generated.read_text())
        self.assertEqual(sorted(p.name for p in self.root.iterdir()), ["run_tcr_sample.sh"])
        subprocess.run(["bash", "-n", str(generated)], check=True)

    def test_existing_edited_job_is_not_overwritten(self):
        generated = self.root / "run_tcr_sample.sh"
        generated.write_text("# my edited sample paths\n")
        result = self.generate()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(generated.read_text(), "# my edited sample paths\n")

    def test_paths_are_literal_and_stage_argument_is_forwarded(self):
        repo, capture = self.stub_runner()
        literal = str(self.root / "reads $(touch INJECTED) `touch INJECTED2` 'quoted'")
        result = self.generate("--output", "custom run.sh", GPS_REPO=str(repo),
                               FASTQ_DIR=literal, WORK_DIR=str(self.root / "output space"),
                               WHITELIST=str(self.root / "whitelist with spaces.txt.gz"), NCORES="12")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(capture.exists())
        subprocess.run(["bash", str(self.root / "custom run.sh"), "export"], cwd=self.root,
                       env=dict(self.environment, CAPTURE=str(capture)), check=True)
        observed = json.loads(capture.read_text())
        self.assertEqual(observed["settings"]["FASTQ_DIR"], literal)
        self.assertEqual(observed["settings"]["WHITELIST"], str(self.root / "whitelist with spaces.txt.gz"))
        self.assertEqual(observed["settings"]["NCORES"], "12")
        self.assertEqual(observed["args"], ["export"])
        self.assertFalse((self.root / "INJECTED").exists())
        self.assertFalse((self.root / "INJECTED2").exists())

    def test_editing_work_path_keeps_whitelist_empty_until_explicitly_set(self):
        repo, capture = self.stub_runner()
        result = self.generate(GPS_REPO=str(repo))
        self.assertEqual(result.returncode, 0, result.stderr)
        generated = self.root / "run_tcr_sample.sh"
        work = str(self.root / "edited working directory")
        lines = generated.read_text().splitlines()
        lines = ["WORK_DIR=" + shlex.quote(work) if line.startswith("WORK_DIR=") else line for line in lines]
        generated.write_text("\n".join(lines) + "\n")
        subprocess.run(["bash", str(generated)], cwd=self.root,
                       env=dict(self.environment, CAPTURE=str(capture)), check=True)
        observed = json.loads(capture.read_text())["settings"]
        self.assertEqual(observed["WORK_DIR"], work)
        self.assertEqual(observed["WHITELIST"], "")

    def test_pipeline_actions_are_not_mistaken_for_generation_options(self):
        result = self.generate("export")
        self.assertEqual(result.returncode, 2)
        self.assertFalse((self.root / "export").exists())


if __name__ == "__main__":
    unittest.main()
