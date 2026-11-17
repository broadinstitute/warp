"""Tests for scripts/image_pins.py. No network: the registry is faked.

Run from the WARP root: python3 -m unittest discover -s scripts/tests -p 'test_*.py'
"""

import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import image_pins  # noqa: E402

REG = "us.gcr.io/broad-gotc-prod"
OLD = "sha256:" + "1" * 64
NEW = "sha256:" + "2" * 64

FILES = {
    # full references; samtools-bwa must not match bwa; comments are skipped
    "tasks/wdl/Align.wdl": (
        'task A {\n'
        '  String docker = "%s/bwa:1.0.0-0.7.17-1660770463"\n'
        "  String other = '%s/samtools-bwa:1.0.0'\n"
        '  # String old = "%s/bwa:0.9"\n'
        '}\n' % (REG, REG, REG)),
    # docker-prefix pipeline: bare references count here
    "pipelines/wdl/p/P.wdl": (
        'workflow P {\n'
        '  String docker_prefix = "%s/"\n'
        '  String bwa_docker = "bwa@%s"\n'
        '  String ubuntu_docker = "ubuntu:20.04"\n'
        '}\n' % (REG, OLD)),
    # no prefix literal: a bare "bwa:..." is not ours
    "pipelines/wdl/q/Q.wdl": 'workflow Q {\n  String x = "bwa:1.0.0"\n}\n',
    "verification/VerifyTasks.wdl": 'task V {\n  String docker = "%s/bwa:1.0.0"\n}\n' % REG,
}


class PinTest(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        for path, text in FILES.items():
            os.makedirs(os.path.join(self.root, os.path.dirname(path)), exist_ok=True)
            with open(os.path.join(self.root, path), "w") as handle:
                handle.write(text)

    def read(self, path):
        with open(os.path.join(self.root, path)) as handle:
            return handle.read()

    def test_pin_rewrites_full_and_bare_references(self):
        changed, skipped = image_pins.pin(self.root, REG, "bwa", "1.1.0-0.7.18-1795000000", NEW)
        align, p, q = self.read("tasks/wdl/Align.wdl"), self.read("pipelines/wdl/p/P.wdl"), self.read("pipelines/wdl/q/Q.wdl")
        self.assertIn('"%s/bwa:1.1.0-0.7.18-1795000000@%s"' % (REG, NEW), align)
        self.assertIn("samtools-bwa:1.0.0'", align)            # other image untouched
        self.assertIn('# String old = "%s/bwa:0.9"' % REG, align)  # comment untouched
        self.assertIn('"bwa:1.1.0-0.7.18-1795000000@%s"' % NEW, p)
        self.assertIn('"ubuntu:20.04"', p)
        self.assertEqual(q, FILES["pipelines/wdl/q/Q.wdl"])
        self.assertEqual(len(changed), 2)
        self.assertEqual(skipped, ["verification/VerifyTasks.wdl:2"])
        self.assertEqual(self.read("verification/VerifyTasks.wdl"), FILES["verification/VerifyTasks.wdl"])

    def test_pin_paths_dry_run_and_shared(self):
        changed, _ = image_pins.pin(self.root, REG, "bwa", "t", NEW, only=["pipelines/wdl/p"])
        self.assertEqual([c.split(":")[0] for c in changed], ["pipelines/wdl/p/P.wdl"])
        self.assertEqual(self.read("tasks/wdl/Align.wdl"), FILES["tasks/wdl/Align.wdl"])
        changed, _ = image_pins.pin(self.root, REG, "bwa", "t2", NEW, include_shared=True, dry_run=True)
        self.assertIn("verification/VerifyTasks.wdl:2", " ".join(changed))
        self.assertEqual(self.read("verification/VerifyTasks.wdl"), FILES["verification/VerifyTasks.wdl"])

    def test_check_all(self):
        def digest_of(registry, image, ref):
            return {"1.0.0-0.7.17-1660770463": OLD, "t": NEW}.get(ref)
        with mock.patch.object(image_pins, "digest_of", digest_of), \
                mock.patch.object(image_pins, "tags_by_digest", lambda r, i: {OLD: ["x"]} if i == "bwa" else {}):
            image_pins.pin(self.root, REG, "bwa", "t", NEW, only=["tasks"])
            problems = {(path, msg.split()[0]) for path, _, msg in image_pins.check(self.root, REG)}
        self.assertEqual(problems, {
            ("tasks/wdl/Align.wdl", "samtools-bwa:1.0.0"),        # no digest
            ("pipelines/wdl/p/P.wdl", "bwa"),                     # digest but no tag
            ("verification/VerifyTasks.wdl", "bwa:1.0.0"),        # no digest
        })

    def test_check_only_added_lines(self):
        git = lambda *a: subprocess.run(["git", "-C", self.root, *a], check=True, capture_output=True)
        git("init", "-q", "-b", "main")
        git("-c", "user.name=t", "-c", "user.email=t@t", "add", ".")
        git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qm", "base")
        image_pins.pin(self.root, REG, "bwa", "t", NEW, only=["tasks"])
        git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-qam", "pin")
        with mock.patch.object(image_pins, "digest_of", lambda r, i, ref: OLD):
            problems = image_pins.check(self.root, REG, base="HEAD~1")
        self.assertEqual([(p, n) for p, n, _ in problems], [("tasks/wdl/Align.wdl", 2)])
        self.assertIn("points at", problems[0][2])

    def test_release_tags_once_and_reuses(self):
        tagged = []
        with mock.patch.object(image_pins, "digest_of", lambda r, i, ref: NEW), \
                mock.patch.object(image_pins, "tags_by_digest", lambda r, i: {NEW: ["branch"]}), \
                mock.patch.object(image_pins, "add_tag", lambda r, i, d, t: tagged.append(t)):
            tag, digest = image_pins.release_tag(REG, "bwa", "branch", "1.1.0-0.7.18", dry_run=False)
        self.assertEqual(digest, NEW)
        self.assertRegex(tag, r"^1\.1\.0-0\.7\.18-[0-9]{10}$")
        self.assertEqual(tagged, [tag])
        with mock.patch.object(image_pins, "digest_of", lambda r, i, ref: NEW), \
                mock.patch.object(image_pins, "tags_by_digest", lambda r, i: {NEW: ["1.1.0-0.7.18-1795000000"]}), \
                mock.patch.object(image_pins, "add_tag", lambda *a: self.fail("must reuse the tag")):
            self.assertEqual(image_pins.release_tag(REG, "bwa", "branch", "1.1.0-0.7.18", False)[0],
                             "1.1.0-0.7.18-1795000000")
        with self.assertRaises(image_pins.Error):
            image_pins.release_tag(REG, "bwa", "branch", "v1", False)


if __name__ == "__main__":
    unittest.main()
