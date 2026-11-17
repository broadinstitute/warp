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


class Repo(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        for path, text in FILES.items():
            os.makedirs(os.path.join(self.root, os.path.dirname(path)), exist_ok=True)
            with open(os.path.join(self.root, path), "w") as handle:
                handle.write(text)

    def read(self, path):
        with open(os.path.join(self.root, path)) as handle:
            return handle.read()


class PinTest(Repo):
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
                mock.patch.object(image_pins, "manifests", lambda r, i: {OLD: {"tag": ["x"]}} if i == "bwa" else {}):
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

    def test_image_version(self):
        script = ('DOCKER_IMAGE_VERSION=1.1.2\nTIMESTAMP=$(date +"%s")\nDRAGMAP_VERSION="1.2.1"  # pinned\n'
                  'SAMTOOLS_VERSION=1.11\nfunction main(){\n    SAMTOOLS_VERSION="$2"\n'
                  '    IMAGE_TAG="$TAG"\n    IMAGE_TAG="${DOCKER_IMAGE_VERSION}-${DRAGMAP_VERSION}-$SAMTOOLS_VERSION-$TIMESTAMP"\n}\n')
        self.assertEqual(image_pins.image_version(script, "x"), "1.1.2-1.2.1-1.11")
        prefixed = 'DOCKER_IMAGE_VERSION=1.0.4\nV="1.0"\n  IMAGE_TAG="atac_barcodes:$DOCKER_IMAGE_VERSION-$V-$TIMESTAMP"\n'
        self.assertEqual(image_pins.image_version(prefixed, "x"), "1.0.4-1.0")
        with self.assertRaises(image_pins.Error):  # computed at build time: needs --version
            image_pins.image_version('DOCKER_IMAGE_VERSION=4.0.1\n IMAGE_TAG="$DOCKER_IMAGE_VERSION-$N-$TIMESTAMP"\n', "x")


PR = {"number": 7, "head": {"ref": "bwa-fix", "sha": "h" * 40}, "base": {"ref": "develop", "sha": "b" * 40}}
WORKFLOW = ("env:\n  DOCKER_REGISTRY: us.gcr.io\n  GCR_PATH: broad-gotc-prod/bwa\n"
            "jobs:\n  build:\n    defaults:\n      run:\n        working-directory: 3rd-party-tools/bwa\n")
BUILD = 'DOCKER_IMAGE_VERSION=%s\nBWA_VERSION="0.7.18"\n    IMAGE_TAG="$DOCKER_IMAGE_VERSION-$BWA_VERSION-$TIMESTAMP"\n'


class ReleaseTest(Repo):
    def release(self, status="completed", base_version="1.0.0", uploaded=1794902700000, tags=(), **kwargs):
        run = {"path": ".github/workflows/build-bwa.yml", "status": status, "conclusion": "success",
               "created_at": "2026-11-17T08:00:00Z", "updated_at": "2026-11-17T08:10:00Z"}  # 1794902400-3000 s
        runs = {"workflow_runs": [run, dict(run, path=".github/workflows/trivy.yml")]}
        files = {(PR["head"]["sha"], ".github/workflows/build-bwa.yml"): WORKFLOW,
                 (PR["head"]["sha"], "3rd-party-tools/bwa/docker_build.sh"): BUILD % "1.1.0",
                 (PR["base"]["sha"], "3rd-party-tools/bwa/docker_build.sh"): BUILD % base_version}
        self.tagged, registry = [], {"bwa-fix": NEW}

        def add_tag(_, image, digest, tag):
            self.tagged.append(tag)
            registry[tag] = digest
        with mock.patch.object(image_pins, "api", lambda path: PR if path.startswith("pulls/") else runs), \
                mock.patch.object(image_pins, "raw", lambda sha, path, missing_ok=False: files.get((sha, path))), \
                mock.patch.object(image_pins, "digest_of", lambda r, i, ref: registry.get(ref)), \
                mock.patch.object(image_pins, "manifests", lambda r, i: {
                    NEW: {"tag": ["bwa-fix", *tags], "timeUploadedMs": str(uploaded)}, OLD: {"tag": []}}), \
                mock.patch.object(image_pins, "add_tag", add_tag):
            image_pins.release(self.root, REG, "https://github.com/broadinstitute/warp-tools/pull/7", **kwargs)

    def test_release_tags_the_pr_build_and_pins_it(self):
        self.release(only=["tasks"])
        self.assertRegex(self.tagged[0], r"^1\.1\.0-0\.7\.18-[0-9]{10}$")
        self.assertIn('"%s/bwa:%s@%s"' % (REG, self.tagged[0], NEW), self.read("tasks/wdl/Align.wdl"))
        self.assertEqual(self.read("pipelines/wdl/p/P.wdl"), FILES["pipelines/wdl/p/P.wdl"])

    def test_release_reuses_a_tag_and_takes_overrides(self):
        self.release(tags=["1.1.0-0.7.18-1795000000"])
        self.assertEqual(self.tagged, [])
        self.assertIn('"bwa:1.1.0-0.7.18-1795000000@%s"' % NEW, self.read("pipelines/wdl/p/P.wdl"))
        self.release(digest=OLD, version="1.1.0-x", uploaded=0)  # an earlier build: no upload-time check
        self.assertRegex(self.tagged[0], r"^1\.1\.0-x-[0-9]{10}$")

    def test_release_refuses(self):
        for kwargs, message in (({"status": "in_progress"}, "still running"),
                                ({"base_version": "1.1.0"}, "bump it"),
                                ({"uploaded": 1794800000000}, "not the build")):
            with self.assertRaisesRegex(image_pins.Error, message):
                self.release(**kwargs)
        self.assertEqual(self.tagged, [])
        self.assertEqual(self.read("tasks/wdl/Align.wdl"), FILES["tasks/wdl/Align.wdl"])


if __name__ == "__main__":
    unittest.main()
