#!/usr/bin/env python3
"""
Static unit tests for embedded Hail code in all_of_us/mitochondria/MitoPostProcessing.wdl.

This test intentionally validates the Python/Hail logic inside the WDL command block
without executing Hail, so it can run quickly in CI-style environments.
"""

import os
import re
import sys

WDL_PATH = os.path.join(
    os.path.dirname(__file__),
    "..", "..", "..",
    "all_of_us", "mitochondria", "MitoPostProcessing.wdl",
)


class TestResult:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.results = []

    def check(self, description, condition, detail=""):
        if condition:
            self.passed += 1
            self.results.append(("PASS", description))
        else:
            self.failed += 1
            message = f"{description}: {detail}" if detail else description
            self.results.append(("FAIL", message))

    def report(self):
        print("=" * 72)
        print("MitoPostProcessing WDL Hail Logic Unit Test")
        print("=" * 72)
        for status, message in self.results:
            marker = "✓" if status == "PASS" else "✗"
            print(f"  {marker} [{status}] {message}")
        print("-" * 72)
        total = self.passed + self.failed
        print(f"  {self.passed} passed, {self.failed} failed, {total} total")
        print("=" * 72)
        return self.failed == 0


def extract_command_block(wdl_text):
    match = re.search(r"command\s*<<<(.*?)>>>", wdl_text, re.DOTALL)
    return match.group(1) if match else None


def extract_python_heredoc(command_block):
    match = re.search(r"python3\s*<<'PYCODE'\n(.*?)\nPYCODE", command_block, re.DOTALL)
    return match.group(1) if match else None


def contains_once(text, pattern):
    return len(re.findall(pattern, text, re.MULTILINE)) == 1


def main():
    wdl_path = os.path.abspath(WDL_PATH)
    if not os.path.exists(wdl_path):
        print(f"ERROR: WDL file not found at {wdl_path}")
        sys.exit(1)

    with open(wdl_path, "r", encoding="utf-8") as handle:
        wdl_text = handle.read()

    command_block = extract_command_block(wdl_text)
    python_block = extract_python_heredoc(command_block or "")

    t = TestResult()

    t.check(
        "Workflow inputs include output_path/input_path/output_base",
        all(token in wdl_text for token in ["String output_path", "String input_path", "String output_base"]),
    )

    t.check(
        "WDL contains an embedded python heredoc",
        python_block is not None,
        detail="Could not find python3 <<'PYCODE' ... PYCODE block",
    )

    if python_block is None:
        ok = t.report()
        sys.exit(0 if ok else 1)

    t.check(
        "Hail initialized once",
        contains_once(python_block, r"^hl\.init\("),
        detail="Expected exactly one hl.init call",
    )

    expected_import_lines = [
        r"^import\s+hail\s+as\s+hl$",
        r"^import\s+matplotlib$",
        r"^import\s+matplotlib\.pyplot\s+as\s+plt$",
        r"^import\s+numpy\s+as\s+np$",
        r"^import\s+pandas\s+as\s+pd$",
        r"^import\s+seaborn\s+as\s+sns$",
    ]
    for import_pattern in expected_import_lines:
        t.check(
            f"Single import line matching `{import_pattern}`",
            contains_once(python_block, import_pattern),
            detail="Import line missing or duplicated",
        )

    t.check(
        "Reads MatrixTable from input_path",
        "filt_annotated_mt_500k = hl.read_matrix_table(input_matrix_path)" in python_block,
    )

    t.check(
        "Uses notebook heteroplasmy threshold upper bound 0.05",
        "(mt.allele_fraction <= 0.05)" in python_block,
    )

    t.check(
        "Includes NUMT warning annotation",
        "numt_fp_warning" in python_block and "low_mtcn_numt_risk" in python_block,
    )

    t.check(
        "Exports VCF from mt_vcf",
        "hl.export_vcf(mt_vcf, vcf_local, tabix=True)" in python_block,
    )

    expected_svg_targets = {
        "variants_per_sample.svg": "variants_per_sample_svg",
        "mito_cn_distribution.svg": "mito_cn_distribution_svg",
        "variant_allele_frequency.svg": "variant_allele_frequency_svg",
        "variant_af_and_allele_fraction.svg": "variant_af_and_allele_fraction_svg",
        "numt_fp_by_mtcn.svg": "numt_fp_by_mtcn_svg",
        "haplogroup_heteroplasmy.svg": "haplogroup_heteroplasmy_svg",
        "haplogroup_homoplasmy.svg": "haplogroup_homoplasmy_svg",
    }

    for svg_name, svg_var in expected_svg_targets.items():
        t.check(
            f"Defines output filename {svg_name}",
            svg_name in python_block,
        )

        t.check(
            f"Saves {svg_name} via savefig",
            re.search(rf"savefig\(\s*{re.escape(svg_var)}\s*,", python_block) is not None,
            detail="savefig call missing",
        )

    t.check(
        "Sample metadata export present",
        "sample_ht.export(sample_metadata_local)" in python_block,
    )

    ok = t.report()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
