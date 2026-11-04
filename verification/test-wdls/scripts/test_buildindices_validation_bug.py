#!/usr/bin/env python3
"""
Test: GTF validation always runs for non-marmoset organisms

Bug Description:
  In the original monolithic BuildStarSingleNucleus task, genome_build
  and genome_source validation against the GTF header ran unconditionally
  for non-marmoset organisms — even when skip_gtf_modification was true
  (the equivalent of run_modify_gtf=false).

  During the refactoring into separate tasks (ModifyGTF, ModifyGTFMarmoset,
  FixGeneNames, etc.), the validation checks were placed inside the
  ModifyGTF task, which is conditional on `run_modify_gtf`.  This meant
  that when `run_modify_gtf = false`, the validation silently didn't run,
  allowing mismatched GTF/genome combinations through to STAR & BWA
  builds without error.

The fix extracts validation into a separate ValidateGTF task that runs
unconditionally for non-marmoset organisms.  This test verifies:

  1. A ValidateGTF task exists.
  2. ValidateGTF checks genome_build AND genome_source.
  3. ValidateGTF is called outside any run_modify_gtf conditional
     (so it runs even when run_modify_gtf = false).
  4. ValidateGTF is called inside the if (!is_marmoset) block
     (matching the original reference behavior).
  5. ModifyGTF does NOT duplicate the validation checks.
"""

import re
import sys
import os

WDL_PATH = os.path.join(
    os.path.dirname(__file__),
    "..", "..", "..",
    "pipelines", "wdl", "build_indices", "BuildIndices.wdl"
)


def read_wdl(path):
    with open(path, "r") as f:
        return f.read()


def extract_task_input_block(wdl_text, task_name):
    """Extract the input { ... } block from a task definition."""
    pattern = rf'task\s+{re.escape(task_name)}\s*\{{(.*?)^\}}'
    match = re.search(pattern, wdl_text, re.DOTALL | re.MULTILINE)
    if not match:
        return None
    task_body = match.group(1)
    input_pattern = r'input\s*\{([^}]*)\}'
    input_match = re.search(input_pattern, task_body)
    if not input_match:
        return None
    return input_match.group(1)


def extract_task_command_block(wdl_text, task_name):
    """Extract the command <<< ... >>> block from a task definition."""
    pattern = rf'task\s+{re.escape(task_name)}\s*\{{(.*?)^\}}'
    match = re.search(pattern, wdl_text, re.DOTALL | re.MULTILINE)
    if not match:
        return None
    task_body = match.group(1)
    cmd_pattern = r'command\s*<<<(.*?)>>>'
    cmd_match = re.search(cmd_pattern, task_body, re.DOTALL)
    if not cmd_match:
        return None
    return cmd_match.group(1)


def extract_workflow_body(wdl_text):
    """Extract the workflow body (between first { and its matching })."""
    start = wdl_text.find("workflow BuildIndices")
    if start == -1:
        return None
    brace_start = wdl_text.index("{", start)
    depth = 0
    for i in range(brace_start, len(wdl_text)):
        if wdl_text[i] == '{':
            depth += 1
        elif wdl_text[i] == '}':
            depth -= 1
            if depth == 0:
                return wdl_text[brace_start + 1:i]
    return None


def find_conditional_blocks(workflow_body):
    """Find all top-level if (...) { ... } blocks and their contents."""
    blocks = []
    i = 0
    while i < len(workflow_body):
        match = re.search(r'\bif\s*\(', workflow_body[i:])
        if not match:
            break
        cond_start = i + match.start()
        paren_start = cond_start + match.end() - match.start() - 1
        depth = 0
        cond_end = paren_start
        for j in range(paren_start, len(workflow_body)):
            if workflow_body[j] == '(':
                depth += 1
            elif workflow_body[j] == ')':
                depth -= 1
                if depth == 0:
                    cond_end = j
                    break
        condition = workflow_body[paren_start + 1:cond_end].strip()

        brace_start = workflow_body.index('{', cond_end)
        depth = 0
        block_end = brace_start
        for j in range(brace_start, len(workflow_body)):
            if workflow_body[j] == '{':
                depth += 1
            elif workflow_body[j] == '}':
                depth -= 1
                if depth == 0:
                    block_end = j
                    break
        block_body = workflow_body[brace_start + 1:block_end]
        blocks.append((condition, block_body))
        i = block_end + 1
    return blocks


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
            msg = f"{description}: {detail}" if detail else description
            self.results.append(("FAIL", msg))

    def report(self):
        print("=" * 70)
        print("BuildIndices GTF Validation Bug Test")
        print("=" * 70)
        for status, desc in self.results:
            marker = "✓" if status == "PASS" else "✗"
            print(f"  {marker} [{status}] {desc}")
        print("-" * 70)
        print(f"  {self.passed} passed, {self.failed} failed, "
              f"{self.passed + self.failed} total")
        print("=" * 70)
        return self.failed == 0


def main():
    wdl_path = os.path.abspath(WDL_PATH)
    if not os.path.exists(wdl_path):
        print(f"ERROR: WDL file not found at {wdl_path}")
        sys.exit(1)

    wdl_text = read_wdl(wdl_path)
    t = TestResult()

    # ──────────────────────────────────────────────────────────────────
    # Test 1: ValidateGTF task exists
    # ──────────────────────────────────────────────────────────────────
    t.check(
        "ValidateGTF task is defined",
        "task ValidateGTF" in wdl_text,
        detail="ValidateGTF task definition not found"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 2: ValidateGTF checks genome_build in its command
    # ──────────────────────────────────────────────────────────────────
    validate_cmd = extract_task_command_block(wdl_text, "ValidateGTF")
    t.check(
        "ValidateGTF command checks genome_build",
        validate_cmd is not None and "genome_build" in validate_cmd,
        detail="genome_build check not found in ValidateGTF command"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 3: ValidateGTF checks genome_source in its command
    # ──────────────────────────────────────────────────────────────────
    t.check(
        "ValidateGTF command checks genome_source",
        validate_cmd is not None and "genome_source" in validate_cmd,
        detail="genome_source check not found in ValidateGTF command"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 4: ValidateGTF is NOT inside any run_modify_gtf conditional
    #         (it must run even when run_modify_gtf = false)
    # ──────────────────────────────────────────────────────────────────
    workflow_body = extract_workflow_body(wdl_text)
    t.check(
        "Workflow body extracted successfully",
        workflow_body is not None,
        detail="Could not parse workflow body"
    )

    if workflow_body:
        cond_blocks = find_conditional_blocks(workflow_body)

        # Find blocks that condition on run_modify_gtf
        modify_gtf_blocks = [
            (cond, body) for cond, body in cond_blocks
            if "run_modify_gtf" in cond
        ]

        validate_in_modify_block = False
        for cond, body in modify_gtf_blocks:
            if "ValidateGTF" in body:
                validate_in_modify_block = True
                break

        t.check(
            "ValidateGTF is NOT inside a run_modify_gtf conditional",
            not validate_in_modify_block,
            detail="ValidateGTF found inside an if (run_modify_gtf ...) block — "
                   "it must run unconditionally"
        )

    # ──────────────────────────────────────────────────────────────────
    # Test 5: ValidateGTF IS called inside the if (!is_marmoset) block
    #         (matching the original reference behavior)
    # ──────────────────────────────────────────────────────────────────
    if workflow_body:
        marmoset_guard_blocks = [
            (cond, body) for cond, body in cond_blocks
            if "is_marmoset" in cond and "run_modify_gtf" not in cond
        ]

        validate_in_marmoset_guard = False
        for cond, body in marmoset_guard_blocks:
            if "ValidateGTF" in body:
                validate_in_marmoset_guard = True
                break

        t.check(
            "ValidateGTF is called inside if (!is_marmoset) block",
            validate_in_marmoset_guard,
            detail="ValidateGTF should be guarded by a marmoset check "
                   "(validation only applies to non-marmoset)"
        )

    # ──────────────────────────────────────────────────────────────────
    # Test 6: ModifyGTF does NOT duplicate the genome_build/source
    #         validation (it should be in ValidateGTF only)
    # ──────────────────────────────────────────────────────────────────
    modify_cmd = extract_task_command_block(wdl_text, "ModifyGTF")
    has_build_check = modify_cmd is not None and (
        "genome_build" in modify_cmd or
        "genome version" in modify_cmd.lower()
    )
    t.check(
        "ModifyGTF does NOT duplicate genome_build validation",
        not has_build_check,
        detail="genome_build validation still found in ModifyGTF command — "
               "it should live in ValidateGTF only"
    )

    has_source_check = modify_cmd is not None and (
        "genome_source" in modify_cmd or
        "source of genome" in modify_cmd.lower()
    )
    t.check(
        "ModifyGTF does NOT duplicate genome_source validation",
        not has_source_check,
        detail="genome_source validation still found in ModifyGTF command — "
               "it should live in ValidateGTF only"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 7: ValidateGTF receives FixGeneNames.fixed_gtf as input
    #         (validation runs on the gene-name-fixed GTF, matching
    #         the original reference behavior)
    # ──────────────────────────────────────────────────────────────────
    validate_call_pattern = re.search(
        r'call\s+ValidateGTF\s*\{[^}]*annotation_gtf\s*=\s*FixGeneNames\.fixed_gtf',
        wdl_text, re.DOTALL
    )
    t.check(
        "ValidateGTF receives FixGeneNames.fixed_gtf",
        validate_call_pattern is not None,
        detail="ValidateGTF should use FixGeneNames.fixed_gtf as input"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 8: FixGeneNames always runs (not inside any conditional)
    # ──────────────────────────────────────────────────────────────────
    if workflow_body:
        fix_in_conditional = False
        for cond, body in cond_blocks:
            if "call FixGeneNames" in body:
                fix_in_conditional = True
                break

        t.check(
            "FixGeneNames is NOT inside any conditional block",
            not fix_in_conditional,
            detail="FixGeneNames should run unconditionally"
        )

    # Check it IS called at the workflow level
    t.check(
        "FixGeneNames is called in the workflow",
        workflow_body is not None and "call FixGeneNames" in workflow_body,
        detail="FixGeneNames call not found in workflow body"
    )

    # ──────────────────────────────────────────────────────────────────
    # Report
    # ──────────────────────────────────────────────────────────────────
    success = t.report()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
