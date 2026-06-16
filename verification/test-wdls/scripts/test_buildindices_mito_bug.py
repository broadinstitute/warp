#!/usr/bin/env python3
"""
Test: Mitochondrial contig deduplication bug guard

Bug Description (develop branch, fixed by PR #1769):
  When mito_accession was provided but run_mitofinder=false, the
  BuildStarSingleNucleus and BuildBWAreference tasks would incorrectly
  remove the mitochondrial contig from the genome FASTA, believing it
  was a duplicate—even though MitoFinder never ran to create one.

The refactored code in rc_buildindices_refactor prevents this bug
structurally:
  1. RemoveDuplicateMitoContig is ONLY called inside the
     `if (run_mitofinder)` conditional block.
  2. BuildStarSingleNucleus and BuildBWAreference no longer accept
     mito_accession as an input, so they cannot perform any contig removal.
  3. Downstream tasks receive `final_genome_fa` which resolves to the
     original genome when mitofinder is not run.

This test parses the WDL to verify these invariants hold.
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
    # Find the task definition
    pattern = rf'task\s+{re.escape(task_name)}\s*\{{(.*?)^\}}'
    match = re.search(pattern, wdl_text, re.DOTALL | re.MULTILINE)
    if not match:
        return None
    task_body = match.group(1)

    # Find the input block within the task
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
    pattern = r'workflow\s+BuildIndices\s*\{(.*?)^}'
    # This simple approach may not work for nested braces;
    # use a brace-counting parser instead.
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
        # Look for 'if ('
        match = re.search(r'\bif\s*\(', workflow_body[i:])
        if not match:
            break
        cond_start = i + match.start()
        # Find the condition expression
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

        # Find the block body
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
        print("BuildIndices Mitochondrial Deduplication Bug Test")
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
    # Test 1: BuildStarSingleNucleus must NOT accept mito_accession
    # ──────────────────────────────────────────────────────────────────
    star_inputs = extract_task_input_block(wdl_text, "BuildStarSingleNucleus")
    t.check(
        "BuildStarSingleNucleus does NOT have mito_accession input",
        star_inputs is not None and "mito_accession" not in star_inputs,
        detail="mito_accession found in BuildStarSingleNucleus inputs"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 2: BuildStarSingleNucleus must NOT have mito contig removal
    #         logic in its command block
    # ──────────────────────────────────────────────────────────────────
    star_cmd = extract_task_command_block(wdl_text, "BuildStarSingleNucleus")
    t.check(
        "BuildStarSingleNucleus command has no mito contig removal",
        star_cmd is not None and "mito_accession" not in star_cmd
        and "duplicate contig" not in star_cmd.lower(),
        detail="mito contig removal logic found in BuildStarSingleNucleus command"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 3: BuildBWAreference must NOT accept mito_accession
    # ──────────────────────────────────────────────────────────────────
    bwa_inputs = extract_task_input_block(wdl_text, "BuildBWAreference")
    t.check(
        "BuildBWAreference does NOT have mito_accession input",
        bwa_inputs is not None and "mito_accession" not in bwa_inputs,
        detail="mito_accession found in BuildBWAreference inputs"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 4: BuildBWAreference must NOT have mito contig removal
    # ──────────────────────────────────────────────────────────────────
    bwa_cmd = extract_task_command_block(wdl_text, "BuildBWAreference")
    t.check(
        "BuildBWAreference command has no mito contig removal",
        bwa_cmd is not None and "mito_accession" not in bwa_cmd
        and "duplicate contig" not in bwa_cmd.lower(),
        detail="mito contig removal logic found in BuildBWAreference command"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 5: RemoveDuplicateMitoContig task exists
    # ──────────────────────────────────────────────────────────────────
    t.check(
        "RemoveDuplicateMitoContig task is defined",
        "task RemoveDuplicateMitoContig" in wdl_text,
        detail="RemoveDuplicateMitoContig task definition not found"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 6: RemoveDuplicateMitoContig is called ONLY inside the
    #         if (run_mitofinder) block
    # ──────────────────────────────────────────────────────────────────
    workflow_body = extract_workflow_body(wdl_text)
    t.check(
        "Workflow body extracted successfully",
        workflow_body is not None,
        detail="Could not parse workflow body"
    )

    if workflow_body:
        cond_blocks = find_conditional_blocks(workflow_body)
        mitofinder_blocks = [
            (cond, body) for cond, body in cond_blocks
            if "run_mitofinder" in cond
        ]
        t.check(
            "Found if (run_mitofinder) conditional block",
            len(mitofinder_blocks) > 0,
            detail="No if (run_mitofinder) block found"
        )

        if mitofinder_blocks:
            mito_block_body = mitofinder_blocks[0][1]
            t.check(
                "RemoveDuplicateMitoContig called inside run_mitofinder block",
                "RemoveDuplicateMitoContig" in mito_block_body,
                detail="RemoveDuplicateMitoContig not found inside "
                       "if (run_mitofinder) block"
            )

        # Check it's NOT called outside
        # Remove the conditional blocks and check the remaining workflow body
        remaining = workflow_body
        for cond, body in cond_blocks:
            # Remove each block from the remaining text
            # Simple approach: replace the body
            remaining = remaining.replace(body, "")
        t.check(
            "RemoveDuplicateMitoContig NOT called outside run_mitofinder block",
            "RemoveDuplicateMitoContig" not in remaining
            or "select_first([RemoveDuplicateMitoContig" in remaining,
            detail="RemoveDuplicateMitoContig is called outside the "
                   "if (run_mitofinder) conditional"
        )

    # ──────────────────────────────────────────────────────────────────
    # Test 7: final_genome_fa falls back to original genome_fa when
    #         mitofinder is not run (select_first pattern)
    # ──────────────────────────────────────────────────────────────────
    fallback_pattern = re.search(
        r'File\s+final_genome_fa\s*=\s*select_first\(\['
        r'RemoveDuplicateMitoContig\.cleaned_fasta,\s*genome_fa\]\)',
        wdl_text
    )
    t.check(
        "final_genome_fa falls back to genome_fa via select_first",
        fallback_pattern is not None,
        detail="Expected: File final_genome_fa = select_first("
               "[RemoveDuplicateMitoContig.cleaned_fasta, genome_fa])"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 8: BuildStarSingleNucleus receives final_genome_fa, not
    #         raw genome_fa
    # ──────────────────────────────────────────────────────────────────
    star_call_pattern = re.search(
        r'call\s+BuildStarSingleNucleus\s*\{[^}]*genome_fa\s*=\s*final_genome_fa',
        wdl_text, re.DOTALL
    )
    t.check(
        "BuildStarSingleNucleus receives final_genome_fa",
        star_call_pattern is not None,
        detail="BuildStarSingleNucleus should use final_genome_fa, "
               "not genome_fa directly"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 9: BuildBWAreference receives final_genome_fa
    # ──────────────────────────────────────────────────────────────────
    bwa_call_pattern = re.search(
        r'call\s+BuildBWAreference\s*\{[^}]*genome_fa\s*=\s*final_genome_fa',
        wdl_text, re.DOTALL
    )
    t.check(
        "BuildBWAreference receives final_genome_fa",
        bwa_call_pattern is not None,
        detail="BuildBWAreference should use final_genome_fa, "
               "not genome_fa directly"
    )

    # ──────────────────────────────────────────────────────────────────
    # Test 10: No run_mitofinder input in BuildStarSingleNucleus
    #          (the develop-branch fix approach leaked this in; the
    #          refactored approach shouldn't need it)
    # ──────────────────────────────────────────────────────────────────
    t.check(
        "BuildStarSingleNucleus does NOT have run_mitofinder input",
        star_inputs is not None and "run_mitofinder" not in star_inputs,
        detail="run_mitofinder found in BuildStarSingleNucleus inputs — "
               "the task should not need to know about mitofinder"
    )

    # ──────────────────────────────────────────────────────────────────
    # Report
    # ──────────────────────────────────────────────────────────────────
    success = t.report()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
