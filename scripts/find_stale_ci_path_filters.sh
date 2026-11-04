#!/usr/bin/env bash
# Report .github/workflows/test_*.yml whose `paths:` filter is out of sync with
# the import closure of the Test WDL that workflow actually runs. Two failure
# modes, both of which silently skip a test on relevant edits:
#   DEAD   - a watched path no longer exists (or is the wrong file entirely)
#   BLIND  - a WDL the test imports is not covered by any watched path/glob
# Read-only; prints findings, exits 0.
#
# ponytail: intentionally NOT wired into CI yet (the test system is flaky and a
# blocking gate would add noise, not signal). Run it by hand before a cleanup PR.
# Future wiring: see AGENTS.md "Stale CI path filters".
set -euo pipefail
cd "$(dirname "$0")/.."

python3 - <<'PY'
import re, os, glob, fnmatch

def imports_of(f):
    d = os.path.dirname(f)
    try: src = open(f).read()
    except OSError: return []
    return [os.path.normpath(os.path.join(d, p))
            for p in re.findall(r'^import\s+"([^"]+)"', src, re.M)
            if not p.startswith("http")]

def closure(start):
    seen, stack = set(), [start]
    while stack:
        c = stack.pop()
        if c in seen: continue
        seen.add(c); stack += imports_of(c)
    return {s for s in seen if s.endswith(".wdl")}

def covered(path, patterns):
    for p in patterns:
        if p == path: return True
        if p.endswith("**") and path.startswith(p.rstrip("*")): return True
        if fnmatch.fnmatch(path, p): return True
    return False

n = 0
for y in sorted(glob.glob(".github/workflows/test_*.yml")):
    txt = open(y).read()
    listed = set(re.findall(r"-\s*'([^']+)'", txt)) | set(re.findall(r'-\s*"([^"]+)"', txt))
    wdl_listed = [p for p in listed if p.endswith(".wdl")]
    tests = sorted(set(re.findall(r'(verification/test-wdls/Test[\w]+\.wdl)', txt)))
    real = set()
    for t in tests:
        if os.path.exists(t): real |= closure(t)
    dead  = sorted(p for p in wdl_listed if not os.path.exists(p))
    blind = sorted(r for r in real if not covered(r, listed))
    if dead or blind or len(tests) != 1:
        print(f"{os.path.basename(y)}: tests={[os.path.basename(t) for t in tests]}")
        for p in dead:  print(f"    DEAD  watched path missing: {p}")
        for p in blind: print(f"    BLIND imported but unwatched: {p}")
        n += 1
print(f"\n{n} workflow(s) with stale path filters.")
PY
