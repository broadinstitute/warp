#!/usr/bin/env bash
# Report WDL `import "..." as NS` lines whose namespace `NS.` is never used.
# Structs are skipped: WDL references struct types by bare name, so a structs/
# import legitimately has no `NS.` usage. Read-only; prints findings, exits 0.
#
# ponytail: intentionally NOT wired into CI yet (the test system is flaky and a
# blocking gate would add noise, not signal). Run it by hand before a cleanup PR.
# Future wiring: see AGENTS.md "Unused imports".
set -euo pipefail
cd "$(dirname "$0")/.."

python3 - <<'PY'
import re, glob, os
files = (glob.glob("pipelines/**/*.wdl", recursive=True)
         + glob.glob("verification/**/*.wdl", recursive=True)
         + glob.glob("tasks/**/*.wdl", recursive=True))
n = 0
for f in sorted(files):
    src = open(f).read()
    body = "\n".join(l for l in src.splitlines() if not l.lstrip().startswith("import "))
    for path, alias in re.findall(r'^import\s+"([^"]+)"(?:\s+as\s+(\w+))?', src, re.M):
        if "structs/" in path.lower():
            continue
        ns = alias or os.path.basename(path)[:-4]
        if not re.search(re.escape(ns) + r'\.', body):
            print(f"{f}: unused import '{ns}' ({path})")
            n += 1
print(f"\n{n} unused non-struct import(s).")
PY
