# Agent Instructions

For all agent-driven work in this repository, follow [.github/copilot-instructions.md](.github/copilot-instructions.md).

## Quick Checklist Before Completing Any WARP Task

1. **Validate all modified WDLs** — and every WDL that imports them (transitively). Use the loop in *WDL Validation (MANDATORY)*.
2. **Version bumps are once per branch/PR.** If the top changelog entry is already an unreleased bump, append bullets to it — do not create a new version.
3. **Cascading bumps.** If a shared task changed, check every dependent pipeline. Either bump+changelog or add a "no functional impact" note.
4. **After merging develop**, grep for duplicate `pipeline_version` lines — keep the higher one.
5. **Update `pipeline_versions.txt`** only when a version actually changes.
6. **Sub-workflow contract.** Removing an input from a shared WDL requires removing it from every caller.
