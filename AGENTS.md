# Agent Instructions

For all agent-driven work in this repository, follow [.github/copilot-instructions.md](.github/copilot-instructions.md).

For all changelog and versioning rules (when to bump, how to write entries, cascading bumps, `pipeline_versions.txt`), see [changelog_style.md](changelog_style.md).

## Quick Checklist Before Completing Any WARP Task

1. **Validate all modified WDLs** — and every WDL that imports them (transitively). Use the loop in *WDL Validation (MANDATORY)*.
2. **Changelog and versioning** — follow [changelog_style.md](changelog_style.md) *Agent / Automation Rules*.
3. **After merging develop**, grep for duplicate `pipeline_version` lines — keep the higher one.
4. **Sub-workflow contract.** Removing an input from a shared WDL requires removing it from every caller.
