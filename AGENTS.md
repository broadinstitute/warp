# Agent Instructions

Single entry point for all agent-driven work in this repository. Start here, follow the links to authoritative references, and consult the *Agent-Specific Notes* below for guidance that is not captured in the human-facing style guides.

## Authoritative References

| Topic | Document | Audience |
| --- | --- | --- |
| Full agent operating instructions (WDL patterns, validation, multi-cloud, CI triggers, merge conflicts) | [.github/copilot-instructions.md](.github/copilot-instructions.md) | Agents |
| Changelog and versioning rules (when to bump, format, cascading bumps, `pipeline_versions.txt`) | [changelog_style.md](changelog_style.md) | Both |
| WDL style guide (formatting, naming, `parameter_meta`, heredoc convention) | [WARP_WDL_Style_Guide.md](WARP_WDL_Style_Guide.md) | Both |
| Docusaurus markdown style (admonitions, tables, tabs, code-block highlighting, cross-refs) | [website/docs/contribution/contribute_to_warp_docs/doc_style.md](website/docs/contribution/contribute_to_warp_docs/doc_style.md) | Both |
| How to build / serve / validate the docs site | [website/docs/contribution/contribute_to_warp_docs/docsite_maintenance.md](website/docs/contribution/contribute_to_warp_docs/docsite_maintenance.md) | Both |

The human-facing guides (last three rows) are canonical for style. Anything *agent-specific* lives here or in [.github/copilot-instructions.md](.github/copilot-instructions.md) — do not duplicate human-facing style content into the agent docs.

## Quick Checklist Before Completing Any WARP Task

1. **Validate all modified WDLs** — and every WDL that imports them (transitively). See *WDL Validation (MANDATORY)* in [.github/copilot-instructions.md](.github/copilot-instructions.md).
2. **Changelog and versioning** — follow [changelog_style.md](changelog_style.md).
3. **Cascading version bumps** — when modifying a shared task (`tasks/wdl/`), every pipeline that imports it (directly or transitively) needs either a patch bump + changelog note **or** an explicit "no functional impact" entry. See the dependency chains in [.github/copilot-instructions.md](.github/copilot-instructions.md#cascading-version-bumps). Validation passing is not sufficient — changelogs must also be updated.
4. **After merging develop**, grep for duplicate `pipeline_version` lines — keep the higher one.
5. **Sub-workflow contract** — removing an input from a shared WDL requires removing it from every caller.
6. **Stale example/test inputs** — when you rename a workflow or remove/rename inputs, audit `pipelines/wdl/<name>/example_inputs/*.json` and `test_inputs/**/*.json`; they break silently because they are not checked by womtool.
7. **Touching a pipeline's interface** — also update the pipeline's docs page under `website/docs/Pipelines/<Name>_Pipeline/README.md` and run `yarn --cwd=website build` to catch broken links.

## Agent-Specific Notes

Conventions that have emerged from agent-driven sessions and that are not (or not yet) covered by the human style guides.

### Pipeline documentation layout

Each pipeline has **two** doc artifacts:

- **In-repo `README.md`** (`pipelines/wdl/<name>/README.md`) — short. Should contain: one-sentence summary, link to the full docs page on the WARP site, minimal "Running the pipeline" snippet, required inputs at a glance, link to `<Pipeline>.changelog.md`. Do **not** mirror the full docs page here.
- **Docs-site page** (`website/docs/Pipelines/<Name>_Pipeline/README.md`) — full documentation. Required Docusaurus frontmatter:

  ```yaml
  ---
  sidebar_position: 1
  slug: /Pipelines/<Name>_Pipeline/README
  ---
  ```

  Plus a sibling `_category_.json` when creating a new category.

When you change a pipeline's interface, update **both** artifacts (or, at minimum, ensure the slim README still points correctly and the docs page reflects the new inputs/outputs).

### Validating documentation changes

See *Documentation → Validate docs changes* in [.github/copilot-instructions.md](.github/copilot-instructions.md). The key distinction: `yarn --cwd=website start` does not fail on broken links — always use `yarn --cwd=website build` before marking docs work complete.

### `input_id` output prefix pattern

When a pipeline emits per-sample artifacts, accept a `String input_id` workflow input and prefix every output filename with `~{input_id}_`. This matches the convention used by scANVI, Optimus, Multiome, and other Skylab-origin pipelines.

### Optimus/Skylab chemistry and reference assets

When adding support for a new 10x chemistry in Optimus:
- Use `tenx_chemistry_version` (Integer, e.g. `4`) for the major chemistry version and `tenx_chemistry_subversion` (optional String, e.g. `"v4_TRU"`) for whitelist variant selection within that version.
- Optimus whitelist files live at `gs://gcp-public-data--broad-references/optimus_whitelists/`. Do **not** use the old `RNA/resources/` path (which contained a `febrary` typo and is incorrect).
- Validate that inputs required by a given chemistry (e.g. `i1_fastq` for v4) are enforced via `ErrorWithMessage` when the chemistry version demands them.

### When you discover a recurring agent mistake

Record it here (under *Agent-Specific Notes*) or in [.github/copilot-instructions.md](.github/copilot-instructions.md) — not in the human-facing style guides. Keep entries short and link out to the canonical reference rather than duplicating it. **Before adding a new note, check whether the information is already in copilot-instructions.md** — prefer a one-line reference over a duplicate section.
