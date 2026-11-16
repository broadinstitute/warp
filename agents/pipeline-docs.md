# Pipeline docs: agent guide

Read this when changing a pipeline's documentation or the Docusaurus site. The general rules are in [AGENTS.md](../AGENTS.md). For Docusaurus markdown conventions follow [doc_style.md](../website/docs/contribution/contribute_to_warp_docs/doc_style.md); for site build/serve commands follow [docsite_maintenance.md](../website/docs/contribution/contribute_to_warp_docs/docsite_maintenance.md).

## Where a pipeline's docs live

Expect a pipeline's documentation to be split across two places. This isn't an enforced policy — it's just how pipelines evolve: a README starts in the pipeline directory, and the full write-up later lands on the docs site.

- **In-repo `README.md`** (`pipelines/wdl/<name>/README.md`) — usually slim and mostly a pointer to the docs-site page: a one-sentence summary, a link to the full docs page on the WARP site, a minimal "Running the pipeline" snippet, required inputs at a glance, and a link to `<Pipeline>.changelog.md`. Keep it that way — don't mirror the full docs page into it.
- **Docs-site page** (`website/docs/Pipelines/<Name>_Pipeline/README.md`) — the full documentation. Required Docusaurus frontmatter:

  ```yaml
  ---
  sidebar_position: 1
  slug: /Pipelines/<Name>_Pipeline/README
  ---
  ```

  Plus a sibling `_category_.json` when creating a new category.

**Do not hardcode version information.** Pipeline documentation, especially the header table on the docs-site page, should not contain hardcoded version numbers or dates. This information becomes stale. Instead, link to the pipeline's `.changelog.md` file as the single source of truth for versioning.

When you change a pipeline's interface, update the docs-site page so its inputs/outputs stay accurate, and make sure the slim in-repo README still points correctly.

## Validating docs changes

Validate with a full build — it catches broken links and bad frontmatter:

```bash
yarn --cwd=website install   # first time only
yarn --cwd=website build
```

`yarn --cwd=website start` is fine for previewing but does **not** fail on broken links — always use `build` before marking docs work complete.

## Cross-page links

Prefer relative paths to other pipeline pages (e.g. `[Multiome](../Multiome_Pipeline/README.md)`). When the target lacks a docs page, link to the GitHub source rather than a non-resolving relative path.
