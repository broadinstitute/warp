# Nextflow Workflows

This directory is reserved for future Nextflow workflows as part of the WARP repository reorganization.

## Structure

When Nextflow workflows are added, they should follow the same organizational pattern as the WDL workflows:

```
pipelines/nextflow/
├── category1/
│   └── workflow1/
├── category2/
│   └── workflow2/
└── ...
```

## Guidelines

- Follow the same naming conventions as WDL pipelines
- Include comprehensive README.md files for each workflow
- Maintain changelog files for version tracking
- Include example input configuration files

## Migration from WDL

When migrating WDL workflows to Nextflow:
1. Preserve the same directory structure pattern
2. Update import paths to reference appropriate task libraries
3. Maintain semantic versioning
4. Update changelog with migration notes

---

For questions about Nextflow workflow development in WARP, please refer to the main [WARP documentation](../../README.md) and [AI coding guidelines](../../.github/copilot-instructions.md).
