# WARP WDL Style Guide Outline

## 1. Introduction

   - Purpose and scope of the style guide
   - Relationship to other WARP documentation
  - How this style guide differs from other WDL style guides (BioWDL, etc.)

## 2. Repository Structure

- Flattened directory organization
  - `pipelines/wdl/` for workflow definitions by category
  - `tasks/wdl/` for reusable task definitions
  - `structs/` for struct definitions
  - `verification/` for test workflows
  - `verification/test-wdls/` for test implementation
- Naming conventions for directories and files
- Import path management with the unified structure

## 3. WDL Basics

- Version declaration (version 1.0)
- File organization hierarchy
- Import statements with relative paths
- Meta parameters for Terra compatibility

## 4. Code Formatting

- Indentation: 2 spaces (WARP-specific, differs from BioWDL)
- Blank lines: Use to separate logical sections
- Line length: No strict limit (differs from BioWDL's 100 character limit)
- Variable naming: 
  - Tasks: Use UpperCamelCase
  - Variables: Use lowercase_underscore (python style, differs from BioWDL)
  - Alias calls: Use UpperCamelCase

## 5. Workflow Structure

   - Required elements:
     - Version string declaration: `String pipeline_version = "x.y.z"`
     - Meta parameters: `allowNestedInputs: true`
   - Input block organization
     - Required inputs first
     - Optional inputs with defaults second
     - Runtime configuration parameters last
   - Cloud provider support inputs and validation
   - Output block organization
   - Error handling with ErrorWithMessage task

## 6. Task Structure

- Standard order of sections:
  - Input block
  - Command block
  - Output block
  - Runtime block
- Input parameters organization
- Command formatting
  - One input argument per line for clarity (WARP-specific)
- Runtime parameter specification
  - Docker image handling for multi-cloud
  - Memory, CPU, disk sizing

## 7. Docker Handling

- Multi-cloud support pattern

## 8. WARP-Tools Integration

- Purpose of the warp-tools repository
  - Centralized location for Docker images and scripts used by WARP workflows
  - Separation of concerns between workflow definitions and execution environments
- When to use warp-tools vs inline WDL code
  - Use warp-tools for:
    - Complex processing scripts (Python, R, etc.)
    - Functionality that may be reused across multiple pipelines
    - Scripts requiring specialized dependencies
    - Tasks involving data transformations or complex algorithms
  - Use inline WDL for:
    - Simple file operations
    - Basic string manipulation
    - Parameter validation
    - Straightforward conditional logic
- Referencing warp-tools in WDL
  - Docker image referencing pattern
  - Path conventions for scripts within containers
- Contribution workflow for adding new tools
  - Docker build guidelines
  - Script documentation requirements
  - Versioning practices

## 9. Python Code in WDL

- Prefer scripts in Docker images over embedded code
  - Embedding complex scripts in WDL command blocks reduces maintainability
  - Package scripts into Docker images whenever possible
  - Benefits: improved testability, readability, and maintainability
- When embedding Python is necessary
  - Keep Python code minimal and focused
  - Place imports at the module level for better readability
  - Avoid nesting function definitions with their own imports
  - Use clear variable naming and add comments for complex logic
- Avoiding Python/WDL syntax conflicts
  - Beware of mixing Python f-strings with WDL interpolation
  - Both use curly braces `{}` which can cause confusion
  - Example: `f'~{prefix}_output.txt'` - the f-string is unnecessary and potentially problematic
  - Prefer simple string interpolation: `'~{prefix}_output.txt'`
- Python best practices in WDL context
  - Validate inputs early
  - Use meaningful variable names
  - Break complex operations into clear, commented steps
  - Consider readability in the WDL context (triple quotes, indentation)

## 10. Version Management

- Semantic versioning requirements
  - Major version: Breaking changes to inputs/outputs or qualitative scientific changes
  - Minor version: Non-breaking additions or non-qualitative changes
  - Patch version: Internal improvements, optimizations, comments
- Version declaration in WDL files
- Version tracking in `pipeline_versions.txt`
- Version consistency between WDL and changelog

## 11. Changelog Requirements

- File naming: `<PipelineName>.changelog.md`
- Entry format:
  ```markdown
  # 1.2.5
  2025-09-19 (Date of Last Commit)

  * Updated docker image to latest version
  ```
- Use active voice, past tense
- Standard messages for common changes
- Breaking vs non-breaking change documentation

## 12. Error Handling

- Use of ErrorWithMessage task
- Input validation patterns
- Standard error messages
- Recovery strategies

## 13. Testing Framework

- Test WDL naming: `Test<PipelineName>.wdl`
- Verification workflow structure
- GetValidationInputs tasks
- Golden file management with update_truth
- Required test inputs and expected outputs

## 14. WDL Validation (Critical)

   - **Mandatory womtool validation**
   - Validation workflow steps
   - Common validation errors
   - Import path validation
   - Handling validation failures

## 15. Multi-cloud Support

- Cloud provider input parameter
- Cloud-specific docker image selection
- Resource optimization per cloud provider
- Cloud-specific runtime attributes

## 16. Release Process

- Build script usage
- Release tagging
- Master tracking
- Dockstore configuration

## 17. Common WDL Patterns

   - Cloud provider detection
   - Resource allocation
   - File handling
   - Scatter-gather implementation
   - Parameter optimization

## 18. Code Review Checklist

- Validation with womtool
- Version consistency verification
- Import path checking
- Changelog format verification
- Docker image verification
- Version update verification

## 19. Best Practices

   - Task reusability
   - Input/output naming conventions
   - Resource allocation
   - Documentation practices
   - WDL modularization
   - Performance optimization

## 20. Git Workflow and Merge Conflict Resolution

   - Conflict resolution priority
   - Version precedence rules
   - Post-merge validation
   - Version consistency validation

## 21. Common Tools and Scripts

- Scripts for automation
- Verification utilities
- Release tools
- Import path fixing

## 22. Appendix
   - Example workflow templates
   - Example task templates
   - Common utility tasks
   - Reference to WDL specification
   - Glossary of terms