---
sidebar_position: 3
---

# Testing Pipelines

WARP pipelines have accompanying automated tests that run on each Pull Request (PR). These tests compare validated outputs to the new PR outputs. For any changes in code shared between pipelines, the tests also confirm which pipelines could be affected and ensure that the PR makes no unexpected changes to the affected pipelines.

:::tip NOTICE 9/29/2020
We have formatted all testing input files for Cromwell 52 or higher. If you are using WARP test input files either directily or as a model of how to configure inputs, these changes may impact you.
:::

## WARP branches and how they relate to testing

WARP has three main branches that are used for different stages of pipeline testing: **develop**, **staging**, and **master**:

| Branch | Purpose |
| --- | --- |
| develop |  Code that has passed plumbing tests; pre-releases for testing |
| staging | Release candidate code that is frozen for longer scientific testing |
| master | Code that has passed scientific testings; published release artifacts |

## Initiating pipeline tests

A PR will initiate a pipeline test if it contains changes to the main workflow WDL, the WDL dependencies (such as tasks), the options JSON file, the pipeline tests, or the test inputs. Smart plumbing and scientific tests compare PR changes to a specified base branch. They then use the WARP directory structure to dynamically determine which pipelines are affected and run tests for all affected pipelines.


## Contact us with questions about testing

If you have questions about the WARP testing process, please reach out to [the WARP team](mailto:warp-pipelines-help@broadinstitute.org).
