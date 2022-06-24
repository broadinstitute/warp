# Testing in WARP

## :book: Table of Contents

- [Testing in WARP](#testing-in-warp)
  - [:book: Table of Contents](#book-table-of-contents)
  - [:dna: Overview](#dna-overview)
    - [Why test?](#why-test)
  - [:page_with_curl: Terms and Concepts](#page_with_curl-terms-and-concepts)
      - [Main Workflow/Pipeline](#main-workflowpipeline)
      - [Wrapper Workflow](#wrapper-workflow)
      - [Validation Workflow](#validation-workflow)
      - [Plumbing vs Scientific data](#plumbing-vs-scientific-data)
      - [Truth Branches](#truth-branches)
  - [:hammer_and_pick: Testing Process](#hammer_and_pick-testing-process)
    - [Scala Framework](#scala-framework)
      - [Command Format](#command-format)

## :dna: Overview 

This README serves as a summary of the testing process for the pipelines that exists in WARP.

### Why test?

Downstream users of these pipelines expect robustness and consistency. When a change is made to a pipeline it is critical to ensure that the change does not affect the scientific validity of pipeline outputs.



## :page_with_curl: Terms and Concepts

#### Main Workflow/Pipeline

Every pipeline that exists in the [pipelines/](../pipelines/) directory of WARP can be though of as a "main" workflow. These are the components that are actually being tested. 

#### Wrapper Workflow

Every main pipeline has an accompanying "wrapper" workflow in the [verification/test-wdls/](../verification/test-wdls/) directory. These wrapper workflows are designed to run the main workflow and pass its outputs along to the validation workflow.

#### Validation Workflow

The validation workflow takes the output of completed pipeline and compares it against a known truth set for that pipeline. As with the wrapper workflow, each main workflow must also have a [validation workflow](../verification/).

#### Plumbing vs Scientific data

Each main workflow should have two sets of test data, plumbing and scientific. The plumbing test data is typically a set of downsampled inputs designed to test that the pipeline runs to completion without error. The scientific test data is a set of full-sized inputs designed to thoroughly test the pipeline.

These test inputs exist alongside each pipeline in their respective `test_inputs` directory.

#### Truth Branches

The testing framework has a concept of "truth branches". This allows users to define a branch of truth output to run against.

In most cases the branch being run against is "master" which corresponds to a google storage bucket containing the agreed upon output of the pipeline as it exists in the "master" git branch.

## :hammer_and_pick: Testing Process

We use Jenkins to run our test suite with two webhooks to kick off tests: PR to dev branch kicks off plumbing tests, PR to main/master kicks off scientific tests.

Jenkins will run a job called `Smart Tests` which can be viewed in the [dsp-jenkins](https://github.com/broadinstitute/dsp-jenkins/blob/master/jobs/gotc-jenkins/WarpSmartTestJob.groovy) repo.

The job of the smart tests is to run the `get_changed_pipeline_worklow_test_args.sh` bash script. The script performs a git diff against the target branch and recursively searches for any pipelines that have been changed either directly or via their imports.

For each pipeline that has changed it will prepare the arguments for the scala framework that submits the wrapper workflow to cromwell.

Finally, it passes these argument to the downstream [Workflow Test](https://github.com/broadinstitute/dsp-jenkins/blob/master/jobs/gotc-jenkins/WarpWorkflowTestsJob.groovy) job which calls the Scala `CloudWorkflow` command and submits the wrapper workflow.

### Scala Framework

As mentioned above, the test framework leverages Scala as a CLI tool to submit the wrapper workflows to Cromwell. The code for the Scala tool can be found in the [tests/broad/scala_test](../tests/broad/scala_test/) directory.

The Scala tool is responsible for preparing the inputs for the wrapper workflow for the given pipeline and submitting that workflow to the cromwell environment specified.

#### Command Format

Command : CloudWorklow

Parameters: 
    - p : Pipeline to test (any) - *required*
    - e : Cromwell environment to run (staging/dev) - *required*
    - t : Type of test to run (Plumbing/Scientific) - *required*
    - b : Truth branch to run against (master/develop/any) -r *required*
    - u : Run the main workflow uncached - *optional*
    - update_truth: Skip validation and just update the truth set specified - *optional*

```bash 
# Run the arrays plumbing test against the master truth set in the cromwell dev environment
$ CloudWorkflow -p Arrays -e dev -b master -t Plumbing

# Update the master truth branch for the arrays plumbing test in the cromwell dev environment
$ CloudWorkflow -p Arrays -e dev -b master -t Plumbing --update-truth

# Update the develop truth branch for the WGS scientific test in the cromwell staging environment
$ CloudWorkflow -p WholeGenomeGermlineSingleSample -e staging -b develop -t Scientific --update-truth

# Run the WGS scientific test against the develop truth set in the cromwell staging environment
$ CloudWorkflow -p WholeGenomeGermlineSingleSample -e staging -b develop -t Scientific 
```