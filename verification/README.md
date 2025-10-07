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
      - [Results vs Storage buckets](#results-vs-storage-buckets)
  - [:hammer_and_pick: Testing Process](#hammer_and_pick-testing-process)
    - [Scala Framework](#scala-framework)
      - [Command Format](#command-format)
    - [What does the wrapper workflow do?](#what-does-the-wrapper-workflow-do)
  - [:mag_right: Working with tests locally](#mag_right-working-with-tests-locally)
    - [Creating a new test](#creating-a-new-test)
    - [Running your new test](#running-your-new-test)

## :dna: Overview 

This README serves as a summary of the testing process for the pipelines that exist in WARP.

### Why test?

Downstream users of these pipelines expect robustness and consistency. When a change is made to a pipeline it is critical to ensure that the change does not affect the scientific validity of pipeline outputs.

## :page_with_curl: Terms and Concepts

#### Main Workflow/Pipeline

Every pipeline that exists in the [pipelines/](../pipelines/) directory of WARP can be though of as a "main" workflow. These are the components that are actually being tested. 

#### Wrapper Workflow

Every main pipeline has an accompanying "wrapper" workflow in the [verification/test-wdls/](../verification/test-wdls/) directory. These wrapper workflows are designed to run the main workflow and pass its outputs along to the validation workflow.

#### Validation Workflow

The validation workflow takes the output of a completed pipeline and compares it against a known truth set for that pipeline. As with the wrapper workflow, each main workflow must also have a [validation workflow](../verification/).

#### Plumbing vs Scientific data

Each main workflow should have two sets of test data, plumbing and scientific. The plumbing test data is typically a set of downsampled inputs designed to test that the pipeline runs to completion without error. The scientific test data is a set of full-sized inputs designed to thoroughly test the pipeline.

These test inputs exist alongside each pipeline in their respective `test_inputs` directory.

#### Truth Branches

The testing framework has a concept of "truth branches". This allows users to define a branch of truth output to run against.

In most cases the branch being run against is "master" which corresponds to a google storage bucket containing the agreed upon output of the pipeline as it exists in the "master" git branch.

#### Results vs Storage buckets

Every workflow run in Cromwell has its output stored in an "execution" bucket. When running the pipeline tests it's important that we store this ouptut in a location that we can later reference if needed.

All test workflows that are run will have their output copied to a "results" bucket -> `gs://broad-gotc-test-results/`

When we want to update the "truth" branches for our pipelines the test workflow will copy the run output to the respective "truth" bucket -> `gs://broad-gotc-test-storage/` (this is achieved with the `--update-truth` flag).

## :hammer_and_pick: Testing Process

We use Jenkins to run our test suite with two webhooks to kick off tests: PR to dev branch kicks off plumbing tests, PR to main/master kicks off scientific tests.

Jenkins will run a job called `Smart Tests` which can be viewed in the [dsp-jenkins](https://github.com/broadinstitute/dsp-jenkins/blob/master/jobs/gotc-jenkins/WarpSmartTestJob.groovy) repo.

The job of the smart tests is to run the `get_changed_pipeline_worklow_test_args.sh` bash script. The script performs a git diff against the target branch and recursively searches for any pipelines that have been changed either directly or indirectly via an imported sub-workflow.

For each pipeline that has changed it will prepare the arguments for the Scala tool that submits the wrapper workflow to cromwell.

Finally, it passes these arguments to the downstream [Workflow Test](https://github.com/broadinstitute/dsp-jenkins/blob/master/jobs/gotc-jenkins/WarpWorkflowTestsJob.groovy) job which calls the Scala `CloudWorkflow` command and submits the wrapper workflow.

### Scala Framework

As mentioned above, the test framework leverages Scala as a CLI tool to submit the wrapper workflows to Cromwell. The code for the Scala tool can be found in the [tests/broad/scala_test](../tests/broad/scala_test/) directory.

The Scala tool is responsible for preparing the inputs for the wrapper workflow of the given pipeline and submitting that workflow to the Cromwell environment specified.

**NOTE** - When looking at the the scala code you will see that there are multiple 'tester' files. A major rewrite of the framework was done to only use one tester command (CloudWorkflow) and to pass in the desired pipeline as an argument. This change makes it easier for contributors and maintainers of these workflows to write tests by pushing the functionality of the test framework away from Scala and into the wrapper workflows.

#### Command Format

```bash
Command: CloudWorkflow [options]
Test a cloud workflow
  -p <value> | --pipeline <value>
        The pipeline to test
  -t <value> | --test <value>
        The type of test to run
  -b <value> | --branch <value>
        The branch of truth data to test against (Defaults to develop)
  -e <value> | --env <value>
        The environment that this should run in [test|staging|dev|prod]
  --update-truth
        Update the truth data with the results of this run.
  -u | --uncached
        Disable call-caching for the main workflow



# Run the arrays plumbing test against the master truth set in the cromwell dev environment
$ CloudWorkflow -p Arrays -e dev -b master -t Plumbing

# Update the master truth branch for the arrays plumbing test in the cromwell dev environment
$ CloudWorkflow -p Arrays -e dev -b master -t Plumbing --update-truth

# Update the develop truth branch for the WGS scientific test in the cromwell staging environment
$ CloudWorkflow -p WholeGenomeGermlineSingleSample -e staging -b develop -t Scientific --update-truth

# Run the WGS scientific test against the develop truth set in the cromwell staging environment
$ CloudWorkflow -p WholeGenomeGermlineSingleSample -e staging -b develop -t Scientific 
```

### What does the wrapper workflow do?

As mentioned previously, each main workflow will have an accompanying wrapper workflow to act as a test harness. The steps that this workflow takes are the following:

1. Execute the main workflow
2. Collect the outputs from the main workflow, these are typically separated into 'regular' outputs and metrics based outputs
3. Copy these outputs to the *results* bucket
4. If updating the truth for the pipeline then copy the outputs to the *truth* bucket and finish
5. Otherwise, get the location for the test and truth inputs
6. Call the validation workflow to compare the current test outputs against the known truth outputs

## :mag_right: Working with tests locally

### Creating a new test

If you are a contributor and would like to add your own test we have a [script](../verification/test-wdls/scripts/) to generate a wrapper workflow as long as the main and validation workflow already exist (it will generate ~95%, inputs will need to be provided for each call to the GetValidationInputs tasks).

Finally you can add your new test to the scala framework by editing the [PipelineTestType](../tests/broad/scala_test/src/main/scala/org/broadinstitute/dsp/pipelines/commandline/PipelineTestType.scala) file and providing the wrapper workflow name, the main workflow name and the path to that workflow relative to the `pipelines/` directory.

Example for *Arrays* workflow:

```bash

$ python3 gentest.py --workflow Arrays --validation VerifyArrays

$ INFO:root: - Workflow found at path -> broadinstitute/warp/pipelines/wdl/arrays/single_sample/Arrays.wdl
$ INFO:root: - Validation wdl found at path -> broadinstitute/warp/verification/VerifyArrays.wdl
$ INFO:root: - Generating inputs for TestArrays.wdl...
$ INFO:root: - Generating subworkflow inputs for Arrays.wdl...
$ INFO:root: - Collecting and formatting inputs for VerifyArrays.wdl...
$ INFO:root: - Building TestArrays.wdl...
$ INFO:root: - Successfully generated TestArrays.wdl!
```

### Running your new test

If you are connected to the Broad VPN you can run the Scala tool locally to build and submit the job to cromwell.

From the root directory of WARP:

```bash

$ cd tests/broad/scala_test
$ sbt
$ run CloudWorkflow -p Arrays -e dev -t Plumbing -b master

[info] Running org.broadinstitute.dsp.pipelines.WorkflowTest CloudWorkflow -p Arrays -e dev -t Plumbing -b master
[info] 2022-06-27T13:33:17,639Z [INFO] Running the TestArrays workflow using Plumbing data
[info] 2022-06-27T13:33:17,929Z [INFO] Generating WDL inputs for -> SimpleInput.json
[info] 2022-06-27T13:33:18,693Z [INFO] Submitting TestArrays job to Cromwell
[info] 2022-06-27T13:33:21,350Z [INFO] Awaiting completion of 1 workflows
[info] 2022-06-27T13:33:21,472Z [INFO] A timing diagram for workflow dev_SimpleInput will be available at https://cromwell.gotc-dev.broadinstitute.org/api/workflows/1/7d3a91d1-3364-4921-965b-b40bc6842155/timing
[info] 2022-06-27T13:33:26,862Z [INFO] Workflows with status of Submitted: [7d3a91d1-3364-4921-965b-b40bc6842155]
[info] 2022-06-27T13:33:57,296Z [INFO] Workflows with status of Running: [7d3a91d1-3364-4921-965b-b40bc6842155]
```