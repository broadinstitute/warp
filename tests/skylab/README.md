This folder contains testing infrastructure for the pipelines in the skylab repository. The infrastructure allows for test execution via circleci on a cromwell server for multiple pipelines present in this repository. The tests in this directory can be run manually by anyone with an access to a cromwell server. Futhermore, the tests tests not marked as manual are  executed upon creation of a pull request (PR) to the repository. The majority of the datasets shown below are obtained by subsetting of full experimental datasets in order to reduce run time. Manual tests utilize full datasets and can be used to verify the output of the pipelines on well described experimental datasets. Some of these datasets have been used on validation reports of the respective pipelines.

**Description of tests and optional test variants**
- optimus: Tests related to the single-cell RNA-seq version of Optimus
  - 4kpbmc: Runs the full 4kpbmc dataset with the sc variant of Optimus. Checksums are not provided (manual)
  - pr: Run a small dataset derived from the 4k PBMC dataset with Optimus. It performs delta testing on the output matrix
  - prV3: Runs a small dataset derived from the 10k PBMC datasets with Optimus in V3 mode. Performs delta testing on the output matrix
- optimus_mouse
  - pr: Runs a small dataset derived from the 1k heart dataset with optimus V2 in single-cell mode. Checks output checksum
- optimus_snrna
  - pr: Runs a small dataset derived from the 1k heart dataset with optimus V2 in single-cell mode. Checks output checksum
  - 2kBrainNucleiAdultMouse: Runs the snRNA-seq 2k AdultBrain Nuclei dataset with the snRNA-seq version of Optimus (manual)
  - 4kpbmc_manual: Runs the *scRNA-seq* 4k PBMC dataset with the *snRNA-seq* version of Optimus (manual)
- sc_atac
  - pr: runs a small dataset with the sc atac-seq pipeline and verifies the output
- smartseq2_multisample
  - pr: Runs 22 cells throught the multi-sample SS2 pipeline in paired-end mode and checks the output files
  - pr_single_end: Runs 22 cells throught the multi-sample SS2 pipeline in single-end mode and checks the output files
- smartseq2_single_sample
  - pr: Runs a single-cell throught the SS2 pipeline in paired-end mode and checks the output files
  - pr_single_end: Runs a single-cell throught the SS2 pipeline in single-end mode and checks the output files

**Running Tests Manually**
1) Set the environment variable BROAD_CROMWELL_KEY to have the contents of the cromwell credentials key JSON file.
2) cd to the top level directory of the repository
3) run the following command ```./test/trigger_test.sh test_name [test_variant]```, where test_name is a directory name (e.g. ```optimus```) and test_variant is a nested directory name (e.g. ```4kpbmc```). test_variant is optional, if not provided the ```pr``` test will be run.

Example:
```
export BROAD_CROMWELL_KEY=`cat ~/identities/credentials.json`
cd ~/skylab/
./test/trigger_test.sh optimus 4kpmc
```

**PR testing infrastructure in Skylab.**
The tests follow the portability spec; the WDLs in this PR represent a cromwell test case for the [portability test](https://docs.google.com/document/d/1ghLoHMbKOPsndA1WgdSAHm5X82p86ryLBiAt1hz6HuI/edit); they are just missing a docker environment to run the WDLs in. 

The PR test is designed to answer the question: "did my changes result in any change to the outputs of my pipeline?"

**Testing Logic & Files**
All files related to testing are in the test/ repository
- `tests/trigger_test.sh` starts a `quay.io/broadinstitute/cromwell-tools` docker and launches tests on caas cromwell
- `tests/test_cromwell_workflow.sh` is called by `tests/trigger_test.sh` within the cromwell-tools docker. It calls the provided wdl, waits for it to finish, and checks it's success state. 
- The scientific and PR tests both contain an "infrastructure wdl" which follows the form `test_${PIPELINE_FOLDER_NAME}_PR.wdl` and is composed of two parts:

  - The pipeline workflow, which is imported and run to completion with call caching, which allows the test to run quickly when nothing has changed that would modify the pipeline. 
  - The checker workflow, which, given the outputs of the just-run pipeline, tests them either against the expected outputs (PR) or a range of acceptable outputs (scientific) 

The tests are designed so that any failure results in `exit 1` within the workflow, after the reason is logged to stderr. If a workflow fails due to a result not matching expectations, the result will be logged in the stderr produced by the checker task.

If a workflow fails due to not running properly, the user will need to go find the part of the workflow that failed. The jenkins job and cromwell backend will both log the workflow ID so that it is easy to go back and find the failed workflow.

**Updating workflows that result in scientific changes:**
In the event that a PR produces a desirable scientific change, the PR test will fail. In order to ensure that future PRs test against the correct "expected values", the test inputs must be updated to reflect the success state of the incoming pipeline. For optimus, this means updating the md5 hashes corresponding to the new output files, which are found in: 
`test/optimus/pr/test_inputs.json` (PR)

More complex tests may require different inputs, but the expected values should always be parameterized in the `test_inputs.json` files to facilitate easy updating. 
