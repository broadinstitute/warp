# The by_chromosome Joint Genotyping pipeline is now deprecated 2025-03-06
# Joint Genotyping By Chromosome

## Background
Joint calling on large numbers of samples (> 10,000) stretches the limits of what one google project can handle, so we've decided to split the computation up by chromosome. In this directory, you will find a script called `by_chromosome_client.sh` that will help you manage the submission and monitoring of running each chromosome.

## Using `by_chromosome_client.sh`
The script makes submitting and monitoring the chromosome workflows easy. It simplifies the management of running the 23 chromosomes by allowing you to perform operations by chromosome instead of by workflow ID. It keeps track of workflow IDs behind the scenes in a local file.

Some of the commands in the script take a chromosome as an argument, specified in the form of `-c chrN`, where `N` is `1-22`, `X`, `Y`, or `part-two`. `part-two` is used to track part two of the joint genotyping.

Only one workflow per chromosome can be tracked at a time. Chromosomes will likely take more than one run because of transient errors, so the log only keeps track of the most recent workflow for each chromosome.

The available commands are:
  - `setup` - Setup symlinks to the inputs and script in the current directory.
  - `submit-part-one -c chrN [-p P]` - Submit `chrN` to cromwell for part one of joint genotyping. Optionally, a number `P` can be given to specify a google project to run on. If not given, the script will choose a random project.
  - `submit-part-two [-p P]` - Submit the results of part-one to cromwell for part two of joint genotyping. Optionally, a number `P` can be given to specify a google project to run on. If not given, the script will choose a random project.
  - `status -c chrN` - Get the status of the current workflow for `chrN`.
  - `monitor [-x] [-c chrN]` - Monitor currently running workflows for all chromosomes. Providing `-x` will expand sub-workflow monitoring. Providing `-c chrN` will monitor only the provided chromosome instead of all running ones.
  - `metadata -c chrN [-x]` - Get the metadata of the current workflow for `chrN`. Providing `-x` will expand sub-workflow metadata.
  - `abort -c chrN` - Abort the workflow for `chrN` so another workflow for `chrN` can be run.
  - `void -c chrN` - Allow another workflow to be run for `chrN`, but don't abort the current workflow.
  - `outputs -c chrN` - Get the outputs of a successful `chrN` workflow.
  - `copy-part-one-outputs` - Copy the outputs of part one out of the execution directory and into a safe location.
  - `copy-part-two-outputs` - Copy the outputs of part two out of the execution directory and into a safe location.
  - `get-failed` - Get all the chromosomes that have failed.
  - `get-succeeded` - Get all the chromosomes that have succeeded.
  - `get-running` - Get all the chromosomes that are still running.
  - `chromosome-cost -c chrN` - Calculate the cost of running `chrN`.
  - `total-cost` - Calculate the total cost of all workflows run.

  
The `monitor` command will probably be the most-used one. It allows you to get a high-level overview of all the chromosomes running in one command. It breaks the workflows down by task, showing how many shards are running, completed, preempted, and failed.

There might be some failures from PAPI errors. The Joint Calling workflow is fully call-cacheable, which means restarting a workflow should pick up from where it left off pretty quickly.
Simply `abort` or `void` the chromosome and resubmit it.

## Joint Calling by Chromosome
If all goes well, and there are no bugs in the pipeline, joint calling by chromosome should be pretty simple.

  0. Clone `dsde-pipelines` into the PO directory using `git clone --depth 1 git@github.com:broadinstitute/dsde-pipelines.git`
  1. Run `dsde-pipelines/pipelines/wdl/dna_seq/germline/joint_genotyping/by_chromosome/by_chromosome_client.sh setup [Exome|WGS]` to set up your directory.
  2. Fill in the constants at the top of `./by_chromosome_client.sh`, which should now be in your directory.
  3. Make sure that this branch of `dsde-pipelines` is never merged. A new branch should be checked out for every run. 
  4. In `./by_chromosome_client.sh`, fill in the variables in the `## PROJECT-SPECIFIC CONSTANTS` section.
  5. For every chromosome being jointly called, run `by_chromosome_client.sh submit-part-one -c chrN`.
  6. Monitor the chromosomes. Re-submit ones that fail. Consult the Lantern team for repeat failures. Do this until all per-chromosome workflows succeed.
  7. Run `./by_chromosome_client.sh copy-part-one-outputs` to copy the part one outputs to a non-volatile location
  8. Run `./by_chromosome_client.sh submit-part-two`
  9. Monitor part two and ensure it finishes, much like part one. Monitoring part two is easily done with `./by_chromosome_client.sh monitor -c part-two`.
  10. Once part two completes, run `./by_chromosome_client.sh copy-part-two-outputs`.

## Becoming `picard-prod` for Submitting Part Two and Copying Cloud Outputs
To copy things around in `broad-gotc-prod-storage` or `broad-exomes-prod-storage`, the `picard-prod` service account needs to be used.
To authenticate as `picard-prod`, run `gcloud auth activate-service-account --key-file <(vault read -format=json -field=data secret/dsde/gotc/prod/picard/picard-account.pem)`.
After this, you will be authenticated as `picard-prod` and the `submit-part-two`, `copy-part-one-outputs`, and `copy-part-two-outputs` commands can be run.

To authenticate as your personal user again, run `gcloud config set account ${BROAD_EMAIL_ADDRESS}`. 
To switch back to `picard-prod`, run `gcloud config set account picard-prod@broad-gotc-prod.iam.gserviceaccount.com`

## Some Notes
The script logs its workflows in the local file `workflows.csv`, so care should be taken to make sure this file is not in an inconsistent state. Only one person should be running the script, and if `workflows.csv` is in source control, care must be taken to make sure that the most up-do-date version of the file is pushed.

The script uses inputs and options template files, so a change in `JointGenotyping.template.input.json` or `JointGenotyping.template.options.json` will affect all future chromosome workflows

`chromosome_scatter_counts.csv` is a hand-curated table of scatter values. They are balanced based on the contig lengths from the sequence dictionary and ensure we don't scatter too wide in Cromwell. We aim to scatter 1000 wide over the entire genome, so please don't change these.

`split_intervals.sh` was used to split a full interval list by chromosome.