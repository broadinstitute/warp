version 1.0

# Workflow: run_admixture_est_rye
# Description: This workflow generate admixture estimations using Rye tool (https://github.com/healthdisparities/rye).

## Copyright Broad Institute, 2023
##
## This WDL pipeline processes a set of files:
## (file formats and examples can be found at https://github.com/healthdisparities/rye/tree/main/examples)
## - eigenvalues: This is computed from genomic data (VCF/BED/PED) using PCA
## - Eigenvectors: This is computed from genomic data (VCF/BED/PED) using PCA
## - pop2group: A population to group mapping file
##
## The pipeline, then, generates admixture estimations for all samples in these files and generates outputs
## in the format *.Q and *.fam which are analogous to the output files produced by ADMIXTURE.
##
## More details can be found at Rye tool github repo (https://github.com/healthdisparities/rye).
##
## Requirements/expectations:
## - Eigenvalue/eigenvectors computed from genomic data (VCF/BED/PED)
## - A population to group mapping file
##
## LICENSING:
## This script is released under the [Appropriate License] (e.g., MIT, BSD-3, etc.). Users are responsible
## for ensuring they are authorized to run all components of this script. Please consult the relevant
## documentation for licensing details of Hail and other tools used in this pipeline.
##
## For information on tool versions and parameters, refer to the specific Docker containers and
## configuration used in this pipeline.

struct RuntimeAttr {
    Float? mem_gb
    Int? disk_gb
    Int? boot_disk_gb
    Int? max_retries
}

workflow run_admixture_est_rye {
    input {
        # Analysis Parameters

        File eigenvalues_file  # Path to the Eigenvalues file
        File eigenvec_file  # Path to the Eigenvectors file
        File pop2group_file  # Path to the mapping file of populations and sub-populations (If any)
        String prefix # A prefix appended to the outputs as a task identifer (e.g. aou_delta)
        Int pcs = 20 # Number of PCs to use (Default = 20)
        Int rounds = 200 # Number of rounds to use for optimization (higher number = more accurate but slower; Default=200)
        Int iter = 100 # Number of iterations to use for optimization (higher number = more accurate but slower; Default=100)

        # VM Parameters
        Int cpus = 16
        # Docker image with Rye tool installed
        String docker_image = "us-central1-docker.pkg.dev/broad-dsde-methods/aou-auxiliary/rye-admixture-estimation-tool:v1.0"
    }
    String pipeline_version= "aou_9.0.0"


    call run_rye {
        # Task inputs mirror workflow inputs
        input:
            eigenvalues_file = eigenvalues_file,
            eigenvec_file = eigenvec_file,
            pop2group_file = pop2group_file,
            prefix = prefix,
            cpus = cpus,
            pcs = pcs,
            rounds = rounds,
            iter = iter,
            docker_image = docker_image
    }
    output {
        File qFile = run_rye.qFile
        File famFile = run_rye.famFile
    }
}

task run_rye {
    input {
        # Task-specific inputs with descriptions
        File eigenvalues_file
        File eigenvec_file
        File pop2group_file
        String prefix
        Int cpus
        Int pcs
        Int rounds
        Int iter
        RuntimeAttr? runtime_attr_override
        String docker_image
    }

    RuntimeAttr runtime_default = object {
                                      # Default runtime attributes
                                      mem_gb: 6.5,
                                      disk_gb: 100,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    command <<<
        set -euxo pipefail

        echo "#### Start Admixture estimation using Rye tool ###"

        # run Rye tool
        /rye/rye.R --eigenval=~{eigenvalues_file} --eigenvec=~{eigenvec_file} --pop2group=~{pop2group_file} --output=~{prefix} --threads=~{cpus} --pcs=~{pcs} --rounds=~{rounds} --iter=~{iter}

        # rename the Q file (easier for delocalization of output file. Change: removing the number of groups from the file name)
        mv /cromwell_root/*Q /cromwell_root/~{prefix}-~{pcs}.Q

        echo "#### DONE ####"
    >>>


    output {
        File qFile = "~{prefix}-~{pcs}.Q"
        File famFile = "~{prefix}-~{pcs}.fam"
    }

    runtime {
        # Runtime settings for the task
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: cpus
        docker: docker_image

    }
}
