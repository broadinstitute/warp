version 1.0

## WORKFLOW DEFINITION
## This workflow splits a multi-sample VCF file into individual single-sample VCF files
workflow SplitMultiSampleVcfWorkflow {
    input {
        File multiSampleVcf
        Int nSamples

        # Optional parameters with defaults
        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 1
        Int memory_mb = 6000
    }

    # Parameter validation
    parameter_meta {
        multiSampleVcf: "Input multi-sample VCF file to be split"
        nSamples: "Number of samples in the input VCF file"
        bcftools_docker: "Docker image containing bcftools"
        cpu: "Number of CPU cores to allocate"
        memory_mb: "Memory allocation in megabytes"
    }

    # Call the SplitMultiSampleVcf task
    call SplitMultiSampleVcf {
        input:
            multiSampleVcf = multiSampleVcf,
            nSamples = nSamples,
            bcftools_docker = bcftools_docker,
            cpu = cpu,
            memory_mb = memory_mb
    }

    # Workflow outputs
    output {
        Array[File] single_sample_vcfs = SplitMultiSampleVcf.single_sample_vcfs
        #Array[File] single_sample_vcf_indices = SplitMultiSampleVcf.single_sample_vcf_indices
    }
}

## TASK DEFINITION
task SplitMultiSampleVcf {
    input {
        File multiSampleVcf
        Int nSamples

        String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int cpu = 1
        Int memory_mb = 6000

        # This calculation is explained in https://github.com/broadinstitute/warp/pull/937
        Int disk_size_gb = ceil(21*nSamples*size(multiSampleVcf, "GiB")/(nSamples+20)) + 10
    }

    parameter_meta {
        multiSampleVcf: "Input multi-sample VCF file"
        nSamples: "Number of samples in the VCF file"
        bcftools_docker: "Docker image with bcftools installed"
        cpu: "Number of CPU cores"
        memory_mb: "Memory in megabytes"
        disk_size_gb: "Disk space in gigabytes"
    }

    command <<<
        set -e -o pipefail

        mkdir out_dir
        bcftools +split ~{multiSampleVcf} -Oz -o out_dir
        #for vcf in out_dir/*.vcf.gz; do
        #  bcftools index -t $vcf
        #done
    >>>

    runtime {
        docker: bcftools_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        noAddress: true
    }

    output {
        Array[File] single_sample_vcfs = glob("out_dir/*.vcf.gz")
        #Array[File] single_sample_vcf_indices = glob("out_dir/*.vcf.gz.tbi")
    }
}