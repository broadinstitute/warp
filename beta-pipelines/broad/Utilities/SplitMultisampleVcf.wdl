version 1.0


# This workflow splits a multi-sample VCF file into individual single-sample VCF files. This is achieved by creating
# "chunk" of samples, processing each chunk to extract the corresponding samples from the multi-sample VCF,
# and then zipping the outputs for each chunk.

# WORKFLOW DEFINITION
workflow SplitMultiSampleVcfWorkflow {
    input {
        File multiSampleVcf
        Boolean createIndexFiles
        String outputLocation

        # Optional parameters w/ defaults
        Int chunkSize = 1000
        Int cpu = 1
        Int memoryMb = 6000
        String bcftoolsDocker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    }

    parameter_meta {
        multiSampleVcf: "Input multi-sample VCF file to be split"
        createIndexFiles: "Whether index files should be created for each individual VCF"
        outputLocation: "GCP location where output vcfs (and indices, if requested) are written to"
        chunkSize: "Number of samples to process in each chunk (default: 1000)"
        cpu: "Number of CPU cores to allocate (default: 1)"
        memoryMb: "Memory allocation in megabytes (default: 6000)"
        bcftoolsDocker: "Docker image containing bcftools (default: us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889)"
    }

    # Extract samples from the provided multi-sample VCF
    call ExtractSamplesFromMultiSampleVcf {
        input:
            multiSampleVcf = multiSampleVcf,
            docker = bcftoolsDocker,
            cpu = cpu,
            memoryMb = memoryMb
    }

    # Break up the extracted samples by the chunk size (and create an output with samples included in the chunk)
    call ProcessSampleList {
        input:
            sampleListFile = ExtractSamplesFromMultiSampleVcf.sampleIdFile,
            chunkSize = chunkSize,
            docker = bcftoolsDocker,
            cpu = cpu,
            memoryMb = memoryMb
    }

    # Extact single-sample VCFs (using the "chunked" sample list generated in the ProcessSampleList task)
    # This is scattered by the number of "chunks" generated (i.e. if 100 samples in chunks of 20,
    # this will scatter 5 wide). Optionally generates index files, and copies outputs to the specified output location.
    scatter (chunkedSampleFile in ProcessSampleList.sampleChunks) {
        call ProcessSampleChunkAndCopyFiles {
            input:
                chunkSize = chunkSize,
                chunkedSampleFile = chunkedSampleFile,
                multiSampleVcf = multiSampleVcf,
                createIndexFiles = createIndexFiles,
                outputLocation = outputLocation,
                docker = bcftoolsDocker,
                cpu = cpu,
                memoryMb = memoryMb
        }
    }
}

# TASK DEFINITIONS
task ExtractSamplesFromMultiSampleVcf {
  input {
    File multiSampleVcf
    String docker
    Int cpu
    Int memoryMb
  }

  parameter_meta {
    multiSampleVcf: "Input multi-sample VCF file to be split"
    docker: "Docker image containing bcftools (default: us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889)"
    cpu: "Number of CPU cores to allocate (default: 1)"
    memoryMb: "Memory allocation in megabytes (default: 6000)"
  }

  command <<<
    set -euo pipefail
    bcftools query -l ~{multiSampleVcf} > sample_ids.txt
  >>>

  output {
    File sampleIdFile = "sample_ids.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "${memoryMb} MiB"
  }
}

task ProcessSampleList {
  input {
    File sampleListFile
    Int chunkSize
    String docker
    Int cpu
    Int memoryMb
  }

  parameter_meta {
    sampleListFile: "File containing ALL sample IDs extracted from the multi-sample VCF"
    chunkSize: "Number of samples to process in each chunk (default: 1000)"
    cpu: "Number of CPU cores to allocate (default: 1)"
    memoryMb: "Memory allocation in megabytes (default: 6000)"
    docker: "Docker image containing bcftools (default: us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889)"
  }

  command <<<
    set -e
    python3 <<CODE

    import os
    output_location = "output_chunks"
    os.mkdir(output_location)

    with open(~{sampleListFile}, "r") as f:
        samples = [line.strip() for line in f if line.strip()]
        print(f"samples: {samples}")

    chunks = [samples[i:i+~{chunkSize}] for i in range(0, len(samples), ~{chunkSize})]

    # Write each chunk to a file
    for i, chunk in enumerate(chunks):
        with open(os.path.join(output_location, f"chunk_{i}.txt"), "w") as out:
            out.write("\n".join(chunk))
            print(f"Wrote chunk {i} with {len(chunk)} samples to chunk_{i}.txt")
    CODE
  >>>

  output {
    Array[File] sampleChunks = glob("output_chunks/chunk_*.txt")
  }

  runtime {
    cpu: cpu
    memory: "${memoryMb} MiB"
    docker: docker
  }
}

task ProcessSampleChunkAndCopyFiles {
  input {
    Int chunkSize
    File chunkedSampleFile
    File multiSampleVcf
    Boolean createIndexFiles
    String outputLocation
    String docker
    Int cpu
    Int memoryMb

    # This calculation is explained in https://github.com/broadinstitute/warp/pull/937
    Int diskSizeGb = ceil(21*chunkSize*size(multiSampleVcf, "GiB")/(chunkSize+20)) + 10
  }

  parameter_meta {
    chunkSize: "Number of samples to process in each chunk (default: 1000)"
    chunkedSampleFile: "File containing a chunk of sample IDs to extract from the multi-sample VCF"
    multiSampleVcf: "Input multi-sample VCF file to be split"
    createIndexFiles: "Whether index files should be created for each individual VCF"
    outputLocation: "GCP location where output vcfs (and indices, if requested) are written to"
    docker: "Docker image containing bcftools (default: us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889)"
    cpu: "Number of CPU cores to allocate (default: 1)"
    memoryMb: "Memory allocation in megabytes (default: 6000)"
    diskSizeGb:  "Disk space in gigabytes"
  }

  command <<<
    set -euo pipefail

    # ------------------------------------------------------------------------------------------------------

    # Set the output directory
    OUTPUT_DIR="output_vcfs"

    # ------------------------------------------------------------------------------------------------------

    # Create output directory
    mkdir -p $OUTPUT_DIR
    echo "Created output directory: $OUTPUT_DIR"

    # ------------------------------------------------------------------------------------------------------

    # Extract the single-sample VCFs from the multi-sample VCF using bcftools (uses the "chunked" sample file to
    # define which samples to extract)

    echo "Extracting individual level samples from multisample VCF"
    bcftools +split ~{multiSampleVcf} --samples-file ~{chunkedSampleFile} -Oz -o $OUTPUT_DIR
    echo "Finished extracting all sample-level VCFs"

    # ------------------------------------------------------------------------------------------------------

    # Ensure output location ends with a single "/"
    CLEANED_OUTPUT_LOCATION="~{outputLocation}"
    [[ "${CLEANED_OUTPUT_LOCATION}" != */ ]] && CLEANED_OUTPUT_LOCATION="${CLEANED_OUTPUT_LOCATION}/"
    echo "Using cleaned output location: ${CLEANED_OUTPUT_LOCATION}"

    # ------------------------------------------------------------------------------------------------------

    # Copy all VCFs to the final output location (optionally generate index files and copy them as well if
    # requested)

    echo "Copying VCF (and generating indices if requested) now"
    for vcf in $OUTPUT_DIR/*.vcf.gz; do
        gsutil cp $vcf ${CLEANED_OUTPUT_LOCATION}
        if [ "${createIndexFiles}" = "true" ]; then
            bcftools index -t $vcf
            gsutil cp $vcf.csi ${CLEANED_OUTPUT_LOCATION}
        fi
    done
    echo "Finished copying all files to ${CLEANED_OUTPUT_LOCATION}"

    # ------------------------------------------------------------------------------------------------------
>>>

  runtime {
    cpu: cpu
    memory: "${memoryMb} MiB"
    docker: docker
    disks: "local-disk ${diskSizeGb} SSD"
    noAddress: true
  }
}
