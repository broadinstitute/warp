version 1.0

# This workflow splits a multi-sample VCF file into individual single-sample VCF files. This is achieved by creating
# "chunk" of samples, processing each chunk to extract the corresponding samples from the multi-sample VCF,
# and then zipping the outputs for each chunk.

# WORKFLOW DEFINITION
workflow SplitMultiSampleVcfWorkflow {
    input {
        File multiSampleVcf
        String workspaceBucket
        String? outputFileNameSuffix

        # Optional parameters w/ defaults
        Boolean createIndexFiles = true
        Int chunkSize = 1000
        Int cpu = 1
        Int memoryMb = 6000
        String bcftoolsDocker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
        Int? diskSizeGb
    }

    String pipeline_version = "2.0.0"
    String gsutil_docker = "gcr.io/google.com/cloudsdktool/cloud-sdk:525.0.0"
    Int calculated_disk_size = ceil(21*chunkSize*size(multiSampleVcf, "GiB")/(chunkSize+20)) + 10
    Int disk_size = select_first([calculated_disk_size, diskSizeGb]) # Default disk size if not provided

    parameter_meta {
        multiSampleVcf: "Input multi-sample VCF file to be split"
        workspaceBucket: "The bucket associated with the workspace."
        outputFileNameSuffix: "If provided, the string will be appended to each output VCF file name (before the .vcf.gz extension)."
        createIndexFiles: "Whether index files should be created for each individual VCF. (default: true)."
        chunkSize: "Number of samples to process in each chunk (default: 1000)"
        cpu: "Number of CPU cores to allocate for each task (default: 1)"
        memoryMb: "Memory allocation in megabytes (default: 6000)"
        bcftoolsDocker: "Docker image containing bcftools (default: us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889)"
    }

    # ONE: Extract samples from the provided multi-sample VCF
    call ExtractSamplesFromMultiSampleVcf {
        input:
            multiSampleVcf = multiSampleVcf,
            docker = bcftoolsDocker,
            cpu = cpu,
            memoryMb = memoryMb,
            diskSizeGb = disk_size
    }

    # TWO: Break up the extracted samples by the chunk size (and create an output with samples included in each chunk)
    call ProcessSampleList {
        input:
            sampleListFile = ExtractSamplesFromMultiSampleVcf.sampleIdFile,
            chunkSize = chunkSize,
            docker = bcftoolsDocker,
            cpu = cpu,
            memoryMb = memoryMb,
            diskSizeGb = disk_size
    }


    scatter (chunkedSampleFile in ProcessSampleList.sampleChunks) {
        # THREE: Extact single-sample VCFs (using the "chunked" sample list generated in the ProcessSampleList task)
        # This is scattered by the number of "chunks" generated (i.e. if 100 samples in chunks of 20,
        # this will scatter 5 wide). Optionally generates index files
        call ExtractSingleSampleVcfs {
            input:
                chunkSize = chunkSize,
                chunkedSampleFile = chunkedSampleFile,
                multiSampleVcf = multiSampleVcf,
                createIndexFiles = createIndexFiles,
                docker = bcftoolsDocker,
                cpu = cpu,
                memoryMb = memoryMb,
                diskSizeGb = disk_size
        }

        # FOUR: Copy each chunk's output files to the destination
        call RenameAndCopyFilesToDestination {
            input:
                multiSampleVcf = multiSampleVcf,
                workspaceBucket = workspaceBucket,
                outputFileNameSuffix = outputFileNameSuffix,
                vcfTar = ExtractSingleSampleVcfs.vcfTar,
                vcfIndexTar = ExtractSingleSampleVcfs.vcfIndexTar,
                createIndexFiles = createIndexFiles,
                docker = gsutil_docker,
                cpu = cpu,
                memoryMb = memoryMb,
                diskSizeGb = disk_size
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
    Int diskSizeGb
  }

  parameter_meta {
    multiSampleVcf: "Input multi-sample VCF file to be split"
    docker: "Docker image containing bcftools (default: us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889)"
    cpu: "Number of CPU cores to allocate for each task (default: 1)"
    memoryMb: "Memory allocation in megabytes (default: 6000)"
    diskSizeGb:  "Disk space in gigabytes"
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
    disks: "local-disk ${diskSizeGb} SSD"
  }
}

task ProcessSampleList {
  input {
    File sampleListFile
    Int chunkSize
    String docker
    Int cpu
    Int memoryMb
    Int diskSizeGb
  }

  parameter_meta {
    sampleListFile: "File containing ALL sample IDs extracted from the multi-sample VCF"
    chunkSize: "Number of samples to process in each chunk (default: 1000)"
    cpu: "Number of CPU cores to allocate for each task (default: 1)"
    memoryMb: "Memory allocation in megabytes (default: 6000)"
    docker: "Docker image containing bcftools (default: us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889)"
    diskSizeGb:  "Disk space in gigabytes"
  }

  command <<<
    set -e
    python3 <<CODE

    import os
    output_location = "output_chunks"
    os.mkdir(output_location)

    with open("~{sampleListFile}", "r") as f:
        samples = [line.strip() for line in f if line.strip()]

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
    disks: "local-disk ${diskSizeGb} SSD"
  }
}

task ExtractSingleSampleVcfs {
  input {
    Int chunkSize
    File chunkedSampleFile
    File multiSampleVcf
    Boolean createIndexFiles
    String docker
    Int cpu
    Int memoryMb
    Int diskSizeGb
  }

  parameter_meta {
    chunkSize: "Number of samples to process in each chunk (default: 1000)"
    chunkedSampleFile: "File containing a chunk of sample IDs to extract from the multi-sample VCF"
    multiSampleVcf: "Input multi-sample VCF file to be split"
    createIndexFiles: "Whether index files should be created for each individual VCF"
    docker: "Docker image containing bcftools (default: us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889)"
    cpu: "Number of CPU cores to allocate for each task (default: 1)"
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

    # Generate index files if requested

    if [ "~{createIndexFiles}" = "true" ]; then
      echo "Index files were requested. Generating index files for each VCF now"
      for vcf in $OUTPUT_DIR/*.vcf.gz; do
        bcftools index -tf $vcf
      done
    fi

    # ------------------------------------------------------------------------------------------------------
    # Create tarballs for output files

    echo "Creating tarball of VCF files"
    tar -czf vcfs.tar.gz -C "$OUTPUT_DIR" $(basename -a $OUTPUT_DIR/*.vcf.gz)

    if [ "~{createIndexFiles}" = "true" ]; then
      echo "Creating tarball of index files"
      tar -czf index_files.tar.gz -C "$OUTPUT_DIR" $(basename -a $OUTPUT_DIR/*.tbi)
    else
      echo "Skipping index tarball – no index files generated."
      touch index_files.tar.gz  # Empty file to avoid downstream failure
    fi

    >>>

  output {
    File vcfTar = "vcfs.tar.gz"
    File vcfIndexTar = "index_files.tar.gz"
  }

  runtime {
    cpu: cpu
    memory: "${memoryMb} MiB"
    docker: docker
    disks: "local-disk ${diskSizeGb} SSD"
    noAddress: true
  }
}

task RenameAndCopyFilesToDestination {

    input {
        File multiSampleVcf
        String workspaceBucket
        String? outputFileNameSuffix
        File vcfTar
        File vcfIndexTar
        Boolean createIndexFiles
        String docker
        Int cpu
        Int memoryMb
        Int diskSizeGb
    }

    parameter_meta {
        multiSampleVcf: "Input multi-sample VCF file to be split"
        workspaceBucket: "The bucket associated with the workspace."
        outputFileNameSuffix: "If provided, the string will be appended to each output VCF file name (before the .vcf.gz extension)."
        vcfTar: "Tar file containing all VCF file paths to be copied to output destinatino"
        vcfIndexTar: "Tar file of all VCF index file paths to be copied to output destinatino"
        createIndexFiles: "Whether index files should be created for each individual VCF"
        docker: "Docker image containing gsutil"
        cpu: "Number of CPU cores to allocate for each task (default: 1)"
        memoryMb: "Memory allocation in megabytes (default: 6000)"
        diskSizeGb:  "Disk space in gigabytes"
    }

    command <<<

    # ------------------------------------------------------------------------------------------------------

    # Generate the output directory where VCFs and index files will be copied
    BASENAME=$(basename "~{multiSampleVcf}")
    SAFE_BASENAME=$(echo "$BASENAME" | sed 's/[^A-Za-z0-9_-]/_/g')
    echo "Subdirectory where files will be written: $SAFE_BASENAME"

    # Normalize workspace bucket to ensure it starts with gs://
    WORKSPACE_BUCKET="~{workspaceBucket}"
    if [[ "$WORKSPACE_BUCKET" != gs://* ]]; then
        WORKSPACE_BUCKET="gs://$WORKSPACE_BUCKET"
    fi

    OUTPUT_DIR="${WORKSPACE_BUCKET}/${SAFE_BASENAME}_output_vcfs/"
    echo "Top-level output directory: $OUTPUT_DIR"


    # ------------------------------------------------------------------------------------------------------
    # Function to optionally rename files with a suffix and update FOFN
    rename_files_with_suffix() {
        local fofn="$1"
        local suffix="$2"
        local ext_pattern="$3"

        if [[ -n "$suffix" ]]; then
            echo "Renaming files in $fofn with suffix '${suffix}' for pattern $ext_pattern"
            local tmp_fofn="${fofn}.tmp"
            > "$tmp_fofn"

            while IFS= read -r filepath; do
                dir=$(dirname "$filepath")
                filename=$(basename "$filepath")
                # Remove the extension pattern from the end
                base="${filename%$ext_pattern}"
                new_name="${base}_${suffix}${ext_pattern#$base}"  # Keep the extension intact
                mv "$filepath" "${dir}/${new_name}"
                echo "${dir}/${new_name}" >> "$tmp_fofn"
            done < "$fofn"

            mv "$tmp_fofn" "$fofn"
        else
            echo "No suffix provided, skipping renaming for $fofn"
        fi
    }

    # ------------------------------------------------------------------------------------------------------

    # Extract VCFs
    echo "Extracting VCF tarball"
    mkdir -p vcfs
    tar -xzf ~{vcfTar} -C vcfs

    $VCF_FOFN = vcf_fofn.txt

    echo "Writing VCF files to FOFN"
    find vcfs -type f -name '*.vcf.gz' -exec realpath {} \; > $VCF_FOFN

    # Call the function to rename the output files
    rename_files_with_suffix $VCF_FOFN "~{outputFileNameSuffix}" ".vcf.gz"

    # ------------------------------------------------------------------------------------------------------
    # Optionally extract index files and write to FOFN
    if [ "~{createIndexFiles}" = "true" ]; then
        echo "Extracting index tarball"
        mkdir -p index_files
        tar -xzf ~{vcfIndexTar} -C index_files

        $VCF_INDEX_FOFN = index_fofn.txt

        echo "Writing index files to FOFN"
        find index_files -type f -name '*.tbi' -exec realpath {} \; > $VCF_INDEX_FOFN

        # Call the function to rename the output files
        rename_files_with_suffix $VCF_INDEX_FOFN "~{outputFileNameSuffix}" ".vcf.gz.tbi"
      else
        echo "No index files provided, creating empty index_fofn.txt"
        touch $VCF_INDEX_FOFN
      fi

    # ------------------------------------------------------------------------------------------------------
    # Copy VCFs
    echo "Copying VCF files to output location"
    gsutil -m cp -I "${OUTPUT_DIR}" < vcf_fofn.txt

    # ------------------------------------------------------------------------------------------------------

    # Copy index files if applicable
    if [ "~{createIndexFiles}" = "true" ]; then
        echo "Copying index files to output location"
        gsutil -m cp -I "${OUTPUT_DIR}" < index_fofn.txt
    fi

    echo "Finished copying all files to ${OUTPUT_DIR}"

    >>>

  runtime {
    cpu: cpu
    memory: "${memoryMb} MiB"
    docker: docker
    noAddress: true
    disks: "local-disk ${diskSizeGb} SSD"
  }
}