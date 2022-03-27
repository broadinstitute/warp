version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines utility tasks used for processing of sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  input {
    File ref_dict
    Int preemptible_tries
  }
  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python3 <<CODE
    with open("~{ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    memory: "2 GiB"
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  input {
    File interval_list
    Int scatter_count
    Int break_bands_at_multiples_of
  }

  command <<<
    set -e
    mkdir out
    java -Xms1000m -Xmx1500m -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=~{break_bands_at_multiples_of} \
      INPUT=~{interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int(stdout())
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
    memory: "2000 MiB"
  }
}

# Convert BAM file to CRAM format
# Note that reading CRAMs directly with Picard is not yet supported
task ConvertToCram {
  input {
    File input_bam
    File ref_fasta
    File ref_fasta_index
    String output_basename
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(2 * size(input_bam, "GiB") + ref_size) + 20

  command <<<
    set -e
    set -o pipefail

    samtools view -C -T ~{ref_fasta} ~{input_bam} | \
    tee ~{output_basename}.cram | \
    md5sum | awk '{print $1}' > ~{output_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    seq_cache_populate.pl -root ./ref/cache ~{ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ~{output_basename}.cram
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    preemptible: preemptible_tries
    memory: "3 GiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_cram = "~{output_basename}.cram"
    File output_cram_index = "~{output_basename}.cram.crai"
    File output_cram_md5 = "~{output_basename}.cram.md5"
  }
}

# Convert CRAM file to BAM format
task ConvertToBam {
  input {
    File input_cram
    File ref_fasta
    File ref_fasta_index
    String output_basename
  }

  command <<<
    set -e
    set -o pipefail

    samtools view -b -o ~{output_basename}.bam -T ~{ref_fasta} ~{input_cram}

    samtools index ~{output_basename}.bam
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
    preemptible: 3
    memory: "3 GiB"
    cpu: "1"
    disks: "local-disk 200 HDD"
  }
  output {
    File output_bam = "~{output_basename}.bam"
    File output_bam_index = "~{output_basename}.bam.bai"
  }
}

# Calculates sum of a list of floats
task SumFloats {
  input {
    Array[Float] sizes
    Int preemptible_tries
  }

  command <<<
    python3 -c 'print(~{sep="+" sizes})'
  >>>
  output {
    Float total_size = read_float(stdout())
  }
  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: preemptible_tries
  }
}

# Print given message to stderr and return an error
task ErrorWithMessage{
  input {
    String message
  }
  command <<<
    >&2 echo "Error: ~{message}"
    exit 1
  >>>

  runtime {
    docker: "ubuntu:20.04"
  }
}

task CopyWorkflowOutputsByPath {
  input {
    String output_file_path
    String copy_bucket_path
    String workflow_name
    String cromwell_url
    String vault_token_path
    String google_account_vault_path

    String docker = "us.gcr.io/broad-gotc-prod/dsde-toolbox:stable_06-10-2021"
    Int memory_mb = 2000
    Int cpu = 1
    Int disk_size_gb = 20
  }
  meta {
    description: "For a test wrapper wdl, copy the results of the workflow being tested from the execution bucket to a given location"
  }

  parameter_meta {
    output_file_path: "File path to parse for cromwell ID"
    copy_bucket_path: "gs:// bucket path to copy workflow results to"
    workflow_name: "Name of the workflow for which the results are copied"
    cromwell_url: "Url for the cromwell server that ran the workflow"
  }

  command <<<
    set -e
    apk add python3
    pip3 install requests

    export PATH=$PATH:/usr/local/google-cloud-sdk/bin
    export VAULT_ADDR=https://clotho.broadinstitute.org:8200
    export VAULT_TOKEN=$(gsutil cat ~{vault_token_path})

    vault read -format=json ~{google_account_vault_path} | jq .data > picard-account.pem
    gcloud auth activate-service-account --key-file=picard-account.pem

    gcloud auth application-default print-access-token
    gcloud auth application-default print-access-token > token.txt

    python3 <<CODE
    import os, re, sys, requests, subprocess

    cromwell_url = "~{cromwell_url}"
    results_path = "~{copy_bucket_path}"

    def get_access_token() -> str:
      instance_url = "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token"
      request_headers = { "Metadata-Flavor": "Google" }

      sys.stdout.write(f"Requesting access token from metadata server... \n")

      # Request an access token from metadata server
      r = requests.get(instance_url, headers=request_headers)
      r.raise_for_status()

      # Extract access token from response
      access_token = r.json()["access_token"]

      sys.stdout.write(f"{r.json()} \n")
      sys.stdout.write(f"Successfully recieved access token from metadata server {access_token} \n")

      return access_token

    def parse_cromwell_id() -> str:
      file_path = "~{output_file_path}"
      pattern =  r"([a-z0-9]{8}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{4}-[a-z0-9]{12})"

      sys.stdout.write(f"Attempting to parse cromwell workflow ID from {file_path} \n")

      match = re.findall(pattern, file_path)

      # Get the second capture group which is the workflow we are actually testing
      if match is not None and match[1]:
        sys.stdout.write(f"Cromwell workflow ID found -> {match[1]} \n")
        return match[1]
      else:
        sys.stderr.write(f"ERROR: Unable to parse cromwell workflow ID from given file path -> {file_path} \n")
        sys.exit(1)

    def get_workflow_outputs(cromwell_id: str, access_token: str) -> list:
      outputs_url = f"{cromwell_url}/api/workflows/1/{cromwell_id}/outputs"
      request_headers = { "Authorization": f"Bearer {access_token}", "accept": "application/json" }

      sys.stdout.write(f"{request_headers} \n")
      sys.stdout.write(f"Querying outputs for cromwell workflow ID -> {cromwell_id} \n")
      sys.stdout.write(f"GET -> {outputs_url} \n")

      # Grab the outputs of the workflow from the cromwell api
      r = requests.get(outputs_url, headers=request_headers)
      sys.stdout.write(str(r.headers))

      sys.stdout.write(str(r.headers))

      outputs = r.json()["outputs"]

      outputs_list, task_outputs = [], [o for o in outputs.values()]

      # Flatten the outputs from each task to single list
      # Only want gs:// paths
      for x in task_outputs:
        if isinstance(x, str):
          outputs_list.append(x)
        if isinstance(x, list):
          outputs_list.extend(x)

      return outputs_list

    cromwell_id, token_file = parse_cromwell_id(), open("token.txt", "r")
    access_token = token_file.readline().rstrip("\n")

    workflow_outputs = get_workflow_outputs(cromwell_id, access_token)

    # Copy every file to results bucket
    for file in workflow_outputs:
      sys.stdout.write(f"...Copying {file} to {results_path}/{file} \n")
      subprocess.run(f"gsutil cp {file} {results_path}")

    # Write all the outputs to a file so we can export from the task
    file = open("output.txt", "w+")
    file.write("\n".join(workflow_outputs))

    sys.stdout.write(f"Copy to {results_path} completed successfully \n")

    CODE
  >>>

  runtime {
    docker: docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu : cpu
  }

  output{
   Array[String] workflow_outputs = read_lines("output.txt")
  }



}