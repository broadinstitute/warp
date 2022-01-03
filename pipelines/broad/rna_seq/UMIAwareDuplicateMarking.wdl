version 1.0

## Copyright Broad Institute, 2021
##
## This WDL pipeline implements UMI Aware Duplicate Marking
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

workflow UMIAwareDuplicateMarking {

  String pipeline_version = "0.1.0"

  input {
    File aligned_bam # aligned bam sorted by the query (read) name. 
    String output_basename
  }

  # First sort the aligned bam by coordinate, so we can group duplicate sets using UMIs in the next step.
  call SortSam as SortSamFirst {
    input:
      input_bam = aligned_bam,
      output_bam_basename = output_basename + ".STAR_aligned.coorinate_sorted",
      sort_order = "coordinate"
  }

  # Further divide each duplicate set (a set of reads with the same insert start and end coordinates)
  # into subsets that share the same UMIs i.e. differenciate PCR duplicates from biological duplicates.
  # (biological duplicates are independent DNA molecules that are sheared such that the inserts are indistinguishable.)
  # input: a coordinate sorted bam
  # output: a coordinate sorted bam with UMIs (what are the generated tags?) .
  
  call GroupByUMIs {
    input:
      bam = SortSamFirst.output_bam,
      bam_index = select_first([SortSamFirst.output_bam_index, "bam_index_not_found"]),
      output_bam_basename = output_basename + ".grouped_by_UMI"
  }

  # input:
  # output: 
  call SortSam as SortSamQueryName {
    input:
      input_bam = GroupByUMIs.grouped_bam,
      output_bam_basename = output_basename + ".grouped.queryname_sorted",
      sort_order = "queryname"
  }

  call MarkDuplicates {
    input:
      bam = SortSamQueryName.output_bam,
      output_basename = output_basename
  }

  call SortSam as SortSamSecond {
    input:
      input_bam = MarkDuplicates.duplicate_marked_bam,
      output_bam_basename = output_basename + ".duplicate_marked.coordinate_sorted",
      sort_order = "coordinate"
  }

  output {
    File duplicate_marked_bam = SortSamSecond.output_bam
    File duplicate_marked_bam_index = select_first([SortSamSecond.output_bam_index, "bam_index_not_found"])
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
  }
}

task MarkDuplicates {
  input {
    File bam
    String output_basename
  }

  String output_bam_basename = output_basename + ".duplicate_marked"

  Int disk_size = ceil(3 * size(bam, "GB")) + 128
  command <<<
    gatk MarkDuplicates -I ~{bam} --READ_ONE_BARCODE_TAG BX -O ~{output_bam_basename}.bam --METRICS_FILE ~{output_basename}.duplicate.metrics --ASSUME_SORT_ORDER queryname
  >>>

  output {
    File duplicate_marked_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{output_basename}.duplicate.metrics"
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "16 GB"
  }
}

task SortSam {
  input {
    File input_bam
    String output_bam_basename
    String sort_order # "queryname" or "coordinate"
  }

  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Float sort_sam_disk_multiplier = 4.0
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 256

  String extra_args = if sort_order == "coordinate" then "CREATE_INDEX=true" else ""


  command {
    java -Xms8192m -jar /usr/picard/picard.jar \
    SortSam \
    INPUT=~{input_bam} \
    OUTPUT=~{output_bam_basename}.bam \
    SORT_ORDER=~{sort_order} \
    CREATE_MD5_FILE=true \
    MAX_RECORDS_IN_RAM=300000 \
    ~{extra_args}

  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.8"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "16 GB"
    preemptible: 0
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File? output_bam_index = "~{output_bam_basename}.bai" # Create index for coordinate sort
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

task GroupByUMIs {
  input {
    File bam
    File bam_index
    String output_bam_basename
  }

  Int disk_space = ceil(2.2 * size(bam, "GB")) + 300
  command <<<
    umi_tools group -I ~{bam} --paired --no-sort-output --output-bam --stdout ~{output_bam_basename}.bam --umi-tag-delimiter "-" \
    --extract-umi-method tag --umi-tag RX --unmapped-reads use
  >>>

  output {
    File grouped_bam = "~{output_bam_basename}.bam"
  }

  runtime {
    # Note to Takuto that this updates umi_tools from 1.0.0 to 1.1.1 (no diffs found in test verification).
    docker : "us.gcr.io/broad-gotc-prod/umi_tools:1.0.0-1.1.1-1638821470"
    disks : "local-disk " + disk_space + " HDD"
    preemptible: 0
    cpu: "8"
    memory: "64 GB" # Sato: is this too much?
  }
}