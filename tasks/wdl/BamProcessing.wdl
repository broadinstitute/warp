version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for BAM file processing of human whole-genome or exome sequencing data.
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

# Sort BAM file by coordinate order
task SortSam {
  input {
    File input_bam
    String output_bam_basename
    Int preemptible_tries
    Int compression_level
    Int additional_disk = 20
    Int memory_multiplier = 1
    #Setting default docker value for workflows that haven't yet been azurized. 
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Float sort_sam_disk_multiplier = 3.25
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + additional_disk

  Int machine_mem_mb = ceil(5000 * memory_multiplier)
  Int java_max_memory_mb = machine_mem_mb - 500
  Int java_inital_memory_mb = machine_mem_mb - 1000

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_inital_memory_mb}m -Xmx~{java_max_memory_mb}m -jar /usr/picard/picard.jar \
      SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000

  }
  runtime {
    docker: docker
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "${machine_mem_mb} MiB"
    preemptible: preemptible_tries
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}


# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  input {
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
    Float total_input_size
    Int compression_level
    Int preemptible_tries

    # The program default for READ_NAME_REGEX is appropriate in nearly every case.
    # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
    # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
    String? read_name_regex
    Int memory_multiplier = 1
    Int additional_disk = 20

    Float? sorting_collection_size_ratio
  }

  # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs and the merged output.
  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
  Float md_disk_multiplier = 3
  Int disk_size = ceil(md_disk_multiplier * total_input_size) + additional_disk

  Float memory_size = 7.5 * memory_multiplier
  Int java_memory_size = (ceil(memory_size) - 2)

  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size}g -jar /usr/picard/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ~{"READ_NAME_REGEX=" + read_name_regex} \
      ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    Array[String] sequence_group_interval
    File dbsnp_vcf
    File dbsnp_vcf_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Float dbsnp_size = size(dbsnp_vcf, "GiB")
  Int disk_size = ceil((size(input_bam, "GiB") / bqsr_scatter) + ref_size + dbsnp_size) + 20

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
    gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms5000m -Xmx5500m" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recalibration_report_filename} \
      --known-sites ~{dbsnp_vcf} \
      --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs} \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "6000 MiB"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibration_report = "~{recalibration_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    Array[String] sequence_group_interval
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int compression_level
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    Int memory_multiplier = 1
    Int additional_disk = 20
    Boolean bin_base_qualities = true
    Boolean somatic = false
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil((size(input_bam, "GiB") * 3 / bqsr_scatter) + ref_size) + additional_disk

  Int memory_size = ceil(3500 * memory_multiplier)
  Int java_memory_mb = memory_size - 500

  Boolean bin_somatic_base_qualities = bin_base_qualities && somatic

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
    gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=~{compression_level} -Xms3000m -Xmx~{java_memory_mb}m" \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{output_bam_basename}.bam \
      -bqsr ~{recalibration_report} \
      ~{true='--static-quantized-quals 10' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 20' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 30' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 40' false='' bin_somatic_base_qualities} \
      ~{true='--static-quantized-quals 50' false='' bin_somatic_base_qualities} \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "~{memory_size} MiB"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
    File recalibrated_bam_checksum = "~{output_bam_basename}.bam.md5"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  input {
    Array[File] input_bqsr_reports
    String output_report_filename
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0"
  }

  command {
    gatk --java-options "-Xms3000m -Xmx3000m" \
      GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O ~{output_report_filename}
    }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "3500 MiB"
    bootDiskSizeGb: 15
    disks: "local-disk 20 HDD"
  }
  output {
    File output_bqsr_report = "~{output_report_filename}"
  }
}

# Combine multiple *sorted* BAM files
task GatherSortedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Float total_input_size
    Int compression_level
    Int preemptible_tries
    Int additional_disk = 20
    Int memory_multiplier = 1
  }

  # Multiply the input bam size by two to account for the input and output
  Int disk_size = ceil(2 * total_input_size) + additional_disk
  Int machine_mem_mb = ceil(3000 * memory_multiplier)
  Int java_max_memory_mb = machine_mem_mb - 500
  Int java_inital_memory_mb = machine_mem_mb - 1000

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_inital_memory_mb}m -Xmx~{java_max_memory_mb}m -jar /usr/picard/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "${machine_mem_mb} MiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

# Combine multiple *unsorted* BAM files
# Note that if/when WDL supports optional outputs, we should merge this task with the sorted version
task GatherUnsortedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Float total_input_size
    Int compression_level
    Int preemptible_tries
  }

  # Multiply the input bam size by two to account for the input and output
  Int disk_size = ceil(2 * total_input_size) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

task GenerateSubsettedContaminationResources {
  input {
    String bait_set_name
    File target_interval_list
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    Int preemptible_tries
  }

  String output_ud = bait_set_name + "." + basename(contamination_sites_ud)
  String output_bed = bait_set_name + "." + basename(contamination_sites_bed)
  String output_mu = bait_set_name + "." + basename(contamination_sites_mu)
  String target_overlap_counts = "target_overlap_counts.txt"

  command <<<
    set -e -o pipefail

    grep -vE "^@" ~{target_interval_list} |
       awk -v OFS='\t' '$2=$2-1' |
       /app/bedtools intersect -c -a ~{contamination_sites_bed} -b - |
       cut -f6 > ~{target_overlap_counts}

    function restrict_to_overlaps() {
        # print lines from whole-genome file from loci with non-zero overlap
        # with target intervals
        WGS_FILE=$1
        EXOME_FILE=$2
        paste ~{target_overlap_counts} $WGS_FILE |
            grep -Ev "^0" |
            cut -f 2- > $EXOME_FILE
        echo "Generated $EXOME_FILE"
    }

    restrict_to_overlaps ~{contamination_sites_ud} ~{output_ud}
    restrict_to_overlaps ~{contamination_sites_bed} ~{output_bed}
    restrict_to_overlaps ~{contamination_sites_mu} ~{output_mu}

  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "3.5 GiB"
    disks: "local-disk 10 HDD"
    docker: "us.gcr.io/broad-gotc-prod/bedtools:2.27.1"
  }
  output {
    File subsetted_contamination_ud = output_ud
    File subsetted_contamination_bed = output_bed
    File subsetted_contamination_mu = output_mu
  }
}

# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling
task CheckContamination {
  input {
    File input_bam
    File input_bam_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File ref_fasta
    File ref_fasta_index
    String output_prefix
    Int preemptible_tries
    Float contamination_underestimation_factor
    Boolean disable_sanity_check = false

     String docker = "us.gcr.io/broad-gotc-prod/verify-bam-id:1.0.1-c1cba76e979904eb69c31520a0d7f5be63c72253-1639071840"
  }

  Int disk_size = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB")) + 30

  command <<<
    set -e

    # creates a ~{output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ~{output_prefix} \
    --BamFile ~{input_bam} \
    --Reference ~{ref_fasta} \
    --UDPath ~{contamination_sites_ud} \
    --MeanPath ~{contamination_sites_mu} \
    --BedPath ~{contamination_sites_bed} \
    ~{true="--DisableSanityCheck" false="" disable_sanity_check} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('~{output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"])/~{contamination_underestimation_factor})
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "7.5 GiB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
    cpu: 2
  }
  output {
    File selfSM = "~{output_prefix}.selfSM"
    Float contamination = read_float(stdout())
  }
}
