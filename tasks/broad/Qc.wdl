version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for QC of human whole-genome or exome sequencing data.
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

# Collect sequencing yield quality metrics
task CollectQualityYieldMetrics {
  input {
    File input_bam
    String metrics_filename

    Int preemptible_tries = 3
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 20

  command {
    java -Xms2000m -Xmx3000m -jar /usr/picard/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=~{input_bam} \
      OQ=true \
      OUTPUT=~{metrics_filename}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3500 MiB"
    preemptible: preemptible_tries
  }
  output {
    File quality_yield_metrics = "~{metrics_filename}"
  }
}

# Collect base quality and insert size metrics
task CollectUnsortedReadgroupBamQualityMetrics {
  input {
    File input_bam
    String output_bam_prefix
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 20

  command {
    java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=QualityScoreDistribution \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=ALL_READS

    touch ~{output_bam_prefix}.insert_size_metrics
    touch ~{output_bam_prefix}.insert_size_histogram.pdf
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    memory: "7000 MiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File base_distribution_by_cycle_pdf = "~{output_bam_prefix}.base_distribution_by_cycle.pdf"
    File base_distribution_by_cycle_metrics = "~{output_bam_prefix}.base_distribution_by_cycle_metrics"
    File insert_size_histogram_pdf = "~{output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "~{output_bam_prefix}.insert_size_metrics"
    File quality_by_cycle_pdf = "~{output_bam_prefix}.quality_by_cycle.pdf"
    File quality_by_cycle_metrics = "~{output_bam_prefix}.quality_by_cycle_metrics"
    File quality_distribution_pdf = "~{output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "~{output_bam_prefix}.quality_distribution_metrics"
  }
}

# Collect alignment summary and GC bias quality metrics
task CollectReadgroupBamQualityMetrics {
  input {
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Boolean collect_gc_bias_metrics = true
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  command {
    # These are optionally generated, but need to exist for Cromwell's sake
    touch ~{output_bam_prefix}.gc_bias.detail_metrics \
      ~{output_bam_prefix}.gc_bias.pdf \
      ~{output_bam_prefix}.gc_bias.summary_metrics

    java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      ~{true='PROGRAM="CollectGcBiasMetrics"' false="" collect_gc_bias_metrics} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=READ_GROUP
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    memory: "7000 MiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File alignment_summary_metrics = "~{output_bam_prefix}.alignment_summary_metrics"
    File gc_bias_detail_metrics = "~{output_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "~{output_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "~{output_bam_prefix}.gc_bias.summary_metrics"
  }
}

# Collect quality metrics from the aggregated bam
task CollectAggregationMetrics {
  input {
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Boolean collect_gc_bias_metrics = true
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  command {
    # These are optionally generated, but need to exist for Cromwell's sake
    touch ~{output_bam_prefix}.gc_bias.detail_metrics \
      ~{output_bam_prefix}.gc_bias.pdf \
      ~{output_bam_prefix}.gc_bias.summary_metrics \
      ~{output_bam_prefix}.insert_size_metrics \
      ~{output_bam_prefix}.insert_size_histogram.pdf

    java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=QualityScoreDistribution \
      ~{true='PROGRAM="CollectGcBiasMetrics"' false="" collect_gc_bias_metrics} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    memory: "7000 MiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File alignment_summary_metrics = "~{output_bam_prefix}.alignment_summary_metrics"
    File bait_bias_detail_metrics = "~{output_bam_prefix}.bait_bias_detail_metrics"
    File bait_bias_summary_metrics = "~{output_bam_prefix}.bait_bias_summary_metrics"
    File gc_bias_detail_metrics = "~{output_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "~{output_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "~{output_bam_prefix}.gc_bias.summary_metrics"
    File insert_size_histogram_pdf = "~{output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "~{output_bam_prefix}.insert_size_metrics"
    File pre_adapter_detail_metrics = "~{output_bam_prefix}.pre_adapter_detail_metrics"
    File pre_adapter_summary_metrics = "~{output_bam_prefix}.pre_adapter_summary_metrics"
    File quality_distribution_pdf = "~{output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "~{output_bam_prefix}.quality_distribution_metrics"
    File error_summary_metrics = "~{output_bam_prefix}.error_summary_metrics"
  }
}

task ConvertSequencingArtifactToOxoG {
  input {
    File pre_adapter_detail_metrics
    File bait_bias_detail_metrics
    String base_name
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int preemptible_tries
    Int memory_multiplier = 1
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(pre_adapter_detail_metrics, "GiB") + size(bait_bias_detail_metrics, "GiB") + ref_size) + 20

  Int memory_size = ceil(4000 * memory_multiplier)
  Int java_memory_size = memory_size - 1000
  Int max_heap = memory_size - 500

  command {
    input_base=$(dirname ~{pre_adapter_detail_metrics})/~{base_name}
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m \
      -jar /usr/picard/picard.jar \
      ConvertSequencingArtifactToOxoG \
      --INPUT_BASE $input_base \
      --OUTPUT_BASE ~{base_name} \
      --REFERENCE_SEQUENCE ~{ref_fasta}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    memory: "~{memory_size} MiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File oxog_metrics = "~{base_name}.oxog_metrics"
  }
}

# Check that the fingerprints of separate readgroups all match
task CrossCheckFingerprints {
  input {
    Array[File] input_bams
    Array[File] input_bam_indexes
    File haplotype_database_file
    String metrics_filename
    Float total_input_size
    Int preemptible_tries
    Float lod_threshold
    String cross_check_by
  }

  Int disk_size = ceil(total_input_size) + 20

  command <<<
    java -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Xmx3000m \
      -jar /usr/picard/picard.jar \
      CrosscheckFingerprints \
      OUTPUT=~{metrics_filename} \
      HAPLOTYPE_MAP=~{haplotype_database_file} \
      EXPECT_ALL_GROUPS_TO_MATCH=true \
      INPUT=~{sep=' INPUT=' input_bams} \
      LOD_THRESHOLD=~{lod_threshold} \
      CROSSCHECK_BY=~{cross_check_by}
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "3500 MiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File cross_check_fingerprints_metrics = "~{metrics_filename}"
  }
}

task CheckFingerprintTask {
  input {
    File? input_bam
    File? input_bam_index
    File? input_vcf
    File? input_vcf_index
    String? input_sample_alias

    File genotypes
    File? genotypes_index
    String expected_sample_alias

    String output_basename
    Float genotype_lod_threshold = 5.0

    File haplotype_database_file
    File? ref_fasta
    File? ref_fasta_index

    Int memory_size = 2500
    Int preemptible_tries = 3

    Boolean allow_lod_zero = false
  }

  Int java_memory_size = memory_size - 1000
  Int max_heap = memory_size - 500

  Int disk_size = ceil(size(input_bam, "GiB") + size(input_vcf, "GiB")) + 20
  # Picard has different behavior depending on whether or not the OUTPUT parameter ends with a '.', so we are explicitly
  #   passing in where we want the two metrics files to go to avoid any potential confusion.
  String summary_metrics_location = "~{output_basename}.fingerprinting_summary_metrics"
  String detail_metrics_location = "~{output_basename}.fingerprinting_detail_metrics"

  File input_file = select_first([input_vcf, input_bam])

  command <<<
    set -e
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
    CheckFingerprint \
      --INPUT ~{input_file} \
      ~{if defined(input_vcf) then "--OBSERVED_SAMPLE_ALIAS \"" + input_sample_alias + "\"" else ""} \
      --GENOTYPES ~{genotypes} \
      --EXPECTED_SAMPLE_ALIAS "~{expected_sample_alias}" \
      ~{if defined(input_bam) then "--IGNORE_READ_GROUPS true" else ""} \
      --HAPLOTYPE_MAP ~{haplotype_database_file} \
      --GENOTYPE_LOD_THRESHOLD ~{genotype_lod_threshold} \
      --SUMMARY_OUTPUT ~{summary_metrics_location} \
      --DETAIL_OUTPUT ~{detail_metrics_location} \
      ~{"--REFERENCE_SEQUENCE " + ref_fasta} \
      ~{true='--EXIT_CODE_WHEN_NO_VALID_CHECKS 0' false='' allow_lod_zero}

    CONTENT_LINE=$(cat ~{summary_metrics_location} |
    grep -n "## METRICS CLASS\tpicard.analysis.FingerprintingSummaryMetrics" |
    cut -f1 -d:)
    CONTENT_LINE=$(($CONTENT_LINE+2))
    sed '8q;d' ~{summary_metrics_location} | cut -f5 > lod
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory_size} MiB"
    preemptible: preemptible_tries
  }

  output {
    File summary_metrics = summary_metrics_location
    File detail_metrics = detail_metrics_location
    Float lod = read_float("lod")
  }
}

task CheckPreValidation {
  input {
    File duplication_metrics
    File chimerism_metrics
    Float max_duplication_in_reasonable_sample
    Float max_chimerism_in_reasonable_sample
    
    Int preemptible_tries = 3
  }

  command <<<
    set -o pipefail
    set -e

    grep -A 1 PERCENT_DUPLICATION ~{duplication_metrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ~{chimerism_metrics} | grep -v OF_PAIR > chimerism.csv

    python3 <<CODE

    import csv
    with open('duplication.csv') as dupfile:
      reader = csv.DictReader(dupfile, delimiter='\t')
      for row in reader:
        with open("duplication_value.txt","w") as file:
          file.write(row['PERCENT_DUPLICATION'])
          file.close()

    with open('chimerism.csv') as chimfile:
      reader = csv.DictReader(chimfile, delimiter='\t')
      for row in reader:
        with open("chimerism_value.txt","w") as file:
          file.write(row['PCT_CHIMERAS'])
          file.close()

    CODE

>>>
  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    preemptible: preemptible_tries
    memory: "2 GiB"
  }
  output {
    Float duplication_rate = read_float("duplication_value.txt")
    Float chimerism_rate = read_float("chimerism_value.txt")
    Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
  }
}

task ValidateSamFile {
  input {
    File input_bam
    File? input_bam_index
    String report_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int? max_output
    Array[String]? ignore
    Boolean? is_outlier_data
    Int preemptible_tries = 0
    Int memory_multiplier = 1
    Int additional_disk = 20

    Int disk_size = ceil(size(input_bam, "GiB") 
                    + size(ref_fasta, "GiB") 
                    + size(ref_fasta_index, "GiB")
                    + size(ref_dict, "GiB")) + additional_disk
  }

  Int memory_size = ceil(16000 * memory_multiplier)
  Int java_memory_size = memory_size - 1000
  Int max_heap = memory_size - 500

  command {
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      INPUT=~{input_bam} \
      OUTPUT=~{report_filename} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      ~{"MAX_OUTPUT=" + max_output} \
      IGNORE=~{default="null" sep=" IGNORE=" ignore} \
      MODE=VERBOSE \
      ~{default='SKIP_MATE_VALIDATION=false' true='SKIP_MATE_VALIDATION=true' false='SKIP_MATE_VALIDATION=false' is_outlier_data} \
      IS_BISULFITE_SEQUENCED=false
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "~{memory_size} MiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File report = "~{report_filename}"
  }
}

task CollectWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    String metrics_filename
    File wgs_coverage_interval_list
    File ref_fasta
    File ref_fasta_index
    Int read_length = 250
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  command {
    java -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "3000 MiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File metrics = "~{metrics_filename}"
  }
}

# Collect raw WGS metrics (commonly used QC thresholds)
task CollectRawWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    String metrics_filename
    File wgs_coverage_interval_list
    File ref_fasta
    File ref_fasta_index
    Int read_length = 250
    Int preemptible_tries
    Int memory_multiplier = 1
    Int additional_disk = 20
  }
  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + additional_disk

  Int memory_size = ceil((if (disk_size < 110) then 5 else 7) * memory_multiplier)
  String java_memory_size = (memory_size - 1) * 1000

  command {
    java -Xms~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectRawWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File metrics = "~{metrics_filename}"
  }
}

task CollectHsMetrics {
  input {
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    String metrics_filename
    File target_interval_list
    File bait_interval_list
    Int preemptible_tries
    Int memory_multiplier = 1
    Int additional_disk = 20
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + additional_disk
  # Try to fit the input bam into memory, within reason.
  Int rounded_bam_size = ceil(size(input_bam, "GiB") + 0.5)
  Int rounded_memory_size = ceil((if (rounded_bam_size > 10) then 10 else rounded_bam_size) * memory_multiplier)
  Int memory_size = if rounded_memory_size < 7 then 7000 else (rounded_memory_size * 1000)
  Int java_memory_size = memory_size - 1000
  Int max_heap = memory_size - 500

  # There are probably more metrics we want to generate with this tool
  command {
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
      CollectHsMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      VALIDATION_STRINGENCY=SILENT \
      TARGET_INTERVALS=~{target_interval_list} \
      BAIT_INTERVALS=~{bait_interval_list} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY \
      OUTPUT=~{metrics_filename}
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "~{memory_size} MiB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File metrics = metrics_filename
  }
}

# Generate a checksum per readgroup
task CalculateReadGroupChecksum {
  input {
    File input_bam
    File input_bam_index
    String read_group_md5_filename
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 80

  command {
    java -Xms1000m -Xmx3500m -jar /usr/picard/picard.jar \
      CalculateReadGroupChecksum \
      INPUT=~{input_bam} \
      OUTPUT=~{read_group_md5_filename}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "6000 MiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File md5_file = "~{read_group_md5_filename}"
  }
}

# Validate a (g)VCF with -gvcf specific validation
task ValidateVCF {
  input {
    File input_vcf
    File input_vcf_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? dbsnp_vcf
    File? dbsnp_vcf_index
    File calling_interval_list
    File? calling_interval_list_index  # if the interval list is a VCF, than an index file makes VcfToIntervalList run faster
    Boolean calling_intervals_defined = true
    Int preemptible_tries = 3
    Boolean is_gvcf = true
    String? extra_args
    #Setting default docker value for workflows that haven't yet been azurized.
    String docker_path = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int machine_mem_mb = 7000
  }

  String calling_interval_list_basename = basename(calling_interval_list)
  String calling_interval_list_index_basename = if calling_intervals_defined then "" else basename(select_first([calling_interval_list_index]))

  Int command_mem_mb = machine_mem_mb - 2000
  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_vcf, "GiB") + size(calling_interval_list, "GiB") + size(dbsnp_vcf, "GiB") + ref_size) + 20

  command {
    set -e
    
    if [ ~{calling_intervals_defined} == "false" ]; then
      # We can't always assume the index was located with the vcf, so make a link so that the paths look the same
      ln -s ~{calling_interval_list} ~{calling_interval_list_basename}
      ln -s ~{calling_interval_list_index} ~{calling_interval_list_index_basename}
      gatk --java-options "-Xms~{command_mem_mb}m -Xmx~{command_mem_mb}m" \
        VcfToIntervalList -I ~{calling_interval_list_basename} -O intervals_from_gvcf.interval_list --VARIANT_ID_METHOD USE_FIRST
      INTERVALS="intervals_from_gvcf.interval_list"
    else
      INTERVALS="~{calling_interval_list}"
    fi
    
    gatk --java-options "-Xms~{command_mem_mb}m -Xmx~{command_mem_mb}m" \
      ValidateVariants \
      -V ~{input_vcf} \
      -R ~{ref_fasta} \
      -L $INTERVALS \
      ~{true="-gvcf" false="" is_gvcf} \
      --validation-type-to-exclude ALLELES \
      ~{"--dbsnp " + dbsnp_vcf} \
      ~{extra_args}
  }
  runtime {
    docker: docker_path
    preemptible: preemptible_tries
    memory: machine_mem_mb + " MiB"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }
}

# Collect variant calling metrics from GVCF output
task CollectVariantCallingMetrics {
  input {
    File input_vcf
    File input_vcf_index
    String metrics_basename
    File dbsnp_vcf
    File dbsnp_vcf_index
    File ref_dict
    File evaluation_interval_list
    Boolean is_gvcf = true
    Int preemptible_tries
    #Setting default docker value for workflows that haven't yet been azurized.
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
  }

  Int disk_size = ceil(size(input_vcf, "GiB") + size(dbsnp_vcf, "GiB")) + 20

  command {
    java -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
      CollectVariantCallingMetrics \
      INPUT=~{input_vcf} \
      OUTPUT=~{metrics_basename} \
      DBSNP=~{dbsnp_vcf} \
      SEQUENCE_DICTIONARY=~{ref_dict} \
      TARGET_INTERVALS=~{evaluation_interval_list} \
      ~{true="GVCF_INPUT=true" false="" is_gvcf}
  }
  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "3000 MiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File summary_metrics = "~{metrics_basename}.variant_calling_summary_metrics"
    File detail_metrics = "~{metrics_basename}.variant_calling_detail_metrics"
  }
}
