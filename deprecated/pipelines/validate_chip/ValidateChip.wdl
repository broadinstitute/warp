version 1.0
# The ValidateChip pipeline is deprecated 2025-03-06
import "../../../../tasks/wdl/IlluminaGenotypingArrayTasks.wdl" as GenotypingTasks
import "../../../../tasks/wdl/InternalArraysTasks.wdl" as InternalTasks

## Copyright Broad Institute, 2019
##
## This WDL pipeline implements validation of a single Illumina array chip type
##
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

workflow ValidateChip {

  String pipeline_version = "1.16.8"

  input {
    String sample_alias
    Int analysis_version_number
    Float call_rate_threshold
    String reported_gender

    String chip_well_barcode
    File red_idat_cloud_path
    File green_idat_cloud_path
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbSNP_vcf
    File dbSNP_vcf_index

    File bead_pool_manifest_file
    String chip_type = basename(bead_pool_manifest_file, ".bpm")

    # for createExtendedIlluminaManifest:
    File chip_manifest_csv_file
    File supported_ref_fasta
    File supported_ref_fasta_index
    File supported_ref_dict
    File chain_file

    File cluster_file

    # For GenotypeConcordance Check:
    File control_sample_vcf_file
    File control_sample_vcf_index_file
    File control_sample_intervals_file
    String control_sample_name
    Float indel_genotype_concordance_threshold

    Int disk_size
    Int preemptible_tries
  }

  call GenotypingTasks.CreateExtendedIlluminaManifest {
    input:
      input_csv = chip_manifest_csv_file,
      output_base_name = chip_type + ".2.0",
      cluster_file = cluster_file,
      dbSNP_vcf_file = dbSNP_vcf,
      dbSNP_vcf_index_file = dbSNP_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      supported_ref_fasta = supported_ref_fasta,
      supported_ref_fasta_index = supported_ref_fasta_index,
      supported_ref_dict = supported_ref_dict,
      chain_file = chain_file,
      disk_size = disk_size,
      preemptible_tries = 0
  }


  call GenotypingTasks.AutoCall {
    input:
      chip_well_barcode = chip_well_barcode,
      green_idat_cloud_path = green_idat_cloud_path,
      red_idat_cloud_path = red_idat_cloud_path,
      bead_pool_manifest_file = bead_pool_manifest_file,
      cluster_file = cluster_file,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call GenotypingTasks.GtcToVcf {
    input:
      vcf_filename = chip_well_barcode + ".vcf.gz",
      input_gtc = AutoCall.gtc_file,
      extended_chip_manifest_file = CreateExtendedIlluminaManifest.output_csv,
      cluster_file = cluster_file,
      bead_pool_manifest_file = bead_pool_manifest_file,
      sample_alias = sample_alias,
      analysis_version_number = analysis_version_number,
      reported_gender = reported_gender,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries,
      pipeline_version = "ValidateChip_v" + pipeline_version
  }

  call GenotypingTasks.SelectIndels {
    input:
      input_vcf_file = GtcToVcf.output_vcf,
      input_vcf_index_file = GtcToVcf.output_vcf_index,
      output_vcf_filename = chip_well_barcode + ".indels.vcf",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call GenotypingTasks.CollectArraysVariantCallingMetrics {
    input:
      input_vcf_file = GtcToVcf.output_vcf,
      input_vcf_index_file = GtcToVcf.output_vcf_index,
      dbSNP_vcf_file = dbSNP_vcf,
      dbSNP_vcf_index_file = dbSNP_vcf_index,
      call_rate_threshold = call_rate_threshold,
      output_metrics_basename = chip_well_barcode,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call GenotypingTasks.ValidateVariants {
    input:
      input_vcf_file = GtcToVcf.output_vcf,
      input_vcf_index_file = GtcToVcf.output_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call GenotypingTasks.VcfToIntervalList {
    input:
      vcf_file = GtcToVcf.output_vcf,
      interval_list_filename = chip_well_barcode + ".interval_list",
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call GenotypingTasks.VcfToIntervalList as IndelVcfToIntervalList {
    input:
      vcf_file = SelectIndels.output_vcf,
      interval_list_filename = chip_well_barcode + ".indels.interval_list",
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call GenotypingTasks.SelectVariants as SelectVariantsForGenotypeConcordance {
    input:
      input_vcf_file = GtcToVcf.output_vcf,
      input_vcf_index_file = GtcToVcf.output_vcf_index,
      output_vcf_filename = chip_well_barcode + ".select_variants.vcf",
      excludeFiltered = true,
      excludeNonVariants = true,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      preemptible_tries = preemptible_tries
  }

  call GenotypingTasks.SelectVariants as SelectIndelVariantsForGenotypeConcordance {
    input:
      input_vcf_file = SelectIndels.output_vcf,
      input_vcf_index_file = SelectIndels.output_vcf_index,
      output_vcf_filename = chip_well_barcode + ".indels.select_variants.vcf",
      excludeFiltered = true,
      excludeNonVariants = true,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      preemptible_tries = preemptible_tries
  }

  call GenotypeConcordance {
    input:
      call_vcf_file = SelectVariantsForGenotypeConcordance.output_vcf,
      call_vcf_index_file = SelectVariantsForGenotypeConcordance.output_vcf_index,
      call_intervals_file = VcfToIntervalList.interval_list_file,
      call_sample_name = chip_well_barcode,
      truth_vcf_file = control_sample_vcf_file,
      truth_vcf_index_file = control_sample_vcf_index_file,
      truth_intervals_file = control_sample_intervals_file,
      truth_sample_name = control_sample_name,
      variant_type = "SNP",
      output_metrics_basename = chip_well_barcode + ".gc",
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  call GenotypeConcordance as IndelGenotypeConcordance {
    input:
      call_vcf_file = SelectIndelVariantsForGenotypeConcordance.output_vcf,
      call_vcf_index_file = SelectIndelVariantsForGenotypeConcordance.output_vcf_index,
      call_intervals_file = IndelVcfToIntervalList.interval_list_file,
      call_sample_name = chip_well_barcode,
      truth_vcf_file = control_sample_vcf_file,
      truth_vcf_index_file = control_sample_vcf_index_file,
      truth_intervals_file = control_sample_intervals_file,
      truth_sample_name = control_sample_name,
      variant_type = "INDEL",
      output_metrics_basename = chip_well_barcode + ".indels.gc",
      concordance_threshold = indel_genotype_concordance_threshold,
      disk_size = disk_size,
      preemptible_tries = preemptible_tries
  }

  output {
    File gtc_file = AutoCall.gtc_file
    File output_vcf = GtcToVcf.output_vcf
    File output_vcf_index = GtcToVcf.output_vcf_index
    File arrays_variant_calling_detail_metrics_file = CollectArraysVariantCallingMetrics.detail_metrics
    File arrays_variant_calling_summary_metrics_file = CollectArraysVariantCallingMetrics.summary_metrics
    File arrays_variant_calling_control_metrics_file = CollectArraysVariantCallingMetrics.control_metrics
    Boolean passes_autocall = CollectArraysVariantCallingMetrics.passes_autocall
    File create_extended_illumina_manifest_extended_csv_file = CreateExtendedIlluminaManifest.output_csv
    File create_extended_illumina_manifest_bad_assays_file = CreateExtendedIlluminaManifest.output_bad_assays_file
    File create_extended_illumina_manifest_report_file = CreateExtendedIlluminaManifest.report_file
    File genotype_concordance_summary_metrics_file = GenotypeConcordance.summary_metrics
    File genotype_concordance_detail_metrics_file  = GenotypeConcordance.detail_metrics
    File genotype_concordance_contingency_metrics_file = GenotypeConcordance.contingency_metrics
    File genotype_concordance_vcf = GenotypeConcordance.output_vcf
    File genotype_concordance_txt_file = GenotypeConcordance.output_txt
    File indel_genotype_concordance_summary_metrics_file = IndelGenotypeConcordance.summary_metrics
    File indel_genotype_concordance_detail_metrics_file  = IndelGenotypeConcordance.detail_metrics
    File indel_genotype_concordance_contingency_metrics_file = IndelGenotypeConcordance.contingency_metrics
    File indel_genotype_concordance_vcf = IndelGenotypeConcordance.output_vcf
    File indel_genotype_concordance_txt_file = IndelGenotypeConcordance.output_txt
    File output_bead_pool_manifest_file = bead_pool_manifest_file
  }
  meta {
    allowNestedInputs: true
  }
}


task GenotypeConcordance {
  input {
    File call_vcf_file
    File call_vcf_index_file
    File call_intervals_file
    String call_sample_name
    File truth_vcf_file
    File truth_vcf_index_file
    File truth_intervals_file
    String truth_sample_name
    String output_metrics_basename
    String variant_type
    Float concordance_threshold=0.99

    Int disk_size
    Int preemptible_tries
  }

  String report_name = "~{output_metrics_basename}.txt"

  command <<<
    set -e

    /gatk/gatk --java-options "-Xms6000m -Xmx6500m" \
      GenotypeConcordance \
      --CALL_VCF ~{call_vcf_file} \
      --CALL_SAMPLE ~{call_sample_name} \
      --TRUTH_VCF ~{truth_vcf_file} \
      --TRUTH_SAMPLE ~{truth_sample_name} \
      --INTERVALS ~{call_intervals_file} \
      --INTERVALS ~{truth_intervals_file} \
      --MISSING_HOM true \
      --OUTPUT ~{output_metrics_basename} \
      --OUTPUT_VCF true

      echo -e "Notes" > ~{report_name}

      for i in 4 5 7 8 13
      do
        NAME=`grep '^VARIANT' ~{output_metrics_basename}.genotype_concordance_summary_metrics | cut -f $i`
        VALUE=`grep '^~{variant_type}'   ~{output_metrics_basename}.genotype_concordance_summary_metrics | cut -f $i`
        echo -e "$NAME\t$VALUE" >> ~{report_name}
      done

      for i in FP FN
      do
        VALUE=`zcat ~{output_metrics_basename}.genotype_concordance.vcf.gz | grep $i | wc -l`
        echo -e "$i\t$VALUE" >> ~{report_name}
      done

      # Strip out comments and blank lines and get the genotype concordance
      # Find the column number of GENOTYPE_CONCORDANCE and retrieve the value in the second line of the column
      # CONCORDANCE set to empty if file has more than 3 lines (should only have column headers and two data lines)
      # Samples that have few indel cause genotype concordance to be reported as ‘?’ in the metrics file
      CONCORDANCE=$(grep -e ~{variant_type} -e VARIANT  ~{output_metrics_basename}.genotype_concordance_summary_metrics |
        sed '/#/d' |
        sed '/^\s*$/d' |
        awk -v col=GENOTYPE_CONCORDANCE -F '\t' \
        'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}}} NR==2{print $c} NR>3{exit 1}')

      if [[ "$CONCORDANCE" =~ ^[+-]?[0-9]+\.?[0-9]*$ || "$CONCORDANCE" =~ "?" ]]
      then
          if awk 'BEGIN{exit ARGV[1]>=ARGV[2]}' "$CONCORDANCE" ~{concordance_threshold}
          then
            # Less than threshold, fail.
            echo "Genotype Concordance Failure"
            exit 1
          fi
      else
          echo "GENOTYPE_CONCORDANCE must be a number and the file should only have 2 lines of data"
          exit 1;
      fi
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    disks: "local-disk " + disk_size + " HDD"
    memory: "7000 MiB"
    preemptible: preemptible_tries
  }

  output {
    File summary_metrics = "~{output_metrics_basename}.genotype_concordance_summary_metrics"
    File detail_metrics = "~{output_metrics_basename}.genotype_concordance_detail_metrics"
    File contingency_metrics = "~{output_metrics_basename}.genotype_concordance_contingency_metrics"
    File output_vcf = "~{output_metrics_basename}.genotype_concordance.vcf.gz"
    File output_txt = report_name
  }
}
