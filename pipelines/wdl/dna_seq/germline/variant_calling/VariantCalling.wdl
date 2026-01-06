version 1.0

import "../../../../../tasks/wdl/GermlineVariantDiscovery.wdl" as Calling
import "../../../../../tasks/wdl/Qc.wdl" as QC
import "../../../../../tasks/wdl/Utilities.wdl" as Utils
import "../../../../../tasks/wdl/BamProcessing.wdl" as BamProcessing
import "../../../../../tasks/wdl/DragenTasks.wdl" as DragenTasks

workflow VariantCalling {


  String pipeline_version = "2.2.7"


  input {
    Boolean run_dragen_mode_variant_calling = false
    Boolean use_spanning_event_genotyping = true
    File calling_interval_list
    File evaluation_interval_list
    Int haplotype_scatter_count
    Int break_bands_at_multiples_of
    Float? contamination
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_str
    File dbsnp_vcf
    File dbsnp_vcf_index
    String base_file_name
    String final_vcf_base_name
    Int agg_preemptible_tries
    Boolean make_gvcf = true
    Boolean make_bamout = false
    Boolean use_gatk3_haplotype_caller = false
    Boolean skip_reblocking = false
    Boolean use_dragen_hard_filtering = false
    String cloud_provider
    String? billing_project
  }

  if (false) {
    String? none = "None"
  }

  # docker images
  String gatk_docker_gcp = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
  String gatk_docker_azure = "terrapublic.azurecr.io/gatk:4.6.1.0"
  String gatk_docker = if cloud_provider == "gcp" then gatk_docker_gcp else gatk_docker_azure

  String gatk_1_3_docker_gcp = "us.gcr.io/broad-gotc-prod/gatk:1.3.0-4.2.6.1-1649964384"
  String gatk_1_3_docker_azure = "us.gcr.io/broad-gotc-prod/gatk:1.3.0-4.2.6.1-1649964384"
  String gatk_1_3_docker = if cloud_provider == "gcp" then gatk_1_3_docker_gcp else gatk_1_3_docker_azure

  String picard_python_docker_gcp = "us.gcr.io/broad-gotc-prod/picard-python:1.0.0-2.26.10-1663951039"
  String picard_python_docker_azure = "dsppipelinedev.azurecr.io/picard-python:1.0.0-2.26.10-1663951039"
  String picard_python_docker = if cloud_provider == "gcp" then picard_python_docker_gcp else picard_python_docker_azure

  String picard_cloud_docker_gcp = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
  String picard_cloud_docker_azure = "dsppipelinedev.azurecr.io/picard-cloud:2.26.10"
  String picard_cloud_docker = if cloud_provider == "gcp" then picard_cloud_docker_gcp else picard_cloud_docker_azure


  # make sure either gcp or azr is supplied as cloud_provider input
  if ((cloud_provider != "gcp") && (cloud_provider != "azure")) {
    call Utils.ErrorWithMessage as ErrorMessageIncorrectInput {
      input:
        message = "cloud_provider must be supplied with either 'gcp' or 'azure'."
    }
  }

  parameter_meta {
    make_bamout: "For CNNScoreVariants to run with a 2D model, a bamout must be created by HaplotypeCaller. The bamout is a bam containing information on how HaplotypeCaller remapped reads while it was calling variants. See https://gatkforums.broadinstitute.org/gatk/discussion/5484/howto-generate-a-bamout-file-showing-how-haplotypecaller-has-remapped-sequence-reads for more details."
    run_dragen_mode_variant_calling: "Run variant calling using the DRAGEN-GATK pipeline, false by default."
  }

  if (run_dragen_mode_variant_calling) {
    call DragenTasks.CalibrateDragstrModel as DragstrAutoCalibration {
      input:
        ref_fasta = ref_fasta,
        ref_fasta_idx = ref_fasta_index,
        ref_dict = ref_dict,
        alignment = input_bam,
        alignment_index = input_bam_index,
        str_table_file = select_first([ref_str]),
        docker = gatk_docker
    }
  }


  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call Utils.ScatterIntervalList as ScatterIntervalList {
    input:
      interval_list = calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      docker = picard_python_docker
  }

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over WGS calling intervals
  scatter (scattered_interval_list in ScatterIntervalList.out) {

    if (use_gatk3_haplotype_caller) {
      call Calling.HaplotypeCaller_GATK35_GVCF as HaplotypeCallerGATK3 {
        input:
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          interval_list = scattered_interval_list,
          gvcf_basename = base_file_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          contamination = contamination,
          preemptible_tries = agg_preemptible_tries,
          hc_scatter = hc_divisor,
          docker = gatk_1_3_docker
      }
    }

    if (!use_gatk3_haplotype_caller) {
      # Generate GVCF by interval
      call Calling.HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4 {
        input:
          contamination = if run_dragen_mode_variant_calling then 0 else contamination,
          input_bam = input_bam,
          input_bam_index = input_bam_index,
          interval_list = scattered_interval_list,
          vcf_basename = base_file_name,
          ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          hc_scatter = hc_divisor,
          make_gvcf = make_gvcf,
          make_bamout = make_bamout,
          run_dragen_mode_variant_calling = run_dragen_mode_variant_calling,
          use_dragen_hard_filtering = use_dragen_hard_filtering,
          use_spanning_event_genotyping = use_spanning_event_genotyping,
          dragstr_model = DragstrAutoCalibration.dragstr_model,
          preemptible_tries = agg_preemptible_tries,
          gatk_docker = gatk_docker,
          billing_project = billing_project
       }

      if (use_dragen_hard_filtering) {
        call Calling.DragenHardFilterVcf as DragenHardFilterVcf {
          input:
            input_vcf = HaplotypeCallerGATK4.output_vcf,
            input_vcf_index = HaplotypeCallerGATK4.output_vcf_index,
            make_gvcf = make_gvcf,
            vcf_basename = base_file_name,
            preemptible_tries = agg_preemptible_tries,
            gatk_docker = gatk_docker
        }
      }

      # If bamout files were created, we need to sort and gather them into one bamout
      if (make_bamout) {
        call BamProcessing.SortSam as SortBamout {
          input:
            input_bam = HaplotypeCallerGATK4.bamout,
            output_bam_basename = final_vcf_base_name,
            preemptible_tries = agg_preemptible_tries,
            compression_level = 2,
            docker = picard_cloud_docker
        }
      }
    }

    File vcfs_to_merge = select_first([HaplotypeCallerGATK3.output_gvcf, DragenHardFilterVcf.output_vcf, HaplotypeCallerGATK4.output_vcf])
    File vcf_indices_to_merge = select_first([HaplotypeCallerGATK3.output_gvcf_index, DragenHardFilterVcf.output_vcf_index, HaplotypeCallerGATK4.output_vcf_index])
  }

  # Combine by-interval (g)VCFs into a single sample (g)VCF file
  String hard_filter_suffix = if use_dragen_hard_filtering then ".hard-filtered" else ""
  String merge_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  call Calling.MergeVCFs as MergeVCFs {
    input:
      input_vcfs = vcfs_to_merge,
      input_vcfs_indexes = vcf_indices_to_merge,
      output_vcf_name = final_vcf_base_name + hard_filter_suffix + merge_suffix,
      preemptible_tries = agg_preemptible_tries,
      docker = picard_cloud_docker
  }

  if (make_gvcf && !skip_reblocking) {
    call Calling.Reblock as Reblock {
      input:
        gvcf = MergeVCFs.output_vcf,
        gvcf_index = MergeVCFs.output_vcf_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        output_vcf_filename = basename(MergeVCFs.output_vcf, ".g.vcf.gz") + ".rb.g.vcf.gz",
        docker_path = gatk_docker
    }
  }

  if (make_bamout) {
    call MergeBamouts {
      input:
        bams = select_all(SortBamout.output_bam),
        output_base_name = final_vcf_base_name
    }
  }

  # Validate the (g)VCF output of HaplotypeCaller
  call QC.ValidateVCF as ValidateVCF {
    input:
      input_vcf = select_first([Reblock.output_vcf, MergeVCFs.output_vcf]),
      input_vcf_index = select_first([Reblock.output_vcf_index, MergeVCFs.output_vcf_index]),
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      calling_interval_list = calling_interval_list,
      is_gvcf = make_gvcf,
      extra_args = if (skip_reblocking == false) then "--no-overlaps" else "",
      docker_path = gatk_docker,
      preemptible_tries = agg_preemptible_tries
  }

  # QC the (g)VCF
  call QC.CollectVariantCallingMetrics as CollectVariantCallingMetrics {
    input:
      input_vcf = select_first([Reblock.output_vcf, MergeVCFs.output_vcf]),
      input_vcf_index = select_first([Reblock.output_vcf_index, MergeVCFs.output_vcf_index]),
      metrics_basename = final_vcf_base_name,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      ref_dict = ref_dict,
      evaluation_interval_list = evaluation_interval_list,
      is_gvcf = make_gvcf,
      preemptible_tries = agg_preemptible_tries,
      docker = picard_cloud_docker
  }

  output {
    File vcf_summary_metrics = CollectVariantCallingMetrics.summary_metrics
    File vcf_detail_metrics = CollectVariantCallingMetrics.detail_metrics
    File output_vcf = select_first([Reblock.output_vcf, MergeVCFs.output_vcf])
    File output_vcf_index = select_first([Reblock.output_vcf_index, MergeVCFs.output_vcf_index])
    File? bamout = MergeBamouts.output_bam
    File? bamout_index = MergeBamouts.output_bam_index
  }
  meta {
    allowNestedInputs: true
  }
}

# This task is here because merging bamout files using Picard produces an error.
task MergeBamouts {

  input {
    Array[File] bams
    String output_base_name
  }

  Int disk_size = ceil(size(bams, "GiB") * 2) + 10

  command <<<
    samtools merge ~{output_base_name}.bam ~{sep=" " bams}
    samtools index ~{output_base_name}.bam
    mv ~{output_base_name}.bam.bai ~{output_base_name}.bai
  >>>

  output {
    File output_bam = "~{output_base_name}.bam"
    File output_bam_index = "~{output_base_name}.bai"
  }

  runtime {
    docker: "biocontainers/samtools:1.3.1"
    memory: "4 GiB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: 3
    cpu: 1
  }
}
