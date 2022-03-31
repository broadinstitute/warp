version 1.0

#import "../jukebox_vc.wdl" as JukeboxVCWdl TODO: put this back when not testing
import "https://raw.githubusercontent.com/broadinstitute/warp/CheckFingerprint_v1.0.1/pipelines/broad/qc/CheckFingerprint.wdl" as FP
# previously: import "../../../../pipelines/broad/qc/CheckFingerprint.wdl" as FP
# not used: import "../../../../tasks/broad/Utilities.wdl" as utils
# possibly add in a File monitoring_log into this during troubleshooting?

workflow BroadInternalJbxWrapper {

  String pipeline_version = "1.0.1"

  input {
  
  # INTERNAL BROAD - TDR AND FINGERPRINTING INPUTS
  String environment
  String sample_lsid
  String output_basename
  File haplotype_database_file = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt"
  File vault_token_path
  String staging_bucket
  
  String? tdr_dataset_uuid
  String? tdr_sample_id
  String? tdr_staging_bucket
  
  #  JUKEBOX VC PIPELINE INPUT, CAN BE RECALIBRATED BAM OR CRAM
  # if skip_alignment_and_markduplicates = true it will work just on the first file in the list
  Array[File] input_cram_bam_list
  Array[File]? input_cram_bam_index

  String base_file_name
  Int? preemptible_tries
  Int? evaluation_preemptible_tries

  File cache_populate_script
  Array[File] ref_fastas_cram

  File ref_fasta
  File ref_fasta_index
  File ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_dict
  File ref_pac
  File ref_sa
  File ref_fasta_sdf
  File sv_calls_excl_regions
  File ref_dbsnp
  File ref_dbsnp_index

  Boolean? use_v_aware_alignment
  File? v_aware_vcf

  File coverage_intervals
  File wgs_calling_interval_list
  File wgs_coverage_interval_list
  Int break_bands_at_multiples_of
  Int haplotype_scatter_count
  Int? read_length
  Int? reads_per_split

  # parameters for rsq filtering
  Boolean filter_by_rsq
  Float rsq_threshold

  # parameters for adapter trimming
  String illumina_adapter_5p
  String illumina_adapter_3p
  Int? umi_length_3p
  Int? umi_length_5p
  Float? error_rate_5p
  Float? error_rate_3p
  Int? min_overlap_5p
  Int? min_overlap_3p

  String broad_gatk_docker
  String crammer_docker
  String jb_gatk_docker
  String gatk_markduplicates_docker
  String jukebox_vc_docker
  String gitc_docker
  String ua_docker

  String? gitc_path_override
  String delly_docker

  File? picard_jar_override

  File? create_report_notebook_override
  String? mark_duplicates_extra_args
  String? hc_extra_args
  String? ua_extra_args
  Array[File]? known_ground_truth
  String? gt_left_sample
  String? gt_right_sample

  Boolean? make_gvcf_override
  Boolean? merge_bam_file_override
  # control flow
  Boolean? align_bwa
  Boolean? skip_alignment
  Boolean? skip_alignment_and_markduplicates
  Boolean? skip_haplotypecaller
  Boolean? skip_svs
  Boolean? converted_library

  # Contamination
  String contamination_sites_path
  File contamination_sites_vcf
  File contamination_sites_vcf_index

  # VCF post-processing
  Array[File] annotation_intervals
  Array[File] ground_truth_files
  File haplotype_database_gt
  File? filtering_model_no_gt

  File af_only_gnomad
  File af_only_gnomad_index
  String reference_version
  File funcotator_germline_data_sources_tar_gz

  Boolean filter_cg_insertions
  File? filtering_blacklist_file

  File? training_blacklist_file
  Int? exome_weight
  String? exome_weight_annotation

  File? interval_list_override
  File runs_file
  String? filtering_model_no_gt_name_override
  Int? delly_memory_override
  String? filtering_model_with_gt_name_override

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size
  Int? additional_metrics_disk # will be added to increase_disk_size

  Boolean collect_statistics
  Boolean? no_address_override
  String? override_input_ending # For drs where there is no extension. Should be "is_cram" or "is_bam"
  Float max_duplication_in_reasonable_sample
  Float max_chimerism_in_reasonable_sample

  Boolean featuremap_generate
  File featuremap_interval_list
  Int featuremap_scatter_count
  Int featuremap_min_mapq
  Int featuremap_snv_identical_bases
  Int featuremap_snv_identical_bases_after
  Int featuremap_min_score
  Int featuremap_limit_score
  String featuremap_extra_args

  File reference_gaps_intervals
  File centromere_intervals

  # When running on Terra, use workspace.name as this input to ensure that all tasks will only cache hit to runs in your
  # own workspace. This will prevent call caching from failing with "Cache Miss (10 failed copy attempts)". Outside of
  # Terra this can be left as the default empty String. This dummy input is only needed for tasks that have no inputs
  # specific to the sample being run (such as GetBwaVersion which does not take in any sample data).
  String dummy_input_for_call_caching = ""

  String collab_sample_id_run_id
  }

  # PARAMETER DEFINITIONS
  parameter_meta {
    sample_lsid: "The sample lsid (an identifier used to retrieve fingerrints from Mercury)"
    input_cram_bam_list: "Input bams/crams. Will be realigned."
    base_file_name: "Basename for the output files."
    environment: "The environment (dev or prod) used for determining which service to use to retrieve Mercury fingerprints"
    vault_token_path: "The path to the vault token used for accessing the Mercury Fingerprint Store"
    output_basename: "String used as a prefix in workflow output files"
    tdr_dataset_uuid: "Optional String used to define the Terra Data Repo dataset to which outputs will be ingested, if populated"
    tdr_staging_bucket: "Optional String defining the GCS bucket to use to stage files for loading to TDR. Workspace bucket is recommended"
  }

  # CALL JUKEBOX VC WDL
  
  #TODO actually import real workflow here
  call JukeboxVC {
    input:
      input_cram_bam_list = input_cram_bam_list,
      input_cram_bam_index = input_cram_bam_index,
      base_file_name = base_file_name,
      preemptible_tries = preemptible_tries,
      evaluation_preemptible_tries = evaluation_preemptible_tries,
      cache_populate_script = cache_populate_script,
      ref_fastas_cram = ref_fastas_cram,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_alt = ref_alt,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_bwt = ref_bwt,
      ref_dict = ref_dict,
      ref_pac = ref_pac,
      ref_sa = ref_sa,
      ref_fasta_sdf = ref_fasta_sdf,
      sv_calls_excl_regions = sv_calls_excl_regions,
      ref_dbsnp = ref_dbsnp,
      ref_dbsnp_index = ref_dbsnp_index,
      use_v_aware_alignment = use_v_aware_alignment,
      v_aware_vcf = v_aware_vcf,
      coverage_intervals = coverage_intervals,
      wgs_calling_interval_list = wgs_calling_interval_list,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      haplotype_scatter_count = haplotype_scatter_count,
      read_length = read_length,
      reads_per_split = reads_per_split,
      filter_by_rsq = filter_by_rsq,
      rsq_threshold = rsq_threshold,
      illumina_adapter_5p = illumina_adapter_5p,
      illumina_adapter_3p = illumina_adapter_3p,
      umi_length_3p = umi_length_3p,
      umi_length_5p = umi_length_5p,
      error_rate_5p = error_rate_5p,
      error_rate_3p = error_rate_3p,
      min_overlap_5p = min_overlap_5p,
      min_overlap_3p = min_overlap_3p,
      broad_gatk_docker = broad_gatk_docker,
      crammer_docker = crammer_docker,
      jb_gatk_docker = jb_gatk_docker,
      gatk_markduplicates_docker = gatk_markduplicates_docker,
      jukebox_vc_docker = jukebox_vc_docker,
      gitc_docker = gitc_docker,
      ua_docker = ua_docker,
      gitc_path_override = gitc_path_override,
      delly_docker = delly_docker,
      picard_jar_override = picard_jar_override,
      create_report_notebook_override = create_report_notebook_override,
      mark_duplicates_extra_args = mark_duplicates_extra_args,
      hc_extra_args = hc_extra_args,
      ua_extra_args = ua_extra_args,
      known_ground_truth = known_ground_truth,
      gt_left_sample = gt_left_sample,
      gt_right_sample = gt_right_sample,
      make_gvcf_override = make_gvcf_override,
      merge_bam_file_override = merge_bam_file_override,
      align_bwa = align_bwa,
      skip_alignment = skip_alignment,
      skip_alignment_and_markduplicates = skip_alignment_and_markduplicates,
      skip_haplotypecaller = skip_haplotypecaller,
      skip_svs = skip_svs,
      converted_library = converted_library,
      contamination_sites_path = contamination_sites_path,
      contamination_sites_vcf = contamination_sites_vcf,
      contamination_sites_vcf_index = contamination_sites_vcf_index,
      annotation_intervals = annotation_intervals,
      ground_truth_files = ground_truth_files,
      haplotype_database_gt = haplotype_database_gt,
      filtering_model_no_gt = filtering_model_no_gt,
      af_only_gnomad = af_only_gnomad,
      af_only_gnomad_index = af_only_gnomad_index,
      reference_version = reference_version,
      funcotator_germline_data_sources_tar_gz = funcotator_germline_data_sources_tar_gz,
      filter_cg_insertions = filter_cg_insertions,
      filtering_blacklist_file = filtering_blacklist_file,
      training_blacklist_file = training_blacklist_file,
      exome_weight = exome_weight,
      exome_weight_annotation = exome_weight_annotation,
      interval_list_override = interval_list_override,
      runs_file = runs_file,
      filtering_model_no_gt_name_override = filtering_model_no_gt_name_override,
      delly_memory_override = delly_memory_override,
      filtering_model_with_gt_name_override = filtering_model_with_gt_name_override,
      increase_disk_size = increase_disk_size,
      additional_metrics_disk = additional_metrics_disk,
      collect_statistics = collect_statistics,
      no_address_override = no_address_override,
      override_input_ending = override_input_ending,
      max_duplication_in_reasonable_sample = max_duplication_in_reasonable_sample,
      max_chimerism_in_reasonable_sample = max_chimerism_in_reasonable_sample,
      featuremap_generate = featuremap_generate,
      featuremap_interval_list = featuremap_interval_list,
      featuremap_scatter_count = featuremap_scatter_count,
      featuremap_min_mapq = featuremap_min_mapq,
      featuremap_snv_identical_bases = featuremap_snv_identical_bases,
      featuremap_snv_identical_bases_after = featuremap_snv_identical_bases_after,
      featuremap_min_score = featuremap_min_score,
      featuremap_limit_score = featuremap_limit_score,
      featuremap_extra_args = featuremap_extra_args,
      reference_gaps_intervals = reference_gaps_intervals,
      centromere_intervals = centromere_intervals,
      dummy_input_for_call_caching = dummy_input_for_call_caching
  }

  # CALL FINGERPRINTING 
  
  call FP.CheckFingerprint as CheckFingerprint {
  # CONFIRMED PUBLIC VERSION 1.0.1 has not changed since release
    input:
      input_bam = JukeboxVC.output_cram,
      input_bam_index = JukeboxVC.output_cram_index,
      sample_alias = JukeboxVC.sample_name,
      sample_lsid = sample_lsid,
      output_basename = output_basename,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      read_fingerprint_from_mercury = true,
      haplotype_database_file = haplotype_database_file,
      environment = environment,
      vault_token_path = vault_token_path
  }

  # CONVERT METRICS TO JSON
  call ConvertMetricsToJson as ConvertWgsMetricsToJson {
    input:
      metric_file = select_first([JukeboxVC.wgs_metrics]),
      basename = base_file_name + ".wgs_metrics",
      collab_sample_id_run_id = collab_sample_id_run_id
  }

  call ConvertMetricsToJson as ConvertAlignmentSummaryMetricsToJson {
    input:
      metric_file = select_first([JukeboxVC.agg_alignment_summary_metrics]),
      basename = base_file_name + ".alignment_metrics",
      collab_sample_id_run_id = collab_sample_id_run_id
  }

  call ConvertMetricsToJson as ConvertGcBiasMetricsToJson {
    input:
      metric_file = select_first([JukeboxVC.agg_gc_bias_summary_metrics]),
      basename = base_file_name + ".gc_metrics",
      collab_sample_id_run_id = collab_sample_id_run_id
  }

  call ConvertMetricsToJson as ConvertGcBiasDetailMetricsToJson {
    input:
      metric_file = select_first([JukeboxVC.agg_gc_bias_detail_metrics]),
      basename = base_file_name + ".gc_detail_metrics",
      collab_sample_id_run_id = collab_sample_id_run_id,
      detail_metrics = true
  }

  call ConvertMetricsToJson as ConvertQualityDistributionMetricsToJson {
    input:
      metric_file = select_first([JukeboxVC.agg_quality_distribution_metrics]),
      basename = base_file_name + ".quality_distribution_metrics",
      collab_sample_id_run_id = collab_sample_id_run_id
  }

  call ConvertMetricsToJson as ConvertDuplicateMetricsToJson {
    input:
      metric_file = select_first([JukeboxVC.duplicate_metrics]),
      basename = base_file_name + ".duplicate_metrics",
      collab_sample_id_run_id = collab_sample_id_run_id
  }

  call ConvertMetricsToJson as ConvertQualityYieldMetricsToJson {
    input:
      metric_file = select_first([JukeboxVC.quality_yield_metrics]),
      basename = base_file_name + ".quality_yield_metrics",
      collab_sample_id_run_id = collab_sample_id_run_id
  }

  call ConvertMetricsToJson as ConvertRawWgsMetricsToJson {
    input:
      metric_file = select_first([JukeboxVC.raw_wgs_metrics]),
      basename = base_file_name + ".raw_wgs_metrics",
      collab_sample_id_run_id = collab_sample_id_run_id
}

  call ConvertContaminationMetricsToJson {
    input:
      metric_file = select_first([JukeboxVC.selfSM]),
      basename = base_file_name + ".contamination_metrics",
      collab_sample_id_run_id = collab_sample_id_run_id
  }
  

  # CALL TDR TASKS TO FORMAT JSON
  
  if (defined(tdr_dataset_uuid) && defined(tdr_sample_id) && defined(tdr_staging_bucket)) {
    call formatPipelineOutputs {
      input:
        collab_sample_id_run_id = collab_sample_id_run_id,
        agg_alignment_summary_metrics = ConvertAlignmentSummaryMetricsToJson.metrics_json,
        agg_gc_bias_detail_metrics = ConvertGcBiasDetailMetricsToJson.metrics_json,
        agg_gc_bias_summary_metrics = ConvertGcBiasMetricsToJson.metrics_json,
        agg_quality_distribution_metrics = ConvertQualityDistributionMetricsToJson.metrics_json,
        duplicate_metrics = ConvertDuplicateMetricsToJson.metrics_json,
        quality_yield_metrics = ConvertQualityYieldMetricsToJson.metrics_json,
        raw_wgs_metrics = ConvertRawWgsMetricsToJson.metrics_json,
        selfSM = ConvertContaminationMetricsToJson.metrics_json,
        wgs_metrics = ConvertWgsMetricsToJson.metrics_json,
        agg_alignment_summary_pdf = JukeboxVC.agg_alignment_summary_pdf,
        agg_gc_bias_pdf = JukeboxVC.agg_gc_bias_pdf,
        agg_quality_distribution_pdf = JukeboxVC.agg_quality_distribution_pdf,
        barcode = JukeboxVC.barcode,
        chimerism_rate = JukeboxVC.chimerism_rate,
        contamination = JukeboxVC.contamination,
        duplication_rate = JukeboxVC.duplication_rate,
        filtered_vcf = JukeboxVC.filtered_vcf,
        filtered_vcf_index = JukeboxVC.filtered_vcf_index,
        flow_order = JukeboxVC.flow_order,
        id = JukeboxVC.id,
        is_outlier_data = JukeboxVC.is_outlier_data,
        model_h5_no_gt = JukeboxVC.model_h5_no_gt,
        model_pkl_no_gt = JukeboxVC.model_pkl_no_gt,
        output_cram = JukeboxVC.output_cram,
        output_cram_index = JukeboxVC.output_cram_index,
        output_cram_md5 = JukeboxVC.output_cram_md5,
        output_gvcf = JukeboxVC.output_gvcf,
        output_gvcf_index = JukeboxVC.output_gvcf_index,
        output_vcf = JukeboxVC.output_vcf,
        output_vcf_index = JukeboxVC.output_vcf_index,
        sample_name = JukeboxVC.sample_name,
        aggregated_metrics_h5 = JukeboxVC.aggregated_metrics_h5,
        aggregated_metrics_json = JukeboxVC.aggregated_metrics_json,
        comparison_output = JukeboxVC.comparison_output,
        coverage_per_motif = JukeboxVC.coverage_per_motif,
        featuremap = JukeboxVC.featuremap,
        featuremap_dataframe = JukeboxVC.featuremap_dataframe,
        featuremap_index = JukeboxVC.featuremap_index,
        featuremap_single_substitutions = JukeboxVC.featuremap_single_substitutions,
        featuremap_single_substitutions_dataframe = JukeboxVC.featuremap_single_substitutions_dataframe,
        featuremap_single_substitutions_index = JukeboxVC.featuremap_single_substitutions_index,
        fingerprints_csv = JukeboxVC.fingerprints_csv,
        funcotator_vcf = JukeboxVC.funcotator_vcf,
        funcotator_vcf_index = JukeboxVC.funcotator_vcf_index,
        gt_sample_name = JukeboxVC.gt_sample_name,
        model_h5 = JukeboxVC.model_h5,
        model_pkl = JukeboxVC.model_pkl,
        no_gt_report_html = JukeboxVC.no_gt_report_html,
        orig_vcf = JukeboxVC.orig_vcf,
        orig_vcf_index = JukeboxVC.orig_vcf_index,
        short_report_h5 = JukeboxVC.short_report_h5,
        short_report_html = JukeboxVC.short_report_html,
        sv_calls = JukeboxVC.sv_calls,
        sv_calls_index = JukeboxVC.sv_calls_index,
        varStats = JukeboxVC.varStats,
        evaluate_report_h5 = JukeboxVC.evaluate_report_h5,
        extended_report_h5 = JukeboxVC.extended_report_h5,
        extended_report_html = JukeboxVC.extended_report_html,
        snp_error_rate = JukeboxVC.snp_error_rate,
        snp_error_rate_plot = JukeboxVC.snp_error_rate_plot,

        output_basename = output_basename,
    
        haplotype_bam = JukeboxVC.haplotype_bam,
        haplotype_bam_index = JukeboxVC.haplotype_bam_index,

        fingerprint_summary_metrics_file = CheckFingerprint.fingerprint_summary_metrics_file,
        fingerprint_detail_metrics_file = CheckFingerprint.fingerprint_detail_metrics_file
    }

    call updateOutputsInTDR {
      input:
        tdr_dataset_uuid = select_first([tdr_dataset_uuid, ""]),
        outputs_json = formatPipelineOutputs.pipeline_outputs_json,
        primary_key = collab_sample_id_run_id,
        staging_bucket = select_first([tdr_staging_bucket, ""])
    }
  }

  output {
    File? picard_fingerprint_summary_metrics = CheckFingerprint.fingerprint_summary_metrics_file
    File? picard_fingerprint_detail_metrics = CheckFingerprint.fingerprint_detail_metrics_file

    File? agg_alignment_summary_metrics = JukeboxVC.agg_alignment_summary_metrics
    File? agg_gc_bias_detail_metrics = JukeboxVC.agg_gc_bias_detail_metrics
    File? agg_gc_bias_summary_metrics = JukeboxVC.agg_gc_bias_summary_metrics
    File? agg_quality_distribution_metrics = JukeboxVC.agg_quality_distribution_metrics
    File? duplicate_metrics = JukeboxVC.duplicate_metrics
    File? quality_yield_metrics = JukeboxVC.quality_yield_metrics
    File? raw_wgs_metrics = JukeboxVC.raw_wgs_metrics
    File? selfSM = JukeboxVC.selfSM
    File? wgs_metrics = JukeboxVC.wgs_metrics
    File? agg_alignment_summary_pdf = JukeboxVC.agg_alignment_summary_pdf
    File? agg_gc_bias_pdf = JukeboxVC.agg_gc_bias_pdf
    File? agg_quality_distribution_pdf = JukeboxVC.agg_quality_distribution_pdf
    String barcode = JukeboxVC.barcode
    Float? chimerism_rate = JukeboxVC.chimerism_rate
    Float? contamination = JukeboxVC.contamination
    Float? duplication_rate = JukeboxVC.duplication_rate
    File? filtered_vcf = JukeboxVC.filtered_vcf
    File? filtered_vcf_index = JukeboxVC.filtered_vcf_index
    String flow_order = JukeboxVC.flow_order
    String id = JukeboxVC.id
    Boolean? is_outlier_data = JukeboxVC.is_outlier_data
    File? model_h5_no_gt = JukeboxVC.model_h5_no_gt
    File? model_pkl_no_gt = JukeboxVC.model_pkl_no_gt
    File? output_cram = JukeboxVC.output_cram
    File? output_cram_index = JukeboxVC.output_cram_index
    File? output_cram_md5 = JukeboxVC.output_cram_md5
    File? output_gvcf = JukeboxVC.output_gvcf
    File? output_gvcf_index = JukeboxVC.output_gvcf_index
    File? output_vcf = JukeboxVC.output_vcf
    File? output_vcf_index = JukeboxVC.output_vcf_index
    String sample_name = JukeboxVC.sample_name
    File? aggregated_metrics_h5 = JukeboxVC.aggregated_metrics_h5
    File? aggregated_metrics_json = JukeboxVC.aggregated_metrics_json
    File? comparison_output = JukeboxVC.comparison_output
    Array[File]? coverage_boxplot_all = JukeboxVC.coverage_boxplot_all
    Array[File]? coverage_depth_parquet_vhq = JukeboxVC.coverage_depth_parquet_vhq
    File? coverage_per_motif = JukeboxVC.coverage_per_motif
    File? featuremap = JukeboxVC.featuremap
    File? featuremap_dataframe = JukeboxVC.featuremap_dataframe
    File? featuremap_index = JukeboxVC.featuremap_index
    File? featuremap_single_substitutions = JukeboxVC.featuremap_single_substitutions
    File? featuremap_single_substitutions_dataframe = JukeboxVC.featuremap_single_substitutions_dataframe
    File? featuremap_single_substitutions_index = JukeboxVC.featuremap_single_substitutions_index
    File? fingerprints_csv = JukeboxVC.fingerprints_csv
    File? funcotator_vcf = JukeboxVC.funcotator_vcf
    File? funcotator_vcf_index = JukeboxVC.funcotator_vcf_index
    String? gt_sample_name = JukeboxVC.gt_sample_name
    File? model_h5 = JukeboxVC.model_h5
    File? model_pkl = JukeboxVC.model_pkl
    File? no_gt_report_html = JukeboxVC.no_gt_report_html
    File? orig_vcf = JukeboxVC.orig_vcf
    File? orig_vcf_index = JukeboxVC.orig_vcf_index
    Array[File]? output_bam = JukeboxVC.output_bam
    Array[File]? output_bam_index = JukeboxVC.output_bam_index
    File? short_report_h5 = JukeboxVC.short_report_h5
    File? short_report_html = JukeboxVC.short_report_html
    File? sv_calls = JukeboxVC.sv_calls
    File? sv_calls_index = JukeboxVC.sv_calls_index
    File? varStats = JukeboxVC.varStats

    File? haplotype_bam = JukeboxVC.haplotype_bam
    File? haplotype_bam_index = JukeboxVC.haplotype_bam_index
  }
}

task ConvertContaminationMetricsToJson {
  input {
      File metric_file
      String basename
      String collab_sample_id_run_id
  }
  command <<<
    python3 <<CODE
    import pandas as pd
    import os
    import json

    df = pd.read_csv("~{metric_file}", sep="\t")
    df_no_nan = df.where(pd.notnull(df), None)
    header_dict = {'id': '~{collab_sample_id_run_id}', 'metric_class': "contamination"}
    d = {'header': header_dict, 'metrics': df.to_dict(orient='records')}
    with open('~{basename}.json', 'w') as json_file:
        json.dump(d, json_file, indent=2)
    CODE
  >>>
  output {
    File metrics_json = "~{basename}.json"
  }
  runtime {
    docker: "broadinstitute/horsefish:tdr_import_v1.1"
    disks: "local-disk 10 HDD"
  }
}

task ConvertMetricsToJson {
  input {
    File metric_file
    String basename
    String collab_sample_id_run_id 
    Boolean detail_metrics = false
  }
  #String nrows_string = if detail_metrics then "" else ", nrows=1"
  
  command <<<
    python3 <<CODE
    import pandas as pd
    import os
    import json

    def parse_metrics(metric_file):
      """Parses Picard metrics file
      Parameters
      ----------
      metric_file : str
          Picard metric file
      Returns
      -------
      res0 : str
          Program header line
      res1 : str
          Picard file metrics class
      res2 : pd.DataFrame
          Picard metrics table
      res3 : pd.DataFrame
          Picard Histogram output
      """
      lines_in_hist = 0
      with open(metric_file) as infile:
          out = next(infile)
          while not out.startswith("# "):
              out = next(infile)
          res0 = out.strip()
      try:
          with open(metric_file) as infile:
              out = next(infile)
              while not out.startswith("## HISTOGRAM"):
                  out = next(infile)

              res3 = pd.read_csv(infile, sep="\t")
              lines_in_hist = len(res3.index) + 2 #extra blank line at the end and the header line need to be counted
      except StopIteration:
          res3 = None
      try:
        with open(metric_file) as infile:
            out = next(infile)
            while not out.startswith("## METRICS CLASS"):
                out = next(infile)

            res1 = out.strip().split('\t')[1].split('.')[-1]
            res2 = pd.read_csv(infile, sep="\t", skipfooter=lines_in_hist, comment="#").replace({float("nan"): None})
      except StopIteration:
        res1 = None
        res2 = None
      return res0, res1, res2, res3

    if os.path.getsize("~{metric_file}") > 0:
        header, metric_class, stats, histogram = parse_metrics("~{metric_file}")
        if metric_class is not None:
            metric_class = metric_class[metric_class.find("$")+1:]
            header_dict = {'id': '~{collab_sample_id_run_id}', 'program_line': header, 'metric_class': metric_class}
        else:
            header_dict = {'id': '~{collab_sample_id_run_id}', 'program_line': header}
        if histogram is not None:
          hist_dict = histogram.to_dict(orient='split')
          del hist_dict['index']
          if stats is not None:
            d = {'header': header_dict, 'metrics': stats.to_dict(orient='records'), 'histogram': hist_dict}
          else:
            d = {'header': header_dict, 'histogram': hist_dict}
        else:
          d = {'header': header_dict, 'metrics': stats.to_dict(orient='records')}
        with open('~{basename}.json', 'w') as json_file:
          json.dump(d, json_file, indent=2)
    CODE
  >>>
  output {
    File metrics_json = "~{basename}.json"
  }
  runtime {
    docker: "broadinstitute/horsefish:tdr_import_v1.1"
    disks: "local-disk 10 HDD"
  }
}

# DEFINE TDR-SPECIFIC TASKS
  
  task formatPipelineOutputs {
    input {
      String output_basename
      String collab_sample_id_run_id
      
      # Metrics Output Files
      #TODO: update these "optional outputs" to latest version (most will no longer be optional)
      String? agg_alignment_summary_metrics = ""
      String? agg_gc_bias_detail_metrics = ""
      String? agg_gc_bias_summary_metrics = ""
      String? agg_quality_distribution_metrics = ""
      String? duplicate_metrics = ""
      String? quality_yield_metrics = ""
      String? raw_wgs_metrics = ""
      String? selfSM = ""
      String? wgs_metrics = ""
      
      # PDF Output Files
      String? agg_alignment_summary_pdf = ""
      String? agg_gc_bias_pdf = ""
      String? agg_quality_distribution_pdf = ""
      
      String barcode
      Float? chimerism_rate = ""
      Float? contamination = ""
      Float? duplication_rate = ""
      String? filtered_vcf = ""
      String? filtered_vcf_index = ""
      String flow_order
      String? haplotype_bam = ""
      String? haplotype_bam_index = ""
      String id
      Boolean? is_outlier_data
      String? model_h5_no_gt = ""
      String? model_pkl_no_gt = ""
      String? output_cram = ""
      String? output_cram_index = ""
      String? output_cram_md5 = ""
      String? output_gvcf = ""
      String? output_gvcf_index = ""
      String? output_vcf = ""
      String? output_vcf_index = ""
      String sample_name
      String? aggregated_metrics_h5 = ""
      String? aggregated_metrics_json = ""
      String? comparison_output = ""
      String? coverage_per_motif = ""
      String? featuremap = ""
      String? featuremap_dataframe = ""
      String? featuremap_index = ""
      String? featuremap_single_substitutions = ""
      String? featuremap_single_substitutions_dataframe = ""
      String? featuremap_single_substitutions_index = ""
      String? fingerprints_csv = ""
      String? funcotator_vcf = ""
      String? funcotator_vcf_index = ""
      String? gt_sample_name = ""
      String? model_h5 = ""
      String? model_pkl = ""
      String? no_gt_report_html = ""
      String? orig_vcf = ""
      String? orig_vcf_index = ""
      String? short_report_h5 = ""
      String? short_report_html = ""
      String? sv_calls = ""
      String? sv_calls_index = ""
      String? varStats = ""
      String? evaluate_report_h5 = ""
      String? extended_report_h5 = ""
      String? extended_report_html = ""
      String? snp_error_rate = ""
      String? snp_error_rate_plot = ""

      String? fingerprint_detail_metrics_file = ""
      String? fingerprint_summary_metrics_file = ""

      Int cpu = 1
      Int memory_mb = 2000
      Int disk_size_gb = 10
    }

    String outputs_json_file_name = "outputs_to_TDR_~{output_basename}.json"

    String outlier_data_string = if defined(is_outlier_data) && is_outlier_data then "True" else if defined(is_outlier_data) then "False" else ""

    command <<<
          python3 << CODE
          import json

          outputs_dict = {}

          outputs_dict["collab_sample_id_run_id"]="~{collab_sample_id_run_id}"       
          outputs_dict["agg_alignment_summary_metrics"]="~{agg_alignment_summary_metrics}"
          outputs_dict["agg_alignment_summary_pdf"]="~{agg_alignment_summary_pdf}"
          outputs_dict["agg_gc_bias_detail_metrics"]="~{agg_gc_bias_detail_metrics}"
          outputs_dict["agg_gc_bias_pdf"]="~{agg_gc_bias_pdf}"
          outputs_dict["agg_gc_bias_summary_metrics"]="~{agg_gc_bias_summary_metrics}"
          outputs_dict["agg_quality_distribution_metrics"]="~{agg_quality_distribution_metrics}"
          outputs_dict["agg_quality_distribution_pdf"]="~{agg_quality_distribution_pdf}"
          outputs_dict["barcode"]="~{barcode}"
          outputs_dict["chimerism_rate"]="~{chimerism_rate}"
          outputs_dict["contamination"]="~{contamination}"
          outputs_dict["duplicate_metrics"]="~{duplicate_metrics}"
          outputs_dict["duplication_rate"]="~{duplication_rate}"
          outputs_dict["filtered_vcf"]="~{filtered_vcf}"
          outputs_dict["filtered_vcf_index"]="~{filtered_vcf_index}"
          outputs_dict["flow_order"]="~{flow_order}"
          outputs_dict["haplotype_bam"]="~{haplotype_bam}"
          outputs_dict["haplotype_bam_index"]="~{haplotype_bam_index}"
          outputs_dict["id"]="~{id}"
          outputs_dict["is_outlier_data"]="~{outlier_data_string}"
          outputs_dict["model_h5_no_gt"]="~{model_h5_no_gt}"
          outputs_dict["model_pkl_no_gt"]="~{model_pkl_no_gt}"
          outputs_dict["output_cram"]="~{output_cram}"
          outputs_dict["output_cram_index"]="~{output_cram_index}"
          outputs_dict["output_cram_md5"]="~{output_cram_md5}"
          outputs_dict["output_gvcf"]="~{output_gvcf}"
          outputs_dict["output_gvcf_index"]="~{output_gvcf_index}"
          outputs_dict["output_vcf"]="~{output_vcf}"
          outputs_dict["output_vcf_index"]="~{output_vcf_index}"
          outputs_dict["quality_yield_metrics"]="~{quality_yield_metrics}"
          outputs_dict["raw_wgs_metrics"]="~{raw_wgs_metrics}"
          outputs_dict["sample_name"]="~{sample_name}"
          outputs_dict["selfSM"]="~{selfSM}"
          outputs_dict["wgs_metrics"]="~{wgs_metrics}"
          outputs_dict["aggregated_metrics_h5"]="~{aggregated_metrics_h5}"
          outputs_dict["aggregated_metrics_json"]="~{aggregated_metrics_json}"
          outputs_dict["comparison_output"]="~{comparison_output}"
          outputs_dict["coverage_per_motif"]="~{coverage_per_motif}"
          outputs_dict["featuremap"]="~{featuremap}"
          outputs_dict["featuremap_dataframe"]="~{featuremap_dataframe}"
          outputs_dict["featuremap_index"]="~{featuremap_index}"
          outputs_dict["featuremap_single_substitutions"]="~{featuremap_single_substitutions}"
          outputs_dict["featuremap_single_substitutions_dataframe"]="~{featuremap_single_substitutions_dataframe}"
          outputs_dict["featuremap_single_substitutions_index"]="~{featuremap_single_substitutions_index}"
          outputs_dict["fingerprints_csv"]="~{fingerprints_csv}"
          outputs_dict["funcotator_vcf"]="~{funcotator_vcf}"
          outputs_dict["funcotator_vcf_index"]="~{funcotator_vcf_index}"
          outputs_dict["gt_sample_name"]="~{gt_sample_name}"
          outputs_dict["model_h5"]="~{model_h5}"
          outputs_dict["model_pkl"]="~{model_pkl}"
          outputs_dict["no_gt_report_html"]="~{no_gt_report_html}"
          outputs_dict["orig_vcf"]="~{orig_vcf}"
          outputs_dict["orig_vcf_index"]="~{orig_vcf_index}"
          outputs_dict["short_report_h5"]="~{short_report_h5}"
          outputs_dict["short_report_html"]="~{short_report_html}"
          outputs_dict["sv_calls"]="~{sv_calls}"
          outputs_dict["sv_calls_index"]="~{sv_calls_index}"
          outputs_dict["varStats"]="~{varStats}"
          outputs_dict["fingerprint_summary_metrics_file"]="~{fingerprint_summary_metrics_file}"
          outputs_dict["fingerprint_detail_metrics_file"]="~{fingerprint_detail_metrics_file}"
          outputs_dict["evaluate_report_h5"]="~{evaluate_report_h5}"
          outputs_dict["extended_report_h5"]="~{extended_report_h5}"
          outputs_dict["extended_report_html"]="~{extended_report_html}"
          outputs_dict["snp_error_rate"]="~{snp_error_rate}"
          outputs_dict["snp_error_rate_plot"]="~{snp_error_rate_plot}"

          # Write full outputs to file
          with open("~{outputs_json_file_name}", 'w') as outputs_file:
              for key, value in outputs_dict.items():
                  if value == "None" or value == "":
                      outputs_dict[key] = None
              outputs_file.write(json.dumps(outputs_dict))
              outputs_file.write("\n")
          CODE
      >>>

    runtime {
        docker: "broadinstitute/horsefish:tdr_import_v1.1"
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File pipeline_outputs_json = outputs_json_file_name
    }
  }

  task updateOutputsInTDR {
    input {
      String staging_bucket
      String tdr_dataset_uuid
      File outputs_json
      String primary_key

      Int cpu = 1
      Int memory_mb = 2000
      Int disk_size_gb = 10
    }

    String tdr_target_table = "sample"
    String primary_key_field = "collab_sample_id_run_id"

    command <<<
      python -u /scripts/export_pipeline_outputs_to_tdr.py \
        -d "~{tdr_dataset_uuid}" \
        -b "~{staging_bucket}" \
        -t "~{tdr_target_table}" \
        -o "~{outputs_json}" \
        --primary_key_field "~{primary_key_field}" \
        --primary_key_value "~{primary_key}"
    >>>

    runtime {
      docker: "broadinstitute/horsefish:tdr_import_v1.1"
      cpu: cpu
      memory: "~{memory_mb} MiB"
      disks: "local-disk ~{disk_size_gb} HDD"
    }

    output {
      File ingest_logs = stdout()
    }
  }

#TODO: delete this fake task used for testing
task JukeboxVC {
  input {
    Array[File] input_cram_bam_list
    Array[File]? input_cram_bam_index
    String base_file_name
    Int? preemptible_tries
    Int? evaluation_preemptible_tries
    File cache_populate_script
    Array[File] ref_fastas_cram
    File ref_fasta
    File ref_fasta_index
    File ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_dict
    File ref_pac
    File ref_sa
    File ref_fasta_sdf
    File sv_calls_excl_regions
    File ref_dbsnp
    File ref_dbsnp_index
    Boolean? use_v_aware_alignment
    File? v_aware_vcf
    File coverage_intervals
    File wgs_calling_interval_list
    File wgs_coverage_interval_list
    Int break_bands_at_multiples_of
    Int haplotype_scatter_count
    Int? read_length
    Int? reads_per_split
    Boolean filter_by_rsq
    Float rsq_threshold
    String illumina_adapter_5p
    String illumina_adapter_3p
    Int? umi_length_3p
    Int? umi_length_5p
    Float? error_rate_5p
    Float? error_rate_3p
    Int? min_overlap_5p
    Int? min_overlap_3p
    String broad_gatk_docker
    String crammer_docker
    String jb_gatk_docker
    String gatk_markduplicates_docker
    String jukebox_vc_docker
    String gitc_docker
    String ua_docker
    String? gitc_path_override
    String delly_docker
    File? picard_jar_override
    File? create_report_notebook_override
    String? mark_duplicates_extra_args
    String? hc_extra_args
    String? ua_extra_args
    Array[File]? known_ground_truth
    String? gt_left_sample
    String? gt_right_sample
    Boolean? make_gvcf_override
    Boolean? merge_bam_file_override
    Boolean? align_bwa
    Boolean? skip_alignment
    Boolean? skip_alignment_and_markduplicates
    Boolean? skip_haplotypecaller
    Boolean? skip_svs
    Boolean? converted_library
    String contamination_sites_path
    File contamination_sites_vcf
    File contamination_sites_vcf_index
    Array[File] annotation_intervals
    Array[File] ground_truth_files
    File haplotype_database_gt
    File? filtering_model_no_gt
    File af_only_gnomad
    File af_only_gnomad_index
    String reference_version
    File funcotator_germline_data_sources_tar_gz
    Boolean filter_cg_insertions
    File? filtering_blacklist_file
    File? training_blacklist_file
    Int? exome_weight
    String? exome_weight_annotation
    File? interval_list_override
    File runs_file
    String? filtering_model_no_gt_name_override
    Int? delly_memory_override
    File? filtering_model_with_gt_name_override
    Int? increase_disk_size
    Int? additional_metrics_disk
    Boolean collect_statistics
    Boolean? no_address_override
    String? override_input_ending
    Float max_duplication_in_reasonable_sample
    Float max_chimerism_in_reasonable_sample
    Boolean featuremap_generate
    File featuremap_interval_list
    Int featuremap_scatter_count
    Int featuremap_min_mapq
    Int featuremap_snv_identical_bases
    Int featuremap_snv_identical_bases_after
    Int featuremap_min_score
    Int featuremap_limit_score
    String featuremap_extra_args
    File reference_gaps_intervals
    File centromere_intervals
    String dummy_input_for_call_caching
  }
  command {
    gsutil cp gs://fc-secure-5881a722-0e34-45a0-9eb7-eae48f047042/000_testOutputs/* .
  }
  output {
  
    File? output_gvcf = "004733-X0003.annotated.g.vcf.gz"
    File? output_gvcf_index = "004733-X0003.annotated.g.vcf.gz.tbi"
    File? orig_vcf = "unfiltered_vcf"
    File? orig_vcf_index = "unfiltered_vcf_index"
    File? output_vcf = "004733-X0003.vcf.gz"
    File? output_vcf_index = "004733-X0003.vcf.gz"
    File? sv_calls = "SVcalls.output_bcf"
    File? sv_calls_index = "SVcalls.output_bcf_index"

    File? no_gt_report_html = "CreateNoGTReport.no_gt_report_html"
    File? funcotator_vcf = "Funcotator.funcotator_output"
    File? funcotator_vcf_index = "Funcotator.funcotator_output_index"
    #MERGE bam file
    File? haplotype_bam = "004733-X0003.bam"
    File? haplotype_bam_index = "004733-X0003.bam.bai"

    Array[File]? output_bam = ["004733-X0003.bam"]
    Array[File]? output_bam_index = ["004733-X0003.bam.bai"]

    File? output_cram = "004733-X0003.cram"
    File? output_cram_index = "004733-X0003.cram.crai"
    File? output_cram_md5 = "004733-X0003.cram.md5"


    File? selfSM = "004733-X0003.selfSM"
    Float? contamination = 0.000152612

    File? comparison_output = "CompareToGroundTruthFilteredVCF.compare_h5"
    Array[File]? comparison_beds = []

    # VCF post-processing
    File? filtered_vcf = "004733-X0003.filtered.vcf.gz"
    File? filtered_vcf_index = "004733-X0003.filtered.vcf.gz.tbi"
    File? fingerprints_csv = "CrosscheckFingerprints_gt.output_vcf_fingerprints"
    File? model_h5 = "TrainModelWithGT.model_h5"
    File? model_pkl = "TrainModelWithGT.model_pkl"
    File? model_h5_no_gt = "TrainModelNoGT.model_h5"
    File? model_pkl_no_gt = "TrainModelNoGT.model_pkl"
    File? evaluate_report_h5 = "EvaluateResults.short_report_h5"
    File? short_report_html = "CreateReport.short_report_html"
    File? extended_report_html = "CreateReport.extended_report_html"
    File? short_report_h5 = "CreateReport.short_report_h5"
    File? extended_report_h5 = "CreateReport.extended_report_h5"
    File? varStats = "CreateReport.varStats"
    Array[File]? report_plots_png = []

    # STATISTIC COLLECTION
    File? quality_yield_metrics = "004733-X0003.unmapped.quality_yield_metrics"
    File? wgs_metrics = "004733-X0003.wgs_metrics"
    File? raw_wgs_metrics = "004733-X0003.raw_wgs_metrics"
    File? duplicate_metrics = "004733-X0003.duplicate_metrics"
    File? agg_alignment_summary_metrics = "004733-X0003.alignment_summary_metrics"
    File? agg_alignment_summary_pdf = "004733-X0003.read_length_histogram.pdf"
    File? agg_gc_bias_detail_metrics = "004733-X0003.gc_bias.detail_metrics"
    File? agg_gc_bias_pdf = "004733-X0003.gc_bias.pdf"
    File? agg_gc_bias_summary_metrics = "004733-X0003.gc_bias.summary_metrics"
    File? agg_quality_distribution_pdf = "004733-X0003.quality_distribution.pdf"
    File? agg_quality_distribution_metrics = "004733-X0003.quality_distribution_metrics"
    File? aggregated_metrics_h5 = "AggregateMetrics.aggregated_metrics_h5"
    File? aggregated_metrics_json = "ConvertAggregatedMetricsToJson.aggregated_metrics_json"
    Float? duplication_rate = 0.043528
    Float? chimerism_rate = 0.008962
    Boolean? is_outlier_data = false

    # COVERAGE COLLECTION

    Array[File]? coverage_depth_parquet_vhq = []
    Array[File]? coverage_boxplot_all = []
    Array[File]? coverage_profileplot_all = []
    String sample_name = "HG03476"
    String? gt_sample_name = "right_sample_name"
    String flow_order = "TGCA"
    String barcode = "CAACATACATCAGAT"
    String id = "X0003"

    # FEATUREMAP
    File? featuremap = "FeatureMapMerge.featuremap"
    File? featuremap_index = "FeatureMapMerge.featuremap_index"
    File? featuremap_single_substitutions = "FeatureMapSingleSubstitutionsMerge.featuremap"
    File? featuremap_single_substitutions_index = "FeatureMapSingleSubstitutionsMerge.featuremap_index"
    File? featuremap_dataframe = "FeatureMapMergeDataframes.featuremap_df"
    File? featuremap_single_substitutions_dataframe = "FeatureMapMergeSingleSubstitutionsDataframes.featuremap_df"

    # SNP Rate
    File? coverage_per_motif = "SNPRate.coverage_per_motif"
    File? snp_error_rate = "SNPRate.snp_error_rate"
    File? snp_error_rate_plot = "SNPRate.snp_error_rate_threshold5"
  }
  runtime {
    docker: "google/cloud-sdk:374.0.0"
    disks : "local-disk 500 HDD"
  }
}