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

  call MergeMetrics {
    input:
      wgs_metrics = select_first([JukeboxVC.wgs_metrics]),
      alignment_summary_metrics = select_first([JukeboxVC.agg_alignment_summary_metrics]),
      gc_bias_summary_metrics = select_first([JukeboxVC.agg_gc_bias_summary_metrics]),
      duplicate_metrics = select_first([JukeboxVC.duplicate_metrics]),
      quality_yield_metrics = select_first([JukeboxVC.quality_yield_metrics]),
      raw_wgs_metrics = select_first([JukeboxVC.raw_wgs_metrics]),
      contamination_metrics = select_first([JukeboxVC.selfSM]),
      fingerprint_summary_metrics = CheckFingerprint.fingerprint_summary_metrics_file,
      output_basename = collab_sample_id_run_id
  }

  # CALL TDR TASKS TO FORMAT JSON
  
  if (defined(tdr_dataset_uuid) && defined(tdr_sample_id) && defined(tdr_staging_bucket)) {
    call formatPipelineOutputs {
      input:
        output_basename = output_basename,
        collab_sample_id_run_id = collab_sample_id_run_id,
        output_gvcf = JukeboxVC.output_gvcf,
        output_gvcf_index = JukeboxVC.output_gvcf_index,
        output_vcf = JukeboxVC.output_vcf,
        output_vcf_index = JukeboxVC.output_vcf_index,
        haplotype_bam = JukeboxVC.haplotype_bam,
        haplotype_bam_index = JukeboxVC.haplotype_bam_index,
        output_cram = JukeboxVC.output_cram,
        output_cram_index = JukeboxVC.output_cram_index,
        output_cram_md5 = JukeboxVC.output_cram_md5,
        selfSM = JukeboxVC.selfSM,
        contamination = JukeboxVC.contamination,
        filtered_vcf = JukeboxVC.filtered_vcf,
        filtered_vcf_index = JukeboxVC.filtered_vcf_index,
        quality_yield_metrics = JukeboxVC.quality_yield_metrics,
        wgs_metrics = JukeboxVC.wgs_metrics,
        raw_wgs_metrics = JukeboxVC.raw_wgs_metrics,
        duplicate_metrics = JukeboxVC.duplicate_metrics,
        agg_alignment_summary_metrics = JukeboxVC.agg_alignment_summary_metrics,
        agg_alignment_summary_pdf = JukeboxVC.agg_alignment_summary_pdf,
        agg_gc_bias_detail_metrics = JukeboxVC.agg_gc_bias_detail_metrics,
        agg_gc_bias_pdf = JukeboxVC.agg_gc_bias_pdf,
        agg_gc_bias_summary_metrics = JukeboxVC.agg_gc_bias_summary_metrics,
        agg_quality_distribution_pdf = JukeboxVC.agg_quality_distribution_pdf,
        agg_quality_distribution_metrics = JukeboxVC.agg_quality_distribution_metrics,
        duplication_rate = JukeboxVC.duplication_rate,
        chimerism_rate = JukeboxVC.chimerism_rate,
        is_outlier_data = JukeboxVC.is_outlier_data,
        sample_name = JukeboxVC.sample_name,
        flow_order = JukeboxVC.flow_order,
        barcode = JukeboxVC.barcode,
        id = JukeboxVC.id,
        unified_metrics = MergeMetrics.unified_metrics,
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

    File output_gvcf = JukeboxVC.output_gvcf
    File output_gvcf_index = JukeboxVC.output_gvcf_index
    File output_vcf = JukeboxVC.output_vcf
    File output_vcf_index = JukeboxVC.output_vcf_index

    File? haplotype_bam = JukeboxVC.haplotype_bam
    File? haplotype_bam_index = JukeboxVC.haplotype_bam_index

    File output_cram = JukeboxVC.output_cram
    File output_cram_index = JukeboxVC.output_cram_index
    File output_cram_md5 = JukeboxVC.output_cram_md5

    File selfSM = JukeboxVC.selfSM
    Float contamination = JukeboxVC.contamination

    File filtered_vcf = JukeboxVC.filtered_vcf
    File filtered_vcf_index = JukeboxVC.filtered_vcf_index

    File quality_yield_metrics = JukeboxVC.quality_yield_metrics
    File wgs_metrics = JukeboxVC.wgs_metrics
    File raw_wgs_metrics = JukeboxVC.raw_wgs_metrics
    File duplicate_metrics = JukeboxVC.duplicate_metrics
    File agg_alignment_summary_metrics = JukeboxVC.agg_alignment_summary_metrics
    File? agg_alignment_summary_pdf = JukeboxVC.agg_alignment_summary_pdf
    File agg_gc_bias_detail_metrics = JukeboxVC.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = JukeboxVC.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = JukeboxVC.agg_gc_bias_summary_metrics
    File agg_quality_distribution_pdf = JukeboxVC.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = JukeboxVC.agg_quality_distribution_metrics
    Float duplication_rate = JukeboxVC.duplication_rate
    Float chimerism_rate = JukeboxVC.chimerism_rate
    Boolean is_outlier_data = JukeboxVC.is_outlier_data

    String sample_name = JukeboxVC.sample_name
    String flow_order = JukeboxVC.flow_order
    String barcode = JukeboxVC.barcode
    String id = JukeboxVC.id
  }
}

task MergeMetrics {
  input {
    File wgs_metrics
    File alignment_summary_metrics
    File gc_bias_summary_metrics
    File duplicate_metrics
    File quality_yield_metrics
    File raw_wgs_metrics
    File contamination_metrics
    File? fingerprint_summary_metrics
    String output_basename

    String docker =  "python:3.8-slim"
    Int cpu = 1
    Int memory_mb = 3000
    Int disk_size_gb = 10
  }

  String out_filename = output_basename + ".unified_metrics.txt"

  command <<<

    #
    # Script transpose a two line TSV
    #
    cat <<-'EOF' > transpose.py
    import csv, sys

    rows = list(csv.reader(sys.stdin, delimiter='\t'))

    for col in range(0, len(rows[0])):
      key = rows[0][col].lower()
      print(f"{key}\t{rows[1][col]}")
    EOF

    #
    # Script clean the keys, replacing space, dash and forward-slash with underscores,
    # and removing comma, single quote, periods, and #
    #
    cat <<-'EOF' > clean.py
    import sys

    for line in sys.stdin:
      (k,v) = line.strip().split("\t")
      transtable = k.maketrans({' ':'_', '-':'_', '/':'_', ',':None, '\'':None, '.':None, '#':None})
      print(f"{k.translate(transtable)}\t{v}")
    EOF

    # Process each metric file, transposing and cleaning if necessary, and pre-pending a source to the metric name

    echo "Processing WGS Metrics"
    cat ~{wgs_metrics} | grep -A 1 "GENOME_TERRITORY" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "wgs_metrics_" $0}' >> ~{out_filename}

    echo "Processing Alignment Summary Metrics"
    cat ~{alignment_summary_metrics} | grep -A 1 "CATEGORY" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "alignment_summary_metrics_" $0}' >> ~{out_filename}

    echo "Processing GCBias Summary Metrics"
    cat ~{gc_bias_summary_metrics} | grep -A 1 "AT_DROPOUT" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "gc_bias_summary_metrics_" $0}' >> ~{out_filename}

    echo "Processing Duplicate Metrics: WARNING this only works for samples from one library"
    cat ~{duplicate_metrics} | grep -A 1 "UNPAIRED_READ_DUPLICATES" | python transpose.py | awk '{print "duplicate_metrics_" $0}' >> ~{out_filename}

    echo "Processing Quality Yield Metrics"
    cat ~{quality_yield_metrics} | grep -A 1 "TOTAL_READS" | python transpose.py | awk '{print "quality_yield_metrics_" $0}' >> ~{out_filename}

    echo "Processing Raw WGS Metrics"
    cat ~{raw_wgs_metrics} | grep -A 1 "GENOME_TERRITORY" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "raw_wgs_metrics_" $0}' >> ~{out_filename}

    echo "Processing Contamination Metrics"
    cat ~{contamination_metrics} | python transpose.py | python clean.py | awk '{print "contamination_metrics_" $0}' >> ~{out_filename}

    if [[ -f "~{fingerprint_summary_metrics}" ]];
    then
      echo "Processing Fingerprint Summary Metrics - only extracting LOD_EXPECTED_SAMPLE"
      cat ~{fingerprint_summary_metrics} | grep -A 1 "LOD_EXPECTED_SAMPLE" | python transpose.py | grep -i "LOD_EXPECTED_SAMPLE" | awk '{print "fp_"$0}' >> ~{out_filename}
    else
      echo "No Fingerprint Summary Metrics found."
      echo "fp_lod_expected_sample	" >> ~{out_filename}
    fi    >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File unified_metrics = out_filename
  }
}

# DEFINE TDR-SPECIFIC TASKS
  
  task formatPipelineOutputs {
    input {
      String output_basename
      String collab_sample_id_run_id
      
      String output_gvcf
      String output_gvcf_index
      String output_vcf
      String output_vcf_index

      String? haplotype_bam = ""
      String? haplotype_bam_index = ""

      String output_cram
      String output_cram_index
      String output_cram_md5

      String selfSM
      Float contamination

      String filtered_vcf
      String filtered_vcf_index

      String quality_yield_metrics
      String wgs_metrics
      String raw_wgs_metrics
      String duplicate_metrics
      String agg_alignment_summary_metrics
      String? agg_alignment_summary_pdf = ""
      String agg_gc_bias_detail_metrics
      String agg_gc_bias_pdf
      String agg_gc_bias_summary_metrics
      String agg_quality_distribution_pdf
      String agg_quality_distribution_metrics
      Float duplication_rate
      Float chimerism_rate
      Boolean is_outlier_data

      String sample_name
      String flow_order
      String barcode
      String id

      String? fingerprint_detail_metrics_file = ""
      String? fingerprint_summary_metrics_file = ""

      File unified_metrics

      Int cpu = 1
      Int memory_mb = 2000
      Int disk_size_gb = 10
    }

    String outputs_json_file_name = "outputs_to_TDR_~{output_basename}.json"

    String outlier_data_string = if is_outlier_data then "True" else "False"

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
          outputs_dict["fingerprint_summary_metrics_file"]="~{fingerprint_summary_metrics_file}"
          outputs_dict["fingerprint_detail_metrics_file"]="~{fingerprint_detail_metrics_file}"

          # explode unified metrics file
          with open("~{unified_metrics}", "r") as infile:
            for row in infile:
              key, value = row.rstrip("\n").split("\t")
              if value == "NA" or value == "" or value == "?" or value == "-":
                outputs_dict[key] = None
              else:
                outputs_dict[key] = value

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

    File output_gvcf ="004733-X0003.annotated.g.vcf.gz"
    File output_gvcf_index = "004733-X0003.annotated.g.vcf.gz.tbi"
    File output_vcf = "004733-X0003.vcf.gz"
    File output_vcf_index = "004733-X0003.vcf.gz"

    File? haplotype_bam = "004733-X0003.bam"
    File? haplotype_bam_index = "004733-X0003.bam.bai"

    File output_cram = "004733-X0003.cram"
    File output_cram_index = "004733-X0003.cram.crai"
    File output_cram_md5 = "004733-X0003.cram.md5"

    File selfSM = "004733-X0003.selfSM"
    Float contamination = 0.000152612

    File filtered_vcf = "004733-X0003.filtered.vcf.gz"
    File filtered_vcf_index = "004733-X0003.filtered.vcf.gz.tbi"

    # STATISTIC COLLECTION
    File quality_yield_metrics = "downsampled_NA12878.unmapped.quality_yield_metrics"
    File wgs_metrics = "downsampled_NA12878.wgs_metrics"
    File raw_wgs_metrics = "downsampled_NA12878.raw_wgs_metrics"
    File duplicate_metrics = "004733-X0003.duplicate_metrics"
    File agg_alignment_summary_metrics = "downsampled_NA12878.alignment_summary_metrics"
    File? agg_alignment_summary_pdf = "004733-X0003.read_length_histogram.pdf"
    File agg_gc_bias_detail_metrics = "004733-X0003.gc_bias.detail_metrics"
    File agg_gc_bias_pdf = "004733-X0003.gc_bias.pdf"
    File agg_gc_bias_summary_metrics = "downsampled_NA12878.gc_bias.summary_metrics"
    File agg_quality_distribution_pdf = "004733-X0003.quality_distribution.pdf"
    File agg_quality_distribution_metrics = "downsampled_NA12878.quality_distribution_metrics"
    Float duplication_rate = 0.043528
    Float chimerism_rate = 0.008962
    Boolean is_outlier_data = false

    String sample_name = "HG03476"
    String flow_order = "TGCA"
    String barcode = "CAACATACATCAGAT"
    String id = "X0003"
  }
  runtime {
    docker: "google/cloud-sdk:374.0.0"
    disks : "local-disk 500 HDD"
  }
}