version 1.0

import "../../../../../../tasks/broad/JukeboxTasks.wdl" as Tasks
import "../../../../../../tasks/broad/Utilities.wdl" as Utilities
import "../../../../../../tasks/broad/GermlineVariantDiscovery.wdl" as VariantDiscoverTasks
import "../../../../../../tasks/broad/Qc.wdl" as QC
import "../../../../../../tasks/broad/JukeboxAlignmentMarkDuplicates.wdl" as AlignmentAndMarkDuplicates
import "../../../../../../tasks/broad/JukeboxQC.wdl" as WorkflowQC
import "../../../../../../structs/dna_seq/JukeboxStructs" as Structs

# CHANGELOG
#  1.1.1     get multiple input cram
#  1.1.2     optional gvcf, copy report outputs to Nexus
#  2.1.0     gatk 5.0 (not for real, I forgot to update the docker)
#  2.1.1     gatk 5.0 for real
#  2.1.2     changed mark duplicates disk logic
#  2.1.3     changed CheckContamination docker and updated args (--adjust-MQ 0), updated default HaplotypeCaller flags to create bam files
#  2.1.4     fixed logic in H1->H2 retry mechanism, new gatk with bugfix
#  2.1.5     new gatk/picard docker for MarkDuplicatesSpark and CrosscheckFingerprints
#  2.1.5.1   light output naming refactoring
#  2.1.5.2   MarkDuplicatesSpark disk fix
#  2.1.5.3   reduced cpu number for MarkDuplicatesSpark
#  2.1.6     minor restructuring in preparation of gatk 0.5.1
#  2.1.7     replaced ultima markduplicates docker with broad one
#  2.2.0     gatk 0.5.1, blacklist in filtering
#  2.2.1     gatk 0.5.1-2 in markduplicates without localization
#  2.2.2     added CollectIntervalCoverageStats, removed maxretries from tasks without default preemptible tries
#  2.3.1     supports V4 of the basecalling pipeline, aggregates sequencing metrics
#  2.3.2     BIOIN-73 Flow for converted libraries
#            BIOIN-80 Contig filtering replaced with allele filtering in the GATK (much lower cost for high coverage runs)
#            BIOIN-72 Expanded variant calling report
#            BIOIN-48 MarkDuplicates for single end reads
#            BIOIN-69 SV calls
#  2.3.2.1   increased memory in CreateReport
#  2.4.1     [BIOIN-91, BIOIN-114] Improved report
#            BIOIN-130 Added new ground truth datasets
#            BIOIN-109 Variant filtering now outputs correct false positive rate
#            BIOIN-49 Flow based mark duplicates turned on by default with optimized parameters
#            BIOIN-112 SOR-based allele filtering in GATK
#  2.4.2     Variant calling report without ground truth [BIOIN-103]
#            GATK updated to 4.2.0 [BIOIN-139]
#            T/N support in Mutect [BIOIN-115]
#            GATK is about 40% faster [BIOIN-138]
#            Variant report outputs added to the aggregated metrics [BIOIN-23]
#            inputs/outputs are ingested into the database
#            Funcotator reports [BIOIN-84]
# 2.5.1      Major changes
#              [BIOIN-74] Contamination is calculated from realigned reads and flow based model rather than PairHMM
#              [BIOIN-97] FeatureMapper integrated into jukebox_vc pipeline
#              [BIOIN-107, BIOIN-156] Optimization of the coverage collection code
#              [BIOIN-142] Filtering alleles without ground truth
#              [BIOIN-151] Expanded test set
#            Minor changes
#              SB and AS_HmerLength are now reported in the VCF
#              Another error code to SW crash and optional SV calls
#              Barcode and read group ID are extracted and sent to metadata
#              Changed markDuplicates parameter (expected to decrease significantly duplicate marking)
#              Changed GATK parameters (to use adaptive pruning and dynamic read disqualification)
# 2.6.1      BIOIN-168 Compatible to BARC 4.2
# 2.6.2      Major changes
#              BIOIN-142 Improved allele filtering and performance on exome
#              BIOIN-146 UA is an optional aligner
#              BIOIN-162 SNP Error rate is calculated for all runs
#            Minor changes
#              Various bug fixes and stability improvements
#               * Fixes [BIOIN-185]
#               * Fixes [BIOIN-195]
#               * Fixes [BIOIN-198]
#               * Fixes [BIOIN-202]
#               * Remedy for [BIOIN-181]
#  2.6.3     Disabled featureMap by default
#            Added exome weights
#            Fixed output parameter
#  2.7.1     [BIOIN-153] GATK outputs isolated weak variants
#            GATK rebased to version 4.2.2.0
#            [BIOIN-224] Additional metrics: RLQ30, RLQ25, MEDIAN_READ_LENGTH output by picard for the FAT table
#            [BIOIN-127] Coverage statistics outptus bigwig files, coverage annotation in comparison scripts much faste
#            Stability improvements and bug fix in the somatic pipeline
#            [BIOIN-192] Metrics from no ground truth report added to the database
#            [BIOIN-169] Removed task to remove symbolic allele
#            Minor fixes:
#             [BIOIN-183] Better output of the report metadata
#             [BIOIN-152] Improved blacklist for somatic pipeline
#  2.7.2     Allow empty adapter in converted libraries
#  3.1.1     [BIOIN-189] contamination is calcuated on the FlowBased model
#            [BIOIN-257] Apply filtering model on the gVCF
#            [BIOIN-228] VCF contain flow-based annotations
#            [BIOIN-269] Improved filteirng with no overfitting and new training sets
#            Other minor fixes
#  3.2.1     [BIOIN-288] filter reads by rsq
#            FlowBasedHMM used by default as the error model in haplotypeCaller
#            Picard overrides in metrics collection
#            Contamination uses cycle-skip files
#            Annotations in the GATK from the allele filteirng
#  3.2.2     Fixed an error when the pipeline was running without the variant calling
#            [BIOIN-260] Rewired the flow of the VCFs in the pipeline so that the final output contains interval annotation and funcotator receives the correct vcf
#  3.2.3     Removed the rsq filtering
#            Switched error model back to FlowBased

workflow JukeboxSingleSample {
  input {
    String pipeline_version = "3.2.3"

    SampleInputs sample_inputs
    ContaminationSites contamination_sites
    RuntimeOptions runtime_options
    AlignmentReferences alignment_references
    VariantCallingSettings variant_calling_settings
    EnvironmentVersions environment_versions
    VcfPostProcessing vcf_post_processing
    ExtraArgs extra_args
  }
  References references = alignment_references.references

  Int preemptibles = select_first([runtime_options.preemptible_tries, 1])

  Int additional_disk = select_first([runtime_options.increase_disk_size, 20])

  # tasks that need to download locally the full bam/cram can randomly fail because a Google Cloud inconsistency, taking a disk > 500 GB will overcome the issue
  Float secure_disk_size_threshold = 510.0
  Int local_ssd_size = 375

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  # Germline single sample VCFs shouldn't get bigger even when the input bam is bigger (after a certain size)
  Int VCF_disk_size = select_first([runtime_options.increase_disk_size, 80])
  # Sometimes the output is larger than the input, or a task can spill to disk. In these cases we need to account for the
  # input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  Float bwa_disk_multiplier = 2.5

  Int compression_level = 2

  String bwa_commandline = "/usr/gitc/bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

  String gitc_path = select_first([environment_versions.gitc_path_override, "/usr/gitc/"])

  # Ensure no # charachters are found in base_file_name, MarkDuplicatesSpark can't handle it
  String base_file_name_sub = sub(sample_inputs.base_file_name, "#", "")

  # VCF post-processing default values
  Int eval_preemptible_tries = select_first([runtime_options.evaluation_preemptible_tries, 0])
  File? interval_list = vcf_post_processing.interval_list_override
  String filtering_model_no_gt_name = select_first([vcf_post_processing.filtering_model_no_gt_name_override, 'rf_model_ignore_gt_incl_hpol_runs'])
  Boolean no_address = select_first([runtime_options.no_address_override, true ])
  Boolean parallel_no_address = select_first([runtime_options.no_address_override, true ])

  Boolean make_gvcf = select_first([vcf_post_processing.make_gvcf_override, false ])
  Boolean merge_bam_file = select_first([vcf_post_processing.merge_bam_file_override, true ])

  Float ref_size = size(references.ref_fasta, "GB") + size(references.ref_fasta_index, "GB") + size(references.ref_dict, "GB")

  call AlignmentAndMarkDuplicates.AlignmentAndMarkDuplicates as AlignmentAndMarkDuplicates {
    input:
      sample_inputs = sample_inputs,
      base_file_name_sub = base_file_name_sub,
      reads_per_split = extra_args.reads_per_split,
      rsq_threshold = extra_args.rsq_threshold,
      crammer_docker = environment_versions.crammer_docker,
      compression_level = compression_level,
      gitc_docker = environment_versions.gitc_docker,
      gitc_path = gitc_path,
      gatk_markduplicates_docker = environment_versions.gatk_markduplicates_docker,
      jukebox_vc_docker = environment_versions.jukebox_vc_docker,
      no_address = no_address,
      parallel_no_address = parallel_no_address,
      dummy_input_for_call_caching = runtime_options.dummy_input_for_call_caching,
      preemptibles = preemptibles,
      additional_disk = additional_disk,
      bwa_commandline = bwa_commandline,
      alignment_references = alignment_references,
      mark_duplicates_extra_args = extra_args.mark_duplicates_extra_args,
      bwa_disk_multiplier = bwa_disk_multiplier,
      local_ssd_size = local_ssd_size,
      ref_size = ref_size,
      monitoring_script = runtime_options.monitoring_script
  }

  Float agg_bam_size = size(AlignmentAndMarkDuplicates.output_bam, "GB")

  Float dynamic_convert_to_cram_disk_size = (2 * agg_bam_size) + ref_size + additional_disk
  Float convert_to_cram_disk_size = if dynamic_convert_to_cram_disk_size > secure_disk_size_threshold then dynamic_convert_to_cram_disk_size else secure_disk_size_threshold

  # Convert the final merged recalibrated BAM file to CRAM format
  call Utilities.ConvertToCram {
    input:
      input_bam = AlignmentAndMarkDuplicates.output_bam,
      ref_fasta = references.ref_fasta,
      ref_fasta_index = references.ref_fasta_index,
      output_basename = sample_inputs.base_file_name,
      preemptible_tries = preemptibles
  }

  Float cram_size = size(ConvertToCram.output_cram, "GB")
  Float dynamic_validate_cram_disk_size = cram_size + ref_size + runtime_options.additional_metrics_disk
  Float validate_cram_disk_size = if dynamic_validate_cram_disk_size > secure_disk_size_threshold then dynamic_validate_cram_disk_size else secure_disk_size_threshold

  # Validate the CRAM file
  call QC.ValidateSamFile as ValidateCram {
    input:
      input_bam = ConvertToCram.output_cram,
      input_bam_index = ConvertToCram.output_cram_index,
      report_filename = sample_inputs.base_file_name + ".cram.validation_report",
      ref_dict = references.ref_dict,
      ref_fasta = references.ref_fasta,
      ref_fasta_index = references.ref_fasta_index,
      ignore = ["MISSING_TAG_NM" ,"INVALID_PLATFORM_VALUE"],
      max_output = 1000000000,
      preemptible_tries = 0,
      is_outlier_data = true #sets SKIP_MATE_VALIDATION=true
  }

  Float deduplicated_bam_size_vc = size(AlignmentAndMarkDuplicates.output_bam, "GB")

  call Tasks.ExtractSampleNameFlowOrder {
    input:
      input_bam = AlignmentAndMarkDuplicates.output_bam,
      references = references,
      preemptible_tries = preemptibles,
      monitoring_script = runtime_options.monitoring_script,
      docker = environment_versions.broad_gatk_docker
  }

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call Utilities.ScatterIntervalList {
    input:
      interval_list = variant_calling_settings.wgs_calling_interval_list,
      scatter_count = variant_calling_settings.haplotype_scatter_count,
      break_bands_at_multiples_of = variant_calling_settings.break_bands_at_multiples_of
  }

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by factor 2 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int hc_divisor = ScatterIntervalList.interval_count / 2
  #Int hc_divisor = 1

  # Call variants in parallel over WGS calling intervals
  scatter (index in range(ScatterIntervalList.interval_count)) {
    # Generate VCF by interval
    call Tasks.HaplotypeCaller as HC1 {
      input:
        input_bam_list = [AlignmentAndMarkDuplicates.output_bam],
        interval_list = ScatterIntervalList.out[index],
        vcf_basename = base_file_name_sub,
        references = references,
        # Divide the total output VCF size and the input bam size to account for the smaller scattered input and output.
        disk_size = ceil(((deduplicated_bam_size_vc + VCF_disk_size) / hc_divisor) + ref_size + additional_disk),
        preemptible_tries = preemptibles,
        docker = environment_versions.jb_gatk_docker,
        gitc_path = gitc_path,
        monitoring_script = runtime_options.monitoring_script,
        extra_args = extra_args.hc_extra_args,
        no_address = parallel_no_address,
        make_gvcf = make_gvcf,
        memory_gb = 12,
        make_bamout = false
     }

     # if cromwell implements optional outputs,
     Boolean h1_success = (size(HC1.output_vcf_index)!=0)
     if (!h1_success) {
        call Tasks.HaplotypeCaller as HC2 {
        input:
          input_bam_list = [AlignmentAndMarkDuplicates.output_bam],
          interval_list = ScatterIntervalList.out[index],
          vcf_basename = base_file_name_sub,
          references = references,
          # Divide the total output VCF size and the input bam size to account for the smaller scattered input and output.
          disk_size = ceil(((deduplicated_bam_size_vc + VCF_disk_size) / hc_divisor) + ref_size + additional_disk),
          preemptible_tries = preemptibles,
          docker = environment_versions.jb_gatk_docker,
          gitc_path = gitc_path,
          monitoring_script = runtime_options.monitoring_script,
          extra_args = extra_args.hc_extra_args,
          no_address = parallel_no_address,
          make_gvcf = make_gvcf,
          memory_gb = 40,
          native_sw = true,
          make_bamout = false
       }
     }

     # gate success on h1_success
     File HC_output_vcf_idx = select_first([ if h1_success then HC1.output_vcf_index else HC2.output_vcf_index])
     File HC_output_vcf =     select_first([ if h1_success then HC1.output_vcf else HC2.output_vcf])
     File monitoring_log =    select_first([ if h1_success then HC1.monitoring_log else HC2.monitoring_log])
     File haplotypes_bam =    select_first([ if h1_success then HC1.haplotypes_bam else HC2.haplotypes_bam])

  }

  # Combine by-interval VCFs into a single sample  file
  call VariantDiscoverTasks.MergeVCFs {
    input:
      input_vcfs = HC_output_vcf,
      input_vcfs_indexes = HC_output_vcf_idx,
      output_vcf_name = base_file_name_sub + (if make_gvcf then ".g.vcf.gz" else ".vcf.gz"),
      preemptible_tries = preemptibles
  }

  # Combine by-interval BAMs into a single sample file
  if(merge_bam_file) {
      call Tasks.MergeBams {
        input:
          monitoring_script = runtime_options.monitoring_script,
          input_bams = haplotypes_bam,
          output_bam_name = base_file_name_sub + ".bam",
          preemptible_tries = preemptibles,
          docker = environment_versions.gitc_docker,
          gitc_path = gitc_path,
          no_address = no_address
      }
  }

  Float vcf_size = size(MergeVCFs.output_vcf, "GB")

  call Tasks.ConvertGVCFtoVCF {
    input:
      monitoring_script = runtime_options.monitoring_script,
      input_gvcf = MergeVCFs.output_vcf,
      input_gvcf_index = MergeVCFs.output_vcf_index,
      output_vcf_name = base_file_name_sub + '.vcf.gz',
      references = references,
      disk_size_gb = VCF_disk_size,
      preemptible_tries = preemptibles,
      no_address = no_address,
      docker = environment_versions.jb_gatk_docker,
      gitc_path = gitc_path
  }

  # VCF post-processings
  call Tasks.AnnotateVCF {
    input :
      monitoring_script = runtime_options.monitoring_script,
      input_vcf = ConvertGVCFtoVCF.output_vcf,
      input_vcf_index = ConvertGVCFtoVCF.output_vcf_index,
      references = references,
      reference_dbsnp = vcf_post_processing.ref_dbsnp,
      reference_dbsnp_index = vcf_post_processing.ref_dbsnp_index,
      flow_order = ExtractSampleNameFlowOrder.flow_order,
      final_vcf_base_name = base_file_name_sub,
      additional_disk = additional_disk,
      preemptible_tries = preemptibles,
      docker = environment_versions.jb_gatk_docker,
      no_address = no_address,
      gitc_path = gitc_path
  }

  call Tasks.AddIntervalAnnotationsToVCF {
    input:
      input_vcf = AnnotateVCF.output_vcf_annotated,
      input_vcf_index = AnnotateVCF.output_vcf_annotated_index,
      final_vcf_base_name = base_file_name_sub,
      annotation_intervals = vcf_post_processing.annotation_intervals,
      preemptible_tries = preemptibles,
      docker = environment_versions.jukebox_vc_docker,
      monitoring_script = runtime_options.monitoring_script,
      no_address = no_address
  }

  call Tasks.TrainModel {
      input:
        input_file = AddIntervalAnnotationsToVCF.output_vcf,
        input_file_index = AddIntervalAnnotationsToVCF.output_vcf_index,
        input_vcf_name = base_file_name_sub,
        blacklist_file = vcf_post_processing.training_blacklist_file,
        ref_fasta = references.ref_fasta,
        ref_index = references.ref_fasta_index,
        runs_file = vcf_post_processing.runs_file,
        apply_model = filtering_model_no_gt_name,
        annotation_intervals = vcf_post_processing.annotation_intervals,
        exome_weight = vcf_post_processing.exome_weight,
        exome_weight_annotation = vcf_post_processing.exome_weight_annotation,
        additional_disk = additional_disk,
        preemptible_tries = eval_preemptible_tries,
        docker = environment_versions.jukebox_vc_docker,
        monitoring_script = runtime_options.monitoring_script,
        no_address = no_address
  }


  call Tasks.AnnotateVCF_AF {
      input :
        monitoring_script = runtime_options.monitoring_script,
        input_vcf = AddIntervalAnnotationsToVCF.output_vcf,
        input_vcf_index = AddIntervalAnnotationsToVCF.output_vcf_index,
        af_only_gnomad = vcf_post_processing.af_only_gnomad,
        af_only_gnomad_index = vcf_post_processing.af_only_gnomad_index,
        final_vcf_base_name = base_file_name_sub,
        additional_disk = additional_disk,
        preemptible_tries = preemptibles,
        docker = environment_versions.jukebox_vc_docker,
        no_address = no_address
  }

  call Tasks.FilterVCF {
    input:
      input_vcf = AnnotateVCF_AF.output_vcf_annotated,
      input_model = select_first([vcf_post_processing.filtering_model_no_gt,TrainModel.model_pkl]),
      runs_file = vcf_post_processing.runs_file,
      references = references,
      model_name = filtering_model_no_gt_name,
      filter_cg_insertions = vcf_post_processing.filter_cg_insertions,
      blacklist_file = vcf_post_processing.filtering_blacklist_file,
      final_vcf_base_name = base_file_name_sub,
      flow_order = ExtractSampleNameFlowOrder.flow_order,
      annotation_intervals = vcf_post_processing.annotation_intervals,
      disk_size_gb = VCF_disk_size,
      preemptible_tries = eval_preemptible_tries,
      docker = environment_versions.jukebox_vc_docker,
      monitoring_script = runtime_options.monitoring_script,
      no_address = no_address
  }

  call Tasks.MoveAnnotationsToGvcf {
    input:
      filtered_vcf = FilterVCF.output_vcf_filtered,
      filtered_vcf_index = FilterVCF.output_vcf_filtered_index,
      gvcf = MergeVCFs.output_vcf,
      gvcf_index = MergeVCFs.output_vcf_index
  }

  call WorkflowQC.QC as CollectStatistics {
    input:
      agg_bam = AlignmentAndMarkDuplicates.output_bam,
      agg_bam_index = AlignmentAndMarkDuplicates.output_bam_index,
      base_file_name = sample_inputs.base_file_name,
      base_file_name_sub = base_file_name_sub,
      agg_bam_size = agg_bam_size,
      ref_size = ref_size,
      additional_metrics_disk = runtime_options.additional_metrics_disk,
      secure_disk_size_threshold = secure_disk_size_threshold,
      references = references,
      contamination_sites = contamination_sites,
      wgs_coverage_interval_list = vcf_post_processing.wgs_coverage_interval_list,
      picard_jar_override = environment_versions.picard_jar_override,
      jb_gatk_docker = environment_versions.jb_gatk_docker,
      gatk_markduplicates_docker = environment_versions.gatk_markduplicates_docker,
      gitc_path = gitc_path,
      gitc_docker = environment_versions.gitc_docker,
      no_address = no_address,
      preemptibles = preemptibles,
      eval_preemptible_tries = eval_preemptible_tries,
      monitoring_script = runtime_options.monitoring_script,
      VCF_disk_size = VCF_disk_size,
      additional_disk = additional_disk,
      max_duplication_in_reasonable_sample = vcf_post_processing.max_duplication_in_reasonable_sample,
      max_chimerism_in_reasonable_sample = vcf_post_processing.max_chimerism_in_reasonable_sample,
      flow_order = ExtractSampleNameFlowOrder.flow_order
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_gvcf = MoveAnnotationsToGvcf.output_gvcf
    File output_gvcf_index = MoveAnnotationsToGvcf.output_gvcf_index
    File output_vcf = ConvertGVCFtoVCF.output_vcf
    File output_vcf_index = ConvertGVCFtoVCF.output_vcf

    #MERGE bam file
    File? haplotype_bam = MergeBams.output_bam
    File? haplotype_bam_index = MergeBams.output_bam_index

    File output_cram = ConvertToCram.output_cram
    File output_cram_index = ConvertToCram.output_cram_index
    File output_cram_md5 = ConvertToCram.output_cram_md5

    File selfSM = CollectStatistics.selfSM
    Float contamination = CollectStatistics.contamination

    # VCF post-processing
    File filtered_vcf = FilterVCF.output_vcf_filtered
    File filtered_vcf_index = FilterVCF.output_vcf_filtered_index

    # STATISTIC COLLECTION
    File quality_yield_metrics = CollectStatistics.quality_yield_metrics
    File wgs_metrics = CollectStatistics.wgs_metrics
    File raw_wgs_metrics = CollectStatistics.raw_wgs_metrics
    File duplicate_metrics = CollectStatistics.duplicate_metrics
    File agg_alignment_summary_metrics = CollectStatistics.agg_alignment_summary_metrics
    File? agg_alignment_summary_pdf = CollectStatistics.agg_alignment_summary_pdf
    File agg_gc_bias_detail_metrics = CollectStatistics.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = CollectStatistics.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = CollectStatistics.agg_gc_bias_summary_metrics
    File agg_quality_distribution_pdf = CollectStatistics.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = CollectStatistics.agg_quality_distribution_metrics
    Float duplication_rate = CollectStatistics.duplication_rate
    Float chimerism_rate = CollectStatistics.chimerism_rate
    Boolean is_outlier_data = CollectStatistics.is_outlier_data

    String sample_name = ExtractSampleNameFlowOrder.sample_name
    String flow_order = ExtractSampleNameFlowOrder.flow_order
    String barcode = ExtractSampleNameFlowOrder.barcode_seq
    String id = ExtractSampleNameFlowOrder.readgroup_id
  }

}