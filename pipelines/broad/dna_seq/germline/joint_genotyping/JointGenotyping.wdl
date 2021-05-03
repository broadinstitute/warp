version 1.0

import "../../../../../tasks/broad/JointGenotypingTasks.wdl" as Tasks

workflow JointGenotyping {

  String pipeline_version = "2.0.0"

  input {
    String sample_name_map
    String? override_vcf_header_file
    File yaml_file
    String? override_hail_docker
    String? override_region

    String cluster_name
    String output_dir
    String callset_name
    String google_project
    String autoscaling_policy_name
    String? mt_to_vcf_parallel_arg
    String data_type
    Boolean overwrite_mt
    Boolean? gather_vcfs

    File haplotype_database


    # gnarly inputs
    String unpadded_intervals_file

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf
    File dbsnp_vcf_index

    Int small_disk
    Int medium_disk
    Int large_disk
    Int huge_disk

    Int snps_variant_recalibration_threshold

    Array[String] snp_recalibration_tranche_values
    Array[String] snp_recalibration_annotation_values
    Array[String] indel_recalibration_tranche_values
    Array[String] indel_recalibration_annotation_values

    File eval_interval_list
    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
    File mills_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf = dbsnp_vcf
    File dbsnp_resource_vcf_index = dbsnp_vcf_index

    # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
    # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
    Float excess_het_threshold = 54.69
    Float snp_filter_level
    Float indel_filter_level
    Int SNP_VQSR_downsampleFactor
    Int indel_VQSR_downsampleFactor
    Boolean use_allele_specific_annotations = true

    Int vcf_count
  }

  String hail_docker = select_first([override_hail_docker, "ruchim/hail-testing:57f5df4"])
  String region = select_first([override_region, "us-central1"])
  String gvcf_header = select_first([override_vcf_header_file, "gs://jg-dev-hail-testing/withWdl/header.g.vcf.gz"])

#  if (!endsWith(output_dir, "/")) {
#    output_dir = output_dir + "/"
#  }

  call StartCluster {
    input:
      google_project = google_project,
      cluster_name = cluster_name,
      region = region,
      hail_docker = hail_docker,
      yaml_file = yaml_file,
      autoscaling_policy_name = autoscaling_policy_name,
      data_type = data_type,
      sample_name_map = sample_name_map
  }

  call CreateMatrixTable {
    input:
      google_project = google_project,
      cluster_name  = StartCluster.created_cluster_name,
      region = region,
      hail_docker = hail_docker,
      sample_name_map = sample_name_map,
      output_dir = output_dir,
      callset_name = callset_name,
      overwrite_mt = overwrite_mt,
      gvcf_header_file = gvcf_header
  }

  String parallel = select_first([mt_to_vcf_parallel_arg, ""])
  call ConvertMatrixTableToVcf {
    input:
      matrix_table = CreateMatrixTable.matrix_table,
      hail_docker = hail_docker,
      output_dir = output_dir,
      callset_name = callset_name,
      parallel = parallel,
      cluster_name = cluster_name,
      region = region,
      google_project = google_project,
      data_type = data_type
  }

  call CreateExtractFingereprintSiteVcf {
    input:
      google_project = google_project,
      cluster_name = cluster_name,
      region = region,
      hail_docker = hail_docker,
      output_dir = output_dir,
      matrix_table = CreateMatrixTable.matrix_table,
      output_fingerprint_vcf_name = callset_name + '.fingerprint.vcf.bgz'
  }


  # Make a 2.5:1 interval number to samples in callset ratio interval list.
  # We allow overriding the behavior by specifying the desired number of vcfs
  # to scatter over for testing / special requests.
  # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
  # exome scatterCountPerSample is 0.05, min scatter 10, max 1000
  Int unboundedScatterCount = vcf_count
  Int scatterCount = if unboundedScatterCount > 10 then unboundedScatterCount else 10 #I think weird things happen if scatterCount is 1 -- IntervalListTools is noop?
  call Tasks.SplitIntervalList {
    input:
      interval_list = unpadded_intervals_file,
      scatter_count = scatterCount,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      disk_size = small_disk,
      sample_names_unique_done = true
  }

  Array[File] unpadded_intervals = SplitIntervalList.output_intervals

  scatter (idx in range(length(unpadded_intervals))) {
    call GnarlyGenotyperOnVcf {
      input:
        combined_gvcf = ConvertMatrixTableToVcf.output_vcf,
        combined_gvcf_index = ConvertMatrixTableToVcf.output_vcf_index, #TabixBGzippedFile.tabix_output,
        interval = unpadded_intervals[idx],
        output_vcf_filename = callset_name + "." + idx + ".vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf
    }

    call Tasks.HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = GnarlyGenotyperOnVcf.output_vcf,
        vcf_index = GnarlyGenotyperOnVcf.output_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        disk_size = medium_disk
    }
  }

  call Tasks.GatherVcfs as SitesOnlyGatherVcf {
   input:
     input_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
     output_vcf_name = callset_name + ".sites_only.vcf.gz",
     disk_size = medium_disk
  }

  call Tasks.IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      use_allele_specific_annotations = use_allele_specific_annotations,
      disk_size = small_disk
  }

  if (num_gvcfs > snps_variant_recalibration_threshold) {
    call Tasks.SNPsVariantRecalibratorCreateModel {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = callset_name + ".snps.recal",
        tranches_filename = callset_name + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        downsampleFactor = SNP_VQSR_downsampleFactor,
        model_report_filename = callset_name + ".snps.model.report",
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        disk_size = small_disk,
        use_allele_specific_annotations = use_allele_specific_annotations
    }

    scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.sites_only_vcf))) {
      call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
        input:
          sites_only_variant_filtered_vcf = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[idx],
          sites_only_variant_filtered_vcf_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[idx],
          recalibration_filename = callset_name + ".snps." + idx + ".recal",
          tranches_filename = callset_name + ".snps." + idx + ".tranches",
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotation_values,
          model_report = SNPsVariantRecalibratorCreateModel.model_report,
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_resource_vcf,
          dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
          disk_size = small_disk,
          machine_mem_gb = 60,
          use_allele_specific_annotations = use_allele_specific_annotations
        }
      }

    call Tasks.GatherTranches as SNPGatherTranches {
      input:
        tranches = SNPsVariantRecalibratorScattered.tranches,
        output_filename = callset_name + ".snps.gathered.tranches",
        disk_size = small_disk,
        mode = "SNP"
      }
  }

  if (num_gvcfs <= snps_variant_recalibration_threshold) {
    call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
      input:
        sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
        sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
        recalibration_filename = callset_name + ".snps.recal",
        tranches_filename = callset_name + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        use_allele_specific_annotations = use_allele_specific_annotations,
        disk_size = small_disk
    }
  }

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    call Tasks.ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = if defined(SNPsVariantRecalibratorScattered.recalibration) then select_first([SNPsVariantRecalibratorScattered.recalibration])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration]),
        snps_recalibration_index = if defined(SNPsVariantRecalibratorScattered.recalibration_index) then select_first([SNPsVariantRecalibratorScattered.recalibration_index])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration_index]),
        snps_tranches = select_first([SNPGatherTranches.tranches_file, SNPsVariantRecalibratorClassic.tranches]),
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        disk_size = medium_disk,
        use_allele_specific_annotations = true
    }

    # for large callsets we need to collect metrics from the shards and gather them later
    if (!is_small_callset) {
      call Tasks.CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = ApplyRecalibration.recalibrated_vcf,
          input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          metrics_filename_prefix = callset_name + "." + idx,
          dbsnp_vcf = dbsnp_vcf,
          dbsnp_vcf_index = dbsnp_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict,
          disk_size = small_disk
      }
    }
 }

 if (!is_small_callset) {
   call Tasks.GatherVariantCallingMetrics {
    input:
      input_details = select_all(CollectMetricsSharded.detail_metrics_file),
      input_summaries = select_all(CollectMetricsSharded.summary_metrics_file),
      output_prefix = callset_name,
      disk_size = medium_disk
    }
 }



 #TODO: sample_name_map is a String. How can we read this and get line count?
 Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
 Int num_gvcfs = length(sample_name_map_lines)
 Boolean is_small_callset = select_first([gather_vcfs, num_gvcfs <= 1000])

 # for small callsets we can gather the VCF shards and then collect metrics on it
 if (is_small_callset) {

   call Tasks.GatherVcfs as FinalGatherVcf {
     input:
       input_vcfs = ApplyRecalibration.recalibrated_vcf,
       output_vcf_name = callset_name + ".vcf.gz",
       disk_size = huge_disk
   }

   call Tasks.CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
     input:
       input_vcf = FinalGatherVcf.output_vcf,
       input_vcf_index = FinalGatherVcf.output_vcf_index,
       metrics_filename_prefix = callset_name,
       dbsnp_vcf = dbsnp_vcf,
       dbsnp_vcf_index = dbsnp_vcf_index,
       interval_list = eval_interval_list,
       ref_dict = ref_dict,
       disk_size = huge_disk
   }
 }

#   call GatherMetrics {
#        input:
#         input_details_fofn = write_lines(select_all(CollectMetricsSharded.detail_metrics_file)),
#         input_summaries_fofn = write_lines(select_all(CollectMetricsSharded.summary_metrics_file)),
#         output_prefix = callset_name,
#         disk_size = medium_disk
#    }
#

  call Tasks.PartitionSampleNameMap {
    input:
      sample_name_map = sample_name_map,
      line_limit = 1000
  }

  Array[File] extract_fingerprint_vcf_list = [CreateExtractFingereprintSiteVcf.fingerprinting_vcf]

  scatter (idx in range(length(PartitionSampleNameMap.partitions))) {

    Array[File] files_in_partition = read_lines(PartitionSampleNameMap.partitions[idx])

    call Tasks.CrossCheckFingerprint as CrossCheckFingerprintsScattered {
      input:
        gvcf_paths = files_in_partition,
        vcf_paths = extract_fingerprint_vcf_list,
        sample_name_map = sample_name_map,
        haplotype_database = haplotype_database,
        output_base_name = callset_name + "." + idx,
        scattered = true
    }
  }

  call Tasks.GatherPicardMetrics as GatherFingerprintingMetrics {
    input:
      metrics_files = CrossCheckFingerprintsScattered.crosscheck_metrics,
      output_file_name = callset_name + ".fingerprintcheck",
      disk_size = small_disk
  }

  Array[File?] output_vcf_files = if defined(FinalGatherVcf.output_vcf) then [FinalGatherVcf.output_vcf] else ApplyRecalibration.recalibrated_vcf
  Array[File?] output_vcf_index_files = if defined(FinalGatherVcf.output_vcf_index) then [FinalGatherVcf.output_vcf_index] else ApplyRecalibration.recalibrated_vcf_index

  File output_detail_metrics_file = select_first([CollectMetricsOnFullVcf.detail_metrics_file, GatherVariantCallingMetrics.detail_metrics_file])
  File output_summary_metrics_file = select_first([CollectMetricsOnFullVcf.summary_metrics_file, GatherVariantCallingMetrics.summary_metrics_file])

  output {
    # Hail output IS the annotation DB so we don't need that

    # Metrics from either the small or large callset
    #File detail_metrics_file = GatherMetrics.detail_metrics_file
    #File summary_metrics_file = GatherMetrics.summary_metrics_file

    # Outputs from the small callset path through the wdl.
    Array[File] output_vcfs = select_all(output_vcf_files)
    Array[File] output_vcf_indices = select_all(output_vcf_index_files)

    File detail_metrics_file = output_detail_metrics_file
    File summary_metrics_file = output_summary_metrics_file
    Array[File] output_intervals = SplitIntervalList.output_intervals

    File crosscheck_metrics = GatherFingerprintingMetrics.gathered_metrics
  }
}


task CreateExtractFingereprintSiteVcf {
  input {
    String google_project
    String cluster_name
    String region
    String hail_docker
    String matrix_table
    String output_dir
    String output_fingerprint_vcf_name
  }

  String output_fingerprint_vcf_bgz_path = output_dir + output_fingerprint_vcf_name
  String output_fingerprint_vcf_gz_path = sub(output_fingerprint_vcf_bgz_path, ".vcf.bgz", ".vcf.gz")

  command {
    set -e

    # Set region
    gcloud config set dataproc/region ~{region}

    /submit_extract_fingerprint_sites.sh ~{google_project} ~{cluster_name} ~{matrix_table} \
    ~{output_fingerprint_vcf_bgz_path}

    gsutil mv ~{output_fingerprint_vcf_bgz_path} ~{output_fingerprint_vcf_gz_path}
  }

  runtime {
    docker: hail_docker
  }

  output {
    String fingerprinting_vcf = "~{output_fingerprint_vcf_gz_path}"
  }
}

task StartCluster {
  input {
    String google_project
    String cluster_name
    String region
    String hail_docker
    File yaml_file
    String autoscaling_policy_name
    String data_type
    String sample_name_map
  }

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_samples = length(sample_name_map_lines)
  Boolean large_exome_callset = if (num_samples >= 10,000 && data_type == 'Exome')
  Boolean large_wgs_callset = if (num_samples >= 10,000 && data_type == 'WGS')
  String size_based_yaml = (if (large_exome_callset || large_wgs_callset ) then "largeCallSet.yaml" else "smallCallSet.yaml")
  File yaml_file = yaml_file + size_based_yaml

  meta {
    volatile: true
  }

  command {
    set -e

    # Check if cluster is already running
    ClusterExists=$(gcloud dataproc clusters list --project='~{google_project}' --region=~{region} \
    | cut -f 1 -d ' ' | grep '^~{cluster_name}$') || true

    if [[ -z "$ClusterExists" ]]
    then

      # import yaml
      gcloud --project=~{google_project} --quiet beta dataproc autoscaling-policies import \
      ~{autoscaling_policy_name} --region=~{region} --source=~{default="/test_cluster.yaml" yaml_file}

      # start cluster
      hailctl dataproc --beta start ~{cluster_name} --project=~{google_project} --quiet --region=~{region} \
      --max-idle 1h --autoscaling-policy=~{autoscaling_policy_name} --packages gnomad
    fi
  }

  runtime {
    docker: hail_docker
  }

  output {
    String created_cluster_name = cluster_name
  }

}

task CreateMatrixTable {
  input {
    String google_project
    String cluster_name
    String region
    String hail_docker
    String sample_name_map
    String output_dir
    String callset_name
    String gvcf_header_file
    Boolean overwrite_mt
  }

  command {
    # Set region
    gcloud config set dataproc/region ~{region}

    # Submit hail workflow
    /test_submit.sh ~{sample_name_map} ~{cluster_name} ~{output_dir} ~{callset_name} ~{google_project} ~{gvcf_header_file} ~{overwrite_mt}
  }

  runtime {
    docker: hail_docker
  }

  output {
    String matrix_table = "~{output_dir}~{callset_name}.mt"
  }
}

task ConvertMatrixTableToVcf {
  input {
    String matrix_table
    String hail_docker
    String output_dir
    String callset_name
    String parallel
    String google_project
    String cluster_name
    String region
    String data_type
  }

  command {
    set -e

    # Set region
    gcloud config set dataproc/region ~{region}

    /submit_mt_to_vcf_to_cluster.sh ~{google_project} ~{cluster_name} ~{matrix_table} \
    ~{output_dir} ~{data_type} ~{callset_name} ~{parallel}
  }

  runtime {
    docker: hail_docker
  }

  output {
    String output_vcf = "~{output_dir}~{callset_name}.sites.vcf.bgz"
    String output_vcf_index = "~{output_dir}~{callset_name}.sites.vcf.bgz.tbi"
  }
}


task GnarlyGenotyperOnVcf {

  input {
    File combined_gvcf
    File combined_gvcf_index
    File interval
    String output_vcf_filename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String dbsnp_vcf

    String gatk_docker = "gcr.io/broad-dsde-methods/gnarly_genotyper:hail_ukbb_300K"
  }

  parameter_meta {
    interval: {
      localization_optional: true
    }
    combined_gvcf: {
      localization_optional: true
    }
    combined_gvcf_index: {
      localization_optional: true
    }
  }

  Int disk_size = ceil(size(combined_gvcf, "GiB") + size(ref_fasta, "GiB") + size(dbsnp_vcf, "GiB") * 3)

  command <<<
    set -e

    gatk --java-options -Xms8g \
      GnarlyGenotyper \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename} \
      -D ~{dbsnp_vcf} \
      --only-output-calls-starting-in-intervals \
      --keep-all-sites \
      -V ~{combined_gvcf} \
      -L ~{interval}
  >>>

  runtime {
    memory: "26 GiB"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }
  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}


task IndelsVariantRecalibratorCreateModel {
  input {
    String recalibration_filename
    String tranches_filename
    Int downsampleFactor
    String model_report_filename

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File mills_resource_vcf
    File axiomPoly_resource_vcf
    File dbsnp_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 4

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"
  }

  parameter_meta {
        sites_only_variant_filtered_vcf: {
          localization_optional: true
        }
        sites_only_variant_filtered_vcf_index: {
          localization_optional: true
        }
      }

  command {
    set -euo pipefail

    gatk --java-options -Xms100g \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      -O ~{recalibration_filename} \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      -mode INDEL \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      --sample-every-Nth-variant ~{downsampleFactor} \
      --output-model ~{model_report_filename} \
      --max-gaussians ~{max_gaussians} \
      -resource:mills,known=false,training=true,truth=true,prior=12 ~{mills_resource_vcf} \
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ~{axiomPoly_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 ~{dbsnp_resource_vcf}
  }
  runtime {
    memory: "104 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
  }
  output {
    File model_report = "${model_report_filename}"
  }
}

task GatherMetrics {
  input{
    File input_details_fofn
    File input_summaries_fofn

    String output_prefix

    Int disk_size
  }

  command <<<
    set -e
    set -o pipefail

    # this is here to deal with the JES bug where commands may be run twice
    rm -rf metrics

    mkdir metrics
    RETRY_LIMIT=5

    count=0
    until cat ${input_details_fofn} | /usr/bin/gsutil -m cp -L cp.log -c -I metrics/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    count=0
    until cat ${input_summaries_fofn} | /usr/bin/gsutil -m cp -L cp.log -c -I metrics/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    INPUT=`cat ${input_details_fofn} | rev | cut -d '/' -f 1 | rev | sed s/.variant_calling_detail_metrics//g | awk '{printf("I=metrics/%s ", $1)}'`

    java -Xmx2g -Xms2g -jar /usr/gitc/picard.jar \
    AccumulateVariantCallingMetrics \
    $INPUT \
    O= ${output_prefix}
  >>>
  runtime {
    memory: "3 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"

  }
  output {
    File detail_metrics_file = "${output_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "${output_prefix}.variant_calling_summary_metrics"
  }
}