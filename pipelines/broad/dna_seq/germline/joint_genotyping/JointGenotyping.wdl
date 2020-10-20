version 1.0

import "../../../../tasks/JointGenotypingTasks.wdl" as Tasks

workflow JointGenotyping {
  input {
    String sample_map
    File? yaml_file
    String? override_hail_docker
    String? override_region

    String cluster_name
    String output_dir
    String callset_name
    String google_project
    String autoscaling_policy_name
    String partitions
    String data_type
    Boolean overwrite_mt
    Boolean? gather_vcfs


    # gnarly inputs
    String unpadded_intervals_file

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf
    File dbsnp_vcf_index

    Int small_disk
    Int medium_disk
    Int huge_disk

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

  String hail_docker = select_first([override_hail_docker, "ruchim/hail-testing:latest"])
  String region = select_first([override_region, "us-central1"])

#  if (!endsWith(output_dir, "/")) {
#    output_dir = output_dir + "/"
#  }

  call StartClusterAndSubmitHailJob {
    input:
      cluster_name = cluster_name,
      hail_docker = hail_docker,
      google_project = google_project,
      region = region,
      yaml_file = yaml_file,
      autoscaling_policy_name = autoscaling_policy_name,
      sample_map = sample_map,
      output_dir = output_dir,
      overwrite_mt = overwrite_mt,
      callset_name = callset_name
  }

  call ConvertMatrixTableToVcf {
    input:
      matrix_table = StartClusterAndSubmitHailJob.matrix_table,
      hail_docker = hail_docker,
      output_dir = output_dir,
      callset_name = callset_name,
      partitions = partitions,
      cluster_name = cluster_name,
      region = region,
      google_project = google_project,
      data_type = data_type
  }

  ### Gnarly genotyper
  call TabixBGzippedFile {
        input:
          zipped_vcf = ConvertMatrixTableToVcf.output_vcf,
          zipped_vcf_path = ConvertMatrixTableToVcf.output_vcf
      }

  # Make a 2.5:1 interval number to samples in callset ratio interval list.
  # We allow overriding the behavior by specifying the desired number of vcfs
  # to scatter over for testing / special requests.
  # Zamboni notes say "WGS runs get 30x more scattering than Exome" and
  # exome scatterCountPerSample is 0.05, min scatter 10, max 1000
  Int unboundedScatterCount = vcf_count
  Int scatterCount = if unboundedScatterCount > 10 then unboundedScatterCount else 10 #I think weird things happen if scatterCount is 1 -- IntervalListTools is noop?
  call SplitIntervalList {
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
        combined_gvcf_index = select_first(TabixBGzippedFile.bucket_tabix_output),
        interval = unpadded_intervals[idx],
        output_vcf_filename = callset_name + "." + idx + ".vcf.gz",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        dbsnp_vcf = dbsnp_vcf
    }

    call HardFilterAndMakeSitesOnlyVcf {
      input:
        vcf = GnarlyGenotyperOnVcf.output_vcf,
        vcf_index = GnarlyGenotyperOnVcf.output_vcf_index,
        excess_het_threshold = excess_het_threshold,
        variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
        sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
        disk_size = medium_disk
    }
  }

  call GatherVcfs as SitesOnlyGatherVcf {
   input:
     input_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf,
     output_vcf_name = callset_name + ".sites_only.vcf.gz",
     disk_size = medium_disk
  }

  call IndelsVariantRecalibrator {
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

  call SNPsVariantRecalibratorCreateModel {
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
    call SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
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

  call GatherTranches as SNPGatherTranches {
    input:
      tranches = SNPsVariantRecalibratorScattered.tranches,
      output_filename = callset_name + ".snps.gathered.tranches",
      disk_size = small_disk
  }

  scatter (idx in range(length(HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf))) {
    call ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
        input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = SNPsVariantRecalibratorScattered.recalibration[idx],
        snps_recalibration_index = SNPsVariantRecalibratorScattered.recalibration_index[idx],
        snps_tranches = SNPGatherTranches.tranches,
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        disk_size = medium_disk,
        use_allele_specific_annotations = true
    }

    # for large callsets we need to collect metrics from the shards and gather them later
    if (!is_small_callset) {
      call CollectVariantCallingMetrics as CollectMetricsSharded {
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

    call Tasks.GatherVariantCallingMetrics {
      input:
        input_details = select_all(CollectMetricsSharded.detail_metrics_file),
        input_summaries = select_all(CollectMetricsSharded.summary_metrics_file),
        output_prefix = callset_name,
        disk_size = medium_disk
    }
 }

 #TODO: sample_map is a String. How can we read this and get line count?
 Array[Array[String]] sample_name_map_lines = read_tsv(sample_map)
 Int num_gvcfs = length(sample_name_map_lines)
 Boolean is_small_callset = select_first([gather_vcfs, num_gvcfs <= 1000])

 # for small callsets we can gather the VCF shards and then collect metrics on it
 if (is_small_callset) {

   call GatherVcfs as FinalGatherVcf {
     input:
       input_vcfs = ApplyRecalibration.recalibrated_vcf,
       output_vcf_name = callset_name + ".vcf.gz",
       disk_size = huge_disk
   }

   call CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
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
#  call PartitionSampleNameMap {
#    input:
#      sample_name_map = sample_name_map,
#      line_limit = 1000
#  }

#  scatter (idx in range(length(PartitionSampleNameMap.partitions))) {
#
#    Array[File] files_in_partition = read_lines(PartitionSampleNameMap.partitions[idx])
#
#    call CrossCheckFingerprint as CrossCheckFingerprintsScattered {
#      input:
#        gvcf_paths = files_in_partition,
#        vcf_path = fingerprinting_vcf,
#        sample_name_map = sample_name_map,
#        haplotype_database = haplotype_database,
#        output_base_name = callset_name + "." + idx,
#        scattered = true
#    }
#  }
#
#  call GatherPicardMetrics as GatherFingerprintingMetrics {
#    input:
#      metrics_files = CrossCheckFingerprintsScattered.crosscheck_metrics,
#      output_file_name = callset_name + ".fingerprintcheck",
#      disk_size = small_disk
#  }

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

    #File crosscheck_metrics = GatherFingerprintingMetrics.gathered_metrics
  }
}

task StartClusterAndSubmitHailJob {
  input {
    String google_project
    String cluster_name
    String region
    String hail_docker
    File? yaml_file
    String autoscaling_policy_name
    String sample_map
    String output_dir
    String callset_name
    Boolean overwrite_mt
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

    # Set region
    gcloud config set dataproc/region ~{region}

    # Submit hail workflow
    /test_submit.sh ~{sample_map} ~{cluster_name} ~{output_dir} ~{callset_name} ~{google_project} ~{overwrite_mt}
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
    String partitions
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
    ~{output_dir} ~{data_type} ~{callset_name} ~{partitions}
  }

  runtime {
    docker: hail_docker
  }

  output {
    String output_vcf = "~{output_dir}~{callset_name}.vcf.bgz"
  }
}

task TabixBGzippedFile {
  input {
    File zipped_vcf
    String zipped_vcf_path
    String localized_tabix_path = zipped_vcf + ".tbi"
    String bucket_tabix_path = zipped_vcf_path + ".tbi"
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"
  }

  command <<<
    gatk IndexFeatureFile -F ~{zipped_vcf}
    gsutil cp ~{zipped_vcf}".tbi" ~{bucket_tabix_path}
  >>>

  runtime {
    memory: "1 GB"
    disks: "local-disk 200 HDD"
    preemptible: 3
    docker: gatk_docker
  }

  output {
      String bucket_tabix_output = bucket_tabix_path
  }
}

task SplitIntervalList {

  input {
    File interval_list
    Int scatter_count
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Boolean sample_names_unique_done
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"
  }

  parameter_meta {
    interval_list: {
      localization_optional: true
    }
  }

  command <<<
    gatk --java-options -Xms3g SplitIntervals \
      -L ~{interval_list} -O  scatterDir -scatter ~{scatter_count} -R ~{ref_fasta} \
      -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW
   >>>

  runtime {
    memory: "3.75 GiB"
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
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

task HardFilterAndMakeSitesOnlyVcf {

  input {
    File vcf
    File vcf_index
    Float excess_het_threshold

    String variant_filtered_vcf_filename
    String sites_only_vcf_filename

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options -Xms3g \
      VariantFiltration \
      --filter-expression "ExcessHet > ~{excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ~{variant_filtered_vcf_filename} \
      -V ~{vcf}

    gatk --java-options -Xms3g \
      MakeSitesOnlyVcf \
      -I ~{variant_filtered_vcf_filename} \
      -O ~{sites_only_vcf_filename}
  >>>

  runtime {
    memory: "3.75 GiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File variant_filtered_vcf = "~{variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "~{variant_filtered_vcf_filename}.tbi"
    File sites_only_vcf = "~{sites_only_vcf_filename}"
    File sites_only_vcf_index = "~{sites_only_vcf_filename}.tbi"
  }
}

task IndelsVariantRecalibrator {
  input {
    String recalibration_filename
    String tranches_filename
    File? model_report

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
  }

  parameter_meta {
        sites_only_variant_filtered_vcf: {
          localization_optional: true
        }
        sites_only_variant_filtered_vcf_index: {
          localization_optional: true
        }
      }

  String model_report_arg = if defined(model_report) then "--input-model $MODEL_REPORT --output-tranches-for-scatter" else ""

  command {
    set -euo pipefail

    MODEL_REPORT=~{model_report}

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
      ~{model_report_arg} \
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
    File recalibration = "${recalibration_filename}"
    File recalibration_index = "${recalibration_filename}.idx"
    File tranches = "${tranches_filename}"
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


task SNPsVariantRecalibratorCreateModel {
  input {
    String recalibration_filename
    String tranches_filename
    Int downsampleFactor
    String model_report_filename

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File hapmap_resource_vcf
    File omni_resource_vcf
    File one_thousand_genomes_resource_vcf
    File dbsnp_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf_index
    File dbsnp_resource_vcf_index
    Int max_gaussians = 6
    Int java_mem = 100

    Int disk_size
    Boolean use_allele_specific_annotations
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.4.1"
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

    gatk --java-options -Xms~{java_mem}g \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      -O ~{recalibration_filename} \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      -mode SNP \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      --sample-every-Nth-variant ~{downsampleFactor} \
      --output-model ~{model_report_filename} \
      --max-gaussians ~{max_gaussians} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
  }
  runtime {
    memory: "104 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
  }
  output {
    File model_report = "${model_report_filename}"
  }
}

task SNPsVariantRecalibrator {
input {
    String recalibration_filename
    String tranches_filename
    File? model_report

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File hapmap_resource_vcf
    File omni_resource_vcf
    File one_thousand_genomes_resource_vcf
    File dbsnp_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf_index
    File dbsnp_resource_vcf_index
    Int max_gaussians = 6

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    Int? machine_mem_gb
    Boolean use_allele_specific_annotations

  }

  Int auto_mem = ceil(2 * size([sites_only_variant_filtered_vcf,
                                hapmap_resource_vcf,
                                omni_resource_vcf,
                                one_thousand_genomes_resource_vcf,
                                dbsnp_resource_vcf],
                     "GiB"))
  Int machine_mem = select_first([machine_mem_gb, if auto_mem < 7 then 7 else auto_mem])
  Int java_mem = machine_mem - 1


  String model_report_arg = if defined(model_report) then "--input-model $MODEL_REPORT --output-tranches-for-scatter" else ""

  command <<<
    set -euo pipefail

    MODEL_REPORT=~{model_report}

    gatk --java-options -Xms~{java_mem}g \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      -O ~{recalibration_filename} \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      -mode SNP \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
       ~{model_report_arg} \
      --max-gaussians ~{max_gaussians} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "~{machine_mem} GiB"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
  }
}

task GatherTranches {

  input {
    Array[File] tranches
    String output_filename
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"
  }

  parameter_meta {
    tranches: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail

    tranches_fofn=~{write_lines(tranches)}

    # Jose says:
    # Cromwell will fall over if we have it try to localize tens of thousands of files,
    # so we manually localize files using gsutil.
    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
    # PAPI doesn't do.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf tranches
    mkdir tranches
    RETRY_LIMIT=5

    count=0
    until cat $tranches_fofn | /usr/bin/gsutil -m cp -L cp.log -c -I tranches/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the tranches from the cloud' && exit 1
    fi

    cat $tranches_fofn | rev | cut -d '/' -f 1 | rev | awk '{print "tranches/" $1}' > inputs.list

    gatk --java-options -Xms6g \
      GatherTranches \
      --input inputs.list \
      --output ~{output_filename}
  >>>

  runtime {
    memory: "7 GiB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File tranches = "~{output_filename}"
  }
}

task ApplyRecalibration {
  input {
    String recalibrated_vcf_filename
    File input_vcf
    File input_vcf_index
    File indels_recalibration
    File indels_recalibration_index
    File indels_tranches
    File snps_recalibration
    File snps_recalibration_index
    File snps_tranches

    Float indel_filter_level
    Float snp_filter_level
    Boolean use_allele_specific_annotations

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"
  }

  command {
    set -euo pipefail

    gatk --java-options -Xms5g \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ~{input_vcf} \
      --recal-file ~{indels_recalibration} \
      --tranches-file ~{indels_tranches} \
      --truth-sensitivity-filter-level ~{indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \


    gatk --java-options -Xms5g \
      ApplyVQSR \
      -O ~{recalibrated_vcf_filename} \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ~{snps_recalibration} \
      --tranches-file ~{snps_tranches} \
      --truth-sensitivity-filter-level ~{snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \

  }
  runtime {
    memory: "7 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
  }

  output {
    File recalibrated_vcf = "${recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "${recalibrated_vcf_filename}.tbi"
  }
}

task GatherVcfs {

  input {
    Array[File] input_vcfs
    String output_vcf_name
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"
  }

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input ~{sep=" --input " input_vcfs} \
      --output ~{output_vcf_name}

    tabix ~{output_vcf_name}
  >>>

  runtime {
    memory: "7 GiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task CollectVariantCallingMetrics {

  input {
    File input_vcf
    File input_vcf_index
    String metrics_filename_prefix
    File dbsnp_vcf
    File dbsnp_vcf_index
    File interval_list
    File ref_dict
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options -Xms6g \
      CollectVariantCallingMetrics \
      --INPUT ~{input_vcf} \
      --DBSNP ~{dbsnp_vcf} \
      --SEQUENCE_DICTIONARY ~{ref_dict} \
      --OUTPUT ~{metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ~{interval_list}
  >>>

  output {
    File detail_metrics_file = "~{metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "~{metrics_filename_prefix}.variant_calling_summary_metrics"
  }

  runtime {
    memory: "7.5 GiB"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
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

task CrossCheckFingerprint {

  input {
    Array[File] gvcf_paths
    File vcf_path
    File sample_name_map
    File haplotype_database
    String output_base_name
    Boolean scattered = false
    Array[String] expected_inconclusive_samples = []
    Int disk = 2000
    String picard_docker = "us.gcr.io/broad-gotc-prod/gatk4-joint-genotyping:yf_fire_crosscheck_picard_with_nio_fast_fail_fast_sample_map"
  }

  parameter_meta {
    gvcf_paths: {
      localization_optional: true
    }
    vcf_path: {
      localization_optional: true
    }
  }

  Int num_gvcfs = length(gvcf_paths)
  Int cpu = if num_gvcfs < 32 then num_gvcfs else 32
  # Compute memory to use based on the CPU count, following the pattern of
  # 3.75GiB / cpu used by GCP's pricing: https://cloud.google.com/compute/pricing
  Int memMb = round(cpu * 3.75 * 1024)

  String output_name = output_base_name + ".fingerprintcheck"

  command <<<
    set -eu

    gvcfInputsList=~{write_lines(gvcf_paths)}
    vcfInputsList=~{vcf_path}

    cp $gvcfInputsList gvcf_inputs.list
    cp $vcfInputsList vcf_inputs.list

    java -Dpicard.useLegacyParser=false -Xms~{memMb - 512}m \
      -jar /usr/gitc/PicardPublicWithCrosscheckNIOandSampleMapping.jar \
      CrosscheckFingerprints \
      --INPUT gvcf_inputs.list \
      --SECOND_INPUT vcf_inputs.list \
      --HAPLOTYPE_MAP ~{haplotype_database} \
      --INPUT_SAMPLE_FILE_MAP ~{sample_name_map} \
      --CROSSCHECK_BY SAMPLE \
      --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
      --NUM_THREADS ~{cpu} \
      --SKIP_INPUT_READABLITY_TEST \
      ~{true='--EXIT_CODE_WHEN_MISMATCH 0' false='' scattered} \
      --OUTPUT ~{output_name}

    if ~{scattered}; then
      # UNEXPECTED_MATCH is not possible with CHECK_SAME_SAMPLE
      matches=$(grep "EXPECTED_MATCH" ~{output_name} | wc -l)

      # check inconclusive samples
      expectedInconclusiveSamples=("~{sep='" "' expected_inconclusive_samples}")
      inconclusiveSamplesCount=0
      inconclusiveSamples=($(grep 'INCONCLUSIVE' ~{output_name} | cut -f 1))
      for sample in ${inconclusiveSamples[@]}; do
        if printf '%s\n' ${expectedInconclusiveSamples[@]} | grep -P '^'${sample}'$'; then
          inconclusiveSamplesCount=$((inconclusiveSamplesCount+1))
        fi
      done

      total_matches=$((inconclusiveSamplesCount + matches))
      if [[ ${total_matches} -eq ~{num_gvcfs} ]]; then
        >&2 echo "Found the correct number of matches (~{num_gvcfs}) for this shard"
      else
        >&2 echo "ERROR: Found $total_matches 'EXPECTED_MATCH' records, but expected ~{num_gvcfs}"
        exit 1
      fi
    fi
  >>>

  runtime {
    memory: memMb + " MiB"
    disks: "local-disk " + disk + " HDD"
    preemptible: 0
    docker: picard_docker
  }

  output {
    File crosscheck_metrics = output_name
  }
}

task PartitionSampleNameMap {
  input {
    File sample_name_map
    Int line_limit
  }

  command {

    cut -f 2 ~{sample_name_map} > sample_paths
    split -l ~{line_limit} -d sample_paths partition_

    # Let the OS catch up with creation of files for glob command
    sleep 1
  }

  output {
    Array[File] partitions = glob("partition_*")
  }

  runtime {
    memory: "1 GiB"
    preemptible: 1
    disks: "local-disk 10 HDD"
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
  }
}

task CrossCheckFingerprint {
    input {
      Array[File] gvcf_paths
      File vcf_path
      File sample_name_map
      File haplotype_database
      String output_base_name
      Boolean scattered = false
      Array[String] expected_inconclusive_samples = []
      String picard_docker = "us.gcr.io/broad-gotc-prod/gatk4-joint-genotyping:yf_fire_crosscheck_picard_with_nio_fast_fail_fast_sample_map"
    }

  Int num_gvcfs = length(gvcf_paths)
  Int cpu = if num_gvcfs < 32 then num_gvcfs else 32
  # Compute memory to use based on the CPU count, following the pattern of
  # 3.75GiB / cpu used by GCP's pricing: https://cloud.google.com/compute/pricing
  Int memMb = round(cpu * 3.75 * 1024)
  Int disk = 100

  String output_name = output_base_name + ".fingerprintcheck"

  command <<<
    set -eu

    gvcfInputsList=~{write_lines(gvcf_paths)}
    vcfInputsList=~{vcf_path}

    cp $gvcfInputsList gvcf_inputs.list
    cp $vcfInputsList vcf_inputs.list

    java -Dpicard.useLegacyParser=false -Xms~{memMb - 512}m \
      -jar /usr/gitc/PicardPublicWithCrosscheckNIOandSampleMapping.jar \
      CrosscheckFingerprints \
      --INPUT gvcf_inputs.list \
      --SECOND_INPUT vcf_inputs.list \
      --HAPLOTYPE_MAP ~{haplotype_database} \
      --INPUT_SAMPLE_FILE_MAP ~{sample_name_map} \
      --CROSSCHECK_BY SAMPLE \
      --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
      --NUM_THREADS ~{cpu} \
      --SKIP_INPUT_READABLITY_TEST \
      ~{true='--EXIT_CODE_WHEN_MISMATCH 0' false='' scattered} \
      --OUTPUT ~{output_name}

    if ~{scattered}; then
      # UNEXPECTED_MATCH is not possible with CHECK_SAME_SAMPLE
      matches=$(grep "EXPECTED_MATCH" ~{output_name} | wc -l)

      # check inconclusive samples
      expectedInconclusiveSamples=("~{sep='" "' expected_inconclusive_samples}")
      inconclusiveSamplesCount=0
      inconclusiveSamples=($(grep 'INCONCLUSIVE' ~{output_name} | cut -f 1))
      for sample in ${inconclusiveSamples[@]}; do
        if printf '%s\n' ${expectedInconclusiveSamples[@]} | grep -P '^'${sample}'$'; then
          inconclusiveSamplesCount=$((inconclusiveSamplesCount+1))
        fi
      done

      total_matches=$((inconclusiveSamplesCount + matches))
      if [[ ${total_matches} -eq ~{num_gvcfs} ]]; then
        >&2 echo "Found the correct number of matches (~{num_gvcfs}) for this shard"
      else
        >&2 echo "ERROR: Found $total_matches 'EXPECTED_MATCH' records, but expected ~{num_gvcfs}"
        exit 1
      fi
    fi
  >>>

  runtime {
    memory: memMb + " MiB"
    disks: "local-disk " + disk + " HDD"
    preemptible: 0
    docker: picard_docker
  }

  output {
    File crosscheck_metrics = output_name
  }
}

task GatherPicardMetrics {
  input {
    Array[File] metrics_files
    String output_file_name
    Int disk_size
  }

  command {
    # Don't use this task to gather tens of thousands of files.
    # Cromwell can't handle it.

    # This cannot gather metrics with histograms

    head -n 7 ~{metrics_files[0]} > ~{output_file_name}

    for metrics_file in ~{sep=' ' metrics_files}; do
      sed -n '1,7d;p' $metrics_file | grep -v '^$' >> ~{output_file_name}
    done
  }

  output {
    File gathered_metrics = "~{output_file_name}"
  }

  runtime {
    cpu: 1
    memory: "3.75 GiB"
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
  }
}