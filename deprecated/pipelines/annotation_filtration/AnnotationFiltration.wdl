version 1.0
# AnnotationFiltration is now deprecated 2025-03-06

import "../../../tasks/wdl/Funcotator.wdl" as Funcotator

workflow AnnotationFiltration {

  String pipeline_version = "1.2.9"

  input {
    Array[File] vcfs

    String ref_version
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? funcotator_interval_list

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    File? custom_data_source_tar_gz
  }

  parameter_meta {
    vcfs: "The VCF files to annotate and filter"
    ref_version: "The version of the genome build provided (hg19 or hg38)"
    gatk_docker: "Version of the gatk docker to use. Please note that when updating this parameter, the correct custom_data_source_tar_gz must be provided. For example, for gatk version 4.1.0.0, the corresponding custom_data_source is the funcotator_dataSources.v1.6.20190124g release."
    custom_data_source_tar_gz: "Custom funcotator data source for AoU. Must be kept up to date with the gatk_docker parameter."
  }

  scatter (vcf in vcfs) {
    String base_vcf = basename(vcf)
    String compression_suffix = if basename(base_vcf, "gz") != base_vcf then "gz" else ""
    String vcf_suffix = if basename(base_vcf, ".g.vcf.~{compression_suffix}") != base_vcf then ".g.vcf.~{compression_suffix}" else ".vcf.~{compression_suffix}"
    String vcf_base_name = basename(base_vcf, vcf_suffix)
    String vcf_index_suffix = if compression_suffix == "gz" then ".tbi" else ".idx"
    File vcf_index = vcf + vcf_index_suffix

    call Funcotator.Funcotate as Funcotate {
      input:
        input_vcf = vcf,
        input_vcf_idx = vcf_index,

        data_sources_tar_gz = custom_data_source_tar_gz,
        transcript_selection_mode = "ALL",
        filter_funcotations = false,

        output_format = "VCF",
        output_file_base_name = vcf_base_name + ".funcotated",
        compress = true,

        reference_version = ref_version,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        interval_list = funcotator_interval_list,
        use_gnomad_exome = true,
        use_gnomad_genome = true,

        gatk_docker = gatk_docker,
        memory_mb = 4000
    }

    call FilterFuncotations {
      input:
        funcotated_vcf = Funcotate.funcotated_vcf,
        funcotated_vcf_index = Funcotate.funcotated_vcf_index,
        output_file_base_name = vcf_base_name,

        ref_version = ref_version,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,

        gatk_docker = gatk_docker,
        machine_memory = 4
    }

    Array[String] input_vcf_to_variant_count = [vcf, FilterFuncotations.significant_variants_count]
  }

  call GatherFiltrationReport {
    input:
      vcf_to_variant_count_vcf = write_tsv(input_vcf_to_variant_count)
  }

  output {
    Array[File] significant_variant_vcfs = FilterFuncotations.significant_variants_vcf
    File filtration_report = GatherFiltrationReport.filtration_report
  }
  meta {
    allowNestedInputs: true
  }
}

task FilterFuncotations {
  # Inputs for this task
  input {
    File funcotated_vcf
    File funcotated_vcf_index

    String output_file_base_name

    String filtered_name = output_file_base_name + ".filtered.vcf.gz"
    String output_name = output_file_base_name + ".clinsig-variants.vcf.gz"
    String significant_variants_file = "significant_variants.count"

    String ref_version
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # Runtime parameters
    Int machine_memory = 3
    String gatk_docker
    Int preemptible_attempts = 3
    Int additional_disk = 0
    Int cpu_threads = 1
    Boolean use_ssd = false
  }


  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  # 3x the funcotated sizes because we localize, then run filter, then run select.
  # We run all the tools in one task because (from testing) the time it takes to spin
  # up the GATK docker is larger than the time it takes to run all 3 tools.
  Int disk_size = ceil((size(funcotated_vcf, "GiB") * 3) + ref_size) + 20 + additional_disk
  String disk_type = if use_ssd then "SSD" else "HDD"


  Int command_memory = (machine_memory - 1) * 1024

  command <<<
    set -e

    # Marks clinically-significant variants with PASS in the filter column.
    # Also adds a CLINSIG attribute explaining which rule causes a variant
    # to be marked as significant.
    gatk --java-options "-Xmx~{command_memory}m" \
      FilterFuncotations \
      --variant ~{funcotated_vcf} \
      --output ~{filtered_name} \
      --ref-version ~{ref_version} \
      --allele-frequency-data-source gnomad

    gatk --java-options "-Xmx~{command_memory}m" \
      SelectVariants \
      --reference ~{ref_fasta} \
      --variant ~{filtered_name} \
      --output ~{output_name} \
      --exclude-filtered

    # CountVariants outputs:
    #   Tool returned:
    #   <N>
    # We just want the value of N.
    gatk --java-options "-Xmx~{command_memory}m" \
      CountVariants \
      --variant ~{output_name} | tail -n 1
  >>>

  runtime {
    docker: gatk_docker
    memory: machine_memory + " GiB"
    bootDiskSizeGb: 15
    disks: "local-disk ~{disk_size} ~{disk_type}"
    preemptible: preemptible_attempts
    cpu: cpu_threads
  }

  output {
    File significant_variants_vcf = output_name
    File significant_variants_vcf_index = output_name + ".tbi"
    Int significant_variants_count = read_int(stdout())
  }
}

task GatherFiltrationReport {

  input {
    File vcf_to_variant_count_vcf
  }

  String report_name = "filtration_report.tsv"

  command <<<
    set -eo pipefail

    echo -e "vcf\thas_significant_variants" > ~{report_name}
    while read vcf variant_count || test -n "$vcf$variant_count"
    do
      if [[ $variant_count -gt 0 ]]
      then
        echo -e "$vcf\ttrue" >> ~{report_name}
      else
        echo -e "$vcf\tfalse" >> ~{report_name}
      fi
    done <"~{vcf_to_variant_count_vcf}"
  >>>

  output {
    File filtration_report = report_name
  }

  runtime {
    memory: "4 GiB"
    docker: "us.gcr.io/google-containers/alpine-with-bash:1.0"
    disks: "local-disk 10 HDD"
  }
}
