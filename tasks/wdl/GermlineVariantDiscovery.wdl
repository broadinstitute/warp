version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for germline variant discovery of human whole-genome or exome sequencing data.
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

task HaplotypeCaller_GATK35_GVCF {
  input {
    File input_bam
    File input_bam_index
    File interval_list
    String gvcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Int preemptible_tries
    Int hc_scatter
    #Setting default docker value for workflows that haven't yet been azurized.
    String docker = "us.gcr.io/broad-gotc-prod/gatk:1.3.0-4.2.6.1-1649964384"
  }

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(((size(input_bam, "GiB") + 30) / hc_scatter) + ref_size) + 20

  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them.
  #
  # Using PrintReads is a temporary solution until we update HaploypeCaller to use GATK4. Once that is done,
  # HaplotypeCaller can stream the required intervals directly from the cloud.
  command {
    /usr/gitc/gatk4/gatk --java-options "-Xms2000m -Xmx9000m"\
      PrintReads \
      -I ~{input_bam} \
      --interval-padding 500 \
      -L ~{interval_list} \
      -O local.sharded.bam \
    && \
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m -Xmx9000m\
      -jar /usr/gitc/GATK35.jar \
      -T HaplotypeCaller \
      -R ~{ref_fasta} \
      -o ~{gvcf_basename}.vcf.gz \
      -I local.sharded.bam \
      -L ~{interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination ~{default=0 contamination} \
      --read_filter OverclippedRead
  }
  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "10000 MiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_gvcf = "~{gvcf_basename}.vcf.gz"
    File output_gvcf_index = "~{gvcf_basename}.vcf.gz.tbi"
  }
}

task HaplotypeCaller_GATK4_VCF {
  input {
    File input_bam
    File input_bam_index
    File interval_list
    String vcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Boolean make_bamout
    Int preemptible_tries
    Int hc_scatter
    Boolean run_dragen_mode_variant_calling = false
    Boolean use_dragen_hard_filtering = false
    Boolean use_spanning_event_genotyping = true
    File? dragstr_model
    #Setting default docker value for workflows that haven't yet been azurized.
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int memory_multiplier = 1
    String? billing_project
  }
  
  Int memory_size_mb = ceil(8000 * memory_multiplier) + 2000

  String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_file_name = vcf_basename + output_suffix

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(((size(input_bam, "GiB") + 30) / hc_scatter) + ref_size) + 50

  String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

  String requester_pays_flag = if defined(billing_project) then "--gcs-project-for-requester-pays ${billing_project}" else ""

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command <<<
    set -e

    # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
    # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
    # memory_size_gb because of Cromwell's retry with more memory feature.
    mem_unit=$(echo "${MEM_UNIT:-}" | cut -c 1)

    if [[ $mem_unit == "M" ]]; then
      available_memory_mb=$(awk "BEGIN {print int($MEM_SIZE)}")
    elif [[ $mem_unit == "G" ]]; then
      available_memory_mb=$(awk "BEGIN {print int($MEM_SIZE * 1024)}")
    else
      available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
    fi

    let java_memory_size_mb=available_memory_mb-1024
    echo Total available memory: ${available_memory_mb} MB >&2
    echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

    gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_file_name} \
      -contamination ~{default=0 contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      ~{true="--dragen-mode" false="" run_dragen_mode_variant_calling} \
      ~{false="--disable-spanning-event-genotyping" true="" use_spanning_event_genotyping} \
      ~{if defined(dragstr_model) then "--dragstr-params-path " + dragstr_model else ""} \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{bamout_arg} \
      ~{requester_pays_flag}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam
  >>>

  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "~{memory_size_mb} MiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File output_vcf = "~{output_file_name}"
    File output_vcf_index = "~{output_file_name}.tbi"
    File bamout = "~{vcf_basename}.bamout.bam"
  }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name
    Int preemptible_tries = 3
    #Setting default docker value for workflows that haven't yet been azurized.
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
  }

  Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
      MergeVcfs \
      INPUT=~{sep=' INPUT=' input_vcfs} \
      OUTPUT=~{output_vcf_name}
  }
  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "3000 MiB"
    disks: "local-disk ~{disk_size} HDD"
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task Reblock {

  input {
    File gvcf
    File gvcf_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String output_vcf_filename
    String docker_path
    Int additional_disk = 20
    String? annotations_to_keep_command
    String? annotations_to_remove_command
    Float? tree_score_cutoff
    Boolean move_filters_to_genotypes = false
  }

  Int disk_size = ceil((size(gvcf, "GiB")) * 4) + additional_disk


  command {
    set -e 

    # We can't always assume the index was located with the gvcf, so make a link so that the paths look the same
    # Use bash basename instead of WDL to support DRS: https://support.terra.bio/hc/en-us/community/posts/4405396480027
    ln -s ~{gvcf} $(basename ~{gvcf})
    ln -s ~{gvcf_index} $(basename ~{gvcf_index})

    gatk --java-options "-Xms3000m -Xmx3000m" \
      ReblockGVCF \
      -R ~{ref_fasta} \
      -V $(basename ~{gvcf}) \
      -do-qual-approx \
      --floor-blocks -GQB 20 -GQB 30 -GQB 40 \
      ~{annotations_to_keep_command} \
      ~{annotations_to_remove_command} \
      ~{"--tree-score-threshold-to-no-call " + tree_score_cutoff} \
      ~{if move_filters_to_genotypes then "--add-site-filters-to-genotype" else ""} \
      -O ~{output_vcf_filename}
  }

  runtime {
    memory: "3750 MiB"
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: 15
    preemptible: 3
    docker: docker_path
  }

  output {
    File output_vcf = output_vcf_filename
    File output_vcf_index = output_vcf_filename + ".tbi"
  }
}

task HardFilterVcf {
  input {
    File input_vcf
    File input_vcf_index
    String vcf_basename
    File interval_list
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
  }

  Int disk_size = ceil(2 * size(input_vcf, "GiB")) + 20
  String output_vcf_name = vcf_basename + ".filtered.vcf.gz"

  command {
    gatk --java-options "-Xms2000m -Xmx2500m" \
      VariantFiltration \
      -V ~{input_vcf} \
      -L ~{interval_list} \
      --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
      --filter-name "HardFiltered" \
      -O ~{output_vcf_name}
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "3000 MiB"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }
}

# This hard filtering matches DRAGEN 3.4.12. For later DRAGEN versions, this needs to be updated.
task DragenHardFilterVcf {
  input {
    File input_vcf
    File input_vcf_index
    Boolean make_gvcf
    String vcf_basename
    Int preemptible_tries
    String gatk_docker
  }

  Int disk_size = ceil(2 * size(input_vcf, "GiB")) + 20

  String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_vcf_name = vcf_basename + ".hard-filtered" + output_suffix

  command {
     gatk --java-options "-Xms2000m -Xmx2500m" \
      VariantFiltration \
      -V ~{input_vcf} \
      --filter-expression "QUAL < 10.4139" \
      --filter-name "DRAGENHardQUAL" \
      -O ~{output_vcf_name}
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "3000 MiB"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }
}

task CNNScoreVariants {
  input {
    File? bamout
    File? bamout_index
    File input_vcf
    File input_vcf_index
    String vcf_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
  }

  Int disk_size = ceil(size(bamout, "GiB") + size(ref_fasta, "GiB") + (size(input_vcf, "GiB") * 2))

  String base_vcf = basename(input_vcf)
  Boolean is_compressed = basename(base_vcf, "gz") != base_vcf
  String vcf_suffix = if is_compressed then ".vcf.gz" else ".vcf"
  String vcf_index_suffix = if is_compressed then ".tbi" else ".idx"
  String output_vcf = base_vcf + ".scored" + vcf_suffix
  String output_vcf_index = output_vcf + vcf_index_suffix

  String bamout_param = if defined(bamout) then "-I ~{bamout}" else ""
  String tensor_type = if defined(bamout) then "read-tensor" else "reference"

  command {
     gatk --java-options "-Xmx10000m" CNNScoreVariants \
       -V ~{input_vcf} \
       -R ~{ref_fasta} \
       -O ~{output_vcf} \
       ~{bamout_param} \
       -tensor-type ~{tensor_type}
  }

  output {
    File scored_vcf = "~{output_vcf}"
    File scored_vcf_index = "~{output_vcf_index}"
  }

  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "15000 MiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
  }
}

task FilterVariantTranches {

  input {
    File input_vcf
    File input_vcf_index
    String vcf_basename
    Array[String] snp_tranches
    Array[String] indel_tranches
    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
    File dbsnp_resource_vcf
    File dbsnp_resource_vcf_index
    String info_key
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
  }

  Int disk_size = ceil(size(hapmap_resource_vcf, "GiB") +
                        size(omni_resource_vcf, "GiB") +
                        size(one_thousand_genomes_resource_vcf, "GiB") +
                        size(dbsnp_resource_vcf, "GiB") +
                        (size(input_vcf, "GiB") * 2)
                      ) + 20

  command {

    gatk --java-options "-Xmx6000m" FilterVariantTranches \
      -V ~{input_vcf} \
      -O ~{vcf_basename}.filtered.vcf.gz \
      ~{sep=" " prefix("--snp-tranche ", snp_tranches)} \
      ~{sep=" " prefix("--indel-tranche ", indel_tranches)} \
      --resource ~{hapmap_resource_vcf} \
      --resource ~{omni_resource_vcf} \
      --resource ~{one_thousand_genomes_resource_vcf} \
      --resource ~{dbsnp_resource_vcf} \
      --info-key ~{info_key} \
      --create-output-variant-index true
  }

  output {
    File filtered_vcf = "~{vcf_basename}.filtered.vcf.gz"
    File filtered_vcf_index = "~{vcf_basename}.filtered.vcf.gz.tbi"
  }

  runtime {
    memory: "7000 MiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
    docker: gatk_docker
  }
}
