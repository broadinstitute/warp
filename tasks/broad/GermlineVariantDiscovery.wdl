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
    /usr/gitc/gatk4/gatk --java-options "-Xms2g" \
      PrintReads \
      -I ~{input_bam} \
      --interval-padding 500 \
      -L ~{interval_list} \
      -O local.sharded.bam \
    && \
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m \
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
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: preemptible_tries
    memory: "10 GiB"
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
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
  String output_file_name = vcf_basename + output_suffix

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(((size(input_bam, "GiB") + 30) / hc_scatter) + ref_size) + 20

  String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command <<<
    set -e
    gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_file_name} \
      -contamination ~{default=0 contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{bamout_arg}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam
  >>>

  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "6.5 GiB"
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
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar /usr/picard/picard.jar \
      MergeVcfs \
      INPUT=~{sep=' INPUT=' input_vcfs} \
      OUTPUT=~{output_vcf_name}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "3 GiB"
    disks: "local-disk ~{disk_size} HDD"
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task HardFilterVcf {
  input {
    File input_vcf
    File input_vcf_index
    String vcf_basename
    File interval_list
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  Int disk_size = ceil(2 * size(input_vcf, "GiB")) + 20
  String output_vcf_name = vcf_basename + ".filtered.vcf.gz"

  command {
     gatk --java-options "-Xms3000m" \
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
    memory: "3 GiB"
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
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
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
     gatk --java-options -Xmx10g CNNScoreVariants \
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
    memory: "15 GiB"
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
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  Int disk_size = ceil(size(hapmap_resource_vcf, "GiB") +
                        size(omni_resource_vcf, "GiB") +
                        size(one_thousand_genomes_resource_vcf, "GiB") +
                        size(dbsnp_resource_vcf, "GiB") +
                        (size(input_vcf, "GiB") * 2)
                      ) + 20

  command {

    gatk --java-options -Xmx6g FilterVariantTranches \
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
    memory: "7 GiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
    docker: gatk_docker
  }
}