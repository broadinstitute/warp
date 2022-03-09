version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data processing for Multi-sample Illumina Genotyping Arrays
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

workflow MultiSampleArrays {

  String pipeline_version = "1.6.1"

  input {
    File samples_fofn
    File sample_indices_fofn
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String callset_name

    Int preemptible_tries
  }

  call SplitFoFnToListFoFn as SampleFofn {
    input:
      fofn = samples_fofn
  }

  call SplitFoFnToListFoFn as IndexFofn {
    input:
      fofn = sample_indices_fofn
  }

  scatter (idx in range(length(SampleFofn.array_of_fofns))) {
    call CombineVCFs as FirstCombine {
      input:
        vcfs = read_lines(SampleFofn.array_of_fofns[idx]),
        vcf_indices = read_lines(IndexFofn.array_of_fofns[idx]),
        combined_vcf_name = callset_name + "." + idx + ".vcf.gz",
        preemptible_tries = preemptible_tries
    }
  }

  call CombineVCFs as FinalCombine {
    input:
      vcfs = FirstCombine.combined_vcf,
      vcf_indices = FirstCombine.combined_vcf_index,
      combined_vcf_name = callset_name + ".vcf.gz",
      preemptible_tries = preemptible_tries
  }

  output {
    File combined_vcf = FinalCombine.combined_vcf
    File combined_vcf_index = FinalCombine.combined_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}


task SplitFoFnToListFoFn {
  input {
    File fofn
  }

  command <<<
  mkdir fofn
  cd fofn

  # splits list of files into fofns of 600 files each
  split -l 600 ~{fofn}
  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    memory: "1 GiB"
  }

  output {
      Array[File] array_of_fofns = glob("fofn/*")
  }
}


task CombineVCFs {
  input {
    Array[File]+ vcfs
    Array[File]+ vcf_indices
    String combined_vcf_name

    Int preemptible_tries
  }

  Int disk_size = ceil((size(vcfs, "GiB") + size(vcf_indices, "GiB")) * 2.1) + 20

  command <<<
    java -Xmx24g -Dpicard.useLegacyParser=false \
         -jar /usr/picard/picard.jar CombineGenotypingArrayVcfs \
         -I   ~{sep=' -I ' vcfs} \
         -O   ~{combined_vcf_name}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    disks: "local-disk " + disk_size + " HDD"
    memory: "26 GiB"
    preemptible: preemptible_tries
  }

  output {
    File combined_vcf       = "~{combined_vcf_name}"
    File combined_vcf_index = "~{combined_vcf_name}.tbi"
  }
}
