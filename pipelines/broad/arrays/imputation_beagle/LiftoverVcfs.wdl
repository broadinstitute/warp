version 1.0

# Liftover VCFs from hg19 to hg38
workflow LiftoverVcfs {

  String pipeline_version = "1.0.0"

  input {
    File vcf_path
    File vcf_index_path

    File liftover_chain

    String docker = "us.gcr.io/broad-gatk/gatk:4.2.6.1"
    Int min_disk_size = 100

    File hg38_reference_fasta
    File hg38_reference_fasta_index
    File hg38_reference_dict

    Int max_retries = 3
    Int preemptible_tries = 3
  }

  String vcf_basename = basename(vcf_path, ".vcf.gz")

  # Lift over the array to hg38.
  call LiftOverArrays {
    input:
      input_vcf = vcf_path,
      input_vcf_index = vcf_index_path,
      liftover_chain = liftover_chain,
      reference_fasta = hg38_reference_fasta,
      reference_dict = hg38_reference_dict,
      output_basename = vcf_basename,
      docker = docker,
      max_retries = max_retries,
      preemptible_tries = preemptible_tries,
      min_disk_size = min_disk_size
  }

  output {
    File hg38_vcf = LiftOverArrays.lifted_over_vcf
    File hg38_vcf_index = LiftOverArrays.lifted_over_vcf_index
  }
}

task LiftOverArrays {
  input {
    File input_vcf
    File input_vcf_index
    File liftover_chain
    File reference_fasta
    File reference_dict
    String output_basename
    String docker
    Int max_retries
    Int preemptible_tries
    Int min_disk_size
  }

  Int disk_size_from_file = (ceil(size(input_vcf, "GiB") + size(liftover_chain, "GiB") + size(reference_fasta, "GiB")) * 2) + 20
  Int disk_size = if ( disk_size_from_file > min_disk_size ) then disk_size_from_file else min_disk_size

  command <<<
    set -euo pipefail

    gatk --java-options "-Xms4g -Xmx15g" \
    LiftoverVcf \
    --INPUT ~{input_vcf} \
    --OUTPUT ~{output_basename}.liftedover.vcf \
    --CHAIN ~{liftover_chain} \
    --REJECT ~{output_basename}.rejected_variants.vcf \
    --REFERENCE_SEQUENCE ~{reference_fasta} \
    --MAX_RECORDS_IN_RAM 100000

    # compress vcf - this creates a file with .gz suffix
    bgzip ~{output_basename}.liftedover.vcf

    # generate new index - this creates a file with .tbi suffix
    tabix ~{output_basename}.liftedover.vcf.gz
  >>>

  runtime {
    docker: docker
    memory: "16 GiB"
    cpu: "1"
    disks: "local-disk ~{disk_size} HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
  }

  output {
    File lifted_over_vcf = "~{output_basename}.liftedover.vcf.gz"
    File lifted_over_vcf_index = "~{output_basename}.liftedover.vcf.gz.tbi"
  }
}
