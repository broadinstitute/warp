version 1.0

workflow ReblockGVCF {

  String pipeline_version = "1.0.1"

  input {
    File gvcf
    File gvcf_index
    String docker_image = "us.gcr.io/broad-gatk/gatk:4.1.4.1"
  }

  String gvcf_basename = basename(gvcf, ".g.vcf.gz")

  call Reblock {
    input:
      gvcf = gvcf,
      gvcf_index = gvcf_index,
      output_vcf_filename = gvcf_basename + ".reblocked.g.vcf.gz",
      docker_image = docker_image
  }

  output {
    File output_vcf = Reblock.output_vcf
    File output_vcf_index = Reblock.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}

task Reblock {

  input {
    File gvcf
    File gvcf_index
    String output_vcf_filename
    String docker_image
  }

  Int disk_size = ceil(size(gvcf, "GiB")) * 2

  command {
    gatk --java-options "-Xms3g -Xmx3g" \
      ReblockGVCF \
      -V ~{gvcf} \
      -drop-low-quals \
      -do-qual-approx \
      --floor-blocks -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 \
      -O ~{output_vcf_filename}
  }

  runtime {
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
    docker: docker_image
  }

  output {
    File output_vcf = output_vcf_filename
    File output_vcf_index = output_vcf_filename + ".tbi"
  }
} 
