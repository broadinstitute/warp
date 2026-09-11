version 1.0

workflow CreateVcfIndex {
  meta {
    description: "Create a tabix index (.tbi) for a VCF file"
    allowNestedInputs: true
  }

  input {
    File vcf_input

    Int disk_size_gb = 1000
    Int cpu = 1
    Int memory_gb = 16
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int preemptible = 3
  }

  call CreateVcfIndexTask {
    input:
      vcf_input = vcf_input,
      disk_size_gb = disk_size_gb,
      cpu = cpu,
      memory_gb = memory_gb,
      gatk_docker = gatk_docker,
      preemptible = preemptible
  }

  output {
    File output_vcf_index = CreateVcfIndexTask.output_vcf_index
  }
}

task CreateVcfIndexTask {
  input {
    File vcf_input

    Int disk_size_gb = 1000
    Int cpu = 1
    Int memory_gb = 16
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
    Int preemptible = 3
  }

  String vcf_basename = basename(vcf_input)

  command <<<
    set -e -o pipefail

    ln -sf ~{vcf_input} ~{vcf_basename}

    bcftools index -t ~{vcf_basename}
  >>>

  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size_gb + " SSD"
    memory: memory_gb + " GiB"
    cpu: cpu
    preemptible: preemptible
    maxRetries: 1
    noAddress: true
  }

  output {
    File output_vcf_index = vcf_basename + ".tbi"
  }
}
