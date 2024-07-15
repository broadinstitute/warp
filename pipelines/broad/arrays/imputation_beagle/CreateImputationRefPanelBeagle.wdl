version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow CreateImputationRefPanelBeagle {
  input {
    Array[File] ref_vcf
    Array[File] ref_vcf_index
    Int disk_size

    String? output_basename

    Boolean make_brefs = true
    Boolean make_interval_lists = true
  }

  scatter (idx in range(length(ref_vcf))) {
    Int? chr = idx + 1
    String? custom_basename_with_chr = output_basename + ".chr" + chr

    if (make_brefs) {
      call BuildBref3 {
        input:
          vcf = ref_vcf[idx],
          disk_size = disk_size,
          output_basename = custom_basename_with_chr
        }
      }
        
      if (make_interval_lists) {
        call CreateRefPanelIntervalLists {
          input:
            ref_panel_vcf = ref_vcf[idx],
            ref_panel_vcf_index = ref_vcf_index[idx],
            output_basename = custom_basename_with_chr,
        }
      }
  }

  output {
    Array[File?] bref3s = BuildBref3.out_bref3
    Array[File?] interval_lists = CreateRefPanelIntervalLists.interval_list
  }
}

task BuildBref3 {
  input {
    File vcf
    String? output_basename
    Int disk_size
  }

  String name_from_file = basename(vcf, ".vcf.gz")
  String name = select_first([output_basename, name_from_file])

  command <<<
    java -jar /usr/gitc/bref3.22Jul22.46e.jar ~{vcf} > ~{name}.bref3
  >>>

  runtime {
    docker: "us-central1-docker.pkg.dev/morgan-fieldeng-gcp/imputation-beagle-development/imputation-beagle:0.0.1-22Jul22.46e-wip-temp-20240227"
    memory: "256 GB"
    cpu: 4
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File out_bref3 = "~{name}.bref3"
  }
}

task CreateRefPanelIntervalLists {
  input {
    File ref_panel_vcf
    File ref_panel_vcf_index

    String? output_basename

    Int disk_size_gb = ceil(2*size(ref_panel_vcf, "GiB")) + 50 # not sure how big the disk size needs to be since we aren't downloading the entire VCF here
    Int cpu = 1
    Int memory_mb = 8000
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
  }

  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  String name_from_file = basename(ref_panel_vcf, ".vcf.gz")
  String basename = select_first([output_basename, name_from_file])

  command {
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    VcfToIntervalList \
    -I ~{ref_panel_vcf} \
    -O ~{basename}.interval_list
  }

  output {
    File interval_list = "~{basename}.interval_list"
  }

  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
  }
}
