version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow CreateImputationRefPanelBeagle {
    input {
        Array[File] ref_vcf
        Int disk_size
    }

    scatter (idx in range(length(ref_vcf))) {
        call BuildBref3 {
            input:
                vcf = ref_vcf[idx],
                disk_size = disk_size
        }
    }

    output {
        Array[File] out_bref3 = BuildBref3.out_bref3
    }
}

task BuildBref3 {
    input {
        File vcf
        Int disk_size
    }

    String name = basename(vcf, ".vcf.gz")

    command <<<
        java -jar bref3.22Jul22.46e.jar ~{vcf} > ~{name}.bref3
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