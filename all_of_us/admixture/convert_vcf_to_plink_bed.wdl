version 1.0


workflow convert_vcf_to_plink_bed {
    input {
        String prefix
        File merged_vcf_shards
        File merged_vcf_shards_idx
    }
    String pipeline_version = "aou_9.0.1"

    call convert_vcf_to_plink_bed {
        input:
            prefix=prefix,
            vcf=merged_vcf_shards,
            vcf_idx=merged_vcf_shards_idx
    }

    output {
        File bed = convert_vcf_to_plink_bed.bed
        File bim = convert_vcf_to_plink_bed.bim
        File fam = convert_vcf_to_plink_bed.fam
    }
}
task convert_vcf_to_plink_bed {
    input {
        String prefix
        File vcf
        File vcf_idx
    }
    parameter_meta {
    }
    command <<<
        set -e
        /app/bin/plink --double-id --vcf ~{vcf} --make-bed --allow-extra-chr --out ~{prefix}
    >>>

    output {
        File bed = "~{prefix}.bed"
        File bim = "~{prefix}.bim"
        File fam = "~{prefix}.fam"
    }

    runtime {
        docker: "mussmann/admixpipe:3.0"
        memory: "31 GB"
        cpu: "4"
        disks: "local-disk 500 HDD"
    }
}