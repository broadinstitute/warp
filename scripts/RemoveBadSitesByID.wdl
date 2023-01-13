version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow RemoveBadSites {
    input {
        Array[File] vcfs
        File chrX_vcf
        File bad_sites_list
    }

    call SortVariantIdsList {
        input:
            idsList = bad_sites_list
    }

    scatter (vcf in vcfs) {
        call UpdateVariantIds {
            input:
                vcf = vcf,
                basename = basename(vcf, ".vcf.gz") + "updated_ids"
        }

        call SortVariantIds {
            input:
                vcf = UpdateVariantIds.output_vcf,
                basename = basename(vcf, ".vcf.gz") + "sorted_ids"
        }
        call RemoveBadSitesFromVcf {
            input:
                vcf = SortVariantIds.output_vcf,
                bad_sites_list = SortVariantIdsList.sorted_ids,
                basename = basename(vcf, ".vcf.gz")
        }

        call VcfToBcf {
            input:
                vcf = RemoveBadSitesFromVcf.out_vcf
        }

        call BuildM3VCF {
            input:
                vcf = RemoveBadSitesFromVcf.out_vcf,
                vcf_index = RemoveBadSitesFromVcf.out_vcf_index
        }
    }

    call UpdateVariantIds as UpdateVariantIdsX {
        input:
            vcf = chrX_vcf,
            basename = basename(chrX_vcf, ".vcf.gz") + "updated_ids"
    }

    call SortVariantIds as SortVariantIdsX {
        input:
            vcf = UpdateVariantIdsX.output_vcf,
            basename = basename(chrX_vcf, ".vcf.gz") + "sorted_ids"
    }

    call RemoveBadSitesFromVcf as RemoveBadSitesFromVcfX {
        input:
            vcf = SortVariantIdsX.output_vcf,
            bad_sites_list = SortVariantIdsList.sorted_ids,
            basename = basename(chrX_vcf, ".vcf.gz")
    }

    call VcfToBcf as VcfToBcfX {
        input:
            vcf = RemoveBadSitesFromVcfX.out_vcf
    }

    call SplitX {
        input:
            input_vcf = RemoveBadSitesFromVcfX.out_vcf,
            par1_end = 2699520,
            par2_start = 154931044,
            par2_end = 155270560
    }

    call BuildM3VCF as BuildM3VCFPAR1 {
        input:
            vcf = SplitX.par1_vcf,
            vcf_index = SplitX.par1_vcf_index
    }

    call BuildM3VCF as BuildM3VCFPAR2 {
        input:
            vcf = SplitX.par2_vcf,
            vcf_index = SplitX.par2_vcf_index
    }

    call BuildM3VCF as BuildM3VCFNON_PAR {
        input:
            vcf = SplitX.non_par_vcf,
            vcf_index = SplitX.non_par_vcf_index
    }

    output {
        Array[File] out_vcfs = flatten([RemoveBadSitesFromVcf.out_vcf,[RemoveBadSitesFromVcfX.out_vcf]])
        Array[File] out_vcf_indices = flatten([RemoveBadSitesFromVcf.out_vcf_index,[RemoveBadSitesFromVcfX.out_vcf_index]])
        Array[File] out_bcfs = flatten([VcfToBcf.out_bcf,[VcfToBcfX.out_bcf]])
        Array[File] out_bcf_indices = flatten([VcfToBcf.out_bcf_index,[VcfToBcfX.out_bcf_index]])
        Array[File] out_m3vcf = flatten([BuildM3VCF.out_m3vcf,[BuildM3VCFPAR1.out_m3vcf, BuildM3VCFPAR2.out_m3vcf, BuildM3VCFNON_PAR.out_m3vcf]])
    }
}

task SplitX {
    input {
        File input_vcf
        Int par1_end
        Int par2_start
        Int par2_end

    }

    Int disk_size = ceil(size(input_vcf, "GB")) + 50
    String name = basename(input_vcf, ".vcf.gz")

    parameter_meta {
        input_vcf : {
            localization_optional : true
        }
    }

    command <<<
        set -xeuo pipefail
        gatk SelectVariants -V ~{input_vcf} -L X:1-~{par1_end} -select "vc.getEnd()<~{par1_end}" -O PAR1.~{name}.vcf.gz
        gatk SelectVariants -V ~{input_vcf} -L X:~{par2_start}-~{par2_end} -select "vc.getStart()>=~{par2_start}" -O PAR2.~{name}.vcf.gz
        gatk SelectVariants -V ~{input_vcf} -L X:~{par1_end + 1}-~{par2_start} -select "vc.getStart()>~{par1_end} && vc.getStart()<~{par2_start}" -O NON_PAR.~{name}.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        disks: "local-disk " + disk_size + "  HDD"
        memory: "16 GB"
    }

    output {
        File par1_vcf = "PAR1.~{name}.vcf.gz"
        File par1_vcf_index = "PAR1.~{name}.vcf.gz.tbi"
        File par2_vcf = "PAR2.~{name}.vcf.gz"
        File par2_vcf_index = "PAR2.~{name}.vcf.gz.tbi"
        File non_par_vcf = "NON_PAR.~{name}.vcf.gz"
        File non_par_vcf_index = "NON_PAR.~{name}.vcf.gz.tbi"
    }
}

task BuildM3VCF {
    input {
        File vcf
        File vcf_index
    }

    String name = basename(vcf, ".vcf.gz")

    command <<<
        /Minimac3Executable/bin/Minimac3 --refHaps ~{vcf} --processReference --prefix ~{name}.cleaned --rsid
    >>>

    runtime {
        docker: "quay.io/ckachuli/minimac3@sha256:9b50682305869b103f2493179c9ea6ee618bd65d96bea17192dc609ce724c586"
        memory: "16 GB"
        disks: "local-disk 100 HDD"
    }

    output {
        File out_m3vcf = "~{name}.cleaned.m3vcf.gz"
    }
}

task VcfToBcf {
    input {
        File vcf
    }

    String name = basename(vcf, ".vcf.gz")

    command <<<

        bcftools view ~{vcf} -O b -o ~{name}.bcf
        bcftools index ~{name}.bcf
    >>>

    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        memory: "3 GiB"
        disks: "local-disk 100 HDD"
    }

    output {
        File out_bcf = "~{name}.bcf"
        File out_bcf_index = "~{name}.bcf.csi"
    }
}

task RemoveBadSitesFromVcf {
    input {
        File vcf
        File bad_sites_list
        String basename

    }

    parameter_meta {
        vcf : {
            localization_optional : true
        }
    }

    command <<<
        set -xeuo pipefail

        cp ~{bad_sites_list} bad_sites.list

        gatk SelectVariants -V ~{vcf} --exclude-ids bad_sites.list -O ~{basename}.cleaned.vcf.gz
    >>>

    runtime {
                docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        disks: "local-disk 100 HDD"
        memory: "16 GB"
    }

    output {
        File out_vcf = "~{basename}.cleaned.vcf.gz"
        File out_vcf_index = "~{basename}.cleaned.vcf.gz.tbi"
    }
}

task UpdateVariantIds {

  input {
    File vcf
    String basename
    Int disk_space =  3*ceil(size(vcf, "GB"))
  }

  command <<<
    bcftools annotate --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' ~{vcf} -O z -o ~{basename}.vcf.gz
  >>>

  output {
    File output_vcf = "~{basename}.vcf.gz"
  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}

task SortVariantIds {
  input {
    File vcf
    String basename
    Int disk_space =  3*ceil(size(vcf, "GB"))
  }

  command <<<
    # what better way is there to do this I really don't know
    zcat ~{vcf} | awk -v OFS='\t' '{split($3, n, ":"); if ( !($1 ~ /^"#"/) && n[4] < n[3])  $3=n[1]":"n[2]":"n[4]":"n[3]; print $0}' | bgzip -c > ~{basename}.vcf.gz
    bcftools index -t ~{basename}.vcf.gz
  >>>

  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"

  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}

task SortVariantIdsList {
  input {
    File idsList
    Int disk_space =  3*ceil(size(idsList, "GB")) + 100
  }

  command <<<
    # what better way is there to do this I really don't know
    awk -v OFS='\t' '{split($1, n, ":"); if ( !($1 ~ /^"#"/) && n[4] < n[3])  $1=n[1]":"n[2]":"n[4]":"n[3]; print $0}' ~{idsList} > sorted.ids.list
  >>>

  output {
    File sorted_ids = "sorted.ids.list"

  }

  runtime {
    docker: "skwalker/imputation:with_vcftools"
    disks: "local-disk " + disk_space + " HDD"
    memory: "16 GB"
  }
}