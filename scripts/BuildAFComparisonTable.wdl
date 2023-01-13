version 1.0

# This script is under review. It is not actively tested or maintained at this time.
workflow BuildAFComparisonTable {
    input {
        File thousandG_vcf
        File thouandG_vcf_index

        String gnomad_vcf_prefix
        String chr_prefix
        Array[String] gnomad_af_expressions
        Array[String] thousandG_af_expressions

        File reference_fasta
        File reference_index
        File reference_dict
        Int? merge_mem

    }

    call MakeSitesOnlyVcf {
        input:
            vcf = thousandG_vcf,
            vcf_index = thouandG_vcf_index
    }

    call RemoveSymbolicAlleles {
        input:
            original_vcf = MakeSitesOnlyVcf.sites_only_vcf,
            original_vcf_index = MakeSitesOnlyVcf.sites_only_vcf_index,
            output_basename = "symbolic_alleles_removed"
    }

    call AddVariantIDs {
        input:
            vcf = RemoveSymbolicAlleles.output_vcf
    }

    scatter (chrom_i in range(22)) {
      Int chrom = chrom_i + 1 ## range(22) gives you 0..21

      call AnnotateWithAF_t {
        input:
          vcf = AddVariantIDs.with_ids_vcf,
          vcf_index  = AddVariantIDs.with_ids_vcf_index,
          gnomad_vcf = gnomad_vcf_prefix + chrom + ".vcf.bgz",
          gnomad_vcf_index = gnomad_vcf_prefix + chrom + ".vcf.bgz.tbi",
          mem = if chrom == 1 then 64 else 32,
          interval = chr_prefix + chrom,
          ref_fasta = reference_fasta,
          ref_fasta_index = reference_index,
          ref_dict = reference_dict,
          output_basename = "annotated.snps_only." + chrom,
          expressions = gnomad_af_expressions
      }

      call VariantsToTable {
        input:
            vcf = AnnotateWithAF_t.annotated_vcf,
            infoFields = flatten([gnomad_af_expressions, thousandG_af_expressions])
      }
    }

    call AnnotateWithAF_t as AnnotateX {
        input:
            vcf = AddVariantIDs.with_ids_vcf,
            vcf_index  = AddVariantIDs.with_ids_vcf_index,
            gnomad_vcf = gnomad_vcf_prefix + "X.vcf.bgz",
            gnomad_vcf_index = gnomad_vcf_prefix + "X.vcf.bgz.tbi",
            mem = 32,
            interval = chr_prefix + "X",
            ref_fasta = reference_fasta,
            ref_fasta_index = reference_index,
            ref_dict = reference_dict,
            output_basename = "annotated.snps_only.X",
            expressions = gnomad_af_expressions
    }

    call VariantsToTable as VariantsToTableX{
        input:
            vcf = AnnotateX.annotated_vcf,
            infoFields = flatten([gnomad_af_expressions, thousandG_af_expressions])
      }

    call MergeTables {
        input:
            tables = flatten([VariantsToTable.variants_table,[VariantsToTableX.variants_table]]),
            mem = merge_mem
    }

    output {
        File table = MergeTables.table
    }
}

task AnnotateWithAF_t {
  input {
    File vcf
    File vcf_index
    File gnomad_vcf
    File gnomad_vcf_index
    String interval
    File ref_fasta
    File ref_fasta_index
    String output_basename
    File ref_dict
    Int mem = 16
    Array[String] expressions

  }

  Int disk_size = 400

  command <<<
    gatk VariantAnnotator -R ~{ref_fasta} -V ~{vcf} -L ~{interval} -O ~{output_basename}.vcf.gz  \
    --resource:gnomad ~{gnomad_vcf}  --expression ~{sep=" --expression " expressions} -LE
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: mem + " GB" # some of the gnomad vcfs are like 38 gigs so maybe need more ?
  }
  output {
    File annotated_vcf = "~{output_basename}.vcf.gz"
    File annotated_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task GatherVCFsCloud {
    input {
        Array[File] vcfs
    }

    Int disk_size = ceil(2* size(vcfs, "GB")) + 10

    parameter_meta {
        vcfs : {
            localization_optional : true
        }
    }

    command <<<
        gatk GatherVcfsCloud -I ~{sep=" -I " vcfs} --gather-type CONVENTIONAL -O merged.vcf.gz
    >>>

    runtime {
            docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
            disks: "local-disk " + disk_size + " HDD"
            memory: "16 GB"
        }

    output {
        File vcf_out = "merged.vcf.gz"
        File vcf_out_index = "merged.vcf.gz.tbi"
    }
}

task MakeSitesOnlyVcf {
    input {
        File vcf
        File vcf_index
    }

    Int disk_size = ceil(2* size(vcf, "GB")) + 10

    command <<<
        gatk MakeSitesOnlyVcf -I ~{vcf} -O sites_only.vcf.gz
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        disks: "local-disk " + disk_size + " HDD"
        memory: "16 GB"
    }

    output {
        File sites_only_vcf = "sites_only.vcf.gz"
        File sites_only_vcf_index = "sites_only.vcf.gz.tbi"
    }
}

task AddVariantIDs {
    input {
        File vcf
    }

    Int disk_size = ceil(2* size(vcf, "GB")) + 10

    command <<<
        bcftools annotate ~{vcf} --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' -Oz -o variants_with_ids.vcf.gz
        bcftools index -t variants_with_ids.vcf.gz
    >>>

    output {
        File with_ids_vcf = "variants_with_ids.vcf.gz"
        File with_ids_vcf_index = "variants_with_ids.vcf.gz.tbi"
    }
   runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        disks: "local-disk " + disk_size + " HDD"
        memory: "4 GB"
   }
}

task VariantsToTable {
    input {
        File vcf
        Array[String] infoFields

    }

    parameter_meta {
        vcf: {
          description: "vcf",
          localization_optional: true
        }
    }

    Int disk_size = ceil(2 * size(vcf, "GB")) + 10

    command <<<
        gatk VariantsToTable -V ~{vcf} -F CHROM -F POS -F ID -F ~{sep=" -F " infoFields} -O variants.tsv
    >>>

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
        disks: "local-disk " + disk_size + " HDD"
        memory: "16 GB"
    }

    output {
        File variants_table = "variants.tsv"
    }
}

task MergeTables {
    input {
        Array[File] tables
        Int mem = 16
    }

    Int disk_size = 10 + ceil(2.2 * size(tables, "GiB"))

    command <<<
        Rscript -<<"EOF"
            library(dplyr)
            library(readr)
            library(purrr)

            t <- c("~{sep = '","' tables}") %>% map(read_tsv, col_types = cols(CHROM=col_character())) %>% reduce(bind_rows)
            write_tsv(t, "merged_results.tsv")

        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse"
        disks: "local-disk " + disk_size + " HDD"
        memory: mem + " GB"
    }

    output {
        File table = "merged_results.tsv"
    }
}

task RemoveSymbolicAlleles {
  input {
    File original_vcf
    File original_vcf_index
    String output_basename
    Int disk_size = ceil(3*(size(original_vcf, "GB") + size(original_vcf_index, "GB")))
  }
  command {
    gatk SelectVariants -V ~{original_vcf} -xl-select-type SYMBOLIC -O ~{output_basename}.vcf.gz
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4 GB"
  }
}