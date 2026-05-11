version 1.0


workflow PostProcessVCF {
    input {
        File vcf_file
        #String output_prefix
        Int dosage_threads
    }
    call BcftoolsDosage {
        input:
            vcf_file = vcf_file,
            threads = dosage_threads
    }
    output {
        File GenotypeDosage = BcftoolsDosage.dosage
        File GenotypeDosagei = BcftoolsDosage.dosage_index

    }
    #call plink2_afreq {
    #    input:
    #        pgen = plink2.plink_pgen,
    #        psam = plink2.plink_psam,
    #        pvar = plink2.plink_pvar, 
    #        output_prefix = output_prefix
    #}

}

task BcftoolsDosage {
    input {
        File vcf_file
        Int threads
    }

    String vcf_basename = basename(vcf_file, ".vcf.gz")

    command <<<
        printf 'CHROM\nPOS\nREF\nALT\n' > 4_columns.tsv
        bcftools query -l ~{vcf_file} > sample_list.tsv
        cat 4_columns.tsv sample_list.tsv > header.tsv
        csvtk transpose header.tsv -T | gzip > header_row.tsv.gz

        #Extract dosage and merge
        bcftools +dosage --threads ~{threads} ~{vcf_file} -- -t GT | tail -n+2 | gzip > dose_matrix.tsv.gz
        zcat header_row.tsv.gz dose_matrix.tsv.gz | bgzip > ~{vcf_basename}.dose.tsv.gz
        tabix -s1 -b2 -e2 -S1 ~{vcf_basename}.dose.tsv.gz
        >>>
    
    runtime {
        docker: "quay.io/eqtlcatalogue/susie-finemapping:v20.08.1"
        memory: "32G"
        cpu: "${threads}"
        disks: "local-disk 500 SSD"
    }

    output {
        File dosage = "~{basename(vcf_file)}.dose.tsv.gz"
        File dosage_index = "~{basename(vcf_file)}.dose.tsv.gz.tbi"
    }
}
