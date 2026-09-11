version 1.0


import "MTtoVCF/FilterMTAndExportToVCF.wdl" as ProcessMT

workflow PrepareGenotypes{
    input {
        String UriMatrixTable 
        File SampleList
        Int AlleleCountThreshold
        Int max_allele_len 
        String OutputPrefix
        String OutpathVCF
        String OutpathMT
        String CloudTmpdir

        File genotype_rscript

    }
       call ProcessMT.FilterMTAndExportToVCF {
        input: 
            UriMatrixTable =  UriMatrixTable,
            OutputBucketVCF = OutpathVCF ,
            OutputPrefix = OutputPrefix ,
            AlleleCountThreshold = AlleleCountThreshold,
            OutputBucketCheckpointMT = OutpathMT,
            CloudTmpdir = CloudTmpdir,
            SampleList = SampleList

    }

    call plink2 {
        input:
            vcf_file = FilterMTAndExportToVCF.PathVCF,
            output_prefix = OutputPrefix,
            new_id_max_allele_len = max_allele_len,
    }
    
    call ComputeGenotypePCs {
        input:
            vcf_file = FilterMTAndExportToVCF.PathVCF,
            output_prefix = OutputPrefix,
            genotype_rscript = genotype_rscript
    }

    output {
        File VCF =  FilterMTAndExportToVCF.PathVCF
        File GenotypePCs =  ComputeGenotypePCs.output_tsv
        File pgen = plink2.pgen 
        File psam = plink2.psam 
        File pvar = plink2.pvar
    } 
}


task plink2 {
    input {
        File vcf_file
        String output_prefix
        Int new_id_max_allele_len
    }

    Int disk_size = ceil(size(vcf_file, "GB") * 2)

    command <<<
        set -e

       # mkdir -p plink_output

        plink2 --vcf "~{vcf_file}" \
        --make-pgen \
        --out "~{output_prefix}" \
        --set-all-var-ids @:#_\$r_\$a \
        --new-id-max-allele-len "~{new_id_max_allele_len}" \
        --output-chr chrM \
        --chr 1-22
    >>>

    runtime {
        docker: "quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0"
        memory: "16G"
        cpu: 4
        disks: "local-disk ~{disk_size} SSD"
    }

    output {
        File pgen  = "~{output_prefix}.pgen"
        File pvar  = "~{output_prefix}.pvar"
        File psam  = "~{output_prefix}.psam"

    }
}



task ComputeGenotypePCs {
    input {
        File vcf_file
        String output_prefix
        File genotype_rscript
    }

    Int disk_size = ceil(size(vcf_file, "GB") * 2)

    command <<<
        set -e

        Rscript "~{genotype_rscript}" \
            --vcf_path "~{vcf_file}" \
            --prefix "~{output_prefix}"
        >>>
    
        runtime {
            docker: "quay.io/jonnguye/genotype_pcs:micromamba"
            memory: "96G"
            cpu: 2
            disks: "local-disk ~{disk_size} SSD"
        }
    
        output {
            File output_tsv = "~{output_prefix}_genetic_PCs.tsv"
        }
}
