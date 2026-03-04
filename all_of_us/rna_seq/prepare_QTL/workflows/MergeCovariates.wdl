version 1.0


workflow MergeCovariates {
    input {
        File GenotypePCs 
        String OutputPrefix 
        File MolecularPCs 
    }
    
    call MergeCovariatesR {
        input:
            GenotypePCs = GenotypePCs,
            OutputPrefix = OutputPrefix,
            MolecularPCs = MolecularPCs
    }
    
    output {
        File QtlCovariates = MergeCovariatesR.QtlCovariates

    }
}


task MergeCovariatesR {
    input {
        File GenotypePCs 
        String OutputPrefix 
        File MolecularPCs 
    }

    command <<<
        Rscript /tmp/MergeCovariates.R \
            --GenotypePCs ~{GenotypePCs} \
            --MolecularPCs ~{MolecularPCs} \
            --OutputPrefix ~{OutputPrefix}
        >>>
    
        runtime {
            docker: "us.gcr.io/broad-gotc-prod/aou_rna_prepareqtl:0.0.1"
            memory: "96G"
            cpu: 1
            disks: "local-disk 100 SSD"
        }
    
        output {
            File QtlCovariates = "~{OutputPrefix}_QTL_covariates.tsv"
        }



}
