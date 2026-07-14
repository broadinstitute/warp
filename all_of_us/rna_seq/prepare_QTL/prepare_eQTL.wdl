version 1.0

import  "calculate_phenotypePCs.wdl" as ComputePCs



task eqtl_prepare_expression {
    input {
        File CountGCT 
        File AnnotationGTF 
        File SampleList 
        String OutputPrefix 
        

        Int memory 
        Int disk_space 
        Int num_threads 

        }
    command <<<
        Rscript /tmp/PrepareExpression.R \
            --CountGCT ${CountGCT} \
            --AnnotationGTF ${AnnotationGTF} \
            --SampleList ${SampleList} \
            --OutputPrefix ${OutputPrefix}

        >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou_rna_prepareqtl:0.0.1"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: num_threads
    }

    output {
        File ExpressionBed="${OutputPrefix}.expression.bed.gz"
    }
}

workflow eQTLPrepareData {
    input {
        String OutputPrefix 
        File CountGCT 
        File AnnotationGTF 
        File SampleList 

        Int memory 
        Int disk_space 
        Int num_threads 
        
    }
    String pipeline_version = "aou_9.0.0" 
    
    call eqtl_prepare_expression {
        input:
            OutputPrefix = OutputPrefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            CountGCT  = CountGCT,
            AnnotationGTF = AnnotationGTF,
            SampleList = SampleList 
    }

    call ComputePCs.PhenotypePCs {
        input:
            BedFile = eqtl_prepare_expression.ExpressionBed,
            OutputPrefix = OutputPrefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    output {
        File BedFile = eqtl_prepare_expression.ExpressionBed 
        File PhenotypePCsOut = PhenotypePCs.OutPhenotypePCs 
    }
}
