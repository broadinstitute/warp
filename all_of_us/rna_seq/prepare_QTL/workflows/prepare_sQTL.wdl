version 1.0
import  "https://raw.githubusercontent.com/AoU-Multiomics-Analysis/prepare_QTL/refs/heads/main/workflows/calculate_phenotypePCs.wdl" as ComputePCs

task PrepareSpliceData {
    input {
        File SampleList 
        File SpliceData 
        String OutputPrefix 
        
        Int memory 
        Int disk_space 
        Int num_threads
    }
    command {
        Rscript /tmp/PrepareSpliceData.R \
            --SpliceData ${SpliceData} \
            --SampleList ${SampleList} \
            --OutputPrefix ${OutputPrefix}
        }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    output {
        File SplicingBed="${OutputPrefix}.splicing.bed.gz"
        #File PhenotypeGroups = "${OutputPrefix}.phenotype_groups.tsv"
    }
 }

workflow sQTLPrepareData  {
    input { 
        File SampleList 
        File SpliceData
        String OutputPrefix

        Int memory 
        Int disk_space 
        Int num_threads 
    } 
    call PrepareSpliceData {
        input:
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            SampleList = SampleList,
            SpliceData = SpliceData,
            OutputPrefix = OutputPrefix
    }

    call ComputePCs.PhenotypePCs {
        input:
            BedFile = PrepareSpliceData.SplicingBed,
            OutputPrefix = OutputPrefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }
    output {
        File BedFile = PrepareSpliceData.SplicingBed 
        File PhenotypePCsOut = PhenotypePCs.OutPhenotypePCs
        #File PhenotypeGroups = PrepareSpliceData.PhenotypeGroups
    }

}
