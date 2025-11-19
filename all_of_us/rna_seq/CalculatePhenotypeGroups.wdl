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
        
        echo listing files
        ls -lh
        }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:b9986a803ad82aee02945e9131f180fc028c73eadb889184e23aff8dabf573fe"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    output {
        File PhenotypeGroups = "${OutputPrefix}.phenotype_groups.tsv"
    }
 }

workflow CalculatePhenotypeGroups  {
    input { 
        File SampleList 
        File SpliceData
        String OutputPrefix

        Int memory 
        Int disk_space 
        Int num_threads 
    }
    String pipeline_version = "aou_9.0.0" 
    call PrepareSpliceData {
        input:
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            SampleList = SampleList,
            SpliceData = SpliceData,
            OutputPrefix = OutputPrefix
    }

    output {
        File PhenotypeGroups = PrepareSpliceData.PhenotypeGroups
    }

}