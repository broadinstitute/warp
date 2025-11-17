version 1.0

task AggregateSusie{
    input{
        #Array[File] SusieParquets
        File SusieParquetsFOFN
        Int Memory
        String OutputPrefix
        Int NumThreads 
    }

    command <<<
    mkdir parquet_inputs
    # while read path; do
    #     gsutil -m cp "${path}" parquet_inputs/
    # done < ~{SusieParquetsFOFN}

    echo $(date +"[%b %d %H:%M:%S]")
    echo parallel copy start
    # Single invocation, reads entire list from STDIN
    gsutil -m cp -I parquet_inputs/ < ~{SusieParquetsFOFN}
    echo $(date +"[%b %d %H:%M:%S]")
    echo parallel copy done


    parquet_files=(parquet_inputs/*)
    printf "%s\n" ${parquet_files[@]} > parquet_file_list.txt
    # awk '{print $1}' ~{SusieParquetsFOFN} | grep -v '^$' > file_paths.txt 

    # mkdir -p localized
    # gsutil -m cp -I localized/ < file_paths.txt 

    # # Write the new local file paths into filelist.txt
    # echo "Listing files in localized directory:"
    # ls -1 "$(pwd)/localized/*"

    # echo "Creating file_paths.txt with local file paths:"
    # ls -1 "$(pwd)/localized/*" > filelist.txt
    
    echo "Running R script"
    Rscript /tmp/merge_susie.R --FilePaths parquet_file_list.txt  --OutputPrefix ~{OutputPrefix}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/aggregate_susie@sha256:5f742ab9e4d05d945eeb0d1faeef3d035068e62f1f591a75015ae7815e6b54bc"
        disks: "local-disk 500 SSD"
        memory: "~{Memory}GB"
        cpu: "~{NumThreads}"
    }
 


    output {
        File MergedSusieParquet = "${OutputPrefix}_SusieMerged.parquet" 
        File MergedSusieTsv = "${OutputPrefix}_SusieMerged.tsv.gz"
        File FileList = "parquet_file_list.txt" 
    } 

}


task AnnotateSusie {
    input {
        File SusieTSV 
        File GencodeGTF
        File PlinkAfreq
        String OutputPrefix
        Int Memory
        File AnnotationENCODE 
        File AnnotationFANTOM5 
        File AnnotationVEP
        File AnnotationVEPIndex 
        File AnnotationGnomad 
        File AnnotationPhyloP 
      

    }
    command <<<
    Rscript /tmp/annotate_susie_data.R \
        --OutputPrefix ~{OutputPrefix} \
        --GencodeGTF ~{GencodeGTF} \
        --PlinkAfreq ~{PlinkAfreq} \
        --SusieTSV ~{SusieTSV} \
        --phyloPBigWig ~{AnnotationPhyloP} \
        --FANTOM5 ~{AnnotationFANTOM5} \
        --gnomadConstraint ~{AnnotationGnomad} \
        --ENCODEcCRES ~{AnnotationENCODE} \
        --VEPAnnotationsTable ~{AnnotationVEP}
    >>>
   runtime {
        docker: "ghcr.io/aou-multiomics-analysis/aggregate_susie@sha256:5f742ab9e4d05d945eeb0d1faeef3d035068e62f1f591a75015ae7815e6b54bc"
        disks: "local-disk 500 SSD"
        memory: "~{Memory}GB"
        cpu: "1"
    }


    output {
        File AnnotatedSusieParquetOut = "~{OutputPrefix}_SusieMerged.annotated.tsv" 
    }

}



workflow AggregateSusieWorkflow {
    input {
        #Array[File] SusieParquets
        File SusieParquetsFOFN
        Int Memory 
        String OutputPrefix
        Int NumThreads

        File GencodeGTF 
        File PlinkAfreq
        File AnnotationPhyloP 
        File AnnotationENCODE 
        File AnnotationFANTOM5 
        File AnnotationVEP 
        File AnnotationGnomad
        File AnnotationVEPIndex 

    }
    String pipeline_version = "aou_9.0.0"
    
    call AggregateSusie {
        input:
            SusieParquetsFOFN = SusieParquetsFOFN,
            OutputPrefix = OutputPrefix,
            Memory = Memory,
            NumThreads = NumThreads
    }

    call AnnotateSusie {
        input:
            SusieTSV = AggregateSusie.MergedSusieTsv,
            GencodeGTF = GencodeGTF,
            PlinkAfreq = PlinkAfreq,
            OutputPrefix = OutputPrefix,
            Memory = Memory,
            AnnotationPhyloP = AnnotationPhyloP,
            AnnotationENCODE = AnnotationENCODE,
            AnnotationFANTOM5 = AnnotationFANTOM5,
            AnnotationVEP = AnnotationVEP,
            AnnotationVEPIndex = AnnotationVEPIndex,
            AnnotationGnomad = AnnotationGnomad
    } 

    output {
        File AnnotatedMergedSusieParquet = AnnotateSusie.AnnotatedSusieParquetOut
        #File MergedSusieParquet = AggregateSusie.MergedSusieParquet
        #File MergedSusieTsv = AggregateSusie.MergedSusieTsv
    }

}