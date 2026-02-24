version 1.0

# Docker fixed to quay.io/biocontainers/htslib@sha256:ff9d466929dc2d587128afc213fc4516d

#task splitPhenotypeBed {
#    input {
#        File TensorQTLPermutations
#    }

    #String baseName = basename(PhenotypeBed, ".gz")

#    command <<<
#        zcat ~{TensorQTLPermutations} | awk '$18 < 0.05' | head -n 100  > significant_qtls.txt 
#        awk 'NR==1 {header=$0; next} {out=$1".txt"; print header > out; print >> out}' significant_qtls.txt 
#    >>>
#
#    output {
#        Array[File] splitFiles = glob("*.txt")
#    }
#    runtime {
#        docker: "quay.io/biocontainers/htslib:1.22.1--h566b1c6_0"
#        disks: "local-disk 500 SSD"
#        memory: "2GB"
#        cpu: "1"
#    }
#}

task PrepInputs {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        String PhenotypeID
        File PhenotypeBed
        File TensorQTLPermutations
        Int NumPrempt
    }
    command <<<
        echo "Extracting headers from files"
        headerPermutations=$(zcat ~{TensorQTLPermutations} | head -n 1)
        headerBed=$(zcat "~{PhenotypeBed}" | head -n 1)
        
        echo "Bed file header:"
        echo $headerBed

        echo "TensorQTL file header"
        echo $headerPermutations
        
        #zcat ~{GenotypeDosages} |  awk 'NR==1 { if ($0 ~ /^#/) print; else print "#" $0; exit }'  > dosage_header.txt
        zcat ~{GenotypeDosages} |  awk 'NR==1 {print $0; exit }'  > dosage_header.txt

        echo "Subsetting bed file"
        zcat ~{PhenotypeBed} | grep "~{PhenotypeID}" \
            | awk 'BEGIN{OFS="\t"} {$2=$2-1000000; $3=$3+1000000; if($2<1) $2=1; print}' \
            > feature.bed
        
        echo $headerBed > temp_header.txt
        cat temp_header.txt feature.bed | bgzip -c - > ~{PhenotypeID}.bed.bgz
        #tabix ~{PhenotypeID}.bed.bgz

        echo "Subsetting TensorQTL file"
        zcat ~{TensorQTLPermutations} | grep "~{PhenotypeID}" > feature.txt
        echo $headerPermutations > temp_header_perm.txt
        cat temp_header_perm.txt feature.txt > ~{PhenotypeID}.tensorqtl.txt

        echo "Subsetting dose file"
        (cat dosage_header.txt; tabix ~{GenotypeDosages} -R ~{PhenotypeID}.bed.bgz) | bgzip -c > ~{PhenotypeID}.dose.tsv.gz
        #tabix  ~{GenotypeDosages} -R ~{PhenotypeID}.bed.bgz | bgzip -c - > ~{PhenotypeID}.dose.tsv.gz
        tabix -s1 -b2 -e2 -S1 "~{PhenotypeID}.dose.tsv.gz"   
    >>>
    
    runtime {
        docker: "quay.io/biocontainers/htslib@sha256:ff9d466929dc2d587128afc213fc4516d936ccc5e7fa39f39d3769f76b471293"
        disks: "local-disk 500 SSD"
        preemptible: "${NumPrempt}"
        memory: "2GB"
        cpu: "1"
    }
    
    output {
        File SubsetBed = "~{PhenotypeID}.bed.bgz"
        #File SubsetBedIndex = "~{PhenotypeID}.bed.bgz.tbi" 
        File SubsetPermutationPvals = "~{PhenotypeID}.tensorqtl.txt"
        File SubsetDosages = "~{PhenotypeID}.dose.tsv.gz"
        File SubsetDosagesIndex = "~{PhenotypeID}.dose.tsv.gz.tbi"
    }

}


task susieR {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        File QTLCovariates
        File TensorQTLPermutations
        File SampleList
        File PhenotypeBed
        Int CisDistance
        String OutputPrefix
        File susie_rscript
        Int memory
        Int NumPrempt
    }

    command <<<
        set -euo pipefail
        Rscript ~{susie_rscript} \
            --genotype_matrix ~{GenotypeDosages} \
            --sample_meta ~{SampleList} \
            --phenotype_list ~{TensorQTLPermutations} \
            --expression_matrix ~{PhenotypeBed} \
            --covariates ~{QTLCovariates} \
            --out_prefix ~{OutputPrefix} \
            --cisdistance ~{CisDistance} \

    >>>

    runtime {
        docker: 'quay.io/kfkf33/susier@sha256:80434421e53169129fb4a786a4698b838ed63f6ae2c6ad0ee0a419e7610ee113'
        memory: "${memory}GB"
        disks: "local-disk 500 SSD"
        bootDiskSizeGb: 25
        preemptible: "${NumPrempt}"
        cpu: "1"
    }

    output {
        File SusieParquet = "${OutputPrefix}.parquet"
        File lbfParquet = "${OutputPrefix}.lbf_variable.parquet"
        File FullSusieParquet = "${OutputPrefix}.full_susie.parquet"
    }
}


#task MergeSusie {
#    input {
#    Array[File] SusieOutput
#    Int memory
#    String OutputPrefix
    #}
#    
#    command <<<
#    for file in ~{sep='\n' SusieOutput}; do
#    echo $file >> filelist.txt
#    done
#
#    Rscript merge_susie.R \ 
#       --FilePaths filelist.txt \
#       --OutputPrefix ~{OutputPrefix}
#   >>>
#
#runtime {
#        docker: 'quay.io/kfkf33/susier:v24.01.1'
#        memory: "${memory}GB"
#        disks: "local-disk 500 SSD"
#        bootDiskSizeGb: 25
#        cpu: "1"
#    }
#
#
#    output {
#    File MergedSusieParquet = "${OutputPrefix}_SusieMerged.parquet" 
#    File MergedSusieTsv = "${OutputPrefix}_SusieMerged.tsv.gz" 
#
#    }
#
#}

workflow susieR_workflow {
    input {
        File GenotypeDosages
        File GenotypeDosageIndex
        File QTLCovariates
        File TensorQTLPermutations
        File SampleList
        File PhenotypeBed
        Int CisDistance
        File susie_rscript
        Int memory
        Int NumPrempt
        String OutputPrefix
        String PhenotypeID
    }
    String pipeline_version = "aou_9.0.0"

    call PrepInputs {
        input:
            TensorQTLPermutations = TensorQTLPermutations,
            PhenotypeID = PhenotypeID,
            GenotypeDosages = GenotypeDosages,
            GenotypeDosageIndex = GenotypeDosageIndex,
            PhenotypeBed = PhenotypeBed,
            NumPrempt = NumPrempt
    }

    call susieR {
        input:
            GenotypeDosages = PrepInputs.SubsetDosages,
            GenotypeDosageIndex = PrepInputs.SubsetDosagesIndex,
            QTLCovariates = QTLCovariates,
            TensorQTLPermutations = PrepInputs.SubsetPermutationPvals,
            SampleList = SampleList,
            PhenotypeBed = PhenotypeBed ,
            CisDistance = CisDistance,
            OutputPrefix = PhenotypeID,
            susie_rscript = susie_rscript,
            memory = memory,
            NumPrempt = NumPrempt

        }
    
    #call MergeSusie {
    #    input:
    #        SusieOutput = susieR.SusieParquet,
    #        OutputPrefix = OutputPrefix 
    #
    #} 
    output {
        File SusieParquet = susieR.SusieParquet
        File SusielbfParquet = susieR.lbfParquet
        File FullSusieParquet = susieR.FullSusieParquet
        File SubsetBed = PrepInputs.SubsetBed
        #File SubsetBedIndex = PrepInputs.SubsetBedIndex
        File SubsetDosages = PrepInputs.SubsetDosages
        File SubsetDosagesIndex = PrepInputs.SubsetDosagesIndex
    }
}