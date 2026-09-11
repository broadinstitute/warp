version 1.0

workflow FilterMT {
    input {
        String UriMatrixTable
        File SampleList
        Int AlleleCountThreshold = 5
        Int AlleleNumberPercentage = 95
        String OutputBucket 
        String OutputPrefix
        String CloudTmpdir
        String Branch = "main"
    }
    String pipeline_version = "aou_9.0.0"
    
    call TaskFilterMT {
        input:
            UriMatrixTable = UriMatrixTable, 
            SampleList = SampleList,
            AlleleCountThreshold = AlleleCountThreshold,
            AlleleNumberPercentage = AlleleNumberPercentage,
            OutputBucket = OutputBucket,
            OutputPrefix = OutputPrefix,
            CloudTmpdir = CloudTmpdir,
            Branch = Branch
    }

    output {
        String FilteredMT = TaskFilterMT.FilteredMT
    }
}

task TaskFilterMT {
    input {
        String UriMatrixTable
        File SampleList
        Int AlleleCountThreshold
        Int AlleleNumberPercentage
        String OutputBucket 
        String OutputPrefix
        String CloudTmpdir
        String Branch
    }

    command <<<
        export SPARK_LOCAL_DIRS=/cromwell_root

        # writes VCF to bucket path 
        # and also generates outpath.txt upon completion 
        # of writing VCF 
        python3 /filter_and_write_mt.py \
            --MatrixTable ~{UriMatrixTable} \
            --SampleList ~{SampleList} \
            --AlleleCount ~{AlleleCountThreshold} \
            --AlleleNumberPercentage ~{AlleleNumberPercentage} \
            --OutputBucket ~{OutputBucket} \
            --OutputPrefix ~{OutputPrefix} \
            --CloudTmpdir ~{CloudTmpdir}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou_rna_mttovcf:0.0.1"
        memory: "256G"
        cpu: 64
        disks: "local-disk 1000 SSD"
    }
    
    output {
        String FilteredMT = read_string('outpath.txt') 
    }
}
