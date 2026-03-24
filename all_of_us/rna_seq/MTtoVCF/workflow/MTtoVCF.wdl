version 1.0

workflow MTtoVCF {
    input {
        String UriMatrixTable
        String OutputBucket 
        String OutputPrefix
        String CloudTmpdir
        String Branch = "main"
    }
    String pipeline_version = "aou_9.0.0"

    call WriteVCF {
        input: 
            PathMT = UriMatrixTable,
            OutputBucket = OutputBucket,
            OutputPrefix = OutputPrefix,
            CloudTmpdir = CloudTmpdir,
            Branch = Branch
    }

    output {
        String PathVCF = WriteVCF.PathVCF
    }
}

    task WriteVCF {
        input {
            String PathMT 
            String OutputBucket 
            String OutputPrefix
            String CloudTmpdir
            String Branch
        }  
    command <<<
        set -e

        export SPARK_LOCAL_DIRS=/cromwell_root

        # writes VCF to bucket path 
        # and also generates outpath.txt upon completion 
        # of writing VCF 
        python3 /ExportVCF.py \
            --MatrixTable ~{PathMT} \
            --OutputBucket ~{OutputBucket} \
            --OutputPrefix ~{OutputPrefix} \
            --CloudTmpdir ~{CloudTmpdir}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/aou_rna_mttovcf:0.0.1"
        memory: "256G"
        cpu: 64
        disks: "local-disk 2000 SSD"
    }
    
    meta {
        author: "Jonathan Nguyen"
    }
    
    output {
        String PathVCF = read_string('outpath.txt') 
    }
}
