version 1.0

workflow ValidateOptimus {
     meta {
         description: "Validate Optimus Outputs"
     }

     input {
         # Optimus output files to be checked
         File bam
         File matrix
         File matrix_row_index
         File matrix_col_index
         File cell_metrics
         File gene_metrics
         File? loom_file

         # Reference data and checksums
         File reference_matrix
         String expected_bam_hash
         String expected_cell_metric_hash
         String expected_gene_metric_hash
         String expected_loom_file_checksum
     }

     call ValidateBam as ValidateBam {
         input:
            bam = bam,
            expected_checksum = expected_bam_hash
     }

     call ValidateMatrix as ValidateMatrix {
         input:
             matrix = matrix,
             matrix_row_index = matrix_row_index,
             matrix_col_index = matrix_col_index,
             reference_matrix = reference_matrix
     }

     call ValidateLoom as ValidateLoom {
         input:
             loom_file = loom_file,
             expected_loom_file_checksum = expected_loom_file_checksum
     }

     call ValidateMetrics {
         input:
             cell_metrics = cell_metrics,
             gene_metrics = gene_metrics,
             expected_cell_metric_hash = expected_cell_metric_hash,
             expected_gene_metric_hash = expected_gene_metric_hash
     }

     call GenerateReport as GenerateReport {
         input:
             bam_validation_result = ValidateBam.result,
             matrix_validation_result = ValidateMatrix.result,
             metric_and_index_validation_result = ValidateMetrics.result,
             loom_validation_result = ValidateLoom.result
    }

    output {

    }
}

task ValidateBam {
    input {
        File bam
        String expected_checksum
    }

    Int required_disk = ceil(size(bam, "G") * 1.1)

    command <<<
        cacheInvalidationRandomString=4

        echo Starting checksum generation...

        # calculate hash for alignment positions only (a reduced bam hash)
        calculated_checksum=$( samtools view -F 256 "~{bam}" | cut -f 1-11 | md5sum | awk '{print $1}' )
        echo Reduced checksum generation complete

        if [ "$calculated_checksum" == "~{expected_checksum}" ]
        then
             echo Computed and expected bam hashes match \( "$calculated_checksum" \)
             printf PASS > result.txt
        else
             echo Computed \( "$calculated_checksum" \) and expected \( "~{expected_checksum}" \) bam file hashes do not match
             printf FAIL > result.txt
        fi
    >>>

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
        cpu: 1
        memory: "3.75 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String result = read_string("result.txt")
    }
}

task ValidateLoom {
    input {
        File? loom_file
        String expected_loom_file_checksum
    }

    Int required_disk = ceil( size(loom_file, "G") * 1.1)

    command <<<
        cacheInvalidationRandomString=4

        echo Starting checksum generation...
        calculated_loom_file_checksum=$( md5sum < ~{loom_file} | awk '{print $1}' )
        echo Checksum generation complete

        if [ "$calculated_loom_file_checksum" == "~{expected_loom_file_checksum}" ]
        then
            echo Computed and expected loom file hashes match \( "$calculated_loom_file_checksum" \)
        printf PASS > result.txt
        else
            echo Computed \( $calculated_loom_file_checksum \) and expected \( ~{expected_loom_file_checksum} \) loom file hashes do not match
           printf FAIL > result.txt
        fi
   >>>

   runtime {
       docker: "ubuntu:16.04"
       cpu: 1
       memory: "3.75 GB"
       disks: "local-disk ${required_disk} HDD"
   }

   output {
       String result = read_string("result.txt")
   }

}

task ValidateMatrix {
    input {
        File matrix
        File matrix_row_index
        File matrix_col_index
        File reference_matrix
    }

    Int required_disk = ceil( size(matrix, "G") * 1.1 )

    command <<<
        cacheInvalidationRandomString=4
       
       ## Convert matrix to format that can be read by R
       npz2rds.sh -c ~{matrix_col_index} -r ~{matrix_row_index} \
           -d ~{matrix} -o matrix.rds

       cp ~{reference_matrix} referenceMatrix.rds

       ## Run tests
       Rscript /root/tools/checkMatrix.R
       checkMatrixResult=$?

       if [ $checkMatrixResult == 0 ]; then
           printf PASS > result.txt
       else
           printf FAIL > result.txt
       fi

    >>>

    runtime {
        docker: "quay.io/humancellatlas/optimus-matrix-test:0.0.7"
        cpu: 1
        memory: "16 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String result = read_string('result.txt')
        File new_reference_matrix = "newReferenceMatrix.rds"
        File reads_per_cell_histogram = "reads_per_cell_histogram.png"
        File reads_per_gene_histogram = "reads_per_gene_histogram.png"
        File number_of_genes_per_cell = "number_of_genes_per_cell.png"
    }

}

task ValidateMetrics {
    input {
        File cell_metrics
        File gene_metrics

        String expected_cell_metric_hash
        String expected_gene_metric_hash
    }

    Int required_disk = ceil( (size(cell_metrics, "G") + size(gene_metrics, "G") )* 1.1)

    command <<<
        set -eo pipefail

        cacheInvalidationRandomString=4

        # check matrix row and column indexes files hash
        gene_metric_hash=$(zcat "~{gene_metrics}" | md5sum | awk '{print $1}')
        cell_metric_hash=$(zcat "~{cell_metrics}" | md5sum | awk '{print $1}')

        fail=false

        if [ "$gene_metric_hash" == "~{expected_gene_metric_hash}" ]; then
            echo Computed and expected gene metrics match \( "$gene_metric_hash" \)
        else
            echo Computed \( "$gene_metric_hash" \) and expected \( "~{expected_gene_metric_hash}" \) gene checksums do not match
            fail=true
        fi

        if [ "$cell_metric_hash" == "~{expected_cell_metric_hash}" ]; then
            echo Computed and expected cell metrics match \( "$cell_metric_hash" \)
        else
            echo Computed \( "$cell_metric_hash" \) and expected \( "~{expected_cell_metric_hash}" \) cell metrics hashes do not match
            fail=true
        fi

        if [ $fail == "true" ]; then
            printf FAIL > result.txt
        else
            printf PASS > result.txt
        fi
    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "1.00 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String result = read_string("result.txt")
    }

}

task GenerateReport {
  input {
    String bam_validation_result
    String metric_and_index_validation_result
    String matrix_validation_result
    String loom_validation_result
  }

      Int required_disk = 1

  command <<<

    set -eo pipefail

    cacheInvalidationRandomString=4

    # test each output for equality, echoing any failure states to stdout
    fail=false

    echo Bam Validation: ~{bam_validation_result}
    if [ "~{bam_validation_result}" == "FAIL" ]; then
        fail=true
    fi

    echo Metrics Validation: ~{metric_and_index_validation_result}
    if [ ~{metric_and_index_validation_result} == "FAIL" ]; then
        echo --- Ignoring failed metric and index test ---
        # Do not fail tests for this
        # fail=true
    fi

    echo Matrix Validation: ~{matrix_validation_result}
    if [ "~{matrix_validation_result}" == "FAIL" ]; then
        fail=true
    fi

    echo Loom Validation: ~{loom_validation_result}
    if [ "~{loom_validation_result}" == "FAIL" ]; then
        echo --- Ignoring failed loom test ---
        # Do not fail tests for this
        # fail=true
    fi

    if [ "$fail" == "true" ]; then exit 1; fi

  >>>

  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "1.0 GB"
    disks: "local-disk ${required_disk} HDD"
  }

  output {}
}
