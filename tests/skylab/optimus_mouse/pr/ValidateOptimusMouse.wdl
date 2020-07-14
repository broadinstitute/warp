version 1.0

task ValidateOptimusMouse {
  input {
    File bam
    File matrix
    File matrix_row_index
    File matrix_col_index
    File gene_metrics
    File cell_metrics

    String expected_bam_hash
    String expected_matrix_hash
    String expected_gene_metric_hash
    String expected_cell_metric_hash
  }

  Int required_disk = ceil((size(bam, "G") + size(matrix, "G")) * 1.1)

  command <<<

    # catch intermittent failures
    set -eo pipefail

    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed. gzipped files are unzipped to avoid hashing malleable
    # metadata


    unzip "~{matrix}"
    matrix_hash=$(find . -name "*.npy" -type f -exec md5sum {} \; | sort -k 2 | md5sum | awk '{print $1}')
    gene_metric_hash=$(zcat "~{gene_metrics}" | md5sum | awk '{print $1}')
    cell_metric_hash=$(zcat "~{cell_metrics}" | md5sum | awk '{print $1}')

    # calculate hash as above, but ignore run-specific bam headers
    bam_hash=$(samtools view "~{bam}" | md5sum | awk '{print $1}')

    # test each output for equivalence, echoing any failure states to stdout
    fail=false
    if [ "$bam_hash" != "~{expected_bam_hash}" ]; then
      >&2 echo "bam_hash (${bam_hash}) did not match expected hash (~{expected_bam_hash})"
      fail=true
    fi

    if [ "$gene_metric_hash" != "~{expected_gene_metric_hash}" ]; then
      >&2 echo "gene_metric_hash ($gene_metric_hash) did not match expected hash (~{expected_gene_metric_hash})"
      fail=true
    fi

    if [ "$cell_metric_hash" != "~{expected_cell_metric_hash}" ]; then
      >&2 echo "cell_metric_hash ($cell_metric_hash) did not match expected hash (~{expected_cell_metric_hash})"
      fail=true
    fi

    if [ $fail == "true" ]; then exit 1; fi

  >>>
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk ${required_disk} HDD"
  }
}
