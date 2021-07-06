version 1.0

task ValidateOptimusDescriptorAnalysisFiles {
  input {
    File optimus_descriptors_analysis_file_intermediate_bam_json
    String expected_optimus_descriptors_analysis_file_intermediate_bam_json_hash
    File optimus_descriptors_analysis_file_intermediate_loom_json
    String expected_optimus_descriptors_analysis_file_intermediate_loom_json_hash
    File optimus_descriptors_analysis_file_intermediate_reference_json
    String expected_optimus_descriptors_analysis_file_intermediate_reference_json_hash
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail
       # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
       # a hash and the filename that was passed. gzipped files are unzipped to avoid hashing malleable
       # metadata

       #testing descriptors/analsysis_file/.bam
       optimus_descriptors_analysis_file_intermediate_bam_json_hash=$(cat "~{optimus_descriptors_analysis_file_intermediate_bam_json}" | md5sum | awk '{print $1}')

       # test each output for equivalence, echoing any failure states to stdout
       fail=false

       if [ "$optimus_descriptors_analysis_file_intermediate_bam_json_hash" != "~{expected_optimus_descriptors_analysis_file_intermediate_bam_json_hash}" ]; then
         >&2 echo "optimus_descriptors_analysis_file_intermediate_bam_json_hash ($optimus_descriptors_analysis_file_intermediate_bam_json_hash) did not match expected hash (~{expected_optimus_descriptors_analysis_file_intermediate_bam_json_hash})"
         fail=true
       fi

     if [ $fail == "true" ]; then exit 1; fi

       #testing descriptors/analsysis_file/.loom
       optimus_descriptors_analysis_file_intermediate_loom_json_hash=$(cat "~{optimus_descriptors_analysis_file_intermediate_loom_json}" | md5sum | awk '{print $1}')

       fail=false

       if [ "$optimus_descriptors_analysis_file_intermediate_loom_json_hash" != "~{expected_optimus_descriptors_analysis_file_intermediate_loom_json_hash}" ]; then
         >&2 echo "optimus_descriptors_analysis_file_intermediate_loom_json_hash ($optimus_descriptors_analysis_file_intermediate_loom_json_hash) did not match expected hash (~{expected_optimus_descriptors_analysis_file_intermediate_loom_json_hash})"
         fail=true
       fi

     if [ $fail == "true" ]; then exit 1; fi

       #testing descriptors/analsysis_file/.fasta
       optimus_descriptors_analysis_file_intermediate_reference_json_hash=$(cat "~{optimus_descriptors_analysis_file_intermediate_reference_json}" | md5sum | awk '{print $1}')

       fail=false

       if [ "$optimus_descriptors_analysis_file_intermediate_reference_json_hash" != "~{expected_optimus_descriptors_analysis_file_intermediate_reference_json_hash}" ]; then
         >&2 echo "optimus_descriptors_analysis_file_intermediate_reference_json_hash ($optimus_descriptors_analysis_file_intermediate_reference_json_hash) did not match expected hash (~{expected_optimus_descriptors_analysis_file_intermediate_reference_json_hash})"
         fail=true
       fi

     if [ $fail == "true" ]; then exit 1; fi

  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateOptimusLinksFiles {
  input {
    File optimus_links_intermediate_loom_json
    String expected_optimus_links_intermediate_loom_json
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail
       # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
       # a hash and the filename that was passed. gzipped files are unzipped to avoid hashing malleable
       # metadata

       #testing links/loom.json
       optimus_links_intermediate_loom_json_hash=$(cat "~{optimus_links_intermediate_loom_json}" | md5sum | awk '{print $1}')

       # test each output for equivalence, echoing any failure states to stdout
       fail=false

       if [ "$optimus_links_intermediate_loom_json_hash" != "~{expected_optimus_links_intermediate_loom_json}" ]; then
         >&2 echo "optimus_links_intermediate_loom_json_hash ($optimus_links_intermediate_loom_json_hash) did not match expected hash (~{expected_optimus_links_intermediate_loom_json})"
         fail=true
       fi

     if [ $fail == "true" ]; then exit 1; fi

  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateOptimusMetadataAnalysisFiles {
  input {
    File optimus_metadata_analysis_file_intermediate_bam_json
    String expected_optimus_metadata_analysis_file_intermediate_bam_json_hash
    File optimus_metadata_analysis_file_intermediate_loom_json
    String expected_optimus_metadata_analysis_file_intermediate_loom_json_hash
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail
       # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
       # a hash and the filename that was passed. gzipped files are unzipped to avoid hashing malleable
       # metadata

       #testing metadata/analsysis_file/.bam
       optimus_metadata_analysis_file_intermediate_bam_json_hash=$(cat "~{optimus_metadata_analysis_file_intermediate_bam_json}" | md5sum | awk '{print $1}')

       # test each output for equivalence, echoing any failure states to stdout
       fail=false

       if [ "$optimus_metadata_analysis_file_intermediate_bam_json_hash" != "~{expected_optimus_metadata_analysis_file_intermediate_bam_json_hash}" ]; then
         >&2 echo "optimus_metadata_analysis_file_intermediate_bam_json_hash ($optimus_metadata_analysis_file_intermediate_bam_json_hash) did not match expected hash (~{expected_optimus_metadata_analysis_file_intermediate_bam_json_hash})"
         fail=true
       fi

     if [ $fail == "true" ]; then exit 1; fi

       #testing metadata/analsysis_file/.loom
       optimus_metadata_analysis_file_intermediate_loom_json_hash=$(cat "~{optimus_metadata_analysis_file_intermediate_loom_json}" | md5sum | awk '{print $1}')

       fail=false

       if [ "$optimus_metadata_analysis_file_intermediate_loom_json_hash" != "~{expected_optimus_metadata_analysis_file_intermediate_loom_json_hash}" ]; then
         >&2 echo "optimus_metadata_analysis_file_intermediate_loom_json_hash ($optimus_metadata_analysis_file_intermediate_loom_json_hash) did not match expected hash (~{expected_optimus_metadata_analysis_file_intermediate_loom_json_hash})"
         fail=true
       fi

     if [ $fail == "true" ]; then exit 1; fi

  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}