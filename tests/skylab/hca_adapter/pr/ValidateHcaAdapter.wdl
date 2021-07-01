version 1.0

task OptimusValidateDescriptorAnalysisFiles {
  input {
    File descriptors_analysis_file_intermediate_bam_json
    String expected_descriptors_analysis_file_intermediate_bam_json_hash
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail
       # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
       # a hash and the filename that was passed. gzipped files are unzipped to avoid hashing malleable
       # metadata

       descriptors_analysis_file_intermediate_bam_json_hash=$(zcat "~{descriptors_analysis_file_intermediate_bam_json}" | md5sum | awk '{print $1}')

       # test each output for equivalence, echoing any failure states to stdout
       fail=false

       if [ "$descriptors_analysis_file_intermediate_bam_json_hash" != "~{expected_descriptors_analysis_file_intermediate_bam_json_hash}" ]; then
         >&2 echo "descriptors_analysis_file_intermediate_bam_json_hash ($descriptors_analysis_file_intermediate_bam_json_hash) did not match expected hash (~{expected_descriptors_analysis_file_intermediate_bam_json_hash})"
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

#task ValidateDescriptorReferenceFiles {
#  input {
#    File descriptors_reference_file_json
#  }
#
#  command <<<
#  >>>
#
#  runtime {
#    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
#    cpu: 1
#    memory: "3.75 GiB"
#    disks: "local-disk ${required_disk} HDD"
#  }
#}
#
#
#task ValidateLinksFiles {
#  input {
#    File links_intermediate_loom_json
#    File links_project_loom_json
#  }
#
#  command <<<
#  >>>
#
#  runtime {
#    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
#    cpu: 1
#    memory: "3.75 GiB"
#    disks: "local-disk ${required_disk} HDD"
#  }
#}
#
#
#task ValidateMetadataAnalysisFiles {
#  input {
#    File metadata_analysis_file_intermediate_loom_json
#    File metadata_analysis_file_intermediate_bam_json
#    File metadata_analysis_file_project_bam_json
#
#  }
#
#  command <<<
#  >>>
#
#  runtime {
#    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
#    cpu: 1
#    memory: "3.75 GiB"
#    disks: "local-disk ${required_disk} HDD"
#  }
#}
#
#task ValidateMetadataAnalysisProcessFiles {
#  input {
#    File metadata_analysis_process_file_optimus_json
#    File metadata_analysis_process_file_adapter_json
#    File metadata_analysis_process_file_merge_matrices_json
#  }
#
#  command <<<
#  >>>
#
#  runtime {
#    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
#    cpu: 1
#    memory: "3.75 GiB"
#    disks: "local-disk ${required_disk} HDD"
#  }
#}
#
#task ValidateMetadataAnalysisProtocolFiles {
#  input {
#    File metadata_analysis_protocol_file_optimus_json
#    File metadata_analysis_protocol_file_merge_matrices_json
#  }
#
#  command <<<
#  >>>
#
#  runtime {
#    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
#    cpu: 1
#    memory: "3.75 GiB"
#    disks: "local-disk ${required_disk} HDD"
#  }
#}
#
#task ValidateMetadataReferenceFile {
#  input {
#    File metadata_reference_file_json
#  }
#
#  command <<<
#  >>>
#
#  runtime {
#    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
#    cpu: 1
#    memory: "3.75 GiB"
#    disks: "local-disk ${required_disk} HDD"
#  }
#}