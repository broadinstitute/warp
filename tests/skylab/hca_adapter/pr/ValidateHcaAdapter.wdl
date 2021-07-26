version 1.0

task ValidateOptimusDescriptorAnalysisFiles {
  input {
    File optimus_descriptors_analysis_file_intermediate_loom_json
    File optimus_descriptors_analysis_file_intermediate_loom_json_truth
    File optimus_descriptors_analysis_file_intermediate_bam_json
    File optimus_descriptors_analysis_file_intermediate_bam_json_truth
    File optimus_descriptors_analysis_file_intermediate_reference_json
    File optimus_descriptors_analysis_file_intermediate_reference_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing descriptors/analsysis_file/.loom
       diff "~{optimus_descriptors_analysis_file_intermediate_loom_json}" "~{optimus_descriptors_analysis_file_intermediate_loom_json_truth}"

         if [ $? -ne 0 ];
         then
          echo "error"
         fi

       #testing descriptors/analsysis_file/.bam
       diff "~{optimus_descriptors_analysis_file_intermediate_bam_json}" "~{optimus_descriptors_analysis_file_intermediate_bam_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi

       #testing descriptors/reference_file/.fasta
       diff "~{optimus_descriptors_analysis_file_intermediate_reference_json}" "~{optimus_descriptors_analysis_file_intermediate_reference_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
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
    File optimus_links_intermediate_loom_json_truth
    }

 command <<<
  set -eo pipefail
  diff "~{optimus_links_intermediate_loom_json}" "~{optimus_links_intermediate_loom_json_truth}"

  if [ $? -ne 0 ];
  then
   echo "error"
  fi
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
    File optimus_metadata_analysis_file_intermediate_bam_json_truth
    File optimus_metadata_analysis_file_intermediate_loom_json
    File optimus_metadata_analysis_file_intermediate_loom_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing metadata/analysis_file/.bam
       diff "~{optimus_metadata_analysis_file_intermediate_bam_json}" "~{optimus_metadata_analysis_file_intermediate_bam_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi

       #testing metadata/analysis_file/.loom
       diff "~{optimus_metadata_analysis_file_intermediate_loom_json}" "~{optimus_metadata_analysis_file_intermediate_loom_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateOptimusMetadataAnalysisProcessFiles {
  input {
    File optimus_metadata_analysis_process_file_intermediate_json
    File optimus_metadata_analysis_process_file_intermediate_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing metadata/analysis_process/.json
       diff "~{optimus_metadata_analysis_process_file_intermediate_json}" "~{optimus_metadata_analysis_process_file_intermediate_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi

  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateOptimusMetadataAnalysisProtocolFiles {
  input {
    File optimus_metadata_analysis_protocol_file_intermediate_json
    File optimus_metadata_analysis_protocol_file_intermediate_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing metadata/analysis_protocol/.json
       diff "~{optimus_metadata_analysis_protocol_file_intermediate_json}" "~{optimus_metadata_analysis_protocol_file_intermediate_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateOptimusMetadataReferenceFiles {
  input {
    File optimus_metadata_reference_file_intermediate_json
    File optimus_metadata_reference_file_intermediate_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing metadata/reference_file/.json
       diff "~{optimus_metadata_reference_file_intermediate_json}" "~{optimus_metadata_reference_file_intermediate_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateOptimusDescriptorProjectLevelAnalysisFiles {
  input {
    File optimus_descriptors_analysis_file_project_loom_json
    File optimus_descriptors_analysis_file_project_loom_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing descriptors/analysis_file/project.loom.json
       diff "~{optimus_descriptors_analysis_file_project_loom_json}" "~{optimus_descriptors_analysis_file_project_loom_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateOptimusLinksProjectLevel {
  input {
    File optimus_links_project_loom_json
    File optimus_links_project_loom_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing links/project.loom.json
       diff "~{optimus_links_project_loom_json}" "~{optimus_links_project_loom_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}


task ValidateOptimusMetadataProjectLevelAnalysisFiles {
  input {
    File optimus_metadata_analysis_file_project_loom_json
    File optimus_metadata_analysis_file_project_loom_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing links/project.loom.json
       diff "~{optimus_metadata_analysis_file_project_loom_json}" "~{optimus_metadata_analysis_file_project_loom_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateOptimusMetadataProjectLevelAnalysisProcess {
  input {
    File optimus_metadata_analysis_process_project_loom_json
    File optimus_metadata_analysis_process_project_loom_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing links/project.loom.json
       diff "~{optimus_metadata_analysis_process_project_loom_json}" "~{optimus_metadata_analysis_process_project_loom_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}


task ValidateOptimusMetadataProjectLevelAnalysisProtocol {
  input {
    File optimus_metadata_analysis_protocol_project_loom_json
    File optimus_metadata_analysis_protocol_project_loom_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing links/project.loom.json
       diff "~{optimus_metadata_analysis_protocol_project_loom_json}" "~{optimus_metadata_analysis_protocol_project_loom_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}


task ValidateSS2DescriptorAnalysisFiles {
  input {
    File ss2_descriptors_analysis_file_intermediate_bam_json
    File ss2_descriptors_analysis_file_intermediate_bam_json_truth
    File ss2_descriptors_analysis_file_intermediate_bai_json
    File ss2_descriptors_analysis_file_intermediate_bai_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing descriptors/analsysis_file/.bam
       diff "~{ss2_descriptors_analysis_file_intermediate_bam_json}" "~{ss2_descriptors_analysis_file_intermediate_bam_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi

       #testing descriptors/analsysis_file/.bai
       diff "~{ss2_descriptors_analysis_file_intermediate_bai_json}" "~{ss2_descriptors_analysis_file_intermediate_bai_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}


task ValidateSS2LinksFiles {
  input {
    File ss2_links_intermediate_loom_json
    File ss2_links_intermediate_loom_json_truth
    File ss2_links_project_loom_json
    File ss2_links_project_loom_json_truth
    }

 command <<<
  set -eo pipefail
  diff "~{ss2_links_intermediate_loom_json}" "~{ss2_links_intermediate_loom_json_truth}"

  if [ $? -ne 0 ];
  then
   echo "error"
  fi

  diff "~{ss2_links_project_loom_json}" "~{ss2_links_project_loom_json_truth}"

    if [ $? -ne 0 ];
    then
     echo "error"
    fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateSS2MetadataAnalysisFiles {
  input {
    File ss2_metadata_analysis_file_intermediate_bam_json
    File ss2_metadata_analysis_file_intermediate_bam_json_truth
    File ss2_metadata_analysis_file_intermediate_bai_json
    File ss2_metadata_analysis_file_intermediate_bai_json_truth
    File ss2_metadata_analysis_file_project_loom_json
    File ss2_metadata_analysis_file_project_loom_json_truth

    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing metadata/analysis_file/.bam
       diff "~{ss2_metadata_analysis_file_intermediate_bam_json}" "~{ss2_metadata_analysis_file_intermediate_bam_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi

       #testing metadata/analysis_file/.bai
       diff "~{ss2_metadata_analysis_file_intermediate_bai_json}" "~{ss2_metadata_analysis_file_intermediate_bai_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi

       #testing metadata/analysis_file/.loom
        diff "~{ss2_metadata_analysis_file_project_loom_json}" "~{ss2_metadata_analysis_file_project_loom_json_truth}"

                 if [ $? -ne 0 ];
                 then
                  echo "error"
                 fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}


task ValidateSS2MetadataAnalysisProcessFiles {
  input {
    File ss2_metadata_analysis_process_file_intermediate_json
    File ss2_metadata_analysis_process_file_intermediate_json_truth
    File ss2_metadata_analysis_process_file_project_json
    File ss2_metadata_analysis_process_file_project_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing metadata/analysis_process/.json
       diff "~{ss2_metadata_analysis_process_file_intermediate_json}" "~{ss2_metadata_analysis_process_file_intermediate_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi

       #testing metadata/analysis_process/projectloom.json
       diff "~{ss2_metadata_analysis_process_file_project_json}" "~{ss2_metadata_analysis_process_file_project_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi

  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}

task ValidateSS2MetadataAnalysisProtocolFiles {
  input {
    File ss2_metadata_analysis_protocol_file_intermediate_json
    File ss2_metadata_analysis_protocol_file_intermediate_json_truth
    }

 command <<<
       # catch intermittent failures
       set -eo pipefail

       #testing metadata/analysis_protocol/.json
       diff "~{ss2_metadata_analysis_protocol_file_intermediate_json}" "~{ss2_metadata_analysis_protocol_file_intermediate_json_truth}"

                if [ $? -ne 0 ];
                then
                 echo "error"
                fi
  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GiB"
    disks: "local-disk 10 HDD"
  }
}