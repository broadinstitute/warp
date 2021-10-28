version 1.0

import "../tasks/AdapterTasks.wdl" as Tasks

workflow CreateSs2AdapterObjects {
  meta {
    description: "Creates json objects for indexing HCA 10x smartseq2 analysis data"
    allowNestedInputs: true
  }

  input {
    File? bam
    File? bai
    File? loom
    Int? ss2_index
    String input_id
    String version_timestamp
    String pipeline_type
    Boolean is_project_level
    String pipeline_version 
    String reference_file_fasta 
    File metadata
  }

  # If bam is defined then run intermediate level tasks
  if (defined(bam)) {
    call Tasks.GetAnalysisFileMetadata as GetIntermediateAnalysisFileMetadata {
      input:
        input_uuid = input_id,
        pipeline_type = pipeline_type,
        ss2_bam_file = select_first([bam]),
        ss2_bai_file = select_first([bai]),
        version_timestamp = version_timestamp,
        project_level = is_project_level
    }

    # pipeline_type is used for a dockstore URL here, so it needs to fit into this example:
    # intermediate level:
    #   "computational_method": "https://dockstore.org/workflows/github.com/broadinstitute/warp/Smartseq2_Single_Sample:SmartSeq2SingleSample_v5.1.1"
    call Tasks.GetAnalysisProtocolMetadata as GetIntermediateAnalysisProtocolMetadata {
      input:
        input_uuid = input_id,
        pipeline_type = "Smartseq2_Single_Sample",
        version_timestamp = version_timestamp,
        project_level = is_project_level,
        pipeline_version = pipeline_version
    }

    call Tasks.GetCloudFileCreationDate  as GetBamFileCreationDate {
      input:
        file_path = select_first([bam])
    }

    call Tasks.GetFileDescriptor as GetBamFileDescriptor {
      input:
        pipeline_type = pipeline_type,
        file_path = select_first([bam]),
        file_path_string = select_first([bam]),
        input_uuid = input_id,
        creation_time = GetBamFileCreationDate.creation_date,
        version_timestamp = version_timestamp
    }

    call Tasks.GetCloudFileCreationDate  as GetBaiFileCreationDate {
      input:
        file_path = select_first([bai])
    }

    call Tasks.GetFileDescriptor as GetBaiFileDescriptor {
      input:
        pipeline_type = pipeline_type,
        file_path = select_first([bai]),
        file_path_string = select_first([bai]),
        input_uuid = input_id,
        creation_time = GetBaiFileCreationDate.creation_date,
        version_timestamp = version_timestamp
    }
  }

  # If project level then loom is defined
  if (defined(loom)) {
    call Tasks.GetAnalysisFileMetadata as GetProjectAnalysisFileMetadata {
      input:
        input_uuid = input_id,
        pipeline_type = pipeline_type,
        version_timestamp = version_timestamp,
        project_level = is_project_level,
        project_loom = select_first([loom])
    }

    # pipeline_type is used for a dockstore URL here, so it needs to fit into this example:
    # project level:
    #   "computational_method": "https://dockstore.org/workflows/github.com/broadinstitute/warp/Smartseq2_Multisample:MultiSampleSmartSeq2_v2.1.4"
    call Tasks.GetAnalysisProtocolMetadata as GetProjectAnalysisProtocolMetadata {
      input:
        input_uuid = input_id,
        pipeline_type = "Smartseq2_Multisample",
        version_timestamp = version_timestamp,
        project_level = is_project_level,
        pipeline_version = pipeline_version
    }

    call Tasks.GetCloudFileCreationDate as GetLoomFileCreationDate {
      input:
        file_path = select_first([loom])
    }

    call Tasks.GetFileDescriptor as GetLoomFileDescriptor {
      input:
        pipeline_type = pipeline_type,
        file_path = select_first([loom]),
        file_path_string = select_first([loom]),
        input_uuid = input_id,
        creation_time = GetLoomFileCreationDate.creation_date,
        version_timestamp = version_timestamp
    }
  }

  # ss2_index will be undefined for project level
  call Tasks.GetAnalysisProcessMetadata {
    input:
      input_uuid = input_id,
      pipeline_type = pipeline_type,
      version_timestamp = version_timestamp,
      references = reference_file_fasta,
      project_level = is_project_level,
      input_file = metadata,
      ss2_index = ss2_index
  }

  output {
    Array[File] analysis_file_outputs = select_first([GetIntermediateAnalysisFileMetadata.analysis_file_outputs, GetProjectAnalysisFileMetadata.analysis_file_outputs])
    Array[File] analysis_process_outputs = GetAnalysisProcessMetadata.analysis_process_outputs
    Array[File] analysis_protocol_outputs = flatten(select_all([GetIntermediateAnalysisProtocolMetadata.analysis_protocol_outputs, GetProjectAnalysisProtocolMetadata.analysis_protocol_outputs]))
    File analysis_file_outputs_json = select_first([GetIntermediateAnalysisFileMetadata.outputs_json, GetProjectAnalysisFileMetadata.outputs_json])
    Array[File]? loom_file_descriptor_outputs = GetLoomFileDescriptor.file_descriptor_outputs
    Array[File]? bam_file_descriptor_outputs = GetBamFileDescriptor.file_descriptor_outputs
    Array[File]? bai_file_descriptor_outputs = GetBaiFileDescriptor.file_descriptor_outputs
  }
}

