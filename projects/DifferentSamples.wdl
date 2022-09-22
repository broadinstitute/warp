version 1.0

workflow DifferentSamples {

  String pipeline_version = "1.0.0"

  input {
    Array[File] mixed_research_and_clinical_sample
  }

  parameter_meta {
    mixed_research_and_clinical_sample: "A single list of files that contains two (or more) types of files that need to be processed differently. Here we will use a file extension to differentiate."
  }
  # Allows you to specify teask level inputs from the input json file, even if they are not workflow inputs
  meta {
    allowNestedInputs: true
  }

  # How to take a different pipeline path for samples in a single array based on sample name
  # Note this could call a task or a subworkflow
  scatter (sample in mixed_research_and_clinical_sample) {
      if (sample == "*clinical*") {
          call ClinicalTask {
              input:
                  clinical_input = sample
          }
      }
      if (sample != "*clinical*") {
          call ResearchTask {
              input:
                  research_input = sample
          }
      }
   }

   # Note that the results of these subworkflow will be arrays of optional types (for example Array[File?]) If the type works, there is no need to explicitly gather. If you want to explictly gather the results, it would look like this:
   Array[File?] optional_clinical_results = ClinicalTask.results
   Array[File?] optional_research_results = ResearchTask.results

   # gather your samples from the clinical and research sub workflows (assuming both could have run) and convert to non-optional type
   Array[File] clincal_and_research_results = flatten([select_all(ClinicalTask.results), select_all(ResearchTask.results)])

   # if you didn't want to unify the results but need the results to be non-optional
   Array[File] clinical_results = select_all(ClinicalTask.results) # Array[Array[File?]]
   Array[File] research_results = select_all(ResearchTask.results)


   output {
       Array[File] unified_results = clincal_and_research_results
  }

}

task ClinicalTask {
    input {
        File clinical_input
        Int? mem_size_mb
        Int? disk_size_gb

    }
    Int dynamic_memory_size = ceil(size(clinical_input,"MiB")) + 1500
    Int dynamic_disk_size = ceil(size(clinical_input,"GiB")) + 20

    Int mem_size = select_first([mem_size_mb, dynamic_memory_size])
    Int disk_size = select_first([disk_size_gb, dynamic_disk_size])

    String output_filename = basename(clinical_input, ".txt") + ".results.txt"

    command <<<
        cp ~{clinical_input} ~{output_filename}
    >>>

    runtime {
        docker: "ubuntu:20.04"
        memory: mem_size + " MiB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File results = output_filename
    }
}

task ResearchTask {
    input {
        File research_input
        Int? mem_size_mb
        Int? disk_size_gb

    }
    Int dynamic_memory_size = ceil(size(research_input,"MiB")) + 1500
    Int dynamic_disk_size = ceil(size(research_input,"GiB")) + 20

    Int mem_size = select_first([mem_size_mb, dynamic_memory_size])
    Int disk_size = select_first([disk_size_gb, dynamic_disk_size])

    String output_filename = basename(research_input, ".txt") + ".results.txt"

    command <<<
        cp ~{research_input} ~{output_filename}
    >>>

    runtime {
        docker: "ubuntu:20.04"
        memory: mem_size + " MiB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File results = output_filename
    }
}


