version 1.0

workflow create_tarball {
    input {
        Array[File] required_files
        Array[File]? optional_files
	  String input_id
    }

    call tar_files {
        input:
            required_files = required_files,
            optional_files = optional_files,
            output_name = input_id + "_outputs.tar.gz"
    }

    output {
        File tarred_output = tar_files.tarball
    }
}

task tar_files {
    input {
        Array[File] required_files
        Array[File]? optional_files
        String output_name = output_name
    }

    command <<<
        tar -czf ${output_name} \
        ${sep=' ' required_files} \
        ~{sep=' ' optional_files}
    >>>

    output {
        File tarball = output_name
    }

    runtime {
        docker: "ubuntu:latest"
    }
}