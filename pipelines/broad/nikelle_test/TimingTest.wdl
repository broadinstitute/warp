version 1.0

import "../../../tasks/broad/nikelles_tasks.wdl" as NikellesTasks

workflow TimingTest {

    input {
        String name
    }

    call NikellesTasks.WriteHelloGreetingWorkflow {
        input:
            name = name
    }

    call WriteHelloGreetingTask {
        input:
            name = name
    }

    output {
        String greeting = WriteHelloGreetingWorkflow.output_file

    }
    meta {
        allowNestedInputs: true
    }
}


task WriteHelloGreetingTask {
    input {
        String name
    }
    command {
        echo "Hello ${name}!" > output.txt
    }
    output {
        File output_file = "output.txt"
    }
    runtime {
        docker: "ubuntu:latest"
        disks: "local-disk 10 HDD"
        memory: "2 GiB"
        preemptible: 3
    }
}