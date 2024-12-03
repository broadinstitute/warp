version 1.0

workflow WriteHelloGreetingWorkflow {

    input {
        String name
    }

    call WriteHelloGreeting {
        input:
            name = name
    }

    output {
        String output_file = WriteHelloGreeting.output_file

    }
    meta {
        allowNestedInputs: true
    }
}


task WriteHelloGreeting {
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