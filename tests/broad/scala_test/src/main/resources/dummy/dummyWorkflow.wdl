version 1.0

workflow DummyWorkflow {
    input {
        String message
        Int dummy_int
        String dummy_string
        Boolean dummy_boolean
        String? dummy_option
    }

    call PrintMessageToStdout {
        input:
            message = message
    }

    output {
        String echoed = PrintMessageToStdout.out
    }
}

task PrintMessageToStdout {
    input {
        String message
    }

    command {
        echo ~{message}
    }

    output {
        String out = read_string(stdout())
    }

    runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    }
}