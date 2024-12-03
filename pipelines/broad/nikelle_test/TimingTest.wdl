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

    output {
        String greeting = WriteHelloGreetingWorkflow.output_file

    }
    meta {
        allowNestedInputs: true
    }
}
