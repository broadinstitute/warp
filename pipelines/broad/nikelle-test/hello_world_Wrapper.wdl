import "https://raw.githubusercontent.com/aawdeh/CellBender/aa-cbwithoutcuda/wdl/cellbender_remove_background_azure.wdl" as helloworld1
import "https://raw.githubusercontent.com/aawdeh/CellBender/aa-cbwithoutcuda/wdl/cellbender_remove_background_azure.wdl" as helloworld2

workflow HelloWorldWrapper {

    String s = "Hello World!"

    call helloworld1.HelloWorld1 as helloworld1 {
        input: s = s
    }

    call helloworld2.HelloWorld2 as helloworld2 {
        input: s = s
    }
}