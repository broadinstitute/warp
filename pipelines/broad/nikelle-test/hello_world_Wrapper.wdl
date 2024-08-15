import "https://raw.githubusercontent.com/broadinstitute/warp/np_trst_firecloud_Api/pipelines/broad/nikelle-test/helloworld1.wdl" as helloworld1
import "https://raw.githubusercontent.com/broadinstitute/warp/np_trst_firecloud_Api/pipelines/broad/nikelle-test/helloworld2.wdl" as helloworld2

workflow HelloWorldWrapper {

    String s = "Hello World!"

    call helloworld1.HelloWorld1 as helloworld1 {
        input: s = s
    }

    call helloworld2.HelloWorld2 as helloworld2 {
        input: s = s
    }
}