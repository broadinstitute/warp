version 1.0

workflow MakeTable {

    String pipeline_version = "aou_9.0.0"

	input {
		Array[File] consensus_calls
		Array[String] sample_ids
	}

	call Combine { input: consensus_calls = consensus_calls, sample_ids = sample_ids }

	output {
		File result = Combine.result
	}
}

task Combine {
	input {
		Array[File] consensus_calls
		Array[String] sample_ids
		Int disk_space = 10
	}

	command <<<
        # prepare a file containing the sample ids, which will be the first column of the result
        touch sample_id_column.txt
        printf "sample_id\n" >> sample_id_column.txt
        for id in ~{sep=' ' sample_ids}; do
            printf "${id}\n" >> sample_id_column.txt
        done

        touch BODY
        for file in ~{sep=' ' consensus_calls}; do
            touch SAMPLE_HEADER
            # Pipe the result of these substitutions into a loop that reads every gene and two genotypes
            while read gene gt1 gt2 anything_else; do
                printf "${gene}_1\t${gene}_2\t" >> SAMPLE_HEADER
                printf "${gt1}\t${gt2}\t" >> BODY
            done < $file

            printf "\n" >> SAMPLE_HEADER
            printf "\n" >> BODY

            # compare the header from this sample to the header from others.  Hopefully identical
            if [ -e HEADER ]; then
                if cmp --silent -- HEADER SAMPLE_HEADER; then
                    echo "Good, headers match."
                else
                    echo "DANGER! Headers don't match"
                    cat HEADER
                    cat SAMPLE_HEADER
                fi
            else
                cp SAMPLE_HEADER HEADER
            fi
            rm SAMPLE_HEADER
        done

        cat HEADER BODY > genotypes.txt
        paste sample_id_column.txt genotypes.txt > result.txt
	>>>

    runtime {
        docker: "continuumio/anaconda:latest"
        disks: "local-disk " + disk_space + " SSD"
    }

    output {
        File result = "result.txt"
    }
}
