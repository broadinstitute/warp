version 1.0

# Workflow: sc5_2_bgenix_index
# Summary: Given a list of BGEN and sample files, output their corresponding index (.bgi) file.

## All BED files referenced in this python script are UCSC bed files (as opposed to PLINK bed files). A UCSC BED
## file is a tab-delimited text file with at least three columns without header: chromosome, start_position, stop_position.

## Copyright Broad Institute, 2023
##
## This WDL pipeline processes a list of BGEN and their corresponding sample files, to create the corresponding
## index files (.bgi) that enable rapid retrieval of data from the BGEN file based on genetic positions.
## It's designed to be used with human genomic data in the context of variant analysis and annotation.
##
## Requirements/expectations:
## - **This script CAN BE RUN in Terra and vanilla cromwell**
## - Human genomic data in BGEN format
##
## LICENSING:
## This script is released under the [Appropriate License] (e.g., MIT, BSD-3, etc.). Users are responsible
## for ensuring they are authorized to run all components of this script. Please consult the relevant
## documentation for licensing details of Hail and other tools used in this pipeline.
##
## For information on tool versions and parameters, refer to the specific Docker containers.


workflow bgenix_bgen {

    input {
        # Input file, will be gs URL while working in Terra
        File bgen_fofn # Manifest file containing paths to all .bgen files
        File bgen_sample_fofn # Manifest file containing paths to all .sample files, contain sample information
                              # associated with the genotypes stored in the main .begn file.
    }

    String pipeline_version = "aou_beta"

    # Variable assignments and operations on the input files
    String bgen_fofn_basename = basename(bgen_fofn) # bgen_fofn_basename will contain the base name (the file name
                                                    # without the directory path) of the bgen_fofn file
    Array[File] bgens = read_lines(bgen_fofn) # bgen array will contain File objects representing each BGEN file.
    Array[File] bgen_samples = read_lines(bgen_sample_fofn) # bgen_samples array will contain File objects representing
                                                            # each sample file.
    Array[Pair[File,File]] bgen_pair = zip(bgens, bgen_samples) # bgen_pair array will contain pairs of File objects,
                                                                # with each pair representing a BGEN file paired with
                                                                # its corresponding sample file.

    # Parallelize the task across elements of array bgen_pair
    scatter (p in bgen_pair) {
        call bgenix_index {
            input:
                bgen=p.left,
                bgen_sample=p.right
        }
    }

    call create_fofn {
        input:
            file_urls1=bgenix_index.bgi,
            output_prefix=bgen_fofn_basename + ".bgi"
    }

    output {
        Array[File] bgis = bgenix_index.bgi
        File bgenix_index_fofn = create_fofn.fofn1
    }
}

task bgenix_index {

    input {
        File bgen # BGEN file
        File bgen_sample # Sample file associated with the BGEN file
    }
    # Variable Assignment
    String bgen_basename = basename(bgen) # bgen_basename will contain the base name (file name without
                                          # directory path) of the bgen file.
    command <<<
        set -e
        cp ~{bgen} ~{bgen_basename}
        bgenix -index -g ~{bgen_basename}
    >>>

    output {
        File bgi = "~{bgen_basename}.bgi"
    }

    runtime {
        docker: "befh/bgen:v1.1.7"
        memory: "15 GB"
        cpu: "2"
        disks: "local-disk 800 HDD"
    }
}

task create_fofn {
    input {
        Array[File] file_urls1 # An array of File objects containing URLs for file paths.
        String output_prefix # Output prefix for the fofn1 file.
    }
    # Variable Assignment
    File fofn1_in = write_lines(file_urls1) # Temporary file containing the contents of file_urls1

    command <<<
        set -e
        cp ~{fofn1_in} ~{output_prefix}.fofn1.txt
    >>>
    output {
        File fofn1 = "~{output_prefix}.fofn1.txt"
    }
    runtime {
        # Runtime settings
        docker: "us.gcr.io/broad-gatk/gatk:4.2.6.1"
        memory: "3 GB"
        cpu: "1"
        disks: "local-disk 100 HDD"
    }
}