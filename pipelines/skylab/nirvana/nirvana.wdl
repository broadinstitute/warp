version 1.0

workflow AnnotateVCFWorkflow {
    input {
        File input_vcf
        String output_annotated_file_name
        Boolean use_reference_disk
        String cloud_provider
        File omim_annotations
    }

    # Determine docker prefix based on cloud provider
    String gcr_docker_prefix = "us.gcr.io/broad-gotc-prod/"
    String acr_docker_prefix = "dsppipelinedev.azurecr.io/"

    String docker_prefix = if cloud_provider == "gcp" then gcr_docker_prefix else acr_docker_prefix

    # Define docker images
    String nirvana_docker_image = "nirvana:np_add_nirvana_docker"

    call AnnotateVCF {
        input:
            input_vcf = input_vcf,
            output_annotated_file_name = output_annotated_file_name,
            use_reference_disk = use_reference_disk,
            cloud_provider = cloud_provider,
            omim_annotations = omim_annotations,
            docker_path = docker_prefix + nirvana_docker_image
    }

    output {
        File positions_annotation_json = AnnotateVCF.positions_annotation_json
        File genes_annotation_json = AnnotateVCF.genes_annotation_json
    }
}

task AnnotateVCF {
    input {
        File input_vcf
        String output_annotated_file_name
        Boolean use_reference_disk
        File omim_annotations
        String cloud_provider
        String docker_path
    }


    String annotation_json_name = output_annotated_file_name + ".json.gz"
    String gene_annotation_json_name = output_annotated_file_name + ".genes.json.gz"
    String positions_annotation_json_name = output_annotated_file_name + ".positions.json.gz"
    String nirvana_location = "/Nirvana/Nirvana.dll"
    String jasix_location = "/Nirvana/Jasix.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace


        if [[ "~{use_reference_disk}" == "true" ]]
        then
            # There's an issue with how the projects/broad-dsde-cromwell-dev/global/images/nirvana-3-18-1-references-2023-01-03
            # disk image was built: while all the reference files do exist on the image they are not at the expected
            # locations. The following code works around this issue and should continue to work even after a corrected
            # version of the Nirvana reference image is deployed into Terra.

            # Find where the reference disk should have been mounted on this VM.  Note this is referred to as a "candidate
            # mount point" because we do not actually confirm this is a reference disk until the following code block.
            CANDIDATE_MOUNT_POINT=$(lsblk | sed -E -n 's!.*(/mnt/[a-f0-9]+).*!\1!p')
            if [[ -z ${CANDIDATE_MOUNT_POINT} ]]; then
                >&2 echo "Could not find a mounted volume that looks like a reference disk, exiting."
                exit 1
            fi

            # Find one particular reference under the mount path. Note this is not the same reference as was specified in the
            # `inputs` section, so this would only be present if the volume we're looking at is in fact a reference disk.
            REFERENCE_FILE="Homo_sapiens.GRCh38.Nirvana.dat"
            REFERENCE_PATH=$(find ${CANDIDATE_MOUNT_POINT} -name "${REFERENCE_FILE}")
            if [[ -z ${REFERENCE_PATH} ]]; then
                >&2 echo "Could not find reference file '${REFERENCE_FILE}' under candidate reference disk mount point '${CANDIDATE_MOUNT_POINT}', exiting."
                exit 1
            fi

            # Take the parent of the parent directory of this file as root of the locally mounted  references:
            DATA_SOURCES_FOLDER="$(dirname $(dirname ${REFERENCE_PATH}))"
        else
            if [[ "~{cloud_provider}" == "azure" ]]; then
                DATA_SOURCES_FOLDER=/cromwell-executions/nirvana_references
            elif [[ "~{cloud_provider}" == "gcp" ]]; then
                DATA_SOURCES_FOLDER=/cromwell_root/nirvana_references
            else
                >&2 echo "Invalid cloud_provider value. Please specify either 'azure' or 'gcp'."
                exit 1
            fi

            mkdir ${DATA_SOURCES_FOLDER}

            # Download the references
            dotnet /Nirvana/Downloader.dll --ga GRCh38 --out ${DATA_SOURCES_FOLDER}

            # As of 2024-01-24 OMIM is no longer included among the bundle of annotation resources pulled down by the
            # Nirvana downloader. As this annotation set is currently central for our VAT logic, special-case link in
            # the OMIM .nsa bundle we downloaded back when we made the Delta reference disk:
            ln ~{omim_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/

        fi

        # Create Nirvana annotations:

        dotnet ~{nirvana_location} \
            -i ~{input_vcf} \
            -c $DATA_SOURCES_FOLDER~{path} \
            --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
            -r $DATA_SOURCES_FOLDER~{path_reference} \
            -o ~{output_annotated_file_name}

        # https://illumina.github.io/NirvanaDocumentation/introduction/parsing-json#jasix
        # Parse out the Genes section into a separate annotated json
        dotnet  ~{jasix_location} \
            --in ~{annotation_json_name} \
            --section genes \
            --out ~{gene_annotation_json_name}

        # Parse out the Positions section into a separate annotated json
        dotnet  ~{jasix_location} \
        --in ~{annotation_json_name} \
        --section positions \
        --out ~{positions_annotation_json_name}
    >>>

    runtime {
        docker: docker_path
        memory: "64 GB"
        cpu: "4"
        preemptible: 3
        maxRetries: 2
        disks: "local-disk 2000 HDD"
    }

    output {
        File genes_annotation_json = "~{gene_annotation_json_name}"
        File positions_annotation_json = "~{positions_annotation_json_name}"
    }
}
