version 1.0

task VerifyPipelineInputs {
    meta {
        description: "Verify that the pipeline input is either a ubam or pair of fastqs with additional info"
    }

    input {
        File? bam
        File? r1_fastq
        File? r2_fastq

        # fastq specific field
        String? platform
        String? library_name
        String? platform_unit
        String? read_group_name
        String? sequencing_center = "BI"

        String docker = "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
        Int cpu = 1
        Int memory_mb = 2000
        Int disk_size_gb = ceil(size(bam, "GB") + size(r1_fastq,"GB") + size(r2_fastq, "GB")) + 10
    }

    command <<<
        set -e
        python3 <<CODE

        fastq_flag = 0
        bam = "~{bam}"
        r1_fastq = "~{r1_fastq}"
        r2_fastq = "~{r2_fastq}"
        platform = "~{platform}"
        library_name = "~{library_name}"
        platform_unit = "~{platform_unit}"
        read_group_name = "~{read_group_name}"
        sequencing_center = "~{sequencing_center}"

        if bam and not r1_fastq and not r2_fastq:
            pass
        elif r1_fastq and r2_fastq and not bam:
            if platform and library_name and platform_unit and read_group_name and sequencing_center:
                fastq_flag += 1
        else:
            raise ValueError("Invalid Input. Input must be either ubam or pair of fastqs with supplemental data")

        with open("output.txt", "w") as f:
            if fastq_flag == 1:
                f.write("true")
            # Remaining case is that only bam is defined
            else:
                f.write("false")

        CODE
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks: "local-disk ${disk_size_gb} HDD"
    }

    output {
        Boolean fastq_run = read_boolean("output.txt")
    }
}

task ExtractUMIs {
    input {
        File bam
        String read1Structure
        String read2Structure

        String docker = "us.gcr.io/broad-gotc-prod/fgbio:1.0.0-1.4.0-1638817487"
        Int cpu = 4
        Int memory_mb = 8000
        Int disk_size_gb = ceil(2.2 * size(bam, "GB")) + 50
    }

    command <<<
        java -jar /usr/gitc/fgbio.jar ExtractUmisFromBam --input ~{bam} \
            --read-structure ~{read1Structure} \
            --read-structure ~{read2Structure} \
            --molecular-index-tags RX \
            --output extractUMIs.out.bam
    >>>

    runtime {
        docker : docker
        cpu : cpu
        memory : "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
        preemptible: 0
    }

    output {
        File bam_umis_extracted = "extractUMIs.out.bam"
    }
}

task STAR {
    input {
        File bam
        File starIndex

        String docker = "us.gcr.io/broad-gotc-prod/samtools-star:1.0.0-1.11-2.7.10a-1642556627"
        Int cpu = 8
        Int memory_mb = 64000
        Int disk_size_gb = ceil(2.2 * size(bam, "GB") + size(starIndex, "GB")) + 250
    }

    command <<<
        echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
        mkdir star_index
        tar -xvf ~{starIndex} -C star_index --strip-components=1

        STAR \
        --runMode alignReads \
        --runThreadN 8 \
        --genomeDir star_index \
        --outSAMtype BAM Unsorted  \
        --readFilesIn ~{bam} \
        --readFilesType SAM PE \
        --readFilesCommand samtools view -h \
        --limitSjdbInsertNsj 1200000 \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignSoftClipAtReferenceEnds Yes \
        --chimSegmentMin 15 \
        --chimMainSegmentMultNmax 1 \
        --chimOutType WithinBAM SoftClip \
        --chimOutJunctionFormat 0 \
        --twopassMode Basic \
        --quantMode TranscriptomeSAM \
        --quantTranscriptomeBan Singleend
    >>>

    runtime {
        # Note: this is 'us.gcr.io/tag-team-160914/neovax-tag-rnaseq:v1', just pulled into a location visible to warp tests
        docker : docker
        cpu : cpu
        memory : "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
        preemptible: 0
    }

    output {
        File aligned_bam = "Aligned.out.bam"
        File transcriptome_bam = "Aligned.toTranscriptome.out.bam"
    }
}

task FastqToUbam {
    input {
        File r1_fastq
        File r2_fastq
        String bam_filename
        String library_name
        String platform
        String platform_unit
        String read_group_name
        String sequencing_center

        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.6"
        Int cpu = 1
        Int memory_mb = 4000
        Int disk_size_gb = ceil(size(r1_fastq, "GiB")*2.2 + size(r2_fastq, "GiB")*2.2) + 50
    }

    String unmapped_bam_output_name = bam_filename + ".u.bam"

    command <<<
        java -jar /usr/picard/picard.jar FastqToSam \
            SORT_ORDER=unsorted \
            F1=~{r1_fastq}\
            F2=~{r2_fastq} \
            SM="~{bam_filename}" \
            LB="~{library_name}" \
            PL="~{platform}" \
            PU="~{platform_unit}" \
            RG="~{read_group_name}" \
            CN="~{sequencing_center}" \
            O="~{unmapped_bam_output_name}"
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File unmapped_bam = unmapped_bam_output_name
    }
}

task CopyReadGroupsToHeader {
    input {
        File bam_with_readgroups
        File bam_without_readgroups

        String docker = "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1.11-1624651616"
        Int cpu = 1
        Int memory_mb = 8000
        Int disk_size_gb = ceil(2.0 * size([bam_with_readgroups, bam_without_readgroups], "GB")) + 10
    }

    String basename = basename(bam_without_readgroups)

    command <<<
        samtools view -H ~{bam_without_readgroups} > header.sam
        samtools view -H ~{bam_with_readgroups} | grep "@RG" >> header.sam
        samtools reheader header.sam ~{bam_without_readgroups} > ~{basename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File output_bam = basename
    }
}

task GetSampleName {
    input {
        File bam

        String docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0"
        Int cpu = 1
        Int memory_mb = 1000
        Int disk_size_gb = 100

    }

    parameter_meta {
        bam : {
            localization_optional : true
        }
    }

    command <<<
        gatk GetSampleName -I ~{bam} -O sample_name.txt
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        String sample_name = read_string("sample_name.txt")
    }
}

task rnaseqc2 {
    input {
        File bam_file
        File genes_gtf
        String sample_id
        File exon_bed

        String docker =  "us.gcr.io/broad-dsde-methods/ckachulis/rnaseqc:2.4.2"
        Int cpu = 1
        Int memory_mb = 10000
        Int disk_size_gb = ceil(size(bam_file, 'GB') + size(genes_gtf, 'GB')) + 100
    }


    command <<<
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC 2")
        rnaseqc ~{genes_gtf} ~{bam_file} . -s ~{sample_id} -v --bed ~{exon_bed}
        echo "  * compressing outputs"
        gzip *.gct
        echo $(date +"[%b %d %H:%M:%S] done")
    >>>

    output {
        File gene_tpm = "${sample_id}.gene_tpm.gct.gz"
        File gene_counts = "${sample_id}.gene_reads.gct.gz"
        File exon_counts = "${sample_id}.exon_reads.gct.gz"
        File fragment_size_histogram = "${sample_id}.fragmentSizes.txt"
        File metrics = "${sample_id}.metrics.tsv"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }
}

task CollectRNASeqMetrics {
    input {
        File input_bam
        File input_bam_index
        String output_bam_prefix
        File ref_dict
        File ref_fasta
        File ref_fasta_index
        File ref_flat
        File ribosomal_intervals

        String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.6"
        Int cpu = 1
        Int memory_mb = 10000
        Int disk_size_gb = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")) + 20
    }

    command <<<
        java -Xms5000m -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
        REF_FLAT=~{ref_flat} \
        RIBOSOMAL_INTERVALS= ~{ribosomal_intervals} \
        STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
        INPUT=~{input_bam} \
        OUTPUT=~{output_bam_prefix}.rna_metrics
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File rna_metrics = output_bam_prefix + ".rna_metrics"
    }
}

task CollectMultipleMetrics {
    input {
        File input_bam
        File input_bam_index
        String output_bam_prefix
        File ref_dict
        File ref_fasta
        File ref_fasta_index

        String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.6"
        Int cpu = 1
        Int memory_mb = 10000
        Int disk_size_gb = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")) + 20
    }

    command <<<
        java -Xms5000m -jar /usr/picard/picard.jar CollectMultipleMetrics \
        INPUT=~{input_bam} \
        OUTPUT=~{output_bam_prefix} \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=CollectAlignmentSummaryMetrics \
        REFERENCE_SEQUENCE=~{ref_fasta}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File alignment_summary_metrics = output_bam_prefix + ".alignment_summary_metrics"
        File insert_size_metrics = output_bam_prefix + ".insert_size_metrics"
        File insert_size_histogram = output_bam_prefix + ".insert_size_histogram.pdf"
        File base_distribution_by_cycle_metrics = output_bam_prefix + ".base_distribution_by_cycle_metrics"
        File base_distribution_by_cycle_pdf = output_bam_prefix + ".base_distribution_by_cycle.pdf"
        File quality_by_cycle_metrics = output_bam_prefix + ".quality_by_cycle_metrics"
        File quality_by_cycle_pdf = output_bam_prefix + ".quality_by_cycle.pdf"
        File quality_distribution_metrics = output_bam_prefix + ".quality_distribution_metrics"
        File quality_distribution_pdf = output_bam_prefix + ".quality_distribution.pdf"
    }
}

task MergeMetrics {
    input {
        File alignment_summary_metrics
        File insert_size_metrics
        File picard_rna_metrics
        File duplicate_metrics
        File rnaseqc2_metrics
        File? fingerprint_summary_metrics
        String output_basename

        String docker =  "python:3.8-slim"
        Int cpu = 1
        Int memory_mb = 3000
        Int disk_size_gb = 10
    }

    String out_filename = output_basename + ".unified_metrics.txt"

    command <<<

        #
        # Script transpose a two line TSV
        #
        cat <<-'EOF' > transpose.py
        import csv, sys

        rows = list(csv.reader(sys.stdin, delimiter='\t'))

        for col in range(0, len(rows[0])):
            key = rows[0][col].lower()
            print(f"{key}\t{rows[1][col]}")
        EOF

        #
        # Script clean the keys, replacing space, dash and forward-slash with underscores,
        # and removing comma, single quote and periods
        #
        cat <<-'EOF' > clean.py
        import sys

        for line in sys.stdin:
            (k,v) = line.strip().lower().split("\t")
            transtable = k.maketrans({' ':'_', '-':'_', '/':'_', ',':None, '\'':None, '.' : None})
            print(f"{k.translate(transtable)}\t{v}")
        EOF

                # Process each metric file, transposing and cleaning if necessary, and pre-pending a source to the metric name

                echo "Processing Alignment Summary Metrics - Only PAIR line"
                cat ~{alignment_summary_metrics} | egrep "(CATEGORY|^PAIR)" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "picard_" $0}' >> ~{out_filename}

                echo "Processing Insert Size Metrics - removing various WIDTH metrics"
                cat ~{insert_size_metrics} | grep -A 1 "MEDIAN_INSERT_SIZE" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP|WIDTH)" | awk '{print "picard_" $0}' >> ~{out_filename}

                echo "Processing Picard RNA Metrics"
                cat ~{picard_rna_metrics} | grep -A 1 "RIBOSOMAL_BASES" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "picard_rna_metrics_" $0}' >> ~{out_filename}

                echo "Processing Duplicate Metrics"
                cat ~{duplicate_metrics} | grep -A 1 "READ_PAIR_DUPLICATES" | python transpose.py | awk '{print "picard_" $0}' >> ~{out_filename}

                echo "Processing RNASeQC2 Metrics"
                cat ~{rnaseqc2_metrics} | python clean.py | awk '{print "rnaseqc2_" $0}' >> ~{out_filename}

                if [[ -f "~{fingerprint_summary_metrics}" ]];
                then
                    echo "Processing Fingerprint Summary Metrics - only extracting LOD_EXPECTED_SAMPLE"
                    cat ~{fingerprint_summary_metrics} | grep -A 1 "LOD_EXPECTED_SAMPLE" | python transpose.py | grep -i "LOD_EXPECTED_SAMPLE" | awk '{print "fp_"$0}' >> ~{out_filename}
                else
                    echo "No Fingerprint Summary Metrics found."
                    echo "fp_lod_expected_sample	" >> ~{out_filename}
                fi    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File unified_metrics = out_filename
    }
}

task SortSamByCoordinate {
    input {
        File input_bam
        String output_bam_basename

        # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
        # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
        Float sort_sam_disk_multiplier = 4.0

        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.6"
        Int cpu = 1
        Int memory_mb = 7500
        Int disk_size_gb = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20
    }

    command <<<
        java -Xms6500m -Xmx7000m  -jar /usr/picard/picard.jar \
        SortSam \
        INPUT=~{input_bam} \
        OUTPUT=~{output_bam_basename}.bam \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        CREATE_MD5_FILE=true \
        MAX_RECORDS_IN_RAM=300000
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File output_bam = "~{output_bam_basename}.bam"
        File output_bam_index = "~{output_bam_basename}.bai"
        File output_bam_md5 = "~{output_bam_basename}.bam.md5"
    }
}

task SortSamByQueryName {
    input {
        File input_bam
        String output_bam_basename

        # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
        # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
        Float sort_sam_disk_multiplier = 4.0

        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.6"
        Int cpu = 1
        Int memory_mb = 7500
        Int disk_size_gb = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20
    }

    command <<<
        java -Xms6500m -Xmx7000m  -jar /usr/picard/picard.jar \
        SortSam \
        INPUT=~{input_bam} \
        OUTPUT=~{output_bam_basename}.bam \
        SORT_ORDER="queryname" \
        CREATE_MD5_FILE=true \
        MAX_RECORDS_IN_RAM=300000
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File output_bam = "~{output_bam_basename}.bam"
        File output_bam_md5 = "~{output_bam_basename}.bam.md5"
    }
}

task GroupByUMIs {
    input {
        File bam
        File bam_index
        String output_bam_basename

        String docker = "us.gcr.io/broad-gotc-prod/umi_tools:1.0.0-1.1.1-1638821470"
        Int cpu = 8
        Int memory_mb = 64000
        Int disk_size_gb = ceil(2.2 * size(bam, "GB")) + 300
    }

    command <<<
        umi_tools group -I ~{bam} --paired --no-sort-output --output-bam --stdout ~{output_bam_basename}.bam --umi-tag-delimiter "-" \
        --extract-umi-method tag --umi-tag RX --unmapped-reads use
    >>>

    output {
        File grouped_bam = "~{output_bam_basename}.bam"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }
}

task MarkDuplicatesUMIAware {
    input {
        File bam
        String output_basename

        String docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
        Int cpu = 1
        Int memory_mb = 16000
        Int disk_size_gb = ceil(3 * size(bam, "GB")) + 128
    }

    String output_bam_basename = output_basename + ".duplicate_marked"

    command <<<
        gatk MarkDuplicates -I ~{bam} --READ_ONE_BARCODE_TAG BX -O ~{output_bam_basename}.bam --METRICS_FILE ~{output_basename}.duplicate.metrics --ASSUME_SORT_ORDER queryname
    >>>

    output {
        File duplicate_marked_bam = "~{output_bam_basename}.bam"
        File duplicate_metrics = "~{output_basename}.duplicate.metrics"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }
}


task formatPipelineOutputs {
  input {
    String output_basename
    String transcriptome_bam
    String transcriptome_bam_index
    String transcriptome_duplicate_metrics
    String output_bam
    String output_bam_index
    String duplicate_metrics
    String rnaseqc2_gene_tpm
    String rnaseqc2_gene_counts
    String rnaseqc2_exon_counts
    String rnaseqc2_fragment_size_histogram
    String rnaseqc2_metrics
    String picard_rna_metrics
    String picard_alignment_summary_metrics
    String picard_insert_size_metrics
    String picard_insert_size_histogram
    String picard_base_distribution_by_cycle_metrics
    String picard_base_distribution_by_cycle_pdf
    String picard_quality_by_cycle_metrics
    String picard_quality_by_cycle_pdf
    String picard_quality_distribution_metrics
    String picard_quality_distribution_pdf
    String? picard_fingerprint_summary_metrics
    String? picard_fingerprint_detail_metrics
    File unified_metrics
  }

  String outputs_json_file_name = "outputs_to_TDR_~{output_basename}.json"

  command <<<
        python3 << CODE
        import json

        outputs_dict = {}

        # write normal outputs, adjusting field names as needed for TDR schema
        outputs_dict["transcriptome_bam"]="~{transcriptome_bam}"
        outputs_dict["transcriptome_bam_index"]="~{transcriptome_bam_index}"
        outputs_dict["transcriptome_duplicate_metrics_file"]="~{transcriptome_duplicate_metrics}"
        outputs_dict["genome_bam"]="~{output_bam}"
        outputs_dict["genome_bam_index"]="~{output_bam_index}"
        outputs_dict["picard_duplicate_metrics_file"]="~{duplicate_metrics}"
        outputs_dict["rnaseqc2_gene_tpm_file"]="~{rnaseqc2_gene_tpm}"
        outputs_dict["rnaseqc2_gene_counts_file"]="~{rnaseqc2_gene_counts}"
        outputs_dict["rnaseqc2_exon_counts_file"]="~{rnaseqc2_exon_counts}"
        outputs_dict["rnaseqc2_fragment_size_histogram_file"]="~{rnaseqc2_fragment_size_histogram}"
        outputs_dict["rnaseqc2_metrics_file"]="~{rnaseqc2_metrics}"
        outputs_dict["picard_rna_metrics_file"]="~{picard_rna_metrics}"
        outputs_dict["picard_alignment_summary_metrics_file"]="~{picard_alignment_summary_metrics}"
        outputs_dict["picard_insert_size_metrics_file"]="~{picard_insert_size_metrics}"
        outputs_dict["picard_insert_size_histogram_file"]="~{picard_insert_size_histogram}"
        outputs_dict["picard_base_distribution_by_cycle_metrics_file"]="~{picard_base_distribution_by_cycle_metrics}"
        outputs_dict["picard_base_distribution_by_cycle_pdf_file"]="~{picard_base_distribution_by_cycle_pdf}"
        outputs_dict["picard_quality_by_cycle_metrics_file"]="~{picard_quality_by_cycle_metrics}"
        outputs_dict["picard_quality_by_cycle_pdf_file"]="~{picard_quality_by_cycle_pdf}"
        outputs_dict["picard_quality_distribution_metrics_file"]="~{picard_quality_distribution_metrics}"
        outputs_dict["picard_quality_distribution_pdf_file"]="~{picard_quality_distribution_pdf}"
        outputs_dict["fp_summary_metrics_file"]="~{picard_fingerprint_summary_metrics}"
        outputs_dict["fp_detail_metrics_file"]="~{picard_fingerprint_detail_metrics}"

        # explode unified metrics file
        with open("~{unified_metrics}", "r") as infile:
            for row in infile:
                key, value = row.rstrip("\n").split("\t")
                outputs_dict[key] = value

        # write full outputs to file
        with open("~{outputs_json_file_name}", 'w') as outputs_file:
            for key, value in outputs_dict.items():
                if value == "None":
                    outputs_dict[key] = None
            outputs_file.write(json.dumps(outputs_dict))
            outputs_file.write("\n")
        CODE
    >>>

  runtime {
      docker: "broadinstitute/horsefish:twisttcap_scripts"
  }

  output {
      File pipeline_outputs_json = outputs_json_file_name
  }
}

task updateOutputsInTDR {
  input {
    String staging_bucket
    String tdr_dataset_uuid
    String tdr_gcp_project_for_query
    File outputs_json
    String sample_id
  }

  String tdr_target_table = "sample"

  command <<<
        python3 << CODE
        import datetime
        import json
        import requests
        import pandas as pd
        from google.cloud import bigquery
        from google.cloud import storage as gcs
        from oauth2client.client import GoogleCredentials
        from pprint import pprint
        from time import sleep

        # populate variables from inputs
        bucket = "~{staging_bucket}".replace("gs://","")
        dataset_id = "~{tdr_dataset_uuid}"
        target_table = "~{tdr_target_table}"
        outputs_json = "~{outputs_json}"
        sample_id = "~{sample_id}"
        gcp_project_for_query = "~{tdr_gcp_project_for_query}"

        # define some utils functions
        def get_access_token():
            """Get access token."""
            scopes = ["https://www.googleapis.com/auth/userinfo.profile", "https://www.googleapis.com/auth/userinfo.email"]
            credentials = GoogleCredentials.get_application_default()
            credentials = credentials.create_scoped(scopes)
            return credentials.get_access_token().access_token

        def get_headers(request_type='get'):
            headers = {"Authorization": "Bearer " + get_access_token(),
                      "accept": "application/json"}
            if request_type == 'post':
                headers["Content-Type"] = "application/json"
            return headers

        def write_file_to_bucket(filename, bucket):
            dir = "tdr"
            control_file_destination = f"{bucket}/{dir}"
            storage_client = gcs.Client()
            dest_bucket = storage_client.get_bucket(bucket)
            blob = dest_bucket.blob(f"{dir}/{filename}")
            blob.upload_from_filename(filename)
            control_file_full_path = f"gs://{bucket}/{dir}/{filename}"
            print(f"Successfully copied {loading_json_filename} to {control_file_full_path}.")
            return control_file_full_path

        def wait_for_job_status_and_result(job_id, wait_sec=10):
            # first check job status
            uri = f"https://data.terra.bio/api/repository/v1/jobs/{job_id}"

            headers = get_headers()
            response = requests.get(uri, headers=headers)
            status_code = response.status_code

            while status_code == 202:
                print(f"job running. checking again in {wait_sec} seconds")
                sleep(wait_sec)
                response = requests.get(uri, headers=headers)
                status_code = response.status_code

            if status_code != 200:
                print(f"error retrieving status for job_id {job_id}")
                return "internal error", response.text

            job_status = response.json()['job_status']
            print(f'job_id {job_id} has status {job_status}')
            # if job status = done, check job result
            if job_status in ['succeeded', 'failed']:
                result_uri = uri + "/result"
                print(f'retrieving job result from {result_uri}')
                response = requests.get(result_uri, headers=get_headers())

            return job_status, response.json()


        # read workflow outputs from file
        print(f"reading data from outputs_json file {outputs_json}")
        with open(outputs_json, "r") as infile:
            outputs_to_load = json.load(infile)

        # recode any paths (files) for TDR ingest
        print("recoding paths for TDR ingest")
        for k in outputs_to_load.keys():
            v = outputs_to_load[k]
            if v is not None and "gs://" in v:
                outputs_to_load[k] = {
                    "sourcePath": v,
                    "targetPath": v.replace("gs://", "/")
                }

        # get BQ access info for TDR dataset
        print("retrieving BQ access info for TDR dataset")
        uri = f"https://data.terra.bio/api/repository/v1/datasets/{dataset_id}?include=ACCESS_INFORMATION"
        response = requests.get(uri, headers=get_headers())
        tables = response.json()['accessInformation']['bigQuery']['tables']
        dataset_table_fq = None  # fq = fully qualified name, i.e. project.dataset.table
        for table_info in tables:
            if table_info['name'] == target_table:
                dataset_table_fq = table_info['qualifiedName']

        # retrieve data for this sample
        print(f"retrieving data for sample_id {sample_id} from {dataset_table_fq}")
        bq = bigquery.Client(gcp_project_for_query)
        query = f"SELECT * FROM \`{dataset_table_fq}\` WHERE sample_id = '{sample_id}'"
        print("using query:" + query)

        executed_query = bq.query(query)
        results = executed_query.result()

        # this avoids the pyarrow error that arises if we use `df_result = result.to_dataframe()`
        df = results.to_dataframe_iterable()
        reader = next(df)
        df_result = pd.DataFrame(reader)

        # break if there's more than one row in TDR for this sample
        print(f"retrieved {len(df_result)} samples matching sample_id {sample_id}")
        assert(len(df_result) == 1)

        # format to a dictionary
        print("formatting results to dictionary")
        input_data_list = []
        for row_id in df_result.index:
            row_dict = {}
            for col in df_result.columns:
                if isinstance(df_result[col][row_id], pd._libs.tslibs.nattype.NaTType):
                    value = None
                elif isinstance(df_result[col][row_id], pd._libs.tslibs.timestamps.Timestamp):
                    print(f'processing timestamp. value pre-formatting: {df_result[col][row_id]}')
                    formatted_timestamp = df_result[col][row_id].strftime('%Y-%m-%dT%H:%M:%S')
                    print(f'value post-formatting: {formatted_timestamp}')
                    value = formatted_timestamp
                else:
                    value = df_result[col][row_id]
                if value is not None:  # don't include empty values
                    row_dict[col] = value
            input_data_list.append(row_dict)

        sample_data_dict = input_data_list[0]

        # update original sample data with workflow outputs
        sample_data_dict.update(outputs_to_load)
        # remove and store datarepo_row_id
        old_datarepo_row_id = sample_data_dict.pop('datarepo_row_id')
        # update version_timestamp field
        new_version_timestamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%S')
        sample_data_dict['version_timestamp'] = new_version_timestamp

        # write update json to disk and upload to staging bucket
        loading_json_filename = f"{sample_id}_{new_version_timestamp}_recoded_ingestDataset.json"
        with open(loading_json_filename, 'w') as outfile:
            outfile.write(json.dumps(sample_data_dict))
            outfile.write("\n")
        load_file_full_path = write_file_to_bucket(loading_json_filename, bucket)

        # ingest data to TDR
        load_json = json.dumps({"format": "json",
                            "path": load_file_full_path,
                            "table": target_table,
                            "resolve_existing_files": True,
                            })
        uri = f"https://data.terra.bio/api/repository/v1/datasets/{dataset_id}/ingest"
        response = requests.post(uri, headers=get_headers('post'), data=load_json)
        load_job_id = response.json()['id']
        job_status, job_info = wait_for_job_status_and_result(load_job_id)
        if job_status != "succeeded":
            print(f"job status {job_status}:")
            print(job_info)

        # soft delete old row
        print("beginning soft delete")
        soft_delete_data = json.dumps({
              "deleteType": "soft", 
              "specType": "jsonArray",
              "tables": [
                {"jsonArraySpec": {"rowIds": [old_datarepo_row_id]},
                 "tableName": target_table}
              ]})
        uri = f"https://data.terra.bio/api/repository/v1/datasets/{dataset_id}/deletes"
        response = requests.post(uri, headers=get_headers('post'), data=soft_delete_data)

        print("probing soft delete job status")
        if "id" not in response.json():
            pprint(response.text)
        else:
            sd_job_id = response.json()['id']

        job_status, job_info = wait_for_job_status_and_result(sd_job_id)
        if job_status != "succeeded":
            print(f"job status {job_status}:")
            print(job_info)

        CODE
    >>>

    runtime {
        docker: "broadinstitute/horsefish:twisttcap_scripts"
    }

    output {
        File ingest_logs = stdout()
    }
}