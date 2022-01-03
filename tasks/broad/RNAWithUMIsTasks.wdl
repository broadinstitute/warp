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

        String docker = "python:3.7.2"
        Int cpu = 1
        Int memory_mb = 1000
        Int disk_size_gb = 1
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

        String docker = "us.gcr.io/broad-gotc-prod/neovax-tag-rnaseq:v1"
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
        File fingerprint_summary_metrics
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

                echo "Processing Fingerprint Summary Metrics - only extracting LOD_EXPECTED_SAMPLE"
                cat ~{fingerprint_summary_metrics} | grep -A 1 "LOD_EXPECTED_SAMPLE" | python transpose.py | grep -i "LOD_EXPECTED_SAMPLE" | awk '{print "fp_"$0}' >> ~{out_filename}
    >>>

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

task SortSamUMIAware {
    input {
        File input_bam
        String output_bam_basename
        String sort_order # "queryname" or "coordinate"

        # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
        # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
        Float sort_sam_disk_multiplier = 4.0

        String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.6"
        Int cpu = 1
        Int memory_mb = 16000
        Int disk_size_gb = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 256
    }

    String extra_args = if sort_order == "coordinate" then "CREATE_INDEX=true" else ""

    command <<<
        java -Xms8192m -jar /usr/picard/picard.jar \
        SortSam \
        INPUT=~{input_bam} \
        OUTPUT=~{output_bam_basename}.bam \
        SORT_ORDER=~{sort_order} \
        CREATE_MD5_FILE=true \
        MAX_RECORDS_IN_RAM=300000 \
        ~{extra_args}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File output_bam = "~{output_bam_basename}.bam"
        File? output_bam_index = "~{output_bam_basename}.bai" # Create index for coordinate sort
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

