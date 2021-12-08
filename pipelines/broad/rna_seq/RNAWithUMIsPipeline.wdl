version 1.0

import "UMIAwareDuplicateMarking.wdl" as UmiMD

workflow RNAWithUMIsPipeline {
	input {
		File bam
		String read1Structure
		String read2Structure
		File starIndex
		String output_basename
		File gtf

		File ref
		File refIndex
		File refDict
		File refFlat
		File ribosomalIntervals
		File exonBedFile
	}

	call ExtractUMIs {
		input:
			bam = bam,
			read1Structure = read1Structure,
			read2Structure = read2Structure
	}

	call STAR {
		input:
			bam = ExtractUMIs.bam_umis_extracted,
			starIndex = starIndex
	}

	call CopyReadGroupsToHeader {
		input:
			bam_with_readgroups = STAR.aligned_bam,
			bam_without_readgroups = STAR.transcriptome_bam
	}

	call UmiMD.UMIAwareDuplicateMarking {
		input:
			aligned_bam = STAR.aligned_bam,
			output_basename = output_basename
	}

	call UmiMD.UMIAwareDuplicateMarking as UMIAwareDuplicateMarkingTranscriptome {
		input:
			aligned_bam = CopyReadGroupsToHeader.output_bam,
			output_basename = output_basename + ".transcriptome"
	}

	### PLACEHOLDER for CROSSCHECK ###

	call GetSampleName {
		input:
			bam = bam
	}

	call rnaseqc2 {
		input:
			bam_file = UMIAwareDuplicateMarking.duplicate_marked_bam,
			genes_gtf = gtf,
			sample_id = GetSampleName.sample_name,
			exon_bed = exonBedFile
	}

	call CollectRNASeqMetrics {
		input:
			input_bam=UMIAwareDuplicateMarking.duplicate_marked_bam,
			input_bam_index=UMIAwareDuplicateMarking.duplicate_marked_bam_index,
			output_bam_prefix=GetSampleName.sample_name,
			ref_dict=refDict,
			ref_fasta=ref,
			ref_fasta_index=refIndex,
			ref_flat=refFlat,
			ribosomal_intervals=ribosomalIntervals,
			preemptible_tries=0
	}

	call CollectMultipleMetrics {
		input:
			input_bam=UMIAwareDuplicateMarking.duplicate_marked_bam,
			input_bam_index=UMIAwareDuplicateMarking.duplicate_marked_bam_index,
			output_bam_prefix=GetSampleName.sample_name,
			ref_dict=refDict,
			ref_fasta=ref,
			ref_fasta_index=refIndex,
			preemptible_tries=0
	}

	# TODO: wire in fingerprint_summary_metrics once we have it. Using static example for now
	call MergeMetrics {
		input:
			alignment_summary_metrics=CollectMultipleMetrics.alignment_summary_metrics,
			insert_size_metrics=CollectMultipleMetrics.insert_size_metrics,
			picard_rna_metrics=CollectRNASeqMetrics.rna_metrics,
			duplicate_metrics=UMIAwareDuplicateMarking.duplicate_metrics,
			rnaseqc2_metrics=rnaseqc2.metrics,
			fingerprint_summary_metrics="gs://broad-dsp-spec-ops/kcibul/rnaseq/metrics/example.duplicate.metrics"
	}


	output {
		File transcriptome_bam = UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam
		File transcriptome_bam_index = UMIAwareDuplicateMarkingTranscriptome.duplicate_marked_bam_index
		File transcriptome_duplicate_metrics = UMIAwareDuplicateMarkingTranscriptome.duplicate_metrics
		File output_bam = UMIAwareDuplicateMarking.duplicate_marked_bam
		File output_bam_index = UMIAwareDuplicateMarking.duplicate_marked_bam_index
		File duplicate_metrics = UMIAwareDuplicateMarking.duplicate_metrics
		File rnaseqc2_gene_tpm = rnaseqc2.gene_tpm
		File rnaseqc2_gene_counts = rnaseqc2.gene_counts
		File rnaseqc2_exon_counts = rnaseqc2.exon_counts
		File rnaseqc2_fragment_size_histogram = rnaseqc2.fragment_size_histogram
		File rnaseqc2_metrics = rnaseqc2.metrics
		File picard_rna_metrics = CollectRNASeqMetrics.rna_metrics
		File picard_alignment_summary_metrics = CollectMultipleMetrics.alignment_summary_metrics
		File picard_insert_size_metrics = CollectMultipleMetrics.insert_size_metrics
		File picard_insert_size_histogram = CollectMultipleMetrics.insert_size_histogram
		File picard_base_distribution_by_cycle_metrics = CollectMultipleMetrics.base_distribution_by_cycle_metrics
		File picard_base_distribution_by_cycle_pdf = CollectMultipleMetrics.base_distribution_by_cycle_pdf
		File picard_quality_by_cycle_metrics = CollectMultipleMetrics.quality_by_cycle_metrics
		File picard_quality_by_cycle_pdf = CollectMultipleMetrics.quality_by_cycle_pdf
		File picard_quality_distribution_metrics = CollectMultipleMetrics.quality_distribution_metrics
		File picard_quality_distribution_pdf = CollectMultipleMetrics.quality_distribution_pdf
		File unified_metrics = MergeMetrics.unified_metrics
	}
}

task STAR {
	input {
		File bam
		File starIndex
	}
	Int disk_space = ceil(2.2 * size(bam, "GB") + size(starIndex, "GB")) + 250

	command <<<
		echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
		mkdir star_index
		tar -xvf ~{starIndex} -C star_index --strip-components=1

		STAR --readFilesIn ~{bam} --readFilesType SAM PE --readFilesCommand samtools view -h \
			--runMode alignReads --genomeDir star_index --outSAMtype BAM Unsorted --runThreadN 8 \
			--limitSjdbInsertNsj 1200000 --outSAMstrandField intronMotif --outSAMunmapped Within \
			--outFilterType BySJout --outFilterMultimapNmax 20 --outFilterScoreMinOverLread 0.33 \
			--outFilterMatchNminOverLread 0.33 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 \
			--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 --alignSoftClipAtReferenceEnds Yes --chimSegmentMin 15 --chimMainSegmentMultNmax 1 \
			--chimOutType WithinBAM SoftClip --chimOutJunctionFormat 0 --twopassMode Basic --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend

	>>>

	runtime {
		# Copied the docker from tag's private location to a public location.
#		docker : "us.gcr.io/tag-team-160914/neovax-tag-rnaseq:v1"
		docker : "us.gcr.io/broad-gotc-prod/neovax-tag-rnaseq:v1"
		disks : "local-disk " + disk_space + " HDD"
		memory : "64GB"
		cpu : "8"
		preemptible: 0
	}

	output {
		File aligned_bam = "Aligned.out.bam"
		File transcriptome_bam = "Aligned.toTranscriptome.out.bam"
	}
}

task ExtractUMIs {
	input {
		File bam
		String read1Structure
		String read2Structure
	}

	Int disk_space = ceil(2.2 * size(bam, "GB")) + 50

	command <<<
		fgbio ExtractUmisFromBam --input ~{bam} \
			--read-structure ~{read1Structure} \
			--read-structure ~{read2Structure} \
			--molecular-index-tags RX \
			--output extractUMIs.out.bam
	>>>

	runtime {
		docker : "quay.io/biocontainers/fgbio@sha256:a8e5cf58c318bffba3b2b694a3640ecd9e8106cee2e33b75710c0e8215138b6e"
		disks : "local-disk " + disk_space + " HDD"
		preemptible: 0
	}

	output {
		File bam_umis_extracted = "extractUMIs.out.bam"
	}
}

task rnaseqc2 {
	input {
		File bam_file
		File genes_gtf
		String sample_id
		File exon_bed
	}
	
	Int disk_space = ceil(size(bam_file, 'GB') + size(genes_gtf, 'GB')) + 100

	command {
		set -euo pipefail
		echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC 2")
		rnaseqc ~{genes_gtf} ~{bam_file} . -s ~{sample_id} -v --bed ~{exon_bed}
		echo "  * compressing outputs"
		gzip *.gct
		echo $(date +"[%b %d %H:%M:%S] done")
	}

	output {
		File gene_tpm = "${sample_id}.gene_tpm.gct.gz"
		File gene_counts = "${sample_id}.gene_reads.gct.gz"
		File exon_counts = "${sample_id}.exon_reads.gct.gz"
		File fragment_size_histogram = "${sample_id}.fragmentSizes.txt"
		File metrics = "${sample_id}.metrics.tsv"
	}

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/ckachulis/rnaseqc:2.4.2"
		memory: "10GB"
		disks: "local-disk " + disk_space + " HDD"
		preemptible: 0
	}
}

task GetSampleName {
	input {
		File bam
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
		docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
		disks: "local-disk 100 HDD"
	}

	output {
		String sample_name = read_string("sample_name.txt")
	}
}

task CopyReadGroupsToHeader {
	input {
		File bam_with_readgroups
		File bam_without_readgroups
	}

	Int disk_size = ceil(2.0 * size([bam_with_readgroups, bam_without_readgroups], "GB")) + 10
	String basename = basename(bam_without_readgroups)

	command <<<
		samtools view -H ~{bam_without_readgroups} > header.sam
		samtools view -H ~{bam_with_readgroups} | grep "@RG" >> header.sam
		samtools reheader header.sam ~{bam_without_readgroups} > ~{basename}
	>>>

	runtime {
		docker: "us.gcr.io/broad-dsde-methods/samtools@sha256:0e49b0a5d91c203b8c07f5242277c2060b4b8ea54df8b1d123f990a1ad0588b2"
		disks: "local-disk " + disk_size + " HDD"
	}

	output {
		File output_bam = basename
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
		Int preemptible_tries
		File ref_flat
		File ribosomal_intervals
	}

	Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
	Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

	# This jar skips the header check of the ribosomal interval
#	File picard_jar = "gs://broad-dsde-methods-takuto/hydro.gen/picard_ignore_ribosomal_header.jar"
	# copied the jar into a gotc-accessible location. TODO - this presumably won't work in Terra?
	File picard_jar = "gs://broad-gotc-test-storage/rna_seq/picard_ignore_ribosomal_header.jar"

	command {
		java -Xms5000m -jar ~{picard_jar} CollectRnaSeqMetrics \
		REF_FLAT=~{ref_flat} \
		RIBOSOMAL_INTERVALS= ~{ribosomal_intervals} \
		STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
		INPUT=~{input_bam} \
		OUTPUT=~{output_bam_prefix}.rna_metrics
	}

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
		memory: "7 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
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
		Int preemptible_tries
	}

	Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
	Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

	command {
		java -Xms5000m -jar /usr/picard/picard.jar CollectMultipleMetrics \
		INPUT=~{input_bam} \
		OUTPUT=~{output_bam_prefix} \
		PROGRAM=CollectInsertSizeMetrics \
		PROGRAM=CollectAlignmentSummaryMetrics \
		REFERENCE_SEQUENCE=~{ref_fasta}
	}

	runtime {
		docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.6"
		memory: "7 GiB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
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
    }
    
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
        cat ~{alignment_summary_metrics} | egrep "(CATEGORY|^PAIR)" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "picard_" $0}' >> unified_metrics.txt

        echo "Processing Insert Size Metrics - removing various WIDTH metrics"
        cat ~{insert_size_metrics} | grep -A 1 "MEDIAN_INSERT_SIZE" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP|WIDTH)" | awk '{print "picard_" $0}' >> unified_metrics.txt

        echo "Processing Picard RNA Metrics"
        cat ~{picard_rna_metrics} | grep -A 1 "RIBOSOMAL_BASES" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "picard_rna_metrics_" $0}' >> unified_metrics.txt

        echo "Processing Duplicate Metrics"
        cat ~{dupe_metrics} | grep -A 1 "READ_PAIR_DUPLICATES" | python transpose.py | awk '{print "picard_" $0}' >> unified_metrics.txt

        echo "Processing RNASeQC2 Metrics"
        cat ~{rnaseqc2_metrics} | python clean.py | awk '{print "rnaseqc2_" $0}' >> unified_metrics.txt

        echo "Processing Fingerprint Summary Metrics - only extracting LOD_EXPECTED_SAMPLE"
        cat ~{fingerprint_summary_metrics} | grep -A 1 "LOD_EXPECTED_SAMPLE" | python transpose.py | grep -i "LOD_EXPECTED_SAMPLE" | awk '{print "fp_"$0}' >> unified_metrics.txt
    >>>

    runtime {
        docker : "python:3.8-slim"
        disks : "local-disk 10 HDD"
        memory : "3GB"
        cpu : "1"
        preemptible: 0
    }

    output {
        File unified_metrics = "unified_metrics.txt"
    }
}    

