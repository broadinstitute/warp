version 1.0

task VerifyPipelineInputs {
  meta {
    description: "Verify that the pipeline input is either a ubam or pair of fastqs with additional info"
  }

  input {
    File? bam
    File? r1_fastq
    File? r2_fastq

    String docker = "us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian"
    Int cpu = 1
    Int memory_mb = 2000
    Int disk_size_gb = ceil(size(bam, "GiB") + size(r1_fastq,"GiB") + size(r2_fastq, "GiB")) + 10
  }

  command <<<
    set -e
    python3 <<CODE

    fastq_flag = 0
    bam = "~{bam}"
    r1_fastq = "~{r1_fastq}"
    r2_fastq = "~{r2_fastq}"

    if bam and not r1_fastq and not r2_fastq:
      pass
    elif r1_fastq and r2_fastq and not bam:
      fastq_flag += 1
    else:
      raise ValueError("Invalid Input. Input must be either ubam or a pair of fastqs")

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
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
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
    Int memory_mb = 5000
    Int disk_size_gb = ceil(2.2 * size(bam, "GiB")) + 20
  }

  command <<<
    java -jar /usr/gitc/fgbio.jar ExtractUmisFromBam \
      --input ~{bam} \
      --read-structure ~{read1Structure} \
      --read-structure ~{read2Structure} \
      --molecular-index-tags RX \
      --output extractUMIs.out.bam
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
    preemptible: 0
  }

  output {
    File bam_umis_extracted = "extractUMIs.out.bam"
  }
}

# Adapter clipping
# https://github.com/OpenGene/fastp
task Fastp {
  input {
    File fastq1
    File fastq2
    String output_prefix
    File adapter_fasta = "gs://gcp-public-data--broad-references/RNA/resources/Illumina_adapters.fasta"

    String docker = "us.gcr.io/broad-gotc-prod/fastp:1.0.0-0.20.1-1649253500"
    Int memory_mb =  ceil(1.5*size(fastq1, "MiB")) + 8192 # Experimentally determined formula for memory allocation
    Int disk_size_gb = 5*ceil(size(fastq1, "GiB")) + 128
    File monitoring_script = "gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
    Int cpu=4
  }

  command {
    bash ~{monitoring_script} > monitoring.log &

    fastp --in1 ~{fastq1} --in2 ~{fastq2} --out1 ~{output_prefix}_read1.fastq.gz --out2 ~{output_prefix}_read2.fastq.gz \
    --disable_quality_filtering \
    --disable_length_filtering \
    --adapter_fasta ~{adapter_fasta} \
    --thread ~{cpu}
  }
  

  runtime {
    docker: docker
    memory: "~{memory_mb} MiB"
    cpu: cpu
    disks: "local-disk ~{disk_size_gb} HDD"
    preemptible: 0
    maxRetries: 2
  }

  output {
    File monitoring_log = "monitoring.log"
    File fastq1_clipped = output_prefix + "_read1.fastq.gz"
    File fastq2_clipped = output_prefix + "_read2.fastq.gz"
  }

}

task STAR {
  input {
    File bam
    File starIndex

    String docker = "us.gcr.io/broad-gotc-prod/samtools-star:1.0.0-1.11-2.7.10a-1642556627"
    Int cpu = 8
    Int memory_mb = ceil((size(starIndex, "GiB")) + 10) * 1500
    Int disk_size_gb = ceil(2.2 * size(bam, "GiB") + size(starIndex, "GiB")) + 150
  }

  command <<<
    echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
    mkdir star_index
    tar -xvf ~{starIndex} -C star_index --strip-components=1

    STAR \
      --runMode alignReads \
      --runThreadN ~{cpu} \
      --genomeDir star_index \
      --outSAMtype BAM Unsorted \
      --readFilesIn ~{bam} \
      --readFilesType SAM PE \
      --readFilesCommand samtools view -h \
      --runRNGseed 12345 \
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
      --quantTranscriptomeBan IndelSoftclipSingleend \
      --alignEndsProtrude 20 ConcordantPair
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
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

    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:3.0.0"
    Int cpu = 1
    Int memory_mb = 4000
    Int disk_size_gb = ceil(size(r1_fastq, "GiB")*2.2 + size(r2_fastq, "GiB")*2.2) + 50
  }

  String unmapped_bam_output_name = bam_filename + ".u.bam"

  Int java_memory_size = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar FastqToSam \
      SORT_ORDER=unsorted \
      F1=~{r1_fastq}\
      F2=~{r2_fastq} \
      SM="~{bam_filename}" \
      LB="~{library_name}" \
      PL="~{platform}" \
      PU="~{platform_unit}" \
      RG="~{read_group_name}" \
      CN="~{sequencing_center}" \
      O="~{unmapped_bam_output_name}" \
      ALLOW_EMPTY_FASTQ=True
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
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
    Int memory_mb = 2500
    Int disk_size_gb = ceil(2.0 * size([bam_with_readgroups, bam_without_readgroups], "GiB")) + 10
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
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File output_bam = basename
  }
}

task GetSampleName {
  input {
    File bam
    String? billing_project

    String docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int cpu = 1
    Int memory_mb = 1000
    Int disk_size_gb = ceil(2.0 * size(bam, "GiB")) + 10
  }
  String requester_pays_flag = if defined(billing_project) then "--gcs-project-for-requester-pays ${billing_project}" else ""

  parameter_meta {
    bam: {
      localization_optional: true
    }
  }

  command <<<
    gatk GetSampleName -I ~{bam} -O sample_name.txt ~{requester_pays_flag}
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
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

    String docker = "us.gcr.io/broad-dsde-methods/ckachulis/rnaseqc@sha256:a4bec726bb51df5fb8c8f640f7dec22fa28c8f7803ef9994b8ec39619b4514ca"
    Int cpu = 1
    Int memory_mb = 8000
    Int disk_size_gb = ceil(size(bam_file, 'GiB') + size(genes_gtf, 'GiB') + size(exon_bed, 'GiB')) + 50
  }

  command <<<
    set -euo pipefail
    # force fragmentSizes histogram output file to exist (even if empty)
    touch ~{sample_id}.fragmentSizes.txt
    echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC 2")
    rnaseqc ~{genes_gtf} ~{bam_file} . -s ~{sample_id} -v --bed ~{exon_bed}
    echo "  * compressing outputs"
    gzip *.gct
    echo $(date +"[%b %d %H:%M:%S] done")
  >>>

  output {
    File gene_tpm = "~{sample_id}.gene_tpm.gct.gz"
    File gene_counts = "~{sample_id}.gene_reads.gct.gz"
    File exon_counts = "~{sample_id}.exon_reads.gct.gz"
    File fragment_size_histogram = "~{sample_id}.fragmentSizes.txt"
    File metrics = "~{sample_id}.metrics.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
    maxRetries: 2
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


    String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.11"
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")) + 20
  }

  Int java_memory_size = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
      REF_FLAT=~{ref_flat} \
      RIBOSOMAL_INTERVALS= ~{ribosomal_intervals} \
      STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix}.rna_metrics
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
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


    String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:3.0.0"
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")) + 20
  }

  Int java_memory_size = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    #plots will not be produced if there are no reads
    touch ~{output_bam_prefix}.insert_size_histogram.pdf
    touch ~{output_bam_prefix}.insert_size_metrics
    touch ~{output_bam_prefix}.base_distribution_by_cycle.pdf
    touch ~{output_bam_prefix}.quality_by_cycle.pdf
    touch ~{output_bam_prefix}.quality_distribution.pdf

    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar CollectMultipleMetrics \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix} \
      REFERENCE_SEQUENCE=~{ref_fasta}
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
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
      value = rows[1][col]
      if value == "?":
        value = "NaN"
      if key in ["median_insert_size", "median_absolute_deviation", "median_read_length", "mad_read_length", "pf_hq_median_mismatches"]:
        value = str(int(float(value)))
      print(f"{key}\t{value}")
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
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
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


    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.11"
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20
  }

  Int java_memory_size = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar SortSam \
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
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
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
    Float sort_sam_disk_multiplier = 6.0


    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.11"
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20
  }

  Int java_memory_size = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="queryname" \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
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

    String docker = "us.gcr.io/broad-gotc-prod/umi_tools:1.0.0-1.1.1-1690198330"
    Int cpu = 2
    Int memory_mb = 64000
    Int disk_size_gb = ceil(2.2 * size([bam, bam_index], "GiB")) + 100

    File monitoring_script = "gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
  }

  command <<<
    bash ~{monitoring_script} > monitoring.log &

    # umi_tools has a bug which lead to using the order of elements in a set to determine tie-breakers in
    # rare edge-cases.  Sets in python are unordered, so this leads to non-determinism.  Setting PYTHONHASHSEED
    # to 0 means that hashes will be unsalted.  While it is not in any way gauranteed by the language that this
    # will remove the non-determinism, in practice, due to implementation details of sets in cpython, this makes seemingly
    # deterministic behavior much more likely

    export PYTHONHASHSEED=0

    umi_tools group -I ~{bam} --paired --no-sort-output --output-bam --stdout ~{output_bam_basename}.bam --umi-tag-delimiter "-" \
    --extract-umi-method tag --umi-tag RX --unmapped-reads use
  >>>

  output {
    File grouped_bam = "~{output_bam_basename}.bam"
    File monitoring_log = "monitoring.log"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
    maxRetries: 1
  }
}

task MarkDuplicatesUMIAware {
  input {
    File bam
    String output_basename
    Boolean remove_duplicates
    Boolean use_umi

    String docker =  "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.11"
    Int cpu = 1
    Int memory_mb = 16000
    Int disk_size_gb = ceil(3 * size(bam, "GiB")) + 60
  }

  String output_bam_basename = output_basename + ".duplicate_marked"

  command <<<
    java -jar /usr/picard/picard.jar MarkDuplicates \
    INPUT=~{bam} \
    OUTPUT=~{output_bam_basename}.bam \
    METRICS_FILE=~{output_basename}.duplicate.metrics \
    REMOVE_DUPLICATES=~{remove_duplicates} \
    ~{true='READ_ONE_BARCODE_TAG=BX' false='' use_umi} \

  >>>

  output {
    File duplicate_marked_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{output_basename}.duplicate.metrics"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }
}

task formatPipelineOutputs {
  input {
    String sample_id
    String transcriptome_bam
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
    Float contamination
    Float contamination_error
    String fastqc_html_report
    String fastqc_percent_reads_with_adapter

    Int cpu = 1
    Int memory_mb = 2000
    Int disk_size_gb = 10
  }

  String outputs_json_file_name = "outputs_to_TDR_~{sample_id}.json"

  command <<<
    python3 << CODE
    import json

    outputs_dict = {}

    # NOTE: we rename some field names to match the TDR schema
    outputs_dict["sample_id"]="~{sample_id}" # primary key
    outputs_dict["transcriptome_bam"]="~{transcriptome_bam}"
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
    outputs_dict["fastqc_html_report"]="~{fastqc_html_report}"
    outputs_dict["fastqc_percent_reads_with_adapter"]="~{fastqc_percent_reads_with_adapter}"
    outputs_dict["contamination"]="~{contamination}"
    outputs_dict["contamination_error"]="~{contamination_error}"

    # explode unified metrics file
    with open("~{unified_metrics}", "r") as infile:
      for row in infile:
        key, value = row.rstrip("\n").split("\t")
        outputs_dict[key] = value

    # write full outputs to file
    with open("~{outputs_json_file_name}", 'w') as outputs_file:
      for key, value in outputs_dict.items():
        if value == "None" or value == "":
          outputs_dict[key] = None
      outputs_file.write(json.dumps(outputs_dict))
      outputs_file.write("\n")
    CODE
  >>>

  runtime {
    docker: "broadinstitute/horsefish:tdr_import_v1.4"
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File pipeline_outputs_json = outputs_json_file_name
  }
}

task updateOutputsInTDR {
  input {
    String tdr_dataset_uuid
    File outputs_json

    Int cpu = 1
    Int memory_mb = 2000
    Int disk_size_gb = 10
  }

  command <<<
    # input args:
    # -d dataset uuid
    # -t target table in dataset
    # -o json of data to ingest
    # -f field to populate with timestamp at ingest (can have multiple)
    python -u /scripts/export_pipeline_outputs_to_tdr.py \
      -d "~{tdr_dataset_uuid}" \
      -t "sample" \
      -o "~{outputs_json}" \
      -f "version_timestamp" \
      -f "analysis_end_time"
  >>>

  runtime {
    docker: "broadinstitute/horsefish:tdr_import_v1.4"
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File ingest_logs = stdout()
  }
}

# GATK CalculateContamination, adapted for RNA-seq data.
# Specifically, we disable two read filters from the default set of read filters
# for a LocusWalker:
# 1. WellformedReadFilter: This filter removes (among others) reads with N's in the CIGAR string, which is 
# common in RNA data.
# 2. MappingQualityAvailableReadFilter: This filter removes reads with MQ=255, which by SAM spec
# means mapping quality is missing. But STAR uses 255 to mean unique mapping, the equivalent of MQ60
# for other aligners.
task CalculateContamination {
  input {
    File bam
    File bam_index
    String base_name
    File ref_fasta
    File ref_dict
    File ref_fasta_index
    File population_vcf
    File population_vcf_index
    # runtime
    String docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int cpu = 1
    Int memory_mb = 8192
    Int disk_size_gb = 256
  }

  parameter_meta {
    bam: { localization_optional: true }
    bam_index: { localization_optional: true }
    ref_fasta: { localization_optional: true }
    ref_fasta_index: { localization_optional: true }
    ref_dict: { localization_optional: true }
  }

  command <<<
    set -e
    gatk --java-options "-Xmx4096m" GetPileupSummaries \
    -R ~{ref_fasta} \
    -I ~{bam} \
    -V ~{population_vcf} \
    -L ~{population_vcf} \
    -O ~{base_name}_pileups.tsv \
    --disable-read-filter WellformedReadFilter \
    --disable-read-filter MappingQualityAvailableReadFilter

    gatk --java-options "-Xmx4096m" CalculateContamination \
    -I ~{base_name}_pileups.tsv \
    -O ~{base_name}_contamination.tsv
  
    grep -v ^sample ~{base_name}_contamination.tsv | awk 'BEGIN{FS="\t"}{print($2)}' > contamination.txt 
    grep -v ^sample ~{base_name}_contamination.tsv | awk 'BEGIN{FS="\t"}{print($3)}' > contamination_error.txt
  >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
      File pileups = "~{base_name}_pileups.tsv"
      File contamination_table = "~{base_name}_contamination.tsv"
      Float contamination = read_float("contamination.txt")
      Float contamination_error = read_float("contamination_error.txt")
  }
}

task SamToFastq {
  input {
    File bam
    String output_prefix
    String docker = "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.11"
    Int memory_mb = 4096
    Int disk_size_gb = 3*ceil(size(bam, "GiB")) + 128
  }

  command {
    java -jar /usr/picard/picard.jar SamToFastq \
    I=~{bam} \
    FASTQ=~{output_prefix}_1.fastq.gz \
    SECOND_END_FASTQ=~{output_prefix}_2.fastq.gz

  }

  runtime {
    docker: docker
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File fastq1 = output_prefix + "_1.fastq.gz"
    File fastq2 = output_prefix + "_2.fastq.gz"
  }
}

task FastQC {
  input {
    File unmapped_bam
    Float? mem = 4
    String docker = "us.gcr.io/tag-public/tag-tools:1.0.0" # sato: This likely needs to be made public
    Int memory_mb = 4096
    Int disk_size_gb = 3*ceil(size(unmapped_bam, "GiB")) + 128
  }

  String bam_basename = basename(unmapped_bam, ".bam")

  command {
    perl /usr/tag/scripts/FastQC/fastqc ~{unmapped_bam} --extract -o ./
    mv ~{bam_basename}_fastqc/fastqc_data.txt ~{bam_basename}_fastqc_data.txt

    tail -n 2 ~{bam_basename}_fastqc_data.txt | head -n 1 | cut -f 2 > ~{bam_basename}_adapter_content.txt
  }
  
  runtime {
    docker: docker
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File fastqc_data = "~{bam_basename}_fastqc_data.txt"
    File fastqc_html = "~{bam_basename}_fastqc.html"
    Float fastqc_percent_reads_with_adapter = if read_string("~{bam_basename}_adapter_content.txt") == "warn" then -1 else read_float("~{bam_basename}_adapter_content.txt")
  }
}

task TransferReadTags {
  input {
    File aligned_bam
    File ubam
    String output_basename
    String docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int memory_mb = 16000
    Int disk_size_gb = ceil(2 * size(aligned_bam, "GiB")) + ceil(2 * size(ubam, "GiB")) + 128
  }

  command <<<
    gatk TransferReadTags \
    -I ~{aligned_bam} \
    --unmapped-sam ~{ubam} \
    -O ~{output_basename}.bam \
    --read-tags RX
  >>>

  output {
    File output_bam = "~{output_basename}.bam"
  }

  runtime {
    docker: docker
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }
}

task PostprocessTranscriptomeForRSEM {
  input {
    String prefix
    File input_bam # the input must be queryname sorted
    Int disk_size_gb = ceil(3*size(input_bam,"GiB")) + 128
    String docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int memory_mb = 16000
  }

  command {
    gatk PostProcessReadsForRSEM \
    -I ~{input_bam} \
    -O ~{prefix}_RSEM_post_processed.bam \
    --use-jdk-deflater
  }

  output {
    File output_bam = "~{prefix}_RSEM_post_processed.bam"
  }

  runtime {
    docker: docker 
    disks: "local-disk ~{disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
  }
}

task CreateEmptyFile {
  input {
    Int disk_size_gb = 128
    String docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int memory_mb = 4096
  }

  command {
    echo "place holder for a bam index file" > empty.txt
  }

  output {
    File empty_file = "empty.txt"
  }

  runtime {
    docker: docker 
    disks: "local-disk ~{disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
  }
}