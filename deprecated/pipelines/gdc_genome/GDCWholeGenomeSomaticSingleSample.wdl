version 1.0

# GDCWholeGenomeSomaticSingleSample is now deprecated 2025-03-06

import "../../../../../../../pipelines/wdl/reprocessing/cram_to_unmapped_bams/CramToUnmappedBams.wdl" as ToUbams
import "../../../../../../../tasks/wdl/CheckContaminationSomatic.wdl" as CheckContamination

struct FastqPairRecord {
    File forward_fastq
    File reverse_fastq
    String readgroup
    String readgroup_id
}

struct FastqSingleRecord {
    File fastq
    String readgroup
    String readgroup_id
}

task bam_readgroup_to_contents {
    input {
        File bam
        Int preemptible = 10
        Int max_retries = 0
        Int cpu = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }
    Int mem = ceil(size(bam, "MiB")) + 2000 + additional_memory_mb
    Int disk_space = ceil(size(bam, "GiB")) + 10 + additional_disk_gb

    parameter_meta {
        bam: {localization_optional: true}
    }

    command <<<
    set -euo pipefail
    export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
    samtools view -H ~{bam} | \
        awk 'BEGIN {                                                                \
                OFS = FS = "\t";                                                    \
                header = "ID\tBC\tCN\tDS\tDT\tFO\tKS\tLB\tPG\tPI\tPL\tPM\tPU\tSM";  \
                split(header, header_ary, "\t");                                    \
                for (i=1; i in header_ary; i++) {                                   \
                    header_pos[header_ary[i]] = i                                   \
                };                                                                  \
                print header                                                        \
             }                                                                      \
             /^@RG/ {                                                               \
                for (i=2; i<=NF; i++) {                                             \
                    split($i, rg, ":");                                             \
                    row_ary[header_pos[rg[1]]] = rg[2];                             \
                };                                                                  \
                row = row_ary[1];                                                   \
                for (i=2; i in header_ary; i++) {                                   \
                    row = row "\t";                                                 \
                    if (i in row_ary)                                               \
                        row = row row_ary[i];                                       \
                };                                                                  \
                delete row_ary;                                                     \
                print row                                                           \
             }'
    >>>

    output {
        Array[Object] readgroups = read_objects(stdout())
    }
    
    runtime {
        docker: "broadgdac/samtools:1.10"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task biobambam_bamtofastq {
    input {
        Int collate = 1
        String exclude = "QCFAIL,SECONDARY,SUPPLEMENTARY"
        File filename
        Int gz = 1
        String inputformat = "bam"
        Int level = 5
        String outputdir = "."
        Int outputperreadgroup = 1
        String outputperreadgroupsuffixF = "_1.fq.gz"
        String outputperreadgroupsuffixF2 = "_2.fq.gz"
        String outputperreadgroupsuffixO = "_o1.fq.gz"
        String outputperreadgroupsuffixO2 = "_o2.fq.gz"
        String outputperreadgroupsuffixS = "_s.fq.gz"
        Int tryoq = 1
        String T = "tempfq"
        Int cpu = 1
        Int preemptible = 2
        Int max_retries = 0
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }
    Int mem = ceil(size(filename, "MiB")) + 2000 + additional_memory_mb
    Int disk_space = ceil(size(filename, "GiB") * 2) + 10 + additional_disk_gb

    command {
        set -euo pipefail
        /usr/local/bin/bamtofastq \
            T=~{T} \
            collate=~{collate} \
            exclude=~{exclude} \
            filename=~{filename} \
            gz=~{gz} \
            inputformat=~{inputformat} \
            level=~{level} \
            outputdir=~{outputdir} \
            outputperreadgroup=~{outputperreadgroup} \
            outputperreadgroupsuffixF=~{outputperreadgroupsuffixF} \
            outputperreadgroupsuffixF2=~{outputperreadgroupsuffixF2} \
            outputperreadgroupsuffixO=~{outputperreadgroupsuffixO} \
            outputperreadgroupsuffixO2=~{outputperreadgroupsuffixO2} \
            outputperreadgroupsuffixS=~{outputperreadgroupsuffixS} \
            tryoq=~{tryoq}
    }

    output {
        Array[File] output_fastq1 = glob("*~{outputperreadgroupsuffixF}")
        Array[File] output_fastq2 = glob("*~{outputperreadgroupsuffixF2}")
        Array[File] output_fastq_o1 = glob("*~{outputperreadgroupsuffixO}")
        Array[File] output_fastq_o2 = glob("*~{outputperreadgroupsuffixO2}")
        Array[File] output_fastq_s = glob("*~{outputperreadgroupsuffixS}")
    }
    
    runtime {
        docker: "broadgdac/biobambam2:2.0.87-release-20180301132713"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
        
task emit_pe_records {
    input {
        Array[File]+ fastq1_files
        Array[File]+ fastq2_files
        Array[Object]+ readgroups
        String fastq1_suffix = "_1.fq.gz"
        Int preemptible = 10
        Int max_retries = 0
        Int cpu = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }
    Int mem = ceil(size(fastq1_files, "MiB") + size(fastq2_files, "MiB")) + 2000 + additional_memory_mb
    Int disk_space = ceil(size(fastq1_files, "GiB") + size(fastq2_files, "GiB")) + 10 + additional_disk_gb

    File readgroups_tsv = write_objects(readgroups)

    parameter_meta {
        fastq1_files: {localization_optional: true}
        fastq2_files: {localization_optional: true}
    }

    command <<<
    set -euo pipefail

    python <<CODE
from csv import DictReader

def basename(bucket_path):
    return bucket_path.rsplit('/', 1)[-1]

fastq1_files = sorted("~{sep=',' fastq1_files}".split(','), key=basename)
fastq2_files = sorted("~{sep=',' fastq2_files}".split(','), key=basename)

readgroups = dict()
with open("~{readgroups_tsv}") as readgroups_tsv:
    for readgroup in DictReader(readgroups_tsv, dialect="excel-tab"):
        readgroups[readgroup["ID"]] = r"@RG\t" + r"\t".join(
            "{}:{}".format(key, value)
            for key, value in readgroup.items() if value)

print("forward_fastq\treverse_fastq\treadgroup\treadgroup_id")
for fastq1, fastq2 in zip(fastq1_files, fastq2_files):
    rg_id = basename(fastq1).rsplit("~{fastq1_suffix}", 1)[0]
    print("\t".join([fastq1, fastq2, readgroups[rg_id], rg_id]))
CODE
    >>>

    output {
        Array[FastqPairRecord] fastq_pair_records = read_objects(stdout())
    }
    
    runtime {
        docker: "python:3.8-slim"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}
        
task emit_se_records {
    input {
        Array[File] fastq_o1_files
        Array[File] fastq_o2_files
        Array[File] fastq_s_files
        Array[Object]+ readgroups
        String fastq_o1_suffix = "_o1.fq.gz"
        String fastq_o2_suffix = "_o2.fq.gz"
        String fastq_s_suffix = "_s.fq.gz"
        Int preemptible = 10
        Int max_retries = 0
        Int cpu = 1
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }
    Int mem = ceil(size(fastq_o1_files, "MiB") + size(fastq_o2_files, "MiB") + size(fastq_s_files, "MiB")) + 2000 + additional_memory_mb
    Int disk_space = ceil(size(fastq_o1_files, "GiB") + size(fastq_o2_files, "GiB") + size(fastq_s_files, "GiB")) + 10 + additional_disk_gb

    File readgroups_tsv = write_objects(readgroups)

    parameter_meta {
        fastq_o1_files: {localization_optional: true}
        fastq_o2_files: {localization_optional: true}
        fastq_s_files: {localization_optional: true}
    }

    command <<<
    set -euo pipefail

    python <<CODE
from csv import DictReader

def basename(bucket_path):
    return bucket_path.rsplit('/', 1)[-1]

def emit_records(fastqs, suffix, readgroups):
    for fastq in fastqs:
        if not fastq:
            return
        rg_id = basename(fastq).rsplit(suffix, 1)[0]
        print("\t".join([fastq, readgroups[rg_id], rg_id]))

fastq_o1_files = "~{sep=',' fastq_o1_files}".split(',')
fastq_o2_files = "~{sep=',' fastq_o2_files}".split(',')
fastq_s_files = "~{sep=',' fastq_s_files}".split(',')

readgroups = dict()
with open("~{readgroups_tsv}") as readgroups_tsv:
    for readgroup in DictReader(readgroups_tsv, dialect="excel-tab"):
        readgroups[readgroup["ID"]] = r"@RG\t" + r"\t".join(
            "{}:{}".format(key, value)
            for key, value in readgroup.items() if value)

print("fastq\treadgroup\treadgroup_id")
emit_records(fastq_o1_files, "~{fastq_o1_suffix}", readgroups)
emit_records(fastq_o2_files, "~{fastq_o2_suffix}", readgroups)
emit_records(fastq_s_files, "~{fastq_s_suffix}", readgroups)
CODE
    >>>

    output {
        Array[FastqSingleRecord] fastq_single_records = read_objects(stdout())
    }
    
    runtime {
        docker: "python:3.8-slim"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task bwa_pe {
    input {
        FastqPairRecord fastq_record
        File ref_fasta
        File ref_fai
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        Int cpu = 16
        Int preemptible = 1
        Int max_retries = 0
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }
    File fastq1 = fastq_record.forward_fastq
    File fastq2 = fastq_record.reverse_fastq
    String readgroup = fastq_record.readgroup
    String outbam = fastq_record.readgroup_id + ".bam"
    Float ref_size =size([ref_fasta, ref_dict, ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa, ref_fai], "GiB")
    Int mem = ceil(size([fastq1, fastq2], "MiB")) + 10000 + additional_memory_mb
    Int disk_space = ceil((size([fastq1, fastq2], "GiB") * 4) + ref_size) + 10 + additional_disk_gb

    command {
        set -euo pipefail
        bwa mem \
            -t ~{cpu} \
            -T 0 \
            -R "~{readgroup}" \
            ~{ref_fasta} \
            ~{fastq1} \
            ~{fastq2} \
        | samtools view \
            -Shb \
            -o ~{outbam} \
            -
    }
    
    output {
        File bam = "~{outbam}"
    }
    
    runtime {
        docker: "broadgdac/bwa:0.7.15-r1142-dirty"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task bwa_se {
    input {
        FastqSingleRecord fastq_record
        File ref_fasta
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        File ref_fai
        Int cpu = 16
        Int preemptible = 1
        Int max_retries = 0
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }
    File fastq = fastq_record.fastq
    String readgroup = fastq_record.readgroup
    String outbam = fastq_record.readgroup_id + ".bam"
    Float ref_size = size([ref_fasta, ref_dict, ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa, ref_fai], "GiB")
    Int mem = ceil(size(fastq, "MiB")) + 10000 + additional_memory_mb
    Int disk_space = ceil((size(fastq, "GiB") * 4) + ref_size) + 10 + additional_disk_gb

    command {
        set -euo pipefail
        bwa mem \
            -t ~{cpu} \
            -T 0 \
            -R "~{readgroup}" \
            ~{ref_fasta} \
            ~{fastq} \
        | samtools view \
            -Shb \
            -o ~{outbam} \
            -
    }
    
    output {
        File bam = "~{outbam}"
    }
    
    runtime {
        docker: "broadgdac/bwa:0.7.15-r1142-dirty"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task picard_markduplicates {
  input {
    Array[File]+ bams
    String outbam

    Int compression_level = 2
    Int preemptible_tries = 1
    Int max_retries = 0
    String validation_stringency = "SILENT"
    String assume_sort_order = "queryname"

    # The program default for READ_NAME_REGEX is appropriate in nearly every case.
    # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
    # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
    String? read_name_regex
    Int memory_multiplier = 1
    Int additional_disk_gb = 20

    Float? sorting_collection_size_ratio
  }
  Float total_input_size = size(bams, "GiB")
  String metrics_filename = outbam + ".metrics"
  # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs and the merged output.
  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
  Float md_disk_multiplier = 3
  Int disk_size = ceil(md_disk_multiplier * total_input_size) + additional_disk_gb

  Float memory_size = 7.5 * memory_multiplier
  Int java_memory_size = (ceil(memory_size) - 2)

  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size}g -jar /usr/picard/picard.jar \
      MarkDuplicates \
      INPUT=~{sep=' INPUT=' bams} \
      OUTPUT=~{outbam} \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=~{validation_stringency} \
      ASSUME_SORT_ORDER=~{assume_sort_order} \
      ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio} \
      ~{"READ_NAME_REGEX=" + read_name_regex}
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File bam = "~{outbam}"
    File metrics = "~{metrics_filename}"
  }
}

task sort_and_index_markdup_bam {
    input {
        File input_bam
        String tmp_prefix = "tmp_srt"
        Int cpu = 8
        Int preemptible = 2
        Int max_retries = 0
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }
    Int mem = ceil(size(input_bam, "MiB")) + 10000 + additional_memory_mb
    Int disk_space = ceil(size(input_bam, "GiB") * 3.25) + 20 + additional_disk_gb
    Int mem_per_thread = floor(mem / cpu * 0.85)
    Int index_threads = cpu - 1
    String output_bam = basename(input_bam)
    String output_bai = basename(input_bam, ".bam") + ".bai"

    command {
        set -euo pipefail
        samtools sort \
            -@ ~{cpu} \
            -o ~{output_bam} \
            -T ~{tmp_prefix} \
            -m ~{mem_per_thread}M \
            ~{input_bam}
        samtools index \
            -b \
            -@ ~{index_threads} \
            ~{output_bam} \
            ~{output_bai}
    }

    output {
        File bam = "~{output_bam}"
        File bai = "~{output_bai}"
    }
    
    runtime {
        docker: "us.gcr.io/broad-gotc-prod/samtools:1.10"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

        
task gatk_baserecalibrator {
    input {
        File bam
        File dbsnp_vcf
        File dbsnp_vcf_index
        File ref_dict
        File ref_fasta
        File ref_fai
        Int cpu = 2
        Int preemptible = 2
        Int max_retries = 0
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }
    String output_grp = basename(bam, ".bam") + "_bqsr.grp"
    Float ref_size = size([ref_fasta, ref_fai, ref_dict], "GiB")
    Float dbsnp_size = size([dbsnp_vcf, dbsnp_vcf_index], "GiB")
    Int mem = ceil(size(bam, "MiB")) + 6000 + additional_memory_mb
    Int jvm_mem = mem - 1000
    Int max_heap = mem - 500
    Int disk_space = ceil(size(bam, "GiB") + ref_size + dbsnp_size) + 20 + additional_disk_gb

    parameter_meta {
        bam: {localization_optional: true}
        dbsnp_vcf: {localization_optional: true}
        dbsnp_vcf_index: {localization_optional: true}
        ref_dict: {localization_optional: true}
        ref_fasta: {localization_optional: true}
        ref_fai: {localization_optional: true}
    }

    command {
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -Xlog:gc=debug:file=gc_log.log -Xms~{jvm_mem}m -Xmx~{max_heap}m" \
            BaseRecalibrator \
                --input ~{bam} \
                --known-sites ~{dbsnp_vcf} \
                --reference ~{ref_fasta} \
                --tmp-dir . \
                --output ~{output_grp}
    }

    output {
        File bqsr_recal_file = "~{output_grp}"
    }

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task gatk_applybqsr {
    input {
        File input_bam
        File bqsr_recal_file
        Boolean emit_original_quals = true
        Int cpu = 2
        Int preemptible = 2
        Int max_retries = 0
        Int additional_memory_mb = 0
        Int additional_disk_gb = 0
    }
    String output_bam = basename(input_bam)
    String output_bai = basename(input_bam, ".bam") + ".bai"
    Int mem = ceil(size(input_bam, "MiB")) + 4000 + additional_memory_mb
    Int jvm_mem = mem - 1000
    Int max_heap = mem - 500
    Int disk_space = ceil((size(input_bam, "GiB") * 3)) + 20 + additional_disk_gb

    parameter_meta {
        input_bam: {localization_optional: true}
    }

    command {
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -Xlog:gc=debug:file=gc_log.log -Xms~{jvm_mem}m -Xmx~{max_heap}m" \
            ApplyBQSR \
                --input ~{input_bam} \
                --bqsr-recal-file ~{bqsr_recal_file} \
                --emit-original-quals ~{emit_original_quals} \
                --tmp-dir . \
                --output ~{output_bam}
    }

    output {
        File bam = "~{output_bam}"
        File bai = "~{output_bai}"
    }
    
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        memory: mem + " MiB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

# Collect quality metrics from the aggregated bam
task collect_insert_size_metrics {
  input {
    File input_bam
    String output_bam_prefix
    Int additional_memory_mb = 0
    Int additional_disk_gb = 0
  }
  Int mem = ceil(size(input_bam, "GiB")) + 7000 + additional_memory_mb
  Int jvm_mem = mem - 1000
  Int max_heap = mem - 500
  Int disk_size = ceil(size(input_bam, "GiB")) + 20 + additional_disk_gb

  command {
    java -Xms~{jvm_mem}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
      CollectInsertSizeMetrics \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix}.insert_size_metrics \
      HISTOGRAM_FILE=~{output_bam_prefix}.insert_size_histogram.pdf
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    memory: mem + " MiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File insert_size_histogram_pdf = "~{output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "~{output_bam_prefix}.insert_size_metrics"
  }
}


workflow GDCWholeGenomeSomaticSingleSample {

    String pipeline_version = "1.3.5"

    input {
        File? input_cram
        File? input_bam
        File? cram_ref_fasta
        File? cram_ref_fasta_index
        File? output_map
        String? unmapped_bam_suffix
        String base_file_name

        File? ubam

        File contamination_vcf
        File contamination_vcf_index
        File dbsnp_vcf
        File dbsnp_vcf_index

        File ref_fasta
        File ref_fai
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
    }

    String outbam = if (defined(ubam) || defined(input_bam)) then basename(select_first([ubam, input_bam]), ".bam") + ".aln.mrkdp.bam"
                    else basename(select_first([input_cram]), ".cram") + ".aln.mrkdp.bam"

    if (!defined(ubam)) {
        call ToUbams.CramToUnmappedBams {
             input:
                 input_cram = input_cram,
                 input_bam = input_bam,
                 ref_fasta = select_first([cram_ref_fasta, ref_fasta]),
                 ref_fasta_index = select_first([cram_ref_fasta_index, ref_fai]),
                 output_map = output_map,
                 base_file_name = base_file_name,
                 unmapped_bam_suffix = unmapped_bam_suffix
        }
    }

    Array[File] ubams = if defined(ubam) then [select_first([ubam])] else select_first([CramToUnmappedBams.unmapped_bams])

    scatter (ubam in ubams) {
        call bam_readgroup_to_contents {
            input: bam = ubam
        }

        call biobambam_bamtofastq {
             input: filename = ubam
        }
    }

    Array[Object] readgroups = flatten(bam_readgroup_to_contents.readgroups)

    Array[File] fastq1 = flatten(biobambam_bamtofastq.output_fastq1)
    Array[File] fastq2 = flatten(biobambam_bamtofastq.output_fastq2)
    Array[File] fastq_o1 = flatten(biobambam_bamtofastq.output_fastq_o1)
    Array[File] fastq_o2 = flatten(biobambam_bamtofastq.output_fastq_o2)
    Array[File] fastq_s = flatten(biobambam_bamtofastq.output_fastq_s)

    Int pe_count = length(fastq1)
    Int o1_count = length(fastq_o1)
    Int o2_count = length(fastq_o2)
    Int s_count = length(fastq_s)

    if (pe_count > 0) {
        call emit_pe_records {
            input:
                fastq1_files = fastq1,
                fastq2_files = fastq2,
                readgroups = readgroups
        }
        scatter (pe_record in emit_pe_records.fastq_pair_records) {
            call bwa_pe {
                input:
                    fastq_record = pe_record,
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    ref_dict = ref_dict,
                    ref_amb = ref_amb,
                    ref_ann = ref_ann,
                    ref_bwt = ref_bwt,
                    ref_pac = ref_pac,
                    ref_sa = ref_sa
            }
        }
    }

    if (o1_count + o2_count + s_count > 0) {
        call emit_se_records {
            input:
                fastq_o1_files = fastq_o1,
                fastq_o2_files = fastq_o2,
                fastq_s_files = fastq_s,
                readgroups = readgroups
        }
        scatter (se_record in emit_se_records.fastq_single_records) {
            call bwa_se {
                input:
                    fastq_record = se_record,
                    ref_fasta = ref_fasta,
                    ref_fai = ref_fai,
                    ref_dict = ref_dict,
                    ref_amb = ref_amb,
                    ref_ann = ref_ann,
                    ref_bwt = ref_bwt,
                    ref_pac = ref_pac,
                    ref_sa = ref_sa
            }
        }
    }

    Array[File] aligned_bams = flatten([select_first([bwa_pe.bam, []]), select_first([bwa_se.bam, []])])

    call picard_markduplicates {
        input:
            bams = aligned_bams,
            outbam = outbam
    }

    call sort_and_index_markdup_bam {
        input: input_bam = picard_markduplicates.bam
    }

    call CheckContamination.CalculateSomaticContamination as check_contamination {
        input:
            reference = ref_fasta,
            reference_dict = ref_dict,
            reference_index = ref_fai,
            tumor_cram_or_bam = sort_and_index_markdup_bam.bam,
            tumor_crai_or_bai = sort_and_index_markdup_bam.bai,
            contamination_vcf = contamination_vcf,
            contamination_vcf_index = contamination_vcf_index
    }

    call gatk_baserecalibrator {
        input:
            bam = sort_and_index_markdup_bam.bam,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            dbsnp_vcf = dbsnp_vcf,
            dbsnp_vcf_index = dbsnp_vcf_index
    }

    call gatk_applybqsr {
        input:
            input_bam = sort_and_index_markdup_bam.bam,
            bqsr_recal_file = gatk_baserecalibrator.bqsr_recal_file
    }

    String output_bam_prefix = basename(gatk_applybqsr.bam, ".bam")

    call collect_insert_size_metrics {
        input:
            input_bam = gatk_applybqsr.bam,
            output_bam_prefix = output_bam_prefix
    }

    output {
        Array[File]? validation_report = CramToUnmappedBams.validation_report
        Array[File]? unmapped_bams = CramToUnmappedBams.unmapped_bams
        File bam = gatk_applybqsr.bam
        File bai = gatk_applybqsr.bai
        File md_metrics = picard_markduplicates.metrics
        File insert_size_metrics = collect_insert_size_metrics.insert_size_metrics
        File insert_size_histogram_pdf = collect_insert_size_metrics.insert_size_histogram_pdf
        File contamination = check_contamination.contamination
    }
    meta {
        allowNestedInputs: true
    }
}
