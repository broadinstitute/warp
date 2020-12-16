version 1.0

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
        Int mem = 2
        Int disk_space = 10
        Int preemptible = 10
        Int max_retries = 0
        Int cpu = 1
    }

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
        memory: mem + " GB"
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
        Int mem = 2
        Int preemptible = 2
        Int max_retries = 0
    }
    Int disk_space = ceil(size(filename, "G") * 2) + 10

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
        memory: mem + " GB"
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
        Int mem = 2
        Int disk_space = 10
        Int preemptible = 10
        Int max_retries = 0
        Int cpu = 1
    }

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
        memory: mem + " GB"
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
        Int mem = 2
        Int disk_space = 10
        Int preemptible = 10
        Int max_retries = 0
        Int cpu = 1
    }

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
        memory: mem + " GB"
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
        Int mem = 10
        Int preemptible = 1
        Int max_retries = 0
    }
    File fastq1 = fastq_record.forward_fastq
    File fastq2 = fastq_record.reverse_fastq
    String readgroup = fastq_record.readgroup
    String outbam = fastq_record.readgroup_id + ".bam"
    Float ref_size =size([ref_fasta, ref_dict, ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa, ref_fai], "G")
    Int disk_space = ceil((size([fastq1, fastq2], "G") * 2) + ref_size) + 10

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
        memory: mem + " GB"
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
        Int mem = 10
        Int preemptible = 1
        Int max_retries = 0
    }
    File fastq = fastq_record.fastq
    String readgroup = fastq_record.readgroup
    String outbam = fastq_record.readgroup_id + ".bam"
    Float ref_size = size([ref_fasta, ref_dict, ref_amb, ref_ann, ref_bwt, ref_pac, ref_sa, ref_fai], "G")
    Int disk_space = ceil((size(fastq, "G") * 2) + ref_size) + 10

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
        memory: mem + " GB"
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
        String validation_stringency = "SILENT"
        String assume_sort_order = "queryname"
        Int cpu = 1
        Int mem = 16
        Int preemptible = 2
        Int max_retries = 0
    }
    String metrics_file = outbam + ".metrics"
    Int jvm_mem = if mem > 1 then mem - 1  else 1
    Int disk_space = ceil(size(bams, "G") * 2.2)

    command {
        set -euo pipefail
        java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
             -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
             -Xloggc:gc_log.log -Xmx~{jvm_mem}g -jar /usr/picard/picard.jar \
            MarkDuplicates \
                INPUT=~{sep=" INPUT=" bams} \
                TMP_DIR=. \
                VALIDATION_STRINGENCY=~{validation_stringency} \
                ASSUME_SORT_ORDER=~{assume_sort_order} \
                OUTPUT=~{outbam} \
                METRICS_FILE=~{metrics_file}
    }

    output {
        File metrics = "~{metrics_file}"
        File bam = "~{outbam}"
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.18.11"
        memory: mem + " GB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

task sort_and_index_markdup_bam {
    input {
        File input_bam
        String tmp_prefix = "tmp_srt"
        Int cpu = 8
        Int mem = 16
        Int preemptible = 2
        Int max_retries = 0
    }
    Int disk_space = ceil(size(input_bam, "G") * 3.25) + 20
    Int mem_per_thread = floor(mem * 1024 / cpu * 0.85)
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
        docker: "broadgdac/samtools:1.10"
        memory: mem + " GB"
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
        Int mem = 6
        Int preemptible = 2
        Int max_retries = 0
    }
    String output_grp = basename(bam, ".bam") + "_bqsr.grp"
    Float ref_size = size([ref_fasta, ref_fai, ref_dict], "G")
    Float dbsnp_size = size([dbsnp_vcf, dbsnp_vcf_index], "G")
    Int jvm_mem = if mem > 1 then mem - 1 else 1
    Int disk_space = ceil(size(bam, "G") + ref_size + dbsnp_size) + 20

    parameter_meta {
        bam: {localization_optional: true}
    }

    command {
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms~{jvm_mem}g" \
            BaseRecalibrator \
                --input ~{bam} \
                --known-sites ~{dbsnp_vcf} \
                --reference ~{ref_fasta} \
                --TMP_DIR . \
                --output ~{output_grp}
    }

    output {
        File bqsr_recal_file = "~{output_grp}"
    }

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.0.7.0"
        memory: mem + " GB"
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
        Int mem = 4
        Int preemptible = 2
        Int max_retries = 0
    }
    String output_bam = basename(input_bam)
    String output_bai = basename(input_bam, ".bam") + ".bai"
    Int jvm_mem = if mem > 1 then mem - 1 else 1
    Int disk_space = ceil((size(input_bam, "G") * 3)) + 20

    parameter_meta {
        input_bam: {localization_optional: true}
    }

    command {
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms~{jvm_mem}g" \
            ApplyBQSR \
                --input ~{input_bam} \
                --bqsr-recal-file ~{bqsr_recal_file} \
                --emit-original-quals ~{emit_original_quals} \
                --TMP_DIR . \
                --output ~{output_bam}
    }

    output {
        File bam = "~{output_bam}"
        File bai = "~{output_bai}"
    }
    
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.0.7.0"
        memory: mem + " GB"
        disks: "local-disk " + disk_space + " HDD"
        preemptible: preemptible
        maxRetries: max_retries
        cpu: cpu
    }
}

workflow GDCWholeGenomeSomaticSingleSample {

    String pipeline_version = "1.0.0"

    input {
        File ubam
        File ref_fasta
        File ref_fai
        File ref_dict
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
    }
    String outbam = basename(ubam, ".bam") + ".aln.mrkdp.bam"

    call bam_readgroup_to_contents {
        input: bam = ubam
    }

    call biobambam_bamtofastq {
        input: filename = ubam
    }

    Int pe_count = length(biobambam_bamtofastq.output_fastq1)
    Int o1_count = length(biobambam_bamtofastq.output_fastq_o1)
    Int o2_count = length(biobambam_bamtofastq.output_fastq_o2)
    Int s_count = length(biobambam_bamtofastq.output_fastq_s)

    if (pe_count > 0) {
        call emit_pe_records {
            input:
                fastq1_files = biobambam_bamtofastq.output_fastq1,
                fastq2_files = biobambam_bamtofastq.output_fastq2,
                readgroups = bam_readgroup_to_contents.readgroups
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
                fastq_o1_files = biobambam_bamtofastq.output_fastq_o1,
                fastq_o2_files = biobambam_bamtofastq.output_fastq_o2,
                fastq_s_files = biobambam_bamtofastq.output_fastq_s,
                readgroups = bam_readgroup_to_contents.readgroups
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

    call gatk_baserecalibrator {
        input:
            bam = sort_and_index_markdup_bam.bam,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict
    }

    call gatk_applybqsr {
        input:
            input_bam = sort_and_index_markdup_bam.bam,
            bqsr_recal_file = gatk_baserecalibrator.bqsr_recal_file
    }

    output {
        File bam = gatk_applybqsr.bam
        File bai = gatk_applybqsr.bai
        File md_metrics = picard_markduplicates.metrics
    }
}
