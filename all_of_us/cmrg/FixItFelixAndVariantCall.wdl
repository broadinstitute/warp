version 1.0

workflow FixItFelixAndVariantCall {
    input {
        File cram_file
        File cram_file_index

        File original_ref_fasta
        File original_ref_fasta_index
        File original_ref_dict

        File masked_ref_fasta = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2.fasta"
        File masked_ref_fasta_index = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2.fasta.fai"
        File masked_ref_dict = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2.dict"
        File masked_ref_amb = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2.fasta.amb"
        File masked_ref_ann = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2.fasta.ann"
        File masked_ref_bwt = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2.fasta.bwt"
        File masked_ref_pac = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2.fasta.pac"
        File masked_ref_sa = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/GRCh38_masked_v2_decoy_excludes_GPRIN2_DUSP22_FANCD2.fasta.sa"

        File true_locations_intervals = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/real_locations.hg38.bed"
        File false_duplications_intervals = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/false_duplications.hg38.bed"
        File combined_true_false_intervals = "gs://fc-ee5da145-1f24-4104-bf83-3a4e318d9968/combined_real_and_false_locations.hg38.bed"

        Boolean generate_gvcf = false
    }

    Int original_cram_size = ceil(size(cram_file, "GB"))
    String pipeline_version = "aou_9.0.1"

    call subset_cram {
        input:
            input_cram = cram_file,
            input_cram_index = cram_file_index,
            ref_fasta = original_ref_fasta,
            ref_fasta_index = original_ref_fasta_index,
            ref_dict = original_ref_dict,
            intervals = combined_true_false_intervals,
            disk_size = original_cram_size + 20
    }

    Int extracted_cram_size = ceil(size(subset_cram.output_bam, "GB"))

    call FixItFelix {
        input:
            reads = subset_cram.output_bam,
            reads_index = subset_cram.output_bam_index,
            intervals = combined_true_false_intervals,
            output_name = subset_cram.sample_basename,
            masked_ref_dict = masked_ref_dict,
            masked_ref_fasta = masked_ref_fasta,
            masked_ref_fasta_index = masked_ref_fasta_index,
            masked_ref_amb = masked_ref_amb,
            masked_ref_ann = masked_ref_ann,
            masked_ref_bwt = masked_ref_bwt,
            masked_ref_pac = masked_ref_pac,
            masked_ref_sa = masked_ref_sa,
            disk_size = extracted_cram_size + 20
    }

    Int fif_cram_size = ceil(size(FixItFelix.output_bam, "GB"))

    call call_variants {
        input:
            input_bam = FixItFelix.output_bam,
            input_bam_index = FixItFelix.output_bam_index,
            ref_fasta = masked_ref_fasta,
            ref_dict = masked_ref_dict,
            ref_fasta_index = masked_ref_fasta_index,
            location_intervals = true_locations_intervals,
            output_name = subset_cram.sample_basename,
            disk_size = fif_cram_size + 20,
            generate_gvcf = generate_gvcf
    }

    output {
        File output_vcf = call_variants.output_vcf
        File output_vcf_index = call_variants.output_vcf_index
        String output_pipeline_version = pipeline_version
    }
}

task subset_cram {
    input {
        File input_cram
        File input_cram_index
        File ref_fasta
        File ref_dict
        File ref_fasta_index
        File intervals
        Int disk_size
    }

    parameter_meta {
        input_cram: {
            localization_optional: true
        }
    }

    Boolean is_cram = sub(basename(input_cram), ".*\\.", "") == "cram"
    String sample_name = if is_cram then basename(input_cram, ".cram") else basename(input_cram, ".bam")

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.4.0.0"
    String gatk_path = "/gatk/gatk"
    String output_bam_name = sample_name + ".bam"
    String output_index = sample_name + ".bai"

    command {
        set -euo pipefail
        ~{gatk_path} --java-options "-Xmx2G" \
        PrintReads \
        -R ~{ref_fasta} \
        -I ~{input_cram} \
        -L ~{intervals} \
        -O ~{output_bam_name}
    }

    runtime {
        docker: gatk_docker
        memory: "3 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
    }

    output {
        String sample_basename = "~{sample_name}"
        File output_bam = "~{output_bam_name}"
        File output_bam_index = "~{output_index}"
    }
}

task FixItFelix {
    input {
        File reads
        File reads_index
        File intervals
        String output_name
        File masked_ref_dict
        File masked_ref_fasta
        File masked_ref_fasta_index
        File masked_ref_amb
        File masked_ref_ann
        File masked_ref_bwt
        File masked_ref_pac
        File masked_ref_sa
        String old_ref = ""  # Optional input with an empty string as the default value
        Int disk_size
        Int threads = 2
    }

    command <<<
        set -euo pipefail
        bam=~{reads}
        bed=~{intervals}
        ref=~{masked_ref_fasta}
        out=~{output_name}
        old_ref=~{old_ref}
        THREAD=~{threads}

        if [[ $bam =~ \.cram$ &&  -z "$old_ref" ]]; then
            echo "Error: For CRAM files, use -f parameter with old reference"
            exit 1
        fi

        echo "BAM/CRAM file: $bam"
        echo "Bed file: $bed"
        echo "Reference: $ref"
        echo "Output: $out"
        echo "Old Ref: $old_ref"

        export PATH="/usr/gitc:$PATH"

        ls /usr/gitc

        start=$(date +%s.%N)

        if [[ -e $ref && -e $bam && -e $bed ]] ; then
            echo -e "\tExtract reads for remapping"
            if [[ $bam =~ \.cram$ ]]; then
                samtools view -F 2316 -hb --reference ${old_ref} -o ${out}_extracted_reads.bam ${bam} `perl -ane '{print "$F[0]:$F[1]-$F[2] "}' $bed`
            else
                samtools view -F 2316 -hb -o ${out}_extracted_reads.bam ${bam} `perl -ane '{print "$F[0]:$F[1]-$F[2] "}' $bed`
            fi

            samtools view -hb -F 2316 ${out}_extracted_reads.bam | samtools sort -@ ${THREAD} -n -o ${out}_extracted_reads_sorted.bam
            samtools view ${out}_extracted_reads_sorted.bam | cut -f1 | uniq -c | awk '$1!=2 {print $2}' > ${out}_nonpairs_rnames.txt
            samtools view -Sh ${out}_extracted_reads_sorted.bam | fgrep -vf ${out}_nonpairs_rnames.txt | samtools view -hb - > ${out}_pairs_only.bam

            samtools sort -@ ${THREAD} -n -o ${out}_original_sorted_by_read_names.bam ${out}_pairs_only.bam
            samtools sort -@ ${THREAD} -o ${out}_original_sorted.bam ${out}_pairs_only.bam
            samtools index ${out}_original_sorted.bam

            echo -e "\tConvert reads to fastq"
            samtools view -H ${out}_extracted_reads.bam | grep "^@RG" | sed 's/	/\\t/g' > ${out}_bug_regions.RG.txt
            samtools fastq "${out}_original_sorted_by_read_names.bam" -1 "${out}_extract_1.fastq" -2 "${out}_extract_2.fastq" -0 /dev/null -s /dev/null -n

            echo -e "\tStart remapping"
            RG=$(head -n 1 ${out}_bug_regions.RG.txt)
            bwa mem -t "${THREAD}" -R "${RG}" "${ref}" "${out}_extract_1.fastq" "${out}_extract_2.fastq" | samtools view -hb - > "${out}_remapped.bam"
            samtools sort -@ "${THREAD}" -o "${out}_remapped_sorted.bam" "${out}_remapped.bam"
            samtools index ${out}_remapped_sorted.bam

            rm ${out}_remapped.bam
            rm ${out}_extracted_reads.bam
            rm ${out}_extracted_reads_sorted.bam
            rm ${out}_original_sorted_by_read_names.bam
            rm ${out}_pairs_only.bam
            rm ${out}_nonpairs_rnames.txt
            rm ${out}_bug_regions.RG.txt
            rm ${out}_extract_2.fastq
            rm ${out}_extract_1.fastq
            echo "Successfully finished"
        else
            echo "Error one or more files do not exist. Check your paths"
            exit 1
        fi

        echo "<out_prefix>_original_sorted.bam <- Sorted BAM files correspond to the extracted paired-end reads (Original mapping)"
        echo "<out_prefix>_remapped_sorted.bam <- Sorted bam file after re-mapping extracted reads to GRCh38 masked v2 reference (mapping to modified reference)"
        echo "Script ended at: $(date +%s.%N)"
    >>>

    runtime {
        preemptible: 0
        docker: "gcr.io/broad-dsde-methods/fixitfelix:1"
        memory: "10 GiB"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_bam = "~{output_name}_remapped_sorted.bam"
        File output_bam_index = "~{output_name}_remapped_sorted.bam.bai"
    }
}

task call_variants {
    input {
        File input_bam
        File input_bam_index
        File ref_fasta
        File ref_dict
        File ref_fasta_index
        File location_intervals
        Int disk_size
        String output_name
        Boolean generate_gvcf
    }

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.4.0.0"
    String gatk_path = "/gatk/gatk"
    String output_filename = output_name + (if generate_gvcf then ".g.vcf.gz" else ".vcf.gz")

    command <<<
        set -euo pipefail
        ~{gatk_path} --java-options "-Xmx3G" \
        HaplotypeCaller \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        -L ~{location_intervals} \
        --dragen-mode \
        -O ~{output_filename} \
        ~{if generate_gvcf then "--emit-ref-confidence GVCF" else ""}

        if [ "~{generate_gvcf}" == "false" ]; then
            ~{gatk_path} --java-options "-Xmx1G" \
            SelectVariants \
            -V ~{output_filename} \
            --exclude-non-variants \
            -O filtered_~{output_filename}
        fi
    >>>

    runtime {
        docker: gatk_docker
        memory: "4 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
    }

    output {
        File output_vcf = "~{if generate_gvcf then output_filename else 'filtered_' + output_filename}"
        File output_vcf_index = "~{if generate_gvcf then output_filename + '.tbi' else 'filtered_' + output_filename + '.tbi'}"
    }
}