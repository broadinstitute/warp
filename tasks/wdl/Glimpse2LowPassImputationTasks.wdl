version 1.0

task MergeSampleChunksVcfsWithPaste {
    input {
        Array[File] input_vcfs
        String output_vcf_basename

        Int disk_size_gb = ceil(2.2 * size(input_vcfs, "GiB") + 50)
        Int mem_gb = 12
        Int cpu = 4
        Int preemptible = 3
    }

    command <<<
        set -euo pipefail

        vcfs=(~{sep=" " input_vcfs})

        mkfifo fifo_0
        mkfifo fifo_to_paste_0

        i=1

        fifos_to_paste=()
        md5sums=()
        bcftools view -h --no-version ${vcfs[0]} | awk '!/^#CHROM/' > header.vcf
        n_lines=$(wc -l header.vcf | cut -d' ' -f1)

        bgzip -d ${vcfs[0]} -o fifo_0 &

        tail +$((n_lines)) fifo_0 | tee fifo_to_paste_0 | cut -f1-5,9 | md5sum > md5sum_0 &

        for vcf in "${vcfs[@]:1}"; do
        fifo_name="fifo_$i"
        mkfifo "$fifo_name"

        fifo_name_to_md5="fifo_to_md5_$i"
        mkfifo "$fifo_name_to_md5"

        fifo_name_to_paste="fifo_to_paste_$i"
        mkfifo "$fifo_name_to_paste"
        fifos_to_paste+=("$fifo_name_to_paste")

        file_name_md5sum="md5sum_$i"
        md5sums+=("$file_name_md5sum")
        n_lines=$(bcftools view -h --no-version $vcf | awk '!/^#CHROM/' | wc -l | cut -d' ' -f1)

        bgzip -d ${vcf} -o "$fifo_name" &
        tail +$((n_lines)) "$fifo_name" | tee "$fifo_name_to_md5" | cut -f 10- > "$fifo_name_to_paste" &
        cut -f1-5,9 "$fifo_name_to_md5" | md5sum > "$file_name_md5sum" &

        ((i++))
        done

        mkfifo fifo_to_cat

        paste fifo_to_paste_0 "${fifos_to_paste[@]}" | tee fifo_to_cat | awk 'NR % 5000000 == 0' | cut -f 1-5 &

        cat header.vcf fifo_to_cat | bgzip -o ~{output_vcf_basename}.vcf.gz

        for md5sum_file in "${md5sums[@]}"; do
        diff <(cat md5sum_0) <(cat $md5sum_file) >> /dev/null || (echo "Fields 1-5,9 do not match for $md5sum_file" && exit 1)
        done

        for fifo in fifo_*; do
        rm $fifo
        done

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/bcftools_bgzip:beagle_imputation_v1.0.0"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: 1
        noAddress: true
    }

    output {
        File output_vcf = "~{output_vcf_basename}.vcf.gz"
    }
}
