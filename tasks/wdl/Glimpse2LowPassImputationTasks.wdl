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

task ExtractAnnotations {
    input {
        File imputed_vcf
        File imputed_vcf_index
        Int batch_index

        String docker_extract_annotations
        Int disk_size_gb = ceil(2 * size(imputed_vcf, "GiB") + 50)
        Int mem_gb = 2
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        gatk VariantsToTable -V ~{imputed_vcf} -O annotations_batch_~{batch_index}.tsv -F CHROM -F POS -F REF -F ALT -F AF -F INFO
        bgzip annotations_batch_~{batch_index}.tsv
    >>>

    runtime {
        docker: docker_extract_annotations
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File annotations = "annotations_batch_~{batch_index}.tsv.gz"
    }
}

task RecomputeAndAnnotate {
    input {
        File merged_vcf
        Array[File] annotations

        Array[Int] num_samples

        String output_basename

        String docker_merge
        Int disk_size_gb = ceil(2.2 * size(merged_vcf, "GiB") + 50)
        Int mem_gb
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        cat <<EOF > script.py
import pandas as pd
import functools

input_filenames = ['~{sep="', '" annotations}']
num_samples = [~{sep=", " num_samples}]
if len(num_samples) != len(input_filenames):
    raise RuntimeError('The number of input annotations does not match the number of input number of samples.')
num_batches = len(input_filenames)

def calculate_af(row):
    return sum([row[f'AF_{i}'] * num_samples[i] for i in range(num_batches)]) / sum(num_samples)
def calculate_info(row):
    aggregated_af = row['AF']
    return 1 if aggregated_af == 0 or aggregated_af == 1 else \
                     1 - \
                    (sum([(1 - row[f'INFO_{i}']) * 2 * num_samples[i] * row[f'AF_{i}'] * (1 - row[f'AF_{i}']) for i in range(num_batches)])) / \
                    (2 * sum(num_samples) * aggregated_af * (1 - aggregated_af))

annotation_dfs = [pd.read_csv(input_filename, sep='\t').rename(columns={'AF': f'AF_{i}', 'INFO': f'INFO_{i}'}) for i, input_filename in enumerate(input_filenames)]
annotations_merged = functools.reduce(lambda left, right: pd.merge(left, right, on=['CHROM', 'POS', 'REF', 'ALT'], how='inner', validate='one_to_one'), annotation_dfs)
annotations_merged['AF'] = annotations_merged.apply(lambda row: calculate_af(row), axis=1)
annotations_merged['INFO'] = annotations_merged.apply(lambda row: calculate_info(row), axis=1)
annotations_merged.to_csv('aggregated_annotations.tsv', sep='\t', columns=['CHROM', 'POS', 'REF', 'ALT', 'AF', 'INFO'], header=False, index=False)

EOF
        python3 script.py

        bgzip aggregated_annotations.tsv
        tabix -s1 -b2 -e2 aggregated_annotations.tsv.gz

        bcftools annotate -a aggregated_annotations.tsv.gz -c CHROM,POS,REF,ALT,AF,INFO -O z -o ~{output_basename}.vcf.gz ~{merged_vcf}
    >>>

    runtime {
        docker: docker_merge
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File merged_imputed_vcf = "~{output_basename}.vcf.gz"
        File aggregated_annotations = "aggregated_annotations.tsv.gz"
    }
}

task MergeQCMetrics {
    input {
        Array[File] qc_metrics
        String docker = "us.gcr.io/broad-dsde-methods/python-data-slim:1.1"
        String output_basename

        Int disk_size_gb = ceil(2.2 * size(qc_metrics, "GiB") + 50)
        Int mem_gb = 4
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -xeuo pipefail


        python3 <<EOF
        import pandas as pd
        qc_metrics = ['~{sep="', '" qc_metrics}']
        merged_qc_metrics = pd.concat([pd.read_csv(qc_metric, sep='\t', dtype=str) for qc_metric in qc_metrics])
        merged_qc_metrics.to_csv('~{output_basename}.qc_metrics.tsv', sep='\t', index=False)
        EOF
    >>>

    output {
        File merged_qc_metrics = "~{output_basename}.qc_metrics.tsv"
    }

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

task MergeCoverageMetrics {
    input {
        Array[File] coverage_metrics
        String docker = "us.gcr.io/broad-dsde-methods/python-data-slim:1.1"
        String output_basename

        Int disk_size_gb = ceil(2.2 * size(coverage_metrics, "GiB") + 50)
        Int mem_gb = 4
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -xeuo pipefail


        python3 <<EOF
        import pandas as pd
        from collections import defaultdict
        coverage_metrics = ['~{sep="', '" coverage_metrics}']
        all_types = defaultdict(lambda: str)
        all_types.update({"Chunk":int, "ID":int})
        merged_coverage_metrics_array = [pd.read_csv(coverage_metric, sep='\t', dtype=all_types) for coverage_metric in coverage_metrics]
        id_offset = 0
        for cov_metric in merged_coverage_metrics_array:
            cov_metric.ID += id_offset
            id_offset = cov_metric.ID.max() + 1
        merged_coverage_metrics = pd.concat(merged_coverage_metrics_array).sort_values(['Chunk','ID'])
        merged_coverage_metrics.to_csv('~{output_basename}.coverage_metrics.txt', sep='\t', index=False)
        EOF
    >>>

    output {
        File merged_coverage_metrics = "~{output_basename}.coverage_metrics.txt"
    }

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

task SelectVariantRecordsOnly {
    input {
        File vcf
        File? vcf_index
        String basename

        Int disk_size_gb = ceil(2*size(vcf, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 3000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command {
        set -e -o pipefail

        # keep alt sites (i.e. remove hom ref sites)
        bcftools view -i 'GT[*]="alt"' -Oz -o ~{basename}.vcf.gz ~{vcf}
        tabix ~{basename}.vcf.gz
    }

    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        maxRetries: 1
        preemptible: 3
        noAddress: true
    }

    output {
        File output_vcf = "~{basename}.vcf.gz"
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
    }
}

task CreateHomRefSitesOnlyVcf {
    input {
        File vcf
        File? vcf_index
        String basename

        Int disk_size_gb = ceil(2*size(vcf, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command {
        set -e -o pipefail

        # create header with only first 8 columns and store that
        bcftools view -h ~{vcf} | grep "^##" > ~{basename}.vcf
        bcftools view -h ~{vcf} | grep -v "^##" | cut -f1-8 >> ~{basename}.vcf

        # append first 8 columns of hom ref sites to previously stored header
        bcftools query -e 'GT[*]="alt"' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n' ~{vcf} >> ~{basename}.vcf

        bgzip ~{basename}.vcf
    }

    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        maxRetries: 1
        preemptible: 3
        noAddress: true
    }

    output {
        File output_vcf = "~{basename}.vcf.gz"
    }
}
