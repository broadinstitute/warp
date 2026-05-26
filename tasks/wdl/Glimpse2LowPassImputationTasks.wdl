version 1.0

task MergeSampleChunksVcfsWithPaste {
    input {
        Array[File] input_vcfs
        String output_vcf_basename

        Int disk_size_gb = ceil(2.2 * size(input_vcfs, "GiB") + 50)
        Int mem_gb = 8
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
        Int cpu = 1
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
        Int disk_size_gb = ceil(2.2 * size(merged_vcf, "GiB") + size(annotations, "GiB") + 50)
        Int mem_gb = 6
        Int cpu = 1
        Int preemptible = 1
        Int chunk_size = 100000
    }

    command <<<
        cat <<EOF > script.py
import pandas as pd
import numpy as np

input_filenames = ['~{sep="', '" annotations}']
num_samples = [~{sep=", " num_samples}]
if len(num_samples) != len(input_filenames):
    raise RuntimeError('The number of input annotations does not match the number of input number of samples.')

total_samples = sum(num_samples)
num_batches = len(input_filenames)
chunk_size = ~{chunk_size}

# Stream all annotation files in parallel chunks rather than loading everything into memory at once.
# This keeps memory usage proportional to chunk_size * num_batches rather than total_sites * num_batches.
readers = [pd.read_csv(f, sep='\t', chunksize=chunk_size) for f in input_filenames]

with open('aggregated_annotations.tsv', 'w') as out:
    for chunks in zip(*readers):
        # Validate that all batches have identical sites for this chunk
        ref_loci = chunks[0][['CHROM', 'POS', 'REF', 'ALT']].reset_index(drop=True)
        for i, chunk in enumerate(chunks[1:], 1):
            if not ref_loci.equals(chunk[['CHROM', 'POS', 'REF', 'ALT']].reset_index(drop=True)):
                raise RuntimeError(f'Sites in chunk do not match between batch 0 and batch {i}. '
                                   f'First mismatch at: {ref_loci[~ref_loci.eq(chunk[["CHROM","POS","REF","ALT"]].reset_index(drop=True)).all(axis=1)].head(1).to_dict("records")}')

        # Vectorized weighted AF across batches
        agg_af = sum(chunks[i]['AF'].values * num_samples[i] for i in range(num_batches)) / total_samples

        # Vectorized weighted INFO across batches
        numerator = sum(
            (1 - chunks[i]['INFO'].values) * 2 * num_samples[i] * chunks[i]['AF'].values * (1 - chunks[i]['AF'].values)
            for i in range(num_batches)
        )
        denominator = 2 * total_samples * agg_af * (1 - agg_af)
        # INFO is defined as 1 for monomorphic sites (AF == 0 or AF == 1)
        polymorphic = (agg_af != 0) & (agg_af != 1)
        agg_info = np.where(polymorphic, 1 - np.divide(numerator, denominator, where=polymorphic, out=np.zeros_like(denominator)), 1.0)

        def round_to_n_sig_figs(x, n):
            if x == 0:
                return 0.0
            return round(float(x), n - 1 - int(np.floor(np.log10(abs(x)))))

        result = ref_loci.copy()
        # Cap INFO and AF values at 3 sig-figs to avoid blowing up the output file size w/ overprecision
        result['AF'] = np.vectorize(round_to_n_sig_figs)(agg_af, 3)
        result['INFO'] = np.vectorize(round_to_n_sig_figs)(agg_info, 3)
        result.to_csv(out, sep='\t', header=False, index=False)

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

task FilterVcfByInfo {
    input {
        File vcf
        File? vcf_index
        Float info_threshold
        String basename

        Int disk_size_gb = ceil(2*size(vcf, "GiB")) + 10
        Int cpu = 1
        Int mem_gb = 4
        String docker = "us.gcr.io/broad-dsde-methods/bcftools_bgzip:beagle_imputation_v1.0.0"
    }

    command <<<
        set -e -o pipefail

        bcftools filter -i 'INFO/INFO >= ~{info_threshold}' -Oz -o ~{basename}.vcf.gz ~{vcf}
        bcftools index -t ~{basename}.vcf.gz
    >>>

    runtime {
        docker: docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: 3
        maxRetries: 1
        noAddress: true
    }

    output {
        File output_vcf = "~{basename}.vcf.gz"
        File output_vcf_index = "~{basename}.vcf.gz.tbi"
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

task ConvertCramManifestToInputArrays {
    input {
        File cram_manifest
    }

    command <<<
        cat <<EOF > script.py
        import sys
        import pandas as pd

        crams_filename = "crams.txt"
        cram_indices_filename = "cram_indices.txt"
        sample_ids_filename = "sample_ids.txt"

        def write_column(column_data, filename):
            """Write column to file, with each value stripped of leading/trailing whitespace."""
            filtered = column_data.fillna('').astype(str).str.strip()
            with open(filename, 'w') as f:
                for value in filtered:
                    f.write(f"{value}\n")

        # Read the manifest
        df = pd.read_csv("~{cram_manifest}", sep='\t')

        # Check for required columns
        required_cols = ['sample_id', 'cram_path', 'cram_index_path']
        missing_cols = [col for col in required_cols if col not in df.columns]

        if missing_cols:
            print(f"Missing required columns in the CRAM manifest: {', '.join(missing_cols)}.", file=sys.stderr)
            exit(1)
        else:
            # Write to output files, stripping leading/trailing whitespace from each value
            write_column(df['sample_id'], sample_ids_filename)
            write_column(df['cram_path'], crams_filename)
            write_column(df['cram_index_path'], cram_indices_filename)
        EOF
        python3 script.py >&2
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "4 GiB"
        preemptible: 3
        maxRetries: 2
        noAddress: true
    }

    output {
        Array[String] crams = read_lines("crams.txt")
        Array[String] cram_indices = read_lines("cram_indices.txt")
        Array[String] sample_ids = read_lines("sample_ids.txt")
    }
}

task UpdateHeader {
    input {
        File vcf
        File ref_dict
        String? pipeline_header_line
        String output_basename

        Int mem_gb = 2
        Int cpu = 1
        Int disk_size_gb = ceil(2.1 * size(vcf, "GiB") + 5)
        Int max_retries = 1
        String docker
    }

    command <<<
        set -xeuo pipefail

        # Set correct reference dictionary
        bcftools view -h --no-version ~{vcf} > old_header.vcf
        java -jar /picard.jar UpdateVcfSequenceDictionary -I old_header.vcf --SD ~{ref_dict} -O updated_header.vcf

        # Add pipeline_header_line if provided
        if [ -n "~{default="" pipeline_header_line}" ]; then
            TOTAL_LINES=$(wc -l < "updated_header.vcf")
            REMOVED_COMMENT_CHARACTER_HEADER_LINE=$(echo "~{pipeline_header_line}" | sed 's/^#*//')
            sed -i "${TOTAL_LINES}i\##${REMOVED_COMMENT_CHARACTER_HEADER_LINE}" updated_header.vcf
        fi

        # Apply the final header (with ref dict and optionally pipeline_header_line) to the VCF
        bcftools reheader -h updated_header.vcf -o ~{output_basename}.vcf.gz ~{vcf}
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " SSD"
        memory: mem_gb + " GiB"
        cpu: cpu
        maxRetries: max_retries
        preemptible: 3
        noAddress: true
    }

    output {
        File output_vcf = "~{output_basename}.vcf.gz"
    }
}

task GatherVcfsNoIndex {
    input {
        Array[File] input_vcfs
        String output_vcf_basename

        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
        Int cpu = 2
        Int memory_mb = 10000
        Int disk_size_gb = ceil(3*size(input_vcfs, "GiB")) + 10
    }
    Int command_mem = memory_mb - 1500
    Int max_heap = memory_mb - 1000

    command <<<
        set -e -o pipefail

        gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
        GatherVcfs \
        -I ~{sep=' -I ' input_vcfs} \
        --REORDER_INPUT_BY_FIRST_VARIANT \
        -O ~{output_vcf_basename}.vcf.gz
    >>>
    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        maxRetries: 1
        noAddress: true
    }
    output {
        File output_vcf = "~{output_vcf_basename}.vcf.gz"
    }
}

task CreateVcfIndexAndMd5 {
    input {
        File vcf_input

        Int disk_size_gb = ceil(1.1*size(vcf_input, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
        Int preemptible = 3
    }

    String vcf_basename = basename(vcf_input, ".vcf.gz")

    command <<<
        set -e -o pipefail

        ln -sf ~{vcf_input} ~{vcf_basename}.vcf.gz

        bcftools index -t ~{vcf_basename}.vcf.gz

        md5sum ~{vcf_basename}.vcf.gz | awk '{ print $1 }' > ~{vcf_basename}.md5sum
    >>>
    runtime {
        docker: gatk_docker
        disks: "local-disk ${disk_size_gb} SSD"
        memory: "${memory_mb} MiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: 1
        noAddress: true
    }
    output {
        File output_vcf = "~{vcf_basename}.vcf.gz"
        File output_vcf_index = "~{vcf_basename}.vcf.gz.tbi"
        File output_vcf_md5sum = "~{vcf_basename}.md5sum"
    }
}

task CollectQCMetrics {
    input {
        File imputed_vcf
        String output_basename

        Int preemptible = 0
        String docker = "mirror.gcr.io/hailgenetics/hail:0.2.126-py3.11"
        Int cpu = 4
        Int mem_gb = 8
    }

    parameter_meta {
        imputed_vcf: {
                         localization_optional: true
                     }
    }

    Int disk_size_gb = ceil(2*size(imputed_vcf, "GiB") + 50)

    command <<<
        set -euo pipefail

        cat <<'EOF' > script.py
        import hail as hl
        import pandas as pd

        # Calculate metrics
        hl.init(default_reference='GRCh38', idempotent=True)
        vcf = hl.import_vcf('~{imputed_vcf}', force_bgz=True)
        qc = hl.sample_qc(vcf)
        qc_pd = qc.cols().flatten() \
        .rename({'sample_qc.' + col: col for col in list(qc['sample_qc'])}) \
        .rename({'s': 'sample_id'}) \
        .to_pandas()
        qc_pd.to_csv('~{output_basename}.qc_metrics.tsv', sep='\t', index=False, float_format='%.4f')
        EOF
        python3 script.py
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        noAddress: true
    }

    output {
        File qc_metrics = "~{output_basename}.qc_metrics.tsv"
    }
}
