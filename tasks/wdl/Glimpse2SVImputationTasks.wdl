version 1.0

task MergeSampleChunksVcfsWithPaste {
    input {
        Array[File] input_vcfs
        String output_vcf_basename

        Int disk_size_gb = ceil(2.2 * size(input_vcfs, "GiB") + 50)
        Int mem_gb = 8
        Int cpu = 4
        Int preemptible = 0
    }

    command <<<
        set -euo pipefail

        vcfs=(~{sep=" " input_vcfs})

        mkfifo fifo_0
        mkfifo fifo_to_paste_0

        i=1

        fifos_to_paste=()
        md5sums=()
        # Keep only meta header lines here. The #CHROM line is merged in the paste stream.
        bcftools view -h --no-version ${vcfs[0]} | awk '!/^#CHROM/' > header.vcf
        n_lines=$(wc -l header.vcf | cut -d' ' -f1)

        mkfifo fifo_to_md5_0

        # Stream starting at #CHROM so sample-name columns are merged across batches.
        bcftools view --no-version ${vcfs[0]} > fifo_0 &
        tail +$((n_lines)) fifo_0 | tee fifo_to_md5_0 > fifo_to_paste_0 &
        tail -n +2 fifo_to_md5_0 | cut -f1-5,9 | md5sum > md5sum_0 &

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

        bcftools view --no-version $vcf > "$fifo_name" &
        tail +$((n_lines)) "$fifo_name" | tee "$fifo_name_to_md5" | cut -f 10- > "$fifo_name_to_paste" &
        tail -n +2 "$fifo_name_to_md5" | cut -f1-5,9 | md5sum > "$file_name_md5sum" &

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
        set -euo pipefail

        # Ensure index is localized so bcftools can use it for random access if needed
        ls ~{imputed_vcf_index} > /dev/null

        printf 'CHROM\tPOS\tREF\tALT\tAF\tINFO\n' > annotations_batch_~{batch_index}.tsv
        bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/INFO\n' \
        ~{imputed_vcf} >> annotations_batch_~{batch_index}.tsv

        bgzip annotations_batch_~{batch_index}.tsv
    >>>

    runtime {
        docker: docker_extract_annotations
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        noAddress: true
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
        Int preemptible = 0
        Int chunk_size = 100000
    }

    command <<<
        set -euo pipefail

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
        noAddress: true
    }

    output {
        File merged_imputed_vcf = "~{output_basename}.vcf.gz"
        File aggregated_annotations = "aggregated_annotations.tsv.gz"
    }
}

task CreateVcfIndexAndMd5 {
    input {
        File vcf_input
        String output_basename

        Int disk_size_gb = ceil(2.1*size(vcf_input, "GiB")) + 10
        Int cpu = 1
        Int memory_mb = 6000
        String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.5.0.0"
        Int preemptible = 3
    }

    command <<<
        set -euo pipefail

        if [[ "~{vcf_input}" == *.bcf ]]; then
            # Normalize BCF input to a bgzipped VCF for downstream compatibility.
            bcftools view -O z -o ~{output_basename}.vcf.gz ~{vcf_input}
        else
            ln -sf ~{vcf_input} ~{output_basename}.vcf.gz
        fi

        bcftools index -t ~{output_basename}.vcf.gz

        md5sum ~{output_basename}.vcf.gz | awk '{ print $1 }' > ~{output_basename}.md5sum
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
        File output_vcf = "~{output_basename}.vcf.gz"
        File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
        File output_vcf_md5sum = "~{output_basename}.md5sum"
    }
}

task FilterVcfByInfo {
    input {
        File vcf
        Float info_threshold
        String output_prefix

        String docker = "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        Int disk_size_gb = ceil(2.2 * size(vcf, "GiB") + 20)
        Int mem_gb = 4
        Int cpu = 1
        Int preemptible = 3
    }

    command <<<
        set -euo pipefail

        bcftools filter -i 'INFO/INFO >= ~{info_threshold}' -O z -o ~{output_prefix}.vcf.gz ~{vcf}
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
        File output_vcf = "~{output_prefix}.vcf.gz"
    }
}

task SplitVcfManifestIntoBatches {
    input {
        Int batch_size
        File gvcf_manifest
    }

    command <<<
        cat <<EOF > script.py
        import sys
        import pandas as pd

        batch_size = ~{batch_size}

        df = pd.read_csv("~{gvcf_manifest}", sep='\t')

        required_cols = ['sample_id', 'gvcf_path', 'gvcf_index_path']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Missing required columns in the VCF manifest: {', '.join(missing_cols)}.", file=sys.stderr)
            sys.exit(1)

        if df[required_cols].isnull().any().any():
            print("The VCF manifest contains empty values in required columns.", file=sys.stderr)
            sys.exit(1)

        if len(df) == 0:
            print("The VCF manifest must contain at least one row.", file=sys.stderr)
            sys.exit(1)

        chunk_num = 0
        for i in range(0, len(df), batch_size):
            df_chunk = df[i : i + batch_size]
            df_chunk.to_csv(f"chunk_{chunk_num:04d}.tsv", sep='\t', index=False)
            chunk_num += 1

        EOF
        python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        disks: "local-disk 10 HDD"
        memory: "1 GiB"
        preemptible: 3
        noAddress: true
    }

    output {
        Array[File] gvcf_manifest_batches = glob("chunk_*")
    }
}

task ConvertInputArraysToManifest {
    input {
        Array[String] sample_ids
        Array[String] gvcf_paths
        Array[String] gvcf_index_paths
        String output_filename = "manifest.tsv"
    }

    command <<<
        set -euo pipefail

        python3 << 'EOF'
        import sys

        sample_ids = ['~{sep="', '" sample_ids}']
        gvcf_paths = ['~{sep="', '" gvcf_paths}']
        gvcf_index_paths = ['~{sep="', '" gvcf_index_paths}']

        if not (len(sample_ids) == len(gvcf_paths) == len(gvcf_index_paths)):
            print(
                f"ERROR: Input arrays have different lengths: sample_ids={len(sample_ids)}, gvcf_paths={len(gvcf_paths)}, gvcf_index_paths={len(gvcf_index_paths)}",
                file=sys.stderr,
            )
            sys.exit(1)

        with open('~{output_filename}', 'w') as f:
            f.write("sample_id\tgvcf_path\tgvcf_index_path\n")
            for sample_id, gvcf_path, gvcf_index_path in zip(sample_ids, gvcf_paths, gvcf_index_paths):
                f.write(f"{sample_id}\t{gvcf_path}\t{gvcf_index_path}\n")
        EOF
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        noAddress: true
    }

    output {
        File output_manifest = "~{output_filename}"
    }
}

task ParseVcfManifestIntoArrays {
    input {
        File gvcf_manifest
    }

    command <<<
        set -euo pipefail

        cat <<EOF > script.py
        import sys
        import pandas as pd

        df = pd.read_csv("~{gvcf_manifest}", sep='\t')

        required_cols = ['sample_id', 'gvcf_path', 'gvcf_index_path']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Missing required columns in the VCF manifest: {', '.join(missing_cols)}.", file=sys.stderr)
            sys.exit(1)

        if df[required_cols].isnull().any().any():
            print("The VCF manifest contains empty values in required columns.", file=sys.stderr)
            sys.exit(1)

        df['gvcf_path'].to_csv('gvcf_paths.txt', index=False, header=False)
        df['gvcf_index_path'].to_csv('gvcf_index_paths.txt', index=False, header=False)
        df['sample_id'].to_csv('sample_ids.txt', index=False, header=False)
        EOF
        python3 script.py
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/python-data-slim:1.0"
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        noAddress: true
    }

    output {
        Array[File] input_gvcfs = read_lines("gvcf_paths.txt")
        Array[File] input_gvcf_idxs = read_lines("gvcf_index_paths.txt")
        Array[String] sample_ids = read_lines("sample_ids.txt")
    }
}

task ConcatBcfs {
    input{
        Array[File] bcfs
        Array[File] bcf_idxs
        String output_prefix
        String? extra_args
    }

    Int disk_gb = ceil(2.1 * size(bcfs, "GiB")) + 10

    command <<<
        set -euox pipefail

        bcftools concat \
            -f ~{write_lines(bcfs)} \
            ~{extra_args} \
            -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File concatenated_bcf = "~{output_prefix}.bcf"
        File concatenated_bcf_idx = "~{output_prefix}.bcf.csi"
    }

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_gb + " SSD"
        bootDiskSizeGb: 10
        preemptible: 3
        maxRetries: 0
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        noAddress: true
    }
}

