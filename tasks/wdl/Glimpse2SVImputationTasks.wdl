version 1.0

task ExtractAnnotations {
    input {
        File imputed_vcf_or_bcf
        File imputed_vcf_or_bcf_index
        Int batch_index
        String? region

        Int disk_size_gb = ceil(2 * size(imputed_vcf_or_bcf, "GiB") + 50)
        Int mem_gb = 2
        Int cpu = 1
        Int preemptible = 3
    }

    command <<<
        set -euo pipefail

        # Ensure index is localized so bcftools can use it for random access if needed
        ls ~{imputed_vcf_or_bcf_index} > /dev/null

        printf 'CHROM\tPOS\tREF\tALT\tAF\tINFO\n' > annotations_batch_~{batch_index}.tsv
        bcftools query \
        ~{if defined(region) then "--regions-overlap 0 -r " + region else ""} \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/INFO\n' \
        ~{imputed_vcf_or_bcf} >> annotations_batch_~{batch_index}.tsv

        bgzip annotations_batch_~{batch_index}.tsv
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
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
        File merged_vcf_or_bcf
        File merged_vcf_or_bcf_idx
        Array[File] annotations
        Array[Int] num_samples
        String output_basename
        String? region
        Float info_filter_threshold = 0.0

        Int disk_size_gb = ceil(2.2 * size(merged_vcf_or_bcf, "GiB") + size(annotations, "GiB") + 50)
        Int mem_gb = 6
        Int cpu = 1
        Int preemptible = 3
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

# Pre-allocate sample array for broadcasting (shape: num_batches x 1)
samples_array = np.array(num_samples)[:, None]

# Define types to speed up CSV reading
dtypes = {'CHROM': str, 'POS': int, 'REF': str, 'ALT': str, 'AF': float, 'INFO': float}
readers = [pd.read_csv(f, sep='\t', chunksize=chunk_size, dtype=dtypes) for f in input_filenames]

def vectorized_sig_fig_round(x, n):
    """Pure NumPy vectorization for significant figure rounding (no Python loops)"""
    out = np.zeros_like(x)
    mask = (x != 0) & np.isfinite(x)
    x_mask = x[mask]
    
    # Calculate powers of 10 for rounding
    power = (n - 1) - np.floor(np.log10(np.abs(x_mask)))
    factor = 10.0 ** power
    
    # Scale, round, unscale
    out[mask] = np.round(x_mask * factor) / factor
    return out

with open('aggregated_annotations.tsv', 'w') as out:
    for chunks in zip(*readers):
        ref_loci = chunks[0][['CHROM', 'POS', 'REF', 'ALT']]
        for i, chunk in enumerate(chunks[1:], 1):
            if not ref_loci.equals(chunk[['CHROM', 'POS', 'REF', 'ALT']]):
                raise RuntimeError(f'Sites in chunk do not match between batch 0 and batch {i}.')

        # Stack into matrices (shape: num_batches x chunk_size) for vectorized math
        af_matrix = np.vstack([c['AF'].to_numpy() for c in chunks])
        info_matrix = np.vstack([c['INFO'].to_numpy() for c in chunks])

        # Vectorized weighted AF
        agg_af = np.sum(af_matrix * samples_array, axis=0) / total_samples

        # Vectorized weighted INFO
        numerator = np.sum((1 - info_matrix) * 2 * samples_array * af_matrix * (1 - af_matrix), axis=0)
        denominator = 2 * total_samples * agg_af * (1 - agg_af)
        
        polymorphic = (agg_af > 0.0) & (agg_af < 1.0)
        agg_info = np.where(
            polymorphic, 
            1.0 - np.divide(numerator, denominator, where=polymorphic, out=np.zeros_like(denominator)), 
            1.0
        )

        # Assign back using the fast rounding function
        result = ref_loci.copy()
        result['AF'] = vectorized_sig_fig_round(agg_af, 3)
        result['INFO'] = vectorized_sig_fig_round(agg_info, 3)
        
        result.to_csv(out, sep='\t', header=False, index=False)

EOF
        python3 script.py

        bgzip aggregated_annotations.tsv
        tabix -s1 -b2 -e2 aggregated_annotations.tsv.gz

        REGION_ARG=""
        if [ -n "~{region}" ]; then
            REGION_ARG="-r ~{region} --regions-overlap 0"
        fi

        # bcftools 1.18 only allows index-on-the-fly format to be specified via ##idx## notation: https://github.com/samtools/bcftools/issues/2008
        # Use quotes around the output argument to prevent bash from treating ## as a comment
        if awk -v t="~{info_filter_threshold}" 'BEGIN { exit !(t > 0.0) }'; then
            bcftools annotate -a aggregated_annotations.tsv.gz -c CHROM,POS,REF,ALT,AF,INFO ${REGION_ARG} -O u ~{merged_vcf_or_bcf} | \
                bcftools filter -i "INFO/INFO >= ~{info_filter_threshold}" -O z --write-index -o "~{output_basename}.vcf.gz##idx##~{output_basename}.vcf.gz.tbi"
        else
            bcftools annotate -a aggregated_annotations.tsv.gz -c CHROM,POS,REF,ALT,AF,INFO ${REGION_ARG} ~{merged_vcf_or_bcf} -O z --write-index -o "~{output_basename}.vcf.gz##idx##~{output_basename}.vcf.gz.tbi"
        fi
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/samtools-suite:v1.1"
        disks: "local-disk " + disk_size_gb + " SSD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
        noAddress: true
    }

    output {
        File merged_imputed_vcf = "~{output_basename}.vcf.gz"
        File merged_imputed_vcf_idx = "~{output_basename}.vcf.gz.tbi"
        File aggregated_annotations = "aggregated_annotations.tsv.gz"
    }
}

task ConcatAndFinalizeVcfs {
    input{
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_basename
        String? extra_args
    }

    Int disk_gb = ceil(2.1 * size(vcfs, "GiB")) + 10

    command <<<
        set -euox pipefail

        # bcftools 1.24 supports explicitly setting the index format via --write-index=tbi
        bcftools concat \
            -f ~{write_lines(vcfs)} \
            ~{extra_args} \
            -Oz --write-index=tbi -o ~{output_basename}.vcf.gz
            
        md5sum ~{output_basename}.vcf.gz | awk '{ print $1 }' > ~{output_basename}.md5sum
    >>>

    output {
        File concatenated_vcf = "~{output_basename}.vcf.gz"
        File concatenated_vcf_idx = "~{output_basename}.vcf.gz.tbi"
        File concatenated_vcf_md5sum = "~{output_basename}.md5sum"
    }

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: 3
        maxRetries: 0
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        noAddress: true
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

        required_cols = ['gvcf_path', 'gvcf_index_path']
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
        Array[String] gvcf_paths
        Array[String] gvcf_index_paths
        String output_filename = "manifest.tsv"
    }

    command <<<
        set -euo pipefail

        python3 << 'EOF'
        import sys

        gvcf_paths = ['~{sep="', '" gvcf_paths}']
        gvcf_index_paths = ['~{sep="', '" gvcf_index_paths}']

        if not (len(gvcf_paths) == len(gvcf_index_paths)):
            print(
                f"ERROR: Input arrays have different lengths: gvcf_paths={len(gvcf_paths)}, gvcf_index_paths={len(gvcf_index_paths)}",
                file=sys.stderr,
            )
            sys.exit(1)

        with open('~{output_filename}', 'w') as f:
            f.write("gvcf_path\tgvcf_index_path\n")
            for gvcf_path, gvcf_index_path in zip(gvcf_paths, gvcf_index_paths):
                f.write(f"{gvcf_path}\t{gvcf_index_path}\n")
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
        File output_gvcf_manifest = "~{output_filename}"
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

        required_cols = ['gvcf_path', 'gvcf_index_path']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Missing required columns in the VCF manifest: {', '.join(missing_cols)}.", file=sys.stderr)
            sys.exit(1)

        if df[required_cols].isnull().any().any():
            print("The VCF manifest contains empty values in required columns.", file=sys.stderr)
            sys.exit(1)

        df['gvcf_path'].to_csv('gvcf_paths.txt', index=False, header=False)
        df['gvcf_index_path'].to_csv('gvcf_index_paths.txt', index=False, header=False)
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
    }
}

task ConcatBcfs {
    input{
        Array[File] bcfs
        Array[File] bcf_idxs
        String output_basename
        String? extra_args
    }

    Int disk_gb = ceil(2.1 * size(bcfs, "GiB")) + 10

    command <<<
        set -euox pipefail

        bcftools concat \
            -f ~{write_lines(bcfs)} \
            ~{extra_args} \
            -Ob -o ~{output_basename}.bcf
        bcftools index ~{output_basename}.bcf
    >>>

    output {
        File concatenated_bcf = "~{output_basename}.bcf"
        File concatenated_bcf_idx = "~{output_basename}.bcf.csi"
    }

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: 3
        maxRetries: 0
        docker: "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
        noAddress: true
    }
}

