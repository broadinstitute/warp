version 1.0

struct VcfdistOutputs {
    File summary_vcf
    File precision_recall_summary_tsv
    File precision_recall_tsv
    File query_tsv
    File truth_tsv
    File phasing_summary_tsv
    File switchflips_tsv
    File superclusters_tsv
    File phase_blocks_tsv
}

workflow VcfdistEvaluation {
    input {
        Array[String] eval_sample_names
        Array[String] truth_sample_names    # "syndip" in Phase 1+2 dipcall truth VCFs; typically same as eval_sample_names for leave out

        # Phased eval VCF (all of its samples) and its index. The evaluation
        # samples are pulled out of it once, in SubsetEvalSamples.
        File eval_vcf
        File eval_vcf_idx
        # Tab-separated, no header required: sample name as in the eval VCF, then sex.
        # Must cover every name in eval_sample_names; other samples are ignored.
        # The sex values are copied as given into sample_manifest_tsv.
        File sample_sex_tsv

        Array[File] truth_vcfs
        Array[File]? truth_vcf_idxs

        # One BED per sample, or a single BED that is used for every sample
        Array[File] confident_regions_bed_files

        String region
        File reference_fasta
        File reference_fasta_fai

        Boolean do_naively_phase = false    # convert / to | in unphased eval genotypes
        Boolean naive_phase_truth = false   # convert / to | in unphased truth genotypes

        Array[File] vcfdist_bed_files
        Array[String] labels_per_stratification
        String? vcfdist_extra_args

        String summarize_evaluations_docker
    }

    # One pass over the eval VCF for all evaluation samples together
    call SubsetEvalSamples { input:
        eval_vcf = eval_vcf,
        eval_vcf_idx = eval_vcf_idx,
        eval_sample_names = eval_sample_names,
        truth_sample_names = truth_sample_names,
        sample_sex_tsv = sample_sex_tsv,
        region = region,
        reference_fasta_fai = reference_fasta_fai
    }

    scatter (i in range(length(eval_sample_names))) {
        File truth_vcf = truth_vcfs[i]
        File confident_regions_bed = if length(confident_regions_bed_files) == 1 then confident_regions_bed_files[0] else confident_regions_bed_files[i]

        call SubsetSampleFromVcf as SubsetSampleFromVcfEval { input:
            vcf = SubsetEvalSamples.subset_vcf,
            vcf_idx = SubsetEvalSamples.subset_vcf_tbi,
            original_sample_name = eval_sample_names[i],
            sample_name = eval_sample_names[i],
            region = region,
            bed_file = confident_regions_bed,
            reference_fasta_fai = reference_fasta_fai,
            do_naively_phase = do_naively_phase
        }

        call SubsetSampleFromVcf as SubsetSampleFromVcfTruth { input:
            vcf = truth_vcfs[i],
            vcf_idx = if defined(truth_vcf_idxs) then select_first([truth_vcf_idxs])[i]  else truth_vcf_idxs,
            original_sample_name = truth_sample_names[i],
            sample_name = eval_sample_names[i],     # rename truth to match eval
            region = region,
            bed_file = confident_regions_bed,
            reference_fasta_fai = reference_fasta_fai,
            do_naively_phase = naive_phase_truth
        }
    }

    scatter (j in range(length(vcfdist_bed_files))) {
        scatter (i in range(length(eval_sample_names))) {
            call Vcfdist { input:
                sample_name = eval_sample_names[i],
                eval_vcf = SubsetSampleFromVcfEval.single_sample_vcf[i],
                truth_vcf = SubsetSampleFromVcfTruth.single_sample_vcf[i],
                bed_file = vcfdist_bed_files[j],
                reference_fasta = reference_fasta,
                extra_args = vcfdist_extra_args
            }
        }
    }

    call SummarizeEvaluations { input:
        labels_per_vcf = labels_per_stratification,
        vcfdist_outputs_per_vcf_and_sample = Vcfdist.outputs
    }

    output {
        # stratification x sample
        Array[Array[VcfdistOutputs]] vcfdist_summary = Vcfdist.outputs
        File evaluation_summary_tsv = SummarizeEvaluations.evaluation_summary_tsv

        # Evaluation samples only, cut from the eval VCF
        File eval_subset_vcf = SubsetEvalSamples.subset_vcf
        File eval_subset_vcf_tbi = SubsetEvalSamples.subset_vcf_tbi
        # eval_sample, truth_sample, sex: one row per sample, in input order
        File sample_manifest_tsv = SubsetEvalSamples.sample_manifest_tsv
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    String? disk_type
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task SubsetEvalSamples {
    input {
        File eval_vcf
        File eval_vcf_idx
        Array[String] eval_sample_names
        Array[String] truth_sample_names
        File sample_sex_tsv
        String region
        File reference_fasta_fai

        # Number of pieces the region is read in at once; keep at or below cpu_cores
        Int n_parallel = 16

        RuntimeAttr? runtime_attr_override
    }

    # The extra space also raises the persistent disk's read speed
    Int disk_gb = ceil(size(eval_vcf, "GiB") * 1.2) + 200

    command <<<
        set -euxo pipefail

        # ---- Manifest: eval sample, truth sample, sex. Fails if a sample has no sex.
        paste ~{write_lines(eval_sample_names)} ~{write_lines(truth_sample_names)} > pairs.tsv
        awk -F'\t' -v OFS='\t' '
            NR == FNR { if ($0 !~ /^#/ && NF >= 2) sex[$1] = $2; next }
            {
                if (!($1 in sex)) { print "no sex given for sample " $1 > "/dev/stderr"; bad = 1; next }
                print $1, $2, sex[$1]
            }
            END { exit bad }
        ' ~{sample_sex_tsv} pairs.tsv > manifest_body.tsv
        printf 'eval_sample\ttruth_sample\tsex\n' | cat - manifest_body.tsv > sample_manifest.tsv

        # ---- Region bounds (a bare contig name means the whole contig)
        region="~{region}"
        if [[ "$region" == *:* ]]; then
            chrom=${region%%:*}
            range=${region#*:}
            range=${range//,/}
            start=${range%-*}
            end=${range#*-}
        else
            chrom=$region
            start=1
            end=$(awk -F'\t' -v c="$chrom" '$1 == c { print $2 }' ~{reference_fasta_fai})
        fi

        # ---- Subset every evaluation sample in one pass over the eval VCF.
        # The region is cut into equal pieces that are read in parallel through
        # the index. --regions-overlap 0 keeps a record only in the piece that
        # contains its POS, so nothing is duplicated at the piece edges.
        n=~{n_parallel}
        step=$(( (end - start + n) / n ))
        pids=()
        : > pieces.txt
        for (( k = 0; k < n; k++ )); do
            s=$(( start + k * step ))
            e=$(( s + step - 1 ))
            if (( s > end )); then break; fi
            if (( e > end )); then e=$end; fi
            out=$(printf 'piece_%03d.bcf' "$k")
            bcftools view "~{eval_vcf}##idx##~{eval_vcf_idx}" \
                -r "$chrom:$s-$e" --regions-overlap 0 \
                -S ~{write_lines(eval_sample_names)} \
                -Ob -o "$out" &
            pids+=($!)
            echo "$out" >> pieces.txt
        done
        for pid in "${pids[@]}"; do
            wait "$pid"
        done

        bcftools concat --file-list pieces.txt -Oz --write-index=tbi -o eval.subset.vcf.gz
    >>>

    output {
        File subset_vcf = "eval.subset.vcf.gz"
        File subset_vcf_tbi = "eval.subset.vcf.gz.tbi"
        File sample_manifest_tsv = "sample_manifest.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "us.gcr.io/broad-gotc-prod/bcftools-vcftools:2.0.0-1.24-0.1.17-1784569943"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SubsetSampleFromVcf {
    input {
        File vcf
        File? vcf_idx
        String original_sample_name
        String sample_name
        String region
        File? bed_file
        File reference_fasta_fai
        Boolean do_naively_phase = false    # convert / to | in unphased inputs

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size([vcf, bed_file], "GiB"))

    command <<<
        set -euxo pipefail

        if ! ~{defined(vcf_idx)}; then
            bcftools index ~{vcf}
        fi

        # must use -T bed_file to intersect with -r region properly
        bcftools view ~{vcf}~{if defined(vcf_idx) then "##idx##" + vcf_idx else ""} \
            -s ~{original_sample_name} \
            -r ~{region} \
            ~{"-T " + bed_file} \
            ~{if (do_naively_phase) then "-Ou | bcftools +setGT -Oz -o " + sample_name + ".subset.vcf.gz -- -t a -n p" else "-Oz -o " + sample_name + ".subset.vcf.gz"}
        echo ~{sample_name} > sample_name.txt
        bcftools reheader ~{sample_name}.subset.vcf.gz \
            -s sample_name.txt \
            --fai ~{reference_fasta_fai} \
            -o ~{sample_name}.subset.reheadered.vcf.gz
        bcftools index -t ~{sample_name}.subset.reheadered.vcf.gz
    >>>
    
    output {
        File single_sample_vcf = "~{sample_name}.subset.reheadered.vcf.gz"
        File single_sample_vcf_tbi = "~{sample_name}.subset.reheadered.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "LOCAL",
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Vcfdist {
    input {
        String sample_name
        File eval_vcf
        File truth_vcf
        File bed_file
        File reference_fasta
        String? extra_args
        Int verbosity = 1
        
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size([eval_vcf, truth_vcf, reference_fasta], "GiB"))

    command <<<
        set -euxo pipefail

        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{reference_fasta} \
            -b ~{bed_file} \
            -v ~{verbosity} \
            ~{extra_args}

        for tsv in $(ls *.tsv); do mv $tsv ~{sample_name}.$tsv; done
        mv summary.vcf ~{sample_name}.summary.vcf
    >>>

    output {
        VcfdistOutputs outputs = {
            "summary_vcf": "~{sample_name}.summary.vcf",
            "precision_recall_summary_tsv": "~{sample_name}.precision-recall-summary.tsv",
            "precision_recall_tsv": "~{sample_name}.precision-recall.tsv",
            "query_tsv": "~{sample_name}.query.tsv",
            "truth_tsv": "~{sample_name}.truth.tsv",
            "phasing_summary_tsv": "~{sample_name}.phasing-summary.tsv",
            "switchflips_tsv": "~{sample_name}.switchflips.tsv",
            "superclusters_tsv": "~{sample_name}.superclusters.tsv",
            "phase_blocks_tsv": "~{sample_name}.phase-blocks.tsv"
        }
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "timd1/vcfdist:v2.6.4"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SummarizeEvaluations {
    input {
        Array[String] labels_per_vcf
        Array[Array[VcfdistOutputs]] vcfdist_outputs_per_vcf_and_sample

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        pip install pandas

        python - --labels_per_vcf_txt ~{write_lines(labels_per_vcf)} \
                 --vcfdist_outputs_per_vcf_and_sample_json ~{write_json(vcfdist_outputs_per_vcf_and_sample)} \
                 <<-'EOF'
        import argparse
        import json
        import pandas as pd

        def summarize(labels_per_vcf_txt,
                      vcfdist_outputs_per_vcf_and_sample_json):
            with open(labels_per_vcf_txt) as f:
                labels = f.read().splitlines()

            with open(vcfdist_outputs_per_vcf_and_sample_json) as f:
                vcfdist_outputs_per_vcf_and_sample = json.load(f)

            summary_dict = {}
            for i, label in enumerate(labels):
                summary_dict[label] = {}
                summary_dict[label]['NUM_VCFDIST_SAMPLES'] = len(vcfdist_outputs_per_vcf_and_sample[i])
                summary_dict[label].update(summarize_vcfdist_outputs_over_samples(vcfdist_outputs_per_vcf_and_sample[i]))

            pd.DataFrame.from_dict(summary_dict, orient='index').to_csv('evaluation_summary.tsv', sep='\t', float_format="%.4f")

        def summarize_vcfdist_outputs_over_samples(vcfdist_outputs_per_sample):
            precision_recall_metrics_dict = {}
            for s, vcfdist_outputs in enumerate(vcfdist_outputs_per_sample):
                precision_recall_metrics_dict[s] = {}
                pr_metrics_df = pd.read_csv(vcfdist_outputs['precision_recall_summary_tsv'], sep='\t', index_col=[0, 1])
                for var_type in ['SNP', 'INDEL', 'SV']:
                    var_type_metrics_dict = pr_metrics_df.loc[var_type, 'NONE'][['TRUTH_TP', 'QUERY_TP', 'TRUTH_FN', 'QUERY_FP', 'PREC', 'RECALL', 'F1_SCORE']].add_prefix(f'{var_type}_').to_dict()
                    precision_recall_metrics_dict[s].update(var_type_metrics_dict)
            return pd.DataFrame.from_dict(precision_recall_metrics_dict, orient='index').mean(axis=0)


        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--labels_per_vcf_txt',
                                type=str)

            parser.add_argument('--vcfdist_outputs_per_vcf_and_sample_json',
                                type=str)

            args = parser.parse_args()

            summarize(args.labels_per_vcf_txt,
                      args.vcfdist_outputs_per_vcf_and_sample_json)

        if __name__ == '__main__':
            main()
        EOF
    >>>

    output {
        File evaluation_summary_tsv = "evaluation_summary.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        disk_type:          "HDD",
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
