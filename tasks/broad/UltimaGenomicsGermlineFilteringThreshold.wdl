version 1.0

import "../../tasks/broad/JointGenotypingTasks.wdl" as Tasks

# Given a joint callset with a "score_key" INFO level annotation, this pipeline chooses a threshold that maximizes
# the F1 score on a single sample in the callset using known truth data for that one sample.
workflow ExtractOptimizeSingleSample { 
    input {
        Array[File] input_vcf
        Array[File] input_vcf_index
        String base_file_name
        String sample_name_calls
        File gtr_vcf
        File gtr_vcf_index
        File gtr_highconf_intervals
        String sample_name_gtr
        String flow_order
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_fasta_sdf
        File runs_file
        Array[File] annotation_intervals
        String score_key = "SCORE"
        Int medium_disk
    }

    scatter (idx in range(length(input_vcf))) {
        call AnnotateSampleVCF {
            input:
                input_vcf = input_vcf[idx],
                input_vcf_index = input_vcf_index[idx],
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                output_basename = base_file_name
        }

        call ExtractSample {
            input:
                sample_name = sample_name_calls,
                input_vcf = AnnotateSampleVCF.output_vcf_file,
                input_vcf_index = AnnotateSampleVCF.output_vcf_index
        }
    }

    call Tasks.GatherVcfs {
        input:
            input_vcfs = ExtractSample.output_vcf_file,
            output_vcf_name = base_file_name + ".extracted.vcf.gz",
            disk_size = medium_disk
    }

    call FilterSampleVCF{
        input:
            input_vcf = GatherVcfs.output_vcf
    }

    call FilterSymbolicAlleles{
        input:
            base_file_name = base_file_name,
            input_vcf = FilterSampleVCF.output_vcf_file,
            input_vcf_index = FilterSampleVCF.output_vcf_index
         }

    call CompareToGroundTruth{
        input:
            gtr_vcf = gtr_vcf,
            gtr_vcf_index = gtr_vcf_index, 
            gtr_highconf = gtr_highconf_intervals,
            left_sample_name = sample_name_calls,
            right_sample_name = sample_name_gtr,
            input_vcf = FilterSymbolicAlleles.output_vcf,
            input_vcf_index = FilterSymbolicAlleles.output_vcf_index,
            input_vcf_name = base_file_name,
            
            ref_fasta = ref_fasta,
            ref_dict = ref_dict,
            ref_fasta_sdf = ref_fasta_sdf,
            runs_file = runs_file,
            annotation_intervals = annotation_intervals,

            disk_size = 2*(ceil(size(input_vcf, "GB") + size(gtr_vcf, "GB") + size(ref_fasta, "GB") + size(ref_fasta_sdf, "GB")))+10,
            flow_order = flow_order

    }

    call EvaluateResults {
        input:
          h5_input = CompareToGroundTruth.compare_h5,
          input_vcf_name = base_file_name,
          disk_size = 10,
          score_key = score_key
    }

    scatter (idx in range(length(input_vcf))) {
        call HardThresholdVCF {
            input:
              thresholds = EvaluateResults.thresholds_report,
              input_vcf = AnnotateSampleVCF.output_vcf_file[idx],
              input_vcf_index = AnnotateSampleVCF.output_vcf_index[idx],
              score_key = score_key,
              output_basename = base_file_name,
              disk_size = 3*ceil(size(AnnotateSampleVCF.output_vcf_file[idx], "GB")) + 14
        }
    }

    output {
        Array[File] output_vcf = HardThresholdVCF.output_vcf
        Array[File] output_vcf_index = HardThresholdVCF.output_vcf_index
        File eval_report_h5 = EvaluateResults.eval_report_h5
        File thresholds_report = EvaluateResults.thresholds_report
    }
}

task ExtractSample {
    input {
        String sample_name
        File input_vcf
        File input_vcf_index
        String docker = "gcr.io/terra-project-249020/jukebox_vc:test_jc_optimize_3d7509"
    }
    String output_vcf = basename(input_vcf, ".vcf.gz") + ".~{sample_name}.vcf.gz"
    command <<<
        set -eo pipefail
        source ~/.bashrc
        conda activate genomics.py3
        bcftools view -s ~{sample_name} ~{input_vcf} -o  tmp.vcf.gz -O z 
        bcftools index -t tmp.vcf.gz
        bcftools view -h tmp.vcf.gz | sed 's/COMBINED_TREE_SCORE,Number=A,/COMBINED_TREE_SCORE,Number=\.,/g' > hdr.fixed.txt
        bcftools reheader -h hdr.fixed.txt tmp.vcf.gz | bcftools view -Oz -o ~{output_vcf} - 
        bcftools index -t ~{output_vcf}
    >>>
    output {
        File output_vcf_file = "~{output_vcf}"
        File output_vcf_index = "~{output_vcf}.tbi"
    }
    runtime {
        memory: "8GB"
        disks: "local-disk " + (ceil(size(input_vcf, "GB")) * 3 + 10) + " HDD"
        docker: docker
        cpu: 1
    }
}

task FilterSampleVCF{
    input{
        File input_vcf
        String docker = "gcr.io/terra-project-249020/jukebox_vc:test_jc_optimize_3d7509"
    }

    String output_vcf = basename(input_vcf, ".vcf.gz") + ".filter_unused.vcf.gz"

    command <<<
        set -e
        source ~/.bashrc
        conda activate genomics.py3
        bcftools view -a -i 'GT[*]="alt"' ~{input_vcf} -o ~{output_vcf} -O z 
        bcftools index -t ~{output_vcf}
    >>>

    output {
        File output_vcf_file = "~{output_vcf}"
        File output_vcf_index = "~{output_vcf}.tbi"
    }

    runtime {
        memory: "8GB"
        disks: "local-disk " + (ceil(size(input_vcf, "GB")) *2 + 10) + " HDD"
        docker: docker
        cpu: 1
    }
}

task FilterSymbolicAlleles {
    input {
        String base_file_name
        File input_vcf
        File input_vcf_index
        String docker = "gcr.io/terra-project-249020/jukebox_vc:test_jc_optimize_3d7509"
    }

    String output_vcf_name = "~{base_file_name}" + ".vcf.gz"
    command <<<
        set -e
        source ~/.bashrc
        conda activate genomics.py3
        gatk --java-options "-Xmx10g"  SelectVariants \
            -V ~{input_vcf} \
            -O ~{output_vcf_name}.tmp.vcf.gz \
            --remove-unused-alternates
        gatk --java-options "-Xmx10g" SelectVariants \
            -V ~{output_vcf_name}.tmp.vcf.gz \
            -O ~{output_vcf_name} \
            --exclude-non-variants \
            --select-type-to-exclude SYMBOLIC
        >>>
    runtime {
        memory: "12 GB"
        cpu: 1
        disks: "local-disk " + (ceil(size(input_vcf, "GB")) *4 +10) + " HDD"
        docker: docker
    }
    output {
        File output_vcf = "~{output_vcf_name}"
        File output_vcf_index = "~{output_vcf_name}.tbi"
        File monitoring_log = "monitoring.log"
    }
}

task CompareToGroundTruth {
  input {
    String left_sample_name
    String right_sample_name
    File input_vcf
    File input_vcf_index

    String input_vcf_name

    File gtr_vcf
    File gtr_vcf_index
    File gtr_highconf

    File? interval_list
    File runs_file
    File ref_fasta
    File ref_dict
    File ref_fasta_sdf
    String flow_order
    Array[File] annotation_intervals
    Int disk_size
    String docker = "gcr.io/terra-project-249020/jukebox_vc:test_jc_optimize_3d7509"
  }

  String used_flow_order = (if flow_order=="" then "TACG" else flow_order)

  command <<<
    set -e

    source ~/.bashrc
    conda activate genomics.py3

    python -m tarfile -e ~{ref_fasta_sdf} ~{ref_fasta}.sdf

    run_comparison_pipeline.py \
            --n_parts 0 \
            --hpol_filter_length_dist 12 10 \
            --input_prefix $(echo "~{input_vcf}" | sed 's/\(.vcf.gz\|.vcf\)$//') \
            --output_file ~{input_vcf_name}.comp.h5 \
            --gtr_vcf ~{gtr_vcf} \
            --highconf_intervals ~{gtr_highconf} \
            --runs_intervals ~{runs_file} \
            --reference ~{ref_fasta} \
            --reference_dict ~{ref_dict} \
            --call_sample_name ~{left_sample_name} \
            --ignore_filter_status \
            --flow_order ~{used_flow_order} \
            --truth_sample_name ~{right_sample_name}\
            --annotate_intervals ~{sep=" --annotate_intervals " annotation_intervals} \
            --n_jobs 8 \
            --output_suffix '' \
            --output_interval interval_~{input_vcf_name}.comp.bed
    >>>
  runtime {
    memory: "32 GB"
    cpu: 16
    disks: "local-disk " + disk_size + " SSD"
    docker: docker
  }
  output {
    File compare_h5 = "~{input_vcf_name}.comp.h5"
    Array[File] comparison_beds = glob("~{input_vcf_name}.*.bed")
    File output_interval = "interval_~{input_vcf_name}.comp.bed"
  }
}

task EvaluateResults {
  input {
    File h5_input
    String input_vcf_name
    Int disk_size
    String docker = "gcr.io/terra-project-249020/jukebox_vc:test_jc_optimize_3d7509"
    String score_key
  }
  command <<<
    set -e

    source ~/.bashrc
    conda activate genomics.py3

    evaluate_concordance.py \
            --input_file ~{h5_input} \
            --output_prefix ~{input_vcf_name}.report \
            --use_for_group_testing variant_type \
            --score_key ~{score_key}

    >>>
  runtime {
    memory: "32 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
  }
  output {
    File eval_report_h5 = "~{input_vcf_name}.report.h5"
    File thresholds_report = "~{input_vcf_name}.report.thresholds.txt"
  }
}

task HardThresholdVCF {
  input { 
    File thresholds
    File input_vcf
    File input_vcf_index
    String output_basename
    String score_key
    Int disk_size
    String docker = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds:mshand-05a76349aac401cc5d08ad0375c45fa9d4c7c864-4.2.6.1-48-g05a76349a"
  }

  command <<<
    set -eo pipefail
    
    tail -n +2 ~{thresholds} > ~{thresholds}.body
    cat ~{thresholds}.body | awk 'BEGIN { FS=","} \
    {print "-filter \"VARIANT_TYPE == \047" $1 "\047 && ~{score_key} \
     < "$2 "\" --filter-name LOW_SCORE_"$1 " "}' | tr '\n' ' ' > tmpfile
    
    echo "Contents of tmpfile"
    cat tmpfile

    params=$( cat tmpfile )
    echo "Params"
    echo $params

    echo gatk --java-options "-Xmx20g"  VariantFiltration \
    -V ~{input_vcf} \
    -O ~{output_basename}.vcf.gz \
    "$params" >tmpfile
    . ./tmpfile
  >>>
  runtime {
    memory: "32 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: docker
  }

  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task AnnotateSampleVCF {
    input {
        File input_vcf
        File input_vcf_index
        String output_basename
        Int disk_size = ceil(size(input_vcf, "GB") * 2) + 50
        String docker = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds:mshand-05a76349aac401cc5d08ad0375c45fa9d4c7c864-4.2.6.1-48-g05a76349a"
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String flow_order = "TGCA"
    }

    command <<<
        gatk --java-options "-Xmx15g" \
            VariantAnnotator \
            -O ~{output_basename}.vcf.gz \
            -V ~{input_vcf} \
            -R ~{ref_fasta} \
            -G StandardFlowBasedAnnotation \
            --flow-order-for-annotations ~{flow_order}
    >>>
    runtime {
        memory: "16 GB"
        disks: "local-disk " + disk_size + " HDD"
        docker: docker
    }

    output {
        File output_vcf_file = "~{output_basename}.vcf.gz"
        File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
    }
}



