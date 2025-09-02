version 1.0

task CalculateChromosomeLength {
  input {
    File ref_dict
    String chrom

    String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
    Int memory_mb = 2000
    Int cpu = 1
    Int disk_size_gb = ceil(2*size(ref_dict, "GiB")) + 5
  }

  command {
    set -e -o pipefail

    grep -P "SN:~{chrom}\t" ~{ref_dict} | sed 's/.*LN://' | sed 's/\t.*//'
  }
  runtime {
    docker: ubuntu_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    maxRetries: 1
    noAddress: true
  }
  output {
    Int chrom_length = read_int(stdout())
  }
}

task GetMissingContigList {
  input {
    File ref_dict
    File included_contigs

    String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
    Int memory_mb = 2000
    Int cpu = 1
    Int disk_size_gb = ceil(2*size(ref_dict, "GiB")) + 5
  }

  command <<<
    set -e -o pipefail

    grep "@SQ" ~{ref_dict} | sed 's/.*SN://' | sed 's/\t.*//' > contigs.txt

    awk 'NR==FNR{arr[$0];next} !($0 in arr)' ~{included_contigs} contigs.txt > missing_contigs.txt
  >>>

  runtime {
    docker: ubuntu_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    noAddress: true
  }

  output {
    Array[String] missing_contigs = read_lines("missing_contigs.txt")
  }
}

task GenerateChunk {
  input {
    Int start
    Int end
    String chrom
    String basename
    File vcf
    File vcf_index

    Int disk_size_gb = ceil(2*size(vcf, "GiB")) + 10
    Int cpu = 1
    Int memory_mb = 8000
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
  }
  Int command_mem = memory_mb - 1500
  Int max_heap = memory_mb - 1000

  command {
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    SelectVariants \
    -V ~{vcf} \
    --select-type-to-include SNP \
    --max-nocall-fraction 0.1 \
    -xl-select-type SYMBOLIC \
    --select-type-to-exclude MIXED \
    --restrict-alleles-to BIALLELIC \
    -L ~{chrom}:~{start}-~{end} \
    -O ~{basename}.vcf.gz \
    --exclude-filtered true \
    -select 'POS >= ~{start}'
  }
  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    maxRetries: 1
    noAddress: true
  }
  parameter_meta {
    vcf: {
      description: "vcf",
      localization_optional: true
    }
    vcf_index: {
      description: "vcf index",
      localization_optional: true
    }
  }
  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"
  }
}

task CountVariantsInChunks {
  input {
    File vcf
    File vcf_index
    File panel_vcf
    File panel_vcf_index

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int cpu = 1
    Int memory_mb = 6000
    Int disk_size_gb = 2 * ceil(size([vcf, vcf_index, panel_vcf, panel_vcf_index], "GiB")) + 20
  }
  Int command_mem = memory_mb - 1500
  Int max_heap = memory_mb - 1000

  command <<<
    set -e -o pipefail

    echo $(gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" CountVariants -V ~{vcf}  | sed 's/Tool returned://') > var_in_original
    echo $(gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" CountVariants -V ~{vcf} -L ~{panel_vcf} | sed 's/Tool returned://') > var_in_reference
  >>>
  output {
    Int var_in_original = read_int("var_in_original")
    Int var_in_reference = read_int("var_in_reference")
  }
  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    noAddress: true
  }
}

task CheckChunks {
  input {
    File vcf
    File vcf_index
    File panel_vcf
    File panel_vcf_index
    Int var_in_original
    Int var_in_reference

    Int disk_size_gb = ceil(2*size([vcf, vcf_index, panel_vcf, panel_vcf_index], "GiB")) + 10
    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 4000
  }
  command <<<
    set -e -o pipefail

    if [ $(( ~{var_in_reference} * 2 - ~{var_in_original})) -gt 0 ] && [ ~{var_in_reference} -gt 3 ]; then
      echo true > valid_file.txt
    else
      echo false > valid_file.txt
    fi

    bcftools convert -Ob ~{vcf} > valid_variants.bcf
    bcftools index -f valid_variants.bcf
  >>>
  output {
    File valid_chunk_bcf = "valid_variants.bcf"
    File valid_chunk_bcf_index = "valid_variants.bcf.csi"
    Boolean valid = read_boolean("valid_file.txt")
  }
  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    noAddress: true
  }
}

task PhaseVariantsEagle {
  input {
    File dataset_bcf
    File dataset_bcf_index
    File reference_panel_bcf
    File reference_panel_bcf_index
    String chrom
    File genetic_map_file
    Int start
    Int end

    String eagle_docker = "us.gcr.io/broad-gotc-prod/imputation-eagle:1.0.0-2.4-1690199702"
    Int cpu = 8
    Int memory_mb = 32000
    Int disk_size_gb = ceil(3 * size([dataset_bcf, reference_panel_bcf, dataset_bcf_index, reference_panel_bcf_index], "GiB")) + 10
  }
  command <<<
    /usr/gitc/eagle  \
    --vcfTarget ~{dataset_bcf}  \
    --vcfRef ~{reference_panel_bcf} \
    --geneticMapFile ~{genetic_map_file} \
    --outPrefix pre_phased_~{chrom} \
    --vcfOutFormat z \
    --bpStart ~{start} \
    --bpEnd ~{end} \
    --allowRefAltSwap
  >>>
  output {
    File dataset_prephased_vcf="pre_phased_~{chrom}.vcf.gz"
  }
  runtime {
    docker: eagle_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    noAddress: true
  }
}

task Minimac4 {
  input {
    File ref_panel
    File phased_vcf
    String prefix
    String chrom
    Int start
    Int end
    Int window

    String minimac4_docker = "us.gcr.io/broad-gotc-prod/imputation-minimac4:1.0.6-1.0.2-1663948783"
    Int cpu = 1
    Int memory_mb = 4000
    Int disk_size_gb = ceil(size(ref_panel, "GiB") + 2*size(phased_vcf, "GiB")) + 50
  }
  command <<<
    set -e -o pipefail

    /usr/gitc/minimac4 \
    --refHaps ~{ref_panel} \
    --haps ~{phased_vcf} \
    --start ~{start} \
    --end ~{end} \
    --window ~{window} \
    --chr ~{chrom} \
    --noPhoneHome \
    --format GT,DS,GP \
    --allTypedSites \
    --prefix ~{prefix} \
    --minRatio 0.00001

    bcftools index -t ~{prefix}.dose.vcf.gz
  >>>
  output {
    File vcf = "~{prefix}.dose.vcf.gz"
    File vcf_index = "~{prefix}.dose.vcf.gz.tbi"
    File info = "~{prefix}.info"
  }
  runtime {
    docker: minimac4_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    noAddress: true
  }
}

task GatherVcfs {
  input {
    Array[File] input_vcfs
    String output_vcf_basename

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int cpu = 1
    Int memory_mb = 16000
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

    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    IndexFeatureFile -I ~{output_vcf_basename}.vcf.gz
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
    File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}

task ReplaceHeader {
  input {
    File vcf_to_replace_header
    File vcf_with_new_header

    Int cpu = 1
    Int memory_mb = 6000
    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
  }

  String output_name = basename(vcf_to_replace_header,".vcf.gz") + ".new_header.vcf.gz"
  Int disk_size_gb = ceil(4*(size(vcf_to_replace_header, "GiB") + size(vcf_with_new_header, "GiB"))) + 20

  command <<<
    set -e -o pipefail

    bcftools view -h ~{vcf_with_new_header} > header.hr

    bcftools reheader -h header.hr ~{vcf_to_replace_header} -o ~{output_name}
  >>>

  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }

  output {
    File output_vcf = "~{output_name}"
  }
}

task UpdateHeader {
  input {
    File vcf
    File vcf_index
    File ref_dict
    String basename
    Boolean disable_sequence_dictionary_validation = true
    String? pipeline_header_line

    Int disk_size_gb = ceil(4*(size(vcf, "GiB") + size(vcf_index, "GiB"))) + 20
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int cpu = 1
    Int memory_mb = 6000
  }
  Int command_mem = memory_mb - 1500
  Int max_heap = memory_mb - 1000
  String disable_sequence_dict_validation_flag = if disable_sequence_dictionary_validation then "--disable-sequence-dictionary-validation" else ""

  command <<<
    set -e -o pipefail

    ## update the header of the merged vcf with the reference dictionary
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    UpdateVCFSequenceDictionary \
    --source-dictionary ~{ref_dict} \
    --output ~{basename}.vcf.gz \
    --replace -V ~{vcf} \
    ~{disable_sequence_dict_validation_flag}

    ## update header with pipeline name and version if provided
    if [ -n ~{default="" pipeline_header_line} ]; then
      mv ~{basename}.vcf.gz temp.vcf.gz
      mv ~{basename}.vcf.gz.tbi temp.vcf.gz.tbi

      bcftools view -h --no-version temp.vcf.gz > header.txt
      TOTAL_LINES=$(wc -l < "header.txt")
      sed -i "${TOTAL_LINES}i\##~{pipeline_header_line}" header.txt

      bcftools reheader -h header.txt -Oz -o ~{basename}.vcf.gz temp.vcf.gz
      bcftools index -t ~{basename}.vcf.gz
    fi
  >>>
  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
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

task RemoveSymbolicAlleles {
  input {
    File original_vcf
    File original_vcf_index
    String output_basename

    Int disk_size_gb = ceil(3*(size(original_vcf, "GiB") + size(original_vcf_index, "GiB")))
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int cpu = 1
    Int memory_mb = 4000
  }
  Int command_mem = memory_mb - 1500
  Int max_heap = memory_mb - 1000

  command {
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    SelectVariants -V ~{original_vcf} -xl-select-type SYMBOLIC -O ~{output_basename}.vcf.gz
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
}

task SeparateMultiallelics {
  input {
    File original_vcf
    File original_vcf_index
    String output_basename

    Int disk_size_gb =  ceil(2*(size(original_vcf, "GiB") + size(original_vcf_index, "GiB"))) + 10
    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 4000
  }
  command {
    set -e -o pipefail

    bcftools norm -m - ~{original_vcf} -Oz -o ~{output_basename}.vcf.gz
    bcftools index -t ~{output_basename}.vcf.gz
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
}

task OptionalQCSites {
  input {
    File input_vcf
    File input_vcf_index
    String output_vcf_basename
    Float? optional_qc_max_missing
    Float? optional_qc_hwe

    String bcftools_vcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 16000
    Int disk_size_gb = ceil(2*(size(input_vcf, "GiB") + size(input_vcf_index, "GiB"))) + 10

  }
  Float max_missing = select_first([optional_qc_max_missing, 0.05])
  Float hwe = select_first([optional_qc_hwe, 0.000001])
  command <<<
    set -e -o pipefail

    ln -sf ~{input_vcf} input.vcf.gz
    ln -sf ~{input_vcf_index} input.vcf.gz.tbi

    # site missing rate < 5% ; hwe p > 1e-6
    vcftools --gzvcf input.vcf.gz --max-missing ~{max_missing} --hwe ~{hwe} --recode -c | bgzip -c > ~{output_vcf_basename}.vcf.gz
    bcftools index -t ~{output_vcf_basename}.vcf.gz # Note: this is necessary because vcftools doesn't have a way to output a zipped vcf, nor a way to index one (hence needing to use bcf).
  >>>
  runtime {
    docker: bcftools_vcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    noAddress: true
  }
  output {
    File output_vcf = "~{output_vcf_basename}.vcf.gz"
    File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}

task MergeSingleSampleVcfs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcf_indices
    String output_vcf_basename

    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int memory_mb = 2000
    Int cpu = 1
    Int disk_size_gb = 3 * ceil(size(input_vcfs, "GiB") + size(input_vcf_indices, "GiB")) + 20
  }
  command <<<
    set -e -o pipefail

    # Move the index file next to the vcf with the corresponding name

    declare -a VCFS=(~{sep=' ' input_vcfs})
    declare -a VCF_INDICES=(~{sep=' ' input_vcf_indices})

    for i in ${VCF_INDICES[@]}; do
      for v in ${VCFS[@]}; do
        if [[ $(basename $i .vcf.gz.tbi) == $(basename $v .vcf.gz) ]]; then
          mv $i $(dirname $v)
        fi
      done
    done

    bcftools merge ~{sep=' ' input_vcfs} -O z -o ~{output_vcf_basename}.vcf.gz
    bcftools index -t ~{output_vcf_basename}.vcf.gz
  >>>
  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} SSD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    noAddress: true
  }
  output {
    File output_vcf = "~{output_vcf_basename}.vcf.gz"
    File output_vcf_index = "~{output_vcf_basename}.vcf.gz.tbi"
  }
}

task CountSamples {
  input {
    File vcf

    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 3000
    Int disk_size_gb = ceil(size(vcf, "GiB")) + 10
  }

  command <<<
    set -e -o pipefail

    bcftools query -l ~{vcf} | wc -l
  >>>
  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    maxRetries: 1
    noAddress: true
  }
  output {
    Int nSamples = read_int(stdout())
  }
}

task AggregateImputationQCMetrics {
  input {
    File infoFile
    Int nSamples
    String basename

    String rtidyverse_docker = "us.gcr.io/broad-dsde-methods/rocker/tidyverse:4.1.0"
    Int cpu = 1
    Int memory_mb = 2000
    Int disk_size_gb = ceil(size(infoFile, "GiB")) + 10
  }
  command <<<
    Rscript -<< "EOF"
    library(dplyr)
    library(readr)
    library(purrr)
    library(ggplot2)

    sites_info <- read_tsv("~{infoFile}")
    genotyped_sites_info <- sites_info %>% filter(Genotyped=="Genotyped") %>% mutate(EmpRsq=as.numeric(EmpRsq))

    nSites_genotyped <- genotyped_sites_info %>% nrow()
    nSites_genotyped_maf_gt_0.05 <- genotyped_sites_info %>% filter(MAF > 0.05) %>% nrow()
    nSites_genotyped_emp_r2_gt_0.6 <- genotyped_sites_info %>% filter(EmpRsq>0.6) %>% nrow()
    nSites_genotyped_maf_gt_0.05_emp_r2_gt_0.6 <- sites_info %>% filter(EmpRsq>0.6, MAF > 0.05) %>% nrow()

    aggregated_metrics <- tibble(total_sites_evaluated=nSites_genotyped, total_sites_evaluated_maf_gt_0.05=nSites_genotyped_maf_gt_0.05, total_sites_emp_r2_gt_0.6=nSites_genotyped_emp_r2_gt_0.6, total_sites_evaluated_maf_gt_0.05_r2_gt_0.6=nSites_genotyped_maf_gt_0.05_emp_r2_gt_0.6)

    write_tsv(aggregated_metrics, "~{basename}_aggregated_imputation_metrics.tsv")

    EOF
  >>>
  runtime {
    docker: rtidyverse_docker
    disks : "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
  output {
    File aggregated_metrics = "~{basename}_aggregated_imputation_metrics.tsv"
  }
}

task StoreChunksInfo {
  input {
    Array[String] chroms
    Array[Int] starts
    Array[Int] ends
    Array[Int] vars_in_array
    Array[Int] vars_in_panel
    Array[Boolean] valids
    String basename

    String rtidyverse_docker = "us.gcr.io/broad-dsde-methods/rocker/tidyverse:4.1.0"
    Int cpu = 1
    Int memory_mb = 2000
    Int disk_size_gb = 10
  }
  command <<<
    Rscript -<< "EOF"
    library(dplyr)
    library(readr)

    chunk_info <- tibble(chrom = c("~{sep='", "' chroms}"), start = c("~{sep='", "' starts}"), ends = c("~{sep='", "' ends}"), vars_in_array = c("~{sep='", "' vars_in_array}"), vars_in_panel = c("~{sep='", "' vars_in_panel}"), chunk_was_imputed = as.logical(c("~{sep='", "' valids}")))
    failed_chunks <- chunk_info %>% filter(!chunk_was_imputed) %>% select(-chunk_was_imputed)
    n_failed_chunks <- nrow(failed_chunks)
    write_tsv(chunk_info, "~{basename}_chunk_info.tsv")
    write_tsv(failed_chunks, "~{basename}_failed_chunks.tsv")
    write(n_failed_chunks, "n_failed_chunks.txt")
    EOF
  >>>
  runtime {
    docker: rtidyverse_docker
    disks : "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    maxRetries: 1
    noAddress: true
  }
  output {
    File chunks_info = "~{basename}_chunk_info.tsv"
    File failed_chunks = "~{basename}_failed_chunks.tsv"
    File n_failed_chunks = "n_failed_chunks.txt"
  }
}

task MergeImputationQCMetrics {
  input {
    Array[File] metrics
    String basename

    String rtidyverse_docker = "us.gcr.io/broad-dsde-methods/rocker/tidyverse:4.1.0"
    Int cpu = 1
    Int memory_mb = 2000
    Int disk_size_gb = ceil(size(metrics, "GiB")) + 10
  }
  command <<<
    Rscript -<< "EOF"
    library(dplyr)
    library(readr)
    library(purrr)
    library(ggplot2)

    metrics <- list("~{sep='", "' metrics}") %>% map(read_tsv) %>% reduce(`+`) %>% mutate(frac_sites_emp_r2_gt_0.6=total_sites_evaluated_maf_gt_0.05/total_sites_evaluated, frac_sites_maf_gt_0.05_emp_r2_gt_0.6=total_sites_evaluated_maf_gt_0.05_r2_gt_0.6/total_sites_evaluated_maf_gt_0.05)

    write_tsv(metrics, "~{basename}_aggregated_imputation_metrics.tsv")
    write(metrics %>% pull(frac_sites_maf_gt_0.05_emp_r2_gt_0.6), "frac_above_maf_5_percent_well_imputed.txt")

    EOF
  >>>
  runtime {
    docker: rtidyverse_docker
    disks : "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
  output {
    File aggregated_metrics = "~{basename}_aggregated_imputation_metrics.tsv"
    Float frac_above_maf_5_percent_well_imputed = read_float("frac_above_maf_5_percent_well_imputed.txt")
  }
}

task SubsetVcfToRegion {
  input {
    File vcf
    File vcf_index
    String output_basename
    String contig
    Int start
    Int end
    Boolean exclude_filtered = false

    Int disk_size_gb = ceil(2*size(vcf, "GiB")) + 10
    Int cpu = 1
    Int memory_mb = 8000
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
  }
  Int command_mem = memory_mb - 1500
  Int max_heap = memory_mb - 1000

  command {
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    SelectVariants \
    -V ~{vcf} \
    -L ~{contig}:~{start}-~{end} \
    -select 'POS >= ~{start}' ~{if exclude_filtered then "--exclude-filtered" else ""} \
    -O ~{output_basename}.vcf.gz
  }

  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    noAddress: true
  }

  parameter_meta {
    vcf: {
       description: "vcf",
       localization_optional: true
     }
    vcf_index: {
       description: "vcf index",
       localization_optional: true
     }
  }

  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task SetIDs {
  input {
    File vcf
    String output_basename

    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 4000
    Int disk_size_gb = ceil(2.2 * size(vcf, "GiB")) + 10
  }
  command <<<
    set -e -o pipefail

    bcftools annotate ~{vcf} --set-id '%CHROM\:%POS\:%REF\:%FIRST_ALT' -Oz -o ~{output_basename}.vcf.gz
    bcftools index -t ~{output_basename}.vcf.gz
  >>>
  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task ExtractIDs {
  input {
    File vcf
    String output_basename

    Int disk_size_gb = 2*ceil(size(vcf, "GiB")) + 10
    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 4000
  }
  command <<<
    bcftools query -f "%ID\n" ~{vcf} -o ~{output_basename}.ids.txt
  >>>
  output {
    File ids = "~{output_basename}.ids.txt"
  }
  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
}

task SelectVariantsByIds {
  input {
    File vcf
    File vcf_index
    File ids
    String basename

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int cpu = 1
    Int memory_mb = 16000
    Int disk_size_gb = ceil(1.2*size(vcf, "GiB")) + 10
  }
  parameter_meta {
    vcf: {
      description: "vcf",
      localization_optional: true
    }
    vcf_index: {
      description: "vcf",
      localization_optional: true
    }
  }
  Int command_mem = memory_mb - 2000
  Int max_heap = memory_mb - 1500

  command <<<
    set -e -o pipefail

    cp ~{ids} sites.list
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    SelectVariants -V ~{vcf} --exclude-filtered --keep-ids sites.list -O ~{basename}.vcf.gz
  >>>
  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size_gb} SSD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"
  }
}

task RemoveAnnotations {
  input {
    File vcf
    String basename

    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 3000
    Int disk_size_gb = ceil(2.2*size(vcf, "GiB")) + 10
  }
  command <<<
    set -e -o pipefail

    bcftools annotate ~{vcf} -x FORMAT,INFO -Oz -o ~{basename}.vcf.gz
    bcftools index -t ~{basename}.vcf.gz
  >>>
  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"
  }
}

task InterleaveVariants {
  input {
    Array[File] vcfs
    String basename

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int cpu = 1
    Int memory_mb = 16000
    Int disk_size_gb = ceil(3.2*size(vcfs, "GiB")) + 10
  }
  Int command_mem = memory_mb - 1500
  Int max_heap = memory_mb - 1000

  command <<<
    set -e -o pipefail

    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
    MergeVcfs -I ~{sep=" -I " vcfs} -O ~{basename}.vcf.gz
  >>>
  runtime {
    docker: gatk_docker
    disks: "local-disk ${disk_size_gb} SSD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
  output {
    File output_vcf = "~{basename}.vcf.gz"
    File output_vcf_index = "~{basename}.vcf.gz.tbi"
  }
}

task FindSitesUniqueToFileTwoOnly {
  input {
    File file1
    File file2

    String ubuntu_docker = "us.gcr.io/broad-dsde-methods/ubuntu:20.04"
    Int cpu = 1
    Int memory_mb = 4000
    Int disk_size_gb = ceil(size(file1, "GiB") + 2*size(file2, "GiB")) + 10
  }
  command <<<
    set -e -o pipefail

    comm -13 <(sort ~{file1} | uniq) <(sort ~{file2} | uniq) > missing_sites.ids
  >>>
  runtime {
    docker: ubuntu_docker
    disks: "local-disk ${disk_size_gb} HDD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    preemptible: 3
    noAddress: true
  }
  output {
    File missing_sites = "missing_sites.ids"
  }
}

task SplitMultiSampleVcf {
 input {
    File multiSampleVcf
    Int nSamples

    String bcftools_docker = "us.gcr.io/broad-gotc-prod/imputation-bcf-vcf:1.0.7-1.10.2-0.1.16-1669908889"
    Int cpu = 1
    Int memory_mb = 6000

    # This calculation is explained in https://github.com/broadinstitute/warp/pull/937
    Int disk_size_gb = ceil(21*nSamples*size(multiSampleVcf, "GiB")/(nSamples+20)) + 10
  }
  command <<<
    set -e -o pipefail

    mkdir out_dir
    bcftools +split ~{multiSampleVcf} -Oz -o out_dir
    for vcf in out_dir/*.vcf.gz; do
      bcftools index -t $vcf
    done
  >>>
  runtime {
    docker: bcftools_docker
    disks: "local-disk ${disk_size_gb} SSD"
    memory: "${memory_mb} MiB"
    cpu: cpu
    noAddress: true
  }
  output {
    Array[File] single_sample_vcfs = glob("out_dir/*.vcf.gz")
    Array[File] single_sample_vcf_indices = glob("out_dir/*.vcf.gz.tbi")
  }
}
