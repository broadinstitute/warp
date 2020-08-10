version 1.0

task Md5Sum {
  input {
    File input_cloud_path
  }

  String md5_output_name = basename(input_cloud_path) + ".md5sum"

  command <<<
    set -e
    set -o pipefail
    md5sum ~{input_cloud_path} | awk '{ print $1 }' > ~{md5_output_name}
  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 10 HDD"
    memory: "3.5 GiB"
    preemptible: 3
  }

  output {
    File md5_cloud_path = md5_output_name
  }
}

task zCall {
  input {
    String zcall_ped_filename
    String zcall_map_filename
    File input_gtc
    File bead_pool_manifest_csv_file
    File zcall_thresholds_file
    Int disk_size
    Int preemptible_tries
  }

  command <<<
    python /usr/gitc/zcall/zCall.py -B ~{bead_pool_manifest_csv_file} -G ~{input_gtc} -T ~{zcall_thresholds_file} > ~{zcall_ped_filename}
    python /usr/gitc/zcall/makeMAPfile.py -B ~{bead_pool_manifest_csv_file} > ~{zcall_map_filename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/zcall:4.0.1-1572616568"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File ped_file = zcall_ped_filename
    File map_file = zcall_map_filename
  }
}

task BpmToNormalizationManifestCsv {
  input {
    File bead_pool_manifest_file
    File cluster_file
    String bead_pool_manifest_csv_file
    Int disk_size
    Int preemptible_tries
  }
  command {
  java -Xms7g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
    BpmToNormalizationManifestCsv \
    --INPUT ~{bead_pool_manifest_file} \
    --CLUSTER_FILE ~{cluster_file} \
    --OUTPUT ~{bead_pool_manifest_csv_file}
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk 10 HDD"
    memory: "7.5 GiB"
    cpu: 2
    preemptible: preemptible_tries
  }

  output {
    File output_file = bead_pool_manifest_csv_file
  }
}

task GtcToVcf {
  input {
    String vcf_filename
    File input_gtc
    File? gender_gtc
    File extended_chip_manifest_file
    File cluster_file
    File bead_pool_manifest_file
    String sample_alias
    Int analysis_version_number
    String reported_gender
    String pipeline_version

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # Set memory to n1-standard-2 or n1-highmem-2
    Float memory = if size(extended_chip_manifest_file, "MiB") > 800 then 13 else 7.5
    Int disk_size
    Int preemptible_tries

    File? fingerprint_genotypes_vcf_file
    File? fingerprint_genotypes_vcf_index_file
  }

  command {
  java -Xms~{floor(memory - 1)}g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
    GtcToVcf \
    --INPUT ~{input_gtc} \
    ~{"--GENDER_GTC " + gender_gtc} \
    --OUTPUT ~{vcf_filename} \
    --MANIFEST ~{extended_chip_manifest_file} \
    --CLUSTER_FILE ~{cluster_file} \
    --BPM_FILE ~{bead_pool_manifest_file} \
    --DO_NOT_ALLOW_CALLS_ON_ZEROED_OUT_ASSAYS true \
    --SAMPLE_ALIAS "~{sample_alias}" \
    --ANALYSIS_VERSION_NUMBER ~{analysis_version_number} \
    --EXPECTED_GENDER "~{reported_gender}" \
    ~{"--FINGERPRINT_GENOTYPES_VCF_FILE " + fingerprint_genotypes_vcf_file} \
    --REFERENCE_SEQUENCE ~{ref_fasta} \
    --MAX_RECORDS_IN_RAM 100000 \
    --CREATE_INDEX true \
    --PIPELINE_VERSION ~{pipeline_version}
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory} GiB"
    cpu: 2
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = vcf_filename
    File output_vcf_index = "~{vcf_filename}.tbi"
  }
}

task VcfToAdpc {
  input {
    File input_vcf
    File? contamination_controls_vcf
    String samples_filename
    String num_markers_filename
    String output_adpc_filename

    Int disk_size
    Int preemptible_tries
  }

  command {
    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
             VcfToAdpc \
             --VCF ~{input_vcf} \
             ~{"--VCF " + contamination_controls_vcf} \
             --SAMPLES_FILE ~{samples_filename} \
             --NUM_MARKERS_FILE ~{num_markers_filename} \
             --OUTPUT ~{output_adpc_filename}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File output_adpc = output_adpc_filename
    File samples_file = samples_filename
    Int num_markers = read_int(num_markers_filename)
  }
}

task VerifyIDIntensity {
  input {
    File input_vcf
    File input_adpc_file
    Int num_samples
    Int num_markers
    String output_filename

    Int disk_size
    Int preemptible_tries
  }

  command {
    /usr/gitc/verifyIDintensity -m ~{num_markers} -n ~{num_samples} -i ~{input_adpc_file} -v -p > ~{output_filename}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/verify-id-intensity:e6354872834fe4262354a6b27bfe85ecc1323677-1561566044"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File output_file = output_filename
  }
}

task CreateVerifyIDIntensityContaminationMetricsFile {
  input {
    File input_file
    String expected_first_sample_name
    File samples_file
    String output_metrics_basefilename

    Int disk_size
    Int preemptible_tries
  }

  command {
    set -eo pipefail

    # Since VerifyIDIntensity might have been run in multi-sample mode and we only want the contamination
    # of the *first* sample, we create a truncated version of the input_file with ONLY THAT FIRST LINE
    TRUNCATED_INPUT_FILE=~{expected_first_sample_name}.one_record.txt
    head -3 ~{input_file} > $TRUNCATED_INPUT_FILE

    # We also verify that the first line in the samples_file is the chip_well_barcode under analysis
    # or else we have a problem.
    FOUND_FIRST_SAMPLE=$(head -1 ~{samples_file})
    if [[ $FOUND_FIRST_SAMPLE != ~{expected_first_sample_name} ]]
    then
      echo "First sample in ~{samples_file} does not match the chip_well_barcode of this analysis"
      exit 1;
    fi

    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
             CreateVerifyIDIntensityContaminationMetricsFile \
             --INPUT $TRUNCATED_INPUT_FILE \
             --OUTPUT ~{output_metrics_basefilename}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File output_metrics_file = "~{output_metrics_basefilename}.verifyidintensity_metrics"
  }
}

task CollectArraysVariantCallingMetrics {
  input {
    File input_vcf_file
    File input_vcf_index_file
    File dbSNP_vcf_file
    File dbSNP_vcf_index_file
    Float call_rate_threshold
    String output_metrics_basename

    Int disk_size
    Int preemptible_tries
  }

  command <<<
    set -eo pipefail

    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
      CollectArraysVariantCallingMetrics \
      --INPUT ~{input_vcf_file} \
      --DBSNP ~{dbSNP_vcf_file} \
      --CALL_RATE_PF_THRESHOLD ~{call_rate_threshold} \
      --OUTPUT ~{output_metrics_basename}

    # Need to determine the disposition from the metrics.
    # Remove all the comments and blank lines from the file
    # Find the column number of AUTOCALL_PF and retrieve the value in the second line of the AUTOCALL_PF column
    # AUTOCALL_PF set to empty if file has more than 2 lines (should only have column headers and one data line)
    AUTOCALL_PF=$(sed '/#/d' ~{output_metrics_basename}.arrays_variant_calling_detail_metrics |
      sed '/^\s*$/d' |
      awk -v col=AUTOCALL_PF -F '\t' \
      'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}}} NR==2{print $c} NR>2{exit 1}')

    if [[ "$AUTOCALL_PF" == "Y" ]]
    then
      echo true > pass.txt
    elif [[ "$AUTOCALL_PF" == "N" ]]
    then
      echo false > pass.txt
    else
      echo "AUTOCALL_PF should only be Y or N and there should only be one line of data"
      exit 1;
    fi
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File summary_metrics = "~{output_metrics_basename}.arrays_variant_calling_summary_metrics"
    File detail_metrics = "~{output_metrics_basename}.arrays_variant_calling_detail_metrics"
    File control_metrics = "~{output_metrics_basename}.arrays_control_code_summary_metrics"
    Boolean passes_autocall = read_boolean("pass.txt")
  }
}

task VcfToIntervalList {
  input {
    File vcf_file
    String interval_list_filename

    Int disk_size
    Int preemptible_tries
  }

  command {
    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
      VcfToIntervalList \
      --INPUT ~{vcf_file} \
      --OUTPUT ~{interval_list_filename}
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File interval_list_file = interval_list_filename
  }
}

task CheckFingerprint {
  input {
    File input_vcf_file
    File input_vcf_index_file
    File? genotypes_vcf_file
    File? genotypes_vcf_index_file
    File haplotype_database_file
    String observed_sample_alias
    String expected_sample_alias
    String output_metrics_basename

    Int disk_size
    Int preemptible_tries
  }

  # Paraphrased from Yossi:
  # Override the default LOD threshold of 5 because if the PL field
  # is missing from the VCF, CheckFingerprint will default to an error
  # rate equivalent to a LOD score of 2, and we don't want to see
  # confident LOD scores w/ no confident SNPs.
  Float genotype_lod_threshold = 1.9
  String summary_metrics_extension = ".fingerprinting_summary_metrics"

  command <<<
    set -e
    java -Xms2g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
      CheckFingerprint \
      --INPUT ~{input_vcf_file} \
      --OBSERVED_SAMPLE_ALIAS "~{observed_sample_alias}" \
      ~{"--GENOTYPES " + genotypes_vcf_file} \
      --EXPECTED_SAMPLE_ALIAS "~{expected_sample_alias}" \
      --HAPLOTYPE_MAP ~{haplotype_database_file} \
      --GENOTYPE_LOD_THRESHOLD ~{genotype_lod_threshold} \
      --OUTPUT ~{output_metrics_basename}

    CONTENT_LINE=$(cat ~{output_metrics_basename}~{summary_metrics_extension} |
    grep -n "## METRICS CLASS\tpicard.analysis.FingerprintingSummaryMetrics" |
    cut -f1 -d:)
    CONTENT_LINE=$(($CONTENT_LINE+2))
    sed '8q;d' ~{output_metrics_basename}~{summary_metrics_extension} | cut -f5 > lod
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File summary_metrics = "~{output_metrics_basename}~{summary_metrics_extension}"
    File detail_metrics = "~{output_metrics_basename}.fingerprinting_detail_metrics"
    Float lod = read_float("lod")
  }
}

task SelectVariants {
  input {
    File input_vcf_file
    File input_vcf_index_file
    String output_vcf_filename

    Boolean excludeFiltered = false
    Boolean excludeNonVariants = false
    File? variant_rsids_file

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Int disk_size
    Int preemptible_tries
  }

  String base_vcf = basename(output_vcf_filename)
  Boolean is_compressed = basename(base_vcf, "gz") != base_vcf
  String vcf_index_suffix = if is_compressed then ".tbi" else ".idx"

  command <<<
    set -eo pipefail

    /gatk/gatk --java-options -Xms2g \
      SelectVariants \
      -V ~{input_vcf_file} \
      ~{true="--exclude-filtered true" false="" excludeFiltered} \
      ~{true="--exclude-non-variants true" false="" excludeNonVariants} \
      ~{"--keep-ids " + variant_rsids_file} \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.3.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = output_vcf_filename
    File output_vcf_index = "~{output_vcf_filename}~{vcf_index_suffix}"
  }
}

task SelectIndels {
  input {
    File input_vcf_file
    File input_vcf_index_file
    String output_vcf_filename

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Int disk_size
    Int preemptible_tries
  }

  command <<<
    /gatk/gatk --java-options -Xms2g \
      SelectVariants \
      -V ~{input_vcf_file} \
      --select-type-to-include INDEL \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.3.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = output_vcf_filename
    File output_vcf_index = "~{output_vcf_filename}.idx"
  }
}


task AutoCall {
  input {
    String chip_well_barcode
    File green_idat_cloud_path
    File red_idat_cloud_path
    File bead_pool_manifest_file
    File? cluster_file
    Boolean? is_gender_autocall
    Int disk_size
    Int preemptible_tries
  }

  String gtc_filename = "~{chip_well_barcode}.gtc"

  command <<<
    set -e
    rm -rf ~{chip_well_barcode}
    mkdir ~{chip_well_barcode}
    mv ~{green_idat_cloud_path} ~{chip_well_barcode}
    mv ~{red_idat_cloud_path} ~{chip_well_barcode}
    # Make an empty gtc file so that if autocall fails PAPI won't fail on this task (for inability to find the output) and the WF can move on.
    touch ~{gtc_filename}

    /usr/gitc/iaap/iaap-cli/iaap-cli \
      gencall \
      ~{bead_pool_manifest_file} \
      ~{cluster_file} \
      . \
      -f ~{chip_well_barcode} \
      ~{default='' true='--gender-estimate-call-rate-threshold -0.1' false='' is_gender_autocall} \
      -g
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.1-1572616845"
    disks: "local-disk " + disk_size + " HDD"
    memory: "7 GiB"
    preemptible: preemptible_tries
  }

  output {
    File gtc_file = gtc_filename
  }
}

task MergePedIntoVcf {
  input {
    File input_vcf
    File input_vcf_index
    File ped_file
    File map_file
    File zcall_thresholds_file
    String zcall_version

    String output_vcf_filename

    Int disk_size
    Int preemptible_tries
  }

  command {
    java -Xms3g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
      MergePedIntoVcf \
      --ORIGINAL_VCF ~{input_vcf} \
      --PED_FILE ~{ped_file} \
      --MAP_FILE ~{map_file} \
      --ZCALL_THRESHOLDS_FILE ~{zcall_thresholds_file} \
      --ZCALL_VERSION ~{zcall_version} \
      --OUTPUT ~{output_vcf_filename} \
      --CREATE_INDEX true
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    memory: "3.5 GiB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }

  output {
    File output_vcf = output_vcf_filename
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

task SubsetArrayVCF {
  input {
    File intervals
    File input_vcf_file
    File input_vcf_index_file
    File ref_fasta
    File ref_fasta_index
    File ref_dict
  }

  String output_name = basename(input_vcf_file, ".vcf") + "_subset.vcf"
  Int disk_size = ceil(size(input_vcf_file, "GiB") * 2 + size(ref_fasta, "GiB"))

  command <<<
    gatk SelectVariants -V  ~{input_vcf_file} -L ~{intervals} -O ~{output_name} -R ~{ref_fasta}
  >>>

  output {
    File output_vcf = output_name
    File output_vcf_index = output_name + ".idx"
  }

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.3.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
  }

}

task GenotypeConcordance {
  input {
    File call_vcf_file
    File call_vcf_index_file
    File call_intervals_file
    String call_sample_name
    File? truth_vcf_file
    File? truth_vcf_index_file
    File? truth_intervals_file
    String? truth_sample_name
    String output_metrics_basename
    Float genotype_concordance_threshold=0.99

    Int disk_size
    Int preemptible_tries
  }

  command <<<
    set -eo pipefail

    java -Xms3g -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
      GenotypeConcordance \
      --CALL_VCF ~{call_vcf_file} \
      --CALL_SAMPLE ~{call_sample_name} \
      --TRUTH_VCF ~{truth_vcf_file} \
      --TRUTH_SAMPLE ~{truth_sample_name} \
      --INTERVALS ~{call_intervals_file} \
      --INTERVALS ~{truth_intervals_file} \
      --OUTPUT ~{output_metrics_basename}

    # Strip out comments and blank lines and get the genotype concordance
    # Find the column number of GENOTYPE_CONCORDANCE and retrieve the value in the second line of the column
    # CONCORDANCE set to empty if file has more than 3 lines (should only have column headers and two data lines)
    # Samples that have few indel cause genotype concordance to be reported as ‘?’ in the metrics file
    CONCORDANCE=$(sed '/#/d' ~{output_metrics_basename}.genotype_concordance_summary_metrics |
      sed '/^\s*$/d' |
      awk -v col=GENOTYPE_CONCORDANCE -F '\t' \
      'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}}} NR==2{print $c} NR>3{exit 1}')

    if [[ "$CONCORDANCE" =~ ^[+-]?[0-9]+\.?[0-9]*$ || "$CONCORDANCE" =~ "?" ]]
    then
        if awk 'BEGIN{exit ARGV[1]>=ARGV[2]}' "$CONCORDANCE" "~{genotype_concordance_threshold}"
        then
          # Less than threshold. Need to blacklist
          echo true > fails_concordance.txt
        else
          # Passes Genotype Concordance. No need to blacklist
          echo false > fails_concordance.txt
        fi
    else
        echo "GENOTYPE_CONCORDANCE must be a number and the file should only have 2 lines of data"
        exit 1;
    fi
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }

  output {
    File summary_metrics = "~{output_metrics_basename}.genotype_concordance_summary_metrics"
    File detail_metrics = "~{output_metrics_basename}.genotype_concordance_detail_metrics"
    File contingency_metrics = "~{output_metrics_basename}.genotype_concordance_contingency_metrics"
    Boolean fails_concordance = read_boolean("fails_concordance.txt")
  }
}

task ValidateVariants {
  input {
    File input_vcf_file
    File input_vcf_index_file

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    Int disk_size
    Int preemptible_tries
  }
  command <<<
    /gatk/gatk --java-options -Xms2g \
      ValidateVariants \
      -V ~{input_vcf_file} \
      --validation-type-to-exclude ALLELES \
      -R ~{ref_fasta}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.1.3.0"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
  }
}
