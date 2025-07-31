version 1.0

# BroadInternalUltimaGenomics is now deprecated 2025-03-06

import "../../../../../../../pipelines/wdl/dna_seq/germline/single_sample/ugwgs/UltimaGenomicsWholeGenomeGermline.wdl" as UltimaGenomicsWholeGenomeGermline
import "../../../../../../../structs/dna_seq/UltimaGenomicsWholeGenomeGermlineStructs.wdl" as Structs
import "../../../../../../../pipelines/wdl/qc/CheckFingerprint.wdl" as FP

workflow BroadInternalUltimaGenomics {

  String pipeline_version = "1.1.4"

  input {
  
  # INTERNAL BROAD - TDR AND FINGERPRINTING INPUTS
  String environment
  String sample_lsid
  String output_basename
  File haplotype_database_file = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt"
  File vault_token_path
  
  String? tdr_dataset_uuid
  
  # UltimaGenomics Pipeline Inputs
  ContaminationSites contamination_sites
  AlignmentReferences alignment_references
  VariantCallingSettings variant_calling_settings
  VcfPostProcessing vcf_post_processing

  # Sample Information
  Array[File]? input_cram_list
  Array[File]? input_bam_list
  String base_file_name

  Float rsq_threshold = 1.0
  Boolean merge_bam_file = true
  Boolean make_haplotype_bam = false
  Int reads_per_split = 20000000
  String filtering_model_no_gt_name = "rf_model_ignore_gt_incl_hpol_runs"

  String collab_sample_id_run_id
  }

  meta {
    allowNestedInputs: true
  }

  # PARAMETER DEFINITIONS
  parameter_meta {
    sample_lsid: "The sample lsid (an identifier used to retrieve fingerrints from Mercury)"
    input_cram_list: "Input crams. Will be realigned."
    base_file_name: "Basename for the output files."
    environment: "The environment (dev or prod) used for determining which service to use to retrieve Mercury fingerprints"
    vault_token_path: "The path to the vault token used for accessing the Mercury Fingerprint Store"
    output_basename: "String used as a prefix in workflow output files"
    tdr_dataset_uuid: "Optional String used to define the Terra Data Repo dataset to which outputs will be ingested, if populated"
  }

  call UltimaGenomicsWholeGenomeGermline.UltimaGenomicsWholeGenomeGermline {
    input:
      contamination_sites = contamination_sites,
      alignment_references = alignment_references,
      variant_calling_settings = variant_calling_settings,
      vcf_post_processing = vcf_post_processing,
      input_cram_list = input_cram_list,
      input_bam_list = input_bam_list,
      base_file_name = base_file_name,
      rsq_threshold = rsq_threshold,
      merge_bam_file = merge_bam_file,
      make_haplotype_bam = make_haplotype_bam,
      reads_per_split = reads_per_split,
      filtering_model_no_gt_name = filtering_model_no_gt_name
  }
  
  call FP.CheckFingerprint as CheckFingerprint {
    input:
      input_bam = UltimaGenomicsWholeGenomeGermline.output_cram,
      input_bam_index = UltimaGenomicsWholeGenomeGermline.output_cram_index,
      sample_alias = UltimaGenomicsWholeGenomeGermline.sample_name,
      sample_lsid = sample_lsid,
      output_basename = output_basename,
      ref_fasta = alignment_references.references.ref_fasta,
      ref_fasta_index = alignment_references.references.ref_fasta_index,
      ref_dict = alignment_references.references.ref_dict,
      read_fingerprint_from_mercury = true,
      haplotype_database_file = haplotype_database_file,
      environment = environment,
      vault_token_path = vault_token_path
  }

  call MergeMetrics {
    input:
      wgs_metrics = select_first([UltimaGenomicsWholeGenomeGermline.wgs_metrics]),
      alignment_summary_metrics = select_first([UltimaGenomicsWholeGenomeGermline.agg_alignment_summary_metrics]),
      gc_bias_summary_metrics = select_first([UltimaGenomicsWholeGenomeGermline.agg_gc_bias_summary_metrics]),
      duplicate_metrics = select_first([UltimaGenomicsWholeGenomeGermline.duplicate_metrics]),
      quality_yield_metrics = select_first([UltimaGenomicsWholeGenomeGermline.quality_yield_metrics]),
      raw_wgs_metrics = select_first([UltimaGenomicsWholeGenomeGermline.raw_wgs_metrics]),
      contamination_metrics = select_first([UltimaGenomicsWholeGenomeGermline.selfSM]),
      fingerprint_summary_metrics = CheckFingerprint.fingerprint_summary_metrics_file,
      output_basename = collab_sample_id_run_id
  }

  # CALL TDR TASKS TO FORMAT JSON
  
  if (defined(tdr_dataset_uuid)) {
    call formatPipelineOutputs {
      input:
        output_basename = output_basename,
        collab_sample_id_run_id = collab_sample_id_run_id,
        output_gvcf = UltimaGenomicsWholeGenomeGermline.output_gvcf,
        output_gvcf_index = UltimaGenomicsWholeGenomeGermline.output_gvcf_index,
        output_vcf = UltimaGenomicsWholeGenomeGermline.output_vcf,
        output_vcf_index = UltimaGenomicsWholeGenomeGermline.output_vcf_index,
        haplotype_bam = UltimaGenomicsWholeGenomeGermline.haplotype_bam,
        haplotype_bam_index = UltimaGenomicsWholeGenomeGermline.haplotype_bam_index,
        output_cram = UltimaGenomicsWholeGenomeGermline.output_cram,
        output_cram_index = UltimaGenomicsWholeGenomeGermline.output_cram_index,
        output_cram_md5 = UltimaGenomicsWholeGenomeGermline.output_cram_md5,
        selfSM = UltimaGenomicsWholeGenomeGermline.selfSM,
        contamination = UltimaGenomicsWholeGenomeGermline.contamination,
        filtered_vcf = UltimaGenomicsWholeGenomeGermline.filtered_vcf,
        filtered_vcf_index = UltimaGenomicsWholeGenomeGermline.filtered_vcf_index,
        quality_yield_metrics = UltimaGenomicsWholeGenomeGermline.quality_yield_metrics,
        wgs_metrics = UltimaGenomicsWholeGenomeGermline.wgs_metrics,
        raw_wgs_metrics = UltimaGenomicsWholeGenomeGermline.raw_wgs_metrics,
        duplicate_metrics = UltimaGenomicsWholeGenomeGermline.duplicate_metrics,
        agg_alignment_summary_metrics = UltimaGenomicsWholeGenomeGermline.agg_alignment_summary_metrics,
        agg_alignment_summary_pdf = UltimaGenomicsWholeGenomeGermline.agg_alignment_summary_pdf,
        agg_gc_bias_detail_metrics = UltimaGenomicsWholeGenomeGermline.agg_gc_bias_detail_metrics,
        agg_gc_bias_pdf = UltimaGenomicsWholeGenomeGermline.agg_gc_bias_pdf,
        agg_gc_bias_summary_metrics = UltimaGenomicsWholeGenomeGermline.agg_gc_bias_summary_metrics,
        agg_quality_distribution_pdf = UltimaGenomicsWholeGenomeGermline.agg_quality_distribution_pdf,
        agg_quality_distribution_metrics = UltimaGenomicsWholeGenomeGermline.agg_quality_distribution_metrics,
        duplication_rate = UltimaGenomicsWholeGenomeGermline.duplication_rate,
        chimerism_rate = UltimaGenomicsWholeGenomeGermline.chimerism_rate,
        is_outlier_data = UltimaGenomicsWholeGenomeGermline.is_outlier_data,
        sample_name = UltimaGenomicsWholeGenomeGermline.sample_name,
        flow_order = UltimaGenomicsWholeGenomeGermline.flow_order,
        barcode = UltimaGenomicsWholeGenomeGermline.barcode,
        id = UltimaGenomicsWholeGenomeGermline.id,
        unified_metrics = MergeMetrics.unified_metrics,
        fingerprint_summary_metrics_file = CheckFingerprint.fingerprint_summary_metrics_file,
        fingerprint_detail_metrics_file = CheckFingerprint.fingerprint_detail_metrics_file
    }

    call updateOutputsInTDR {
      input:
        tdr_dataset_uuid = select_first([tdr_dataset_uuid, ""]),
        outputs_json = formatPipelineOutputs.pipeline_outputs_json
    }
  }

  output {
    File? picard_fingerprint_summary_metrics = CheckFingerprint.fingerprint_summary_metrics_file
    File? picard_fingerprint_detail_metrics = CheckFingerprint.fingerprint_detail_metrics_file

    File output_gvcf = UltimaGenomicsWholeGenomeGermline.output_gvcf
    File output_gvcf_index = UltimaGenomicsWholeGenomeGermline.output_gvcf_index
    File output_vcf = UltimaGenomicsWholeGenomeGermline.output_vcf
    File output_vcf_index = UltimaGenomicsWholeGenomeGermline.output_vcf_index

    File? haplotype_bam = UltimaGenomicsWholeGenomeGermline.haplotype_bam
    File? haplotype_bam_index = UltimaGenomicsWholeGenomeGermline.haplotype_bam_index

    File output_cram = UltimaGenomicsWholeGenomeGermline.output_cram
    File output_cram_index = UltimaGenomicsWholeGenomeGermline.output_cram_index
    File output_cram_md5 = UltimaGenomicsWholeGenomeGermline.output_cram_md5

    File selfSM = UltimaGenomicsWholeGenomeGermline.selfSM
    Float contamination = UltimaGenomicsWholeGenomeGermline.contamination

    File filtered_vcf = UltimaGenomicsWholeGenomeGermline.filtered_vcf
    File filtered_vcf_index = UltimaGenomicsWholeGenomeGermline.filtered_vcf_index

    File quality_yield_metrics = UltimaGenomicsWholeGenomeGermline.quality_yield_metrics
    File wgs_metrics = UltimaGenomicsWholeGenomeGermline.wgs_metrics
    File raw_wgs_metrics = UltimaGenomicsWholeGenomeGermline.raw_wgs_metrics
    File duplicate_metrics = UltimaGenomicsWholeGenomeGermline.duplicate_metrics
    File agg_alignment_summary_metrics = UltimaGenomicsWholeGenomeGermline.agg_alignment_summary_metrics
    File? agg_alignment_summary_pdf = UltimaGenomicsWholeGenomeGermline.agg_alignment_summary_pdf
    File agg_gc_bias_detail_metrics = UltimaGenomicsWholeGenomeGermline.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = UltimaGenomicsWholeGenomeGermline.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = UltimaGenomicsWholeGenomeGermline.agg_gc_bias_summary_metrics
    File agg_quality_distribution_pdf = UltimaGenomicsWholeGenomeGermline.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = UltimaGenomicsWholeGenomeGermline.agg_quality_distribution_metrics
    Float duplication_rate = UltimaGenomicsWholeGenomeGermline.duplication_rate
    Float chimerism_rate = UltimaGenomicsWholeGenomeGermline.chimerism_rate
    Boolean is_outlier_data = UltimaGenomicsWholeGenomeGermline.is_outlier_data

    String sample_name = UltimaGenomicsWholeGenomeGermline.sample_name
    String flow_order = UltimaGenomicsWholeGenomeGermline.flow_order
    String barcode = UltimaGenomicsWholeGenomeGermline.barcode
    String id = UltimaGenomicsWholeGenomeGermline.id
  }
}

task MergeMetrics {
  input {
    File wgs_metrics
    File alignment_summary_metrics
    File gc_bias_summary_metrics
    File duplicate_metrics
    File quality_yield_metrics
    File raw_wgs_metrics
    File contamination_metrics
    File? fingerprint_summary_metrics
    String output_basename

    String docker =  "python:3.8-slim"
    Int cpu = 1
    Int memory_mb = 3000
    Int disk_size_gb = 10
  }

  String out_filename = output_basename + ".unified_metrics.txt"

  command <<<

    #
    # Script transpose a two line TSV
    #
    cat <<-'EOF' > transpose.py
    import csv, sys

    rows = list(csv.reader(sys.stdin, delimiter='\t'))

    for col in range(0, len(rows[0])):
      key = rows[0][col].lower()
      print(f"{key}\t{rows[1][col]}")
    EOF

    #
    # Script clean the keys, replacing space, dash and forward-slash with underscores,
    # and removing comma, single quote, periods, and #
    #
    cat <<-'EOF' > clean.py
    import sys

    for line in sys.stdin:
      (k,v) = line.strip().split("\t")
      transtable = k.maketrans({' ':'_', '-':'_', '/':'_', ',':None, '\'':None, '.':None, '#':None})
      print(f"{k.translate(transtable)}\t{v}")
    EOF

    # Process each metric file, transposing and cleaning if necessary, and pre-pending a source to the metric name

    echo "Processing WGS Metrics"
    cat ~{wgs_metrics} | grep -A 1 "GENOME_TERRITORY" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "wgs_metrics_" $0}' >> ~{out_filename}

    echo "Processing Alignment Summary Metrics"
    cat ~{alignment_summary_metrics} | grep -A 1 "CATEGORY" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "alignment_summary_metrics_" $0}' >> ~{out_filename}

    echo "Processing GCBias Summary Metrics"
    cat ~{gc_bias_summary_metrics} | grep -A 1 "AT_DROPOUT" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "gc_bias_summary_metrics_" $0}' >> ~{out_filename}

    echo "Processing Duplicate Metrics: WARNING this only works for samples from one library"
    cat ~{duplicate_metrics} | grep -A 1 "UNPAIRED_READ_DUPLICATES" | python transpose.py | awk '{print "duplicate_metrics_" $0}' >> ~{out_filename}

    echo "Processing Quality Yield Metrics"
    cat ~{quality_yield_metrics} | grep -A 1 "TOTAL_READS" | python transpose.py | awk '{print "quality_yield_metrics_" $0}' >> ~{out_filename}

    echo "Processing Raw WGS Metrics"
    cat ~{raw_wgs_metrics} | grep -A 1 "GENOME_TERRITORY" | python transpose.py | grep -Eiv "(SAMPLE|LIBRARY|READ_GROUP)" | awk '{print "raw_wgs_metrics_" $0}' >> ~{out_filename}

    echo "Processing Contamination Metrics"
    cat ~{contamination_metrics} | python transpose.py | python clean.py | awk '{print "contamination_metrics_" $0}' >> ~{out_filename}

    if [[ -f "~{fingerprint_summary_metrics}" ]];
    then
      echo "Processing Fingerprint Summary Metrics - only extracting LOD_EXPECTED_SAMPLE"
      cat ~{fingerprint_summary_metrics} | grep -A 1 "LOD_EXPECTED_SAMPLE" | python transpose.py | grep -i "LOD_EXPECTED_SAMPLE" | awk '{print "fp_"$0}' >> ~{out_filename}
    else
      echo "No Fingerprint Summary Metrics found."
      echo "fp_lod_expected_sample	" >> ~{out_filename}
    fi    >>>

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_mb} MiB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }

  output {
    File unified_metrics = out_filename
  }
}

# DEFINE TDR-SPECIFIC TASKS
  
  task formatPipelineOutputs {
    input {
      String output_basename
      String collab_sample_id_run_id
      
      String output_gvcf
      String output_gvcf_index
      String output_vcf
      String output_vcf_index

      String? haplotype_bam = ""
      String? haplotype_bam_index = ""

      String output_cram
      String output_cram_index
      String output_cram_md5

      String selfSM
      Float contamination

      String filtered_vcf
      String filtered_vcf_index

      String quality_yield_metrics
      String wgs_metrics
      String raw_wgs_metrics
      String duplicate_metrics
      String agg_alignment_summary_metrics
      String? agg_alignment_summary_pdf = ""
      String agg_gc_bias_detail_metrics
      String agg_gc_bias_pdf
      String agg_gc_bias_summary_metrics
      String agg_quality_distribution_pdf
      String agg_quality_distribution_metrics
      Float duplication_rate
      Float chimerism_rate
      Boolean is_outlier_data

      String sample_name
      String flow_order
      String barcode
      String id

      String? fingerprint_detail_metrics_file = ""
      String? fingerprint_summary_metrics_file = ""

      File unified_metrics

      Int cpu = 1
      Int memory_mb = 2000
      Int disk_size_gb = 10
    }

    String outputs_json_file_name = sub("outputs_to_TDR_~{output_basename}.json", " ", "")

    String outlier_data_string = if is_outlier_data then "True" else "False"

    command <<<
          python3 << CODE
          import json

          outputs_dict = {}

          outputs_dict["collab_sample_id_run_id"]="~{collab_sample_id_run_id}"
          outputs_dict["agg_alignment_summary_metrics"]="~{agg_alignment_summary_metrics}"
          outputs_dict["agg_alignment_summary_pdf"]="~{agg_alignment_summary_pdf}"
          outputs_dict["agg_gc_bias_detail_metrics"]="~{agg_gc_bias_detail_metrics}"
          outputs_dict["agg_gc_bias_pdf"]="~{agg_gc_bias_pdf}"
          outputs_dict["agg_gc_bias_summary_metrics"]="~{agg_gc_bias_summary_metrics}"
          outputs_dict["agg_quality_distribution_metrics"]="~{agg_quality_distribution_metrics}"
          outputs_dict["agg_quality_distribution_pdf"]="~{agg_quality_distribution_pdf}"
          outputs_dict["barcode"]="~{barcode}"
          outputs_dict["chimerism_rate"]="~{chimerism_rate}"
          outputs_dict["contamination"]="~{contamination}"
          outputs_dict["duplicate_metrics"]="~{duplicate_metrics}"
          outputs_dict["duplication_rate"]="~{duplication_rate}"
          outputs_dict["filtered_vcf"]="~{filtered_vcf}"
          outputs_dict["filtered_vcf_index"]="~{filtered_vcf_index}"
          outputs_dict["flow_order"]="~{flow_order}"
          outputs_dict["haplotype_bam"]="~{haplotype_bam}"
          outputs_dict["haplotype_bam_index"]="~{haplotype_bam_index}"
          outputs_dict["id"]="~{id}"
          outputs_dict["is_outlier_data"]="~{outlier_data_string}"
          outputs_dict["output_cram"]="~{output_cram}"
          outputs_dict["output_cram_index"]="~{output_cram_index}"
          outputs_dict["output_cram_md5"]="~{output_cram_md5}"
          outputs_dict["output_gvcf"]="~{output_gvcf}"
          outputs_dict["output_gvcf_index"]="~{output_gvcf_index}"
          outputs_dict["output_vcf"]="~{output_vcf}"
          outputs_dict["output_vcf_index"]="~{output_vcf_index}"
          outputs_dict["quality_yield_metrics"]="~{quality_yield_metrics}"
          outputs_dict["raw_wgs_metrics"]="~{raw_wgs_metrics}"
          outputs_dict["sample_name"]="~{sample_name}"
          outputs_dict["selfSM"]="~{selfSM}"
          outputs_dict["wgs_metrics"]="~{wgs_metrics}"
          outputs_dict["fingerprint_summary_metrics_file"]="~{fingerprint_summary_metrics_file}"
          outputs_dict["fingerprint_detail_metrics_file"]="~{fingerprint_detail_metrics_file}"

          # explode unified metrics file
          with open("~{unified_metrics}", "r") as infile:
            for row in infile:
              key, value = row.rstrip("\n").split("\t")
              if value == "NA" or value == "" or value == "?" or value == "-":
                outputs_dict[key] = None
              else:
                outputs_dict[key] = value

          # Write full outputs to file
          with open("~{outputs_json_file_name}", 'w') as outputs_file:
              for key, value in outputs_dict.items():
                  if value == "None" or value == "":
                      outputs_dict[key] = None
              outputs_file.write(json.dumps(outputs_dict))
              outputs_file.write("\n")
          CODE
      >>>

    runtime {
        docker: "broadinstitute/horsefish:tdr_import_v1.4"
        cpu: cpu
        memory: "${memory_mb} MiB"
        disks : "local-disk ${disk_size_gb} HDD"
    }

    output {
        File pipeline_outputs_json = outputs_json_file_name
    }
  }

  task updateOutputsInTDR {
    input {
      String tdr_dataset_uuid
      File outputs_json

      Int cpu = 1
      Int memory_mb = 2000
      Int disk_size_gb = 10
    }

    String tdr_target_table = "sample"

    command <<<
      python -u /scripts/export_pipeline_outputs_to_tdr.py \
        -d "~{tdr_dataset_uuid}" \
        -t "~{tdr_target_table}" \
        -o "~{outputs_json}"
    >>>

    runtime {
      docker: "broadinstitute/horsefish:tdr_import_v1.4"
      cpu: cpu
      memory: "~{memory_mb} MiB"
      disks: "local-disk ~{disk_size_gb} HDD"
    }

    output {
      File ingest_logs = stdout()
    }
  }