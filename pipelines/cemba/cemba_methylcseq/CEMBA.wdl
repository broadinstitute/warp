version 1.0

# see https://docs.google.com/document/d/1mO0nwvrUtHFx_TevX7jQqOZ-SaBkGu1bFCvW4bO3rKM/ for more info on pipeline
workflow CEMBA {
    input {
      # name of outputs and intermediate files of pipeline
      String output_base_sample_name

      # compressed read 1 and read 2 paired inputs
      File fastq_r1_gzipped_input
      File fastq_r2_gzipped_input

      File? monitoring_script

      # a list of cell barcodes that are being used
      File? barcode_white_list

      # location of barcode
      Int barcode_start_pos
      Int barcode_length

      # bowtie2 indexes used for mapping
      # using BuildCembaReference.wdl
      # reference dictionary for building VCF and attaching barcodes
      File reference_fasta
      File reference_fasta_index
      File reference_dictionary
      File fwd_converted_reference_fasta
      File rev_converted_reference_fasta
      Array[File] fwd_bowtie2_index_files
      Array[File] rev_bowtie2_index_files

      # trimming cutadapt options
      Int quality_cutoff
      Int min_length_paired_end_trim
      Int min_length_single_end_trim
      String read1_adapter_seq
      String read2_adapter_seq
      Int cut_length

      # paired end vs single end mapping boolean/option
      Boolean paired_end_run

      # should mark duplicates in picard mark or removes duplicates
      Boolean remove_duplicates

      # default add barcodes to read 1 in outputted bam (compitable with single end alignment only)
      Boolean extract_and_attach_barcodes_in_single_end_run

      # filter map quality option (optional)
      Int? min_map_quality

      # default names for adding a read group
      String read_group_library_name = "Methylation"
      String read_group_platform_name = "Illumina"
      String read_group_platform_unit_name = "snmC-Seq"
    }

    # version of this pipeline
    String pipeline_version = "1.1.5"

  # trim off hardcoded sequence adapters
  call Trim as TrimAdapters {
    input:
      fastq_input_read_a = fastq_r1_gzipped_input,
      fastq_input_read_b = fastq_r2_gzipped_input,
      quality_cutoff = quality_cutoff,
      min_length = min_length_paired_end_trim,
      read1_adapter_seq = read1_adapter_seq,
      read2_adapter_seq = read2_adapter_seq,
      output_base_name = output_base_sample_name,
      monitoring_script = monitoring_script
  }

  if (extract_and_attach_barcodes_in_single_end_run && !paired_end_run) {
    # produce an unmapped bam to tag the barcodes
    call CreateUnmappedBam {
      input:
        fastq_input = TrimAdapters.trimmed_fastqs[0],
        output_base_name = output_base_sample_name + ".R1",
        monitoring_script = monitoring_script
    }

    # extract and tag the barcodes to unmapped bam
    call ExtractCellBarcodes {
      input:
        fastq_input = TrimAdapters.trimmed_fastqs[0],
        unmapped_bam_input = CreateUnmappedBam.unmapped_bam_output,
        barcode_white_list = barcode_white_list,
        barcode_start_pos = barcode_start_pos,
        barcode_length = barcode_length,
        output_base_name = output_base_sample_name + ".R1",
        monitoring_script = monitoring_script
    }
  }

  # trim off Degenerate bases H = [A, T or C]/primer index sequence of Read 1
  call Trim as TrimSingleRead1 {
    input:
      fastq_input_read_a = TrimAdapters.trimmed_fastqs[0],
      min_length = min_length_single_end_trim,
      cut_length = cut_length,
      output_base_name = output_base_sample_name + ".R1",
      monitoring_script = monitoring_script
  }

  # trim off the C/T tail appended by Adaptase of read2
  call Trim as TrimSingleRead2 {
    input:
      fastq_input_read_a = TrimAdapters.trimmed_fastqs[1],
      min_length = min_length_single_end_trim,
      cut_length = cut_length,
      output_base_name = output_base_sample_name + ".R2",
      monitoring_script = monitoring_script
  }

  if (paired_end_run) {
    # map as paired end
    call Align as MapReadsPairedEnd {
      input:
        paired_end_run = true,
        directional = true,
        fastq_input_read_a = TrimSingleRead1.trimmed_fastqs[0],
        fastq_input_read_b = TrimSingleRead2.trimmed_fastqs[0],
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        fwd_converted_reference_fasta = fwd_converted_reference_fasta,
        fwd_bowtie2_index_files = fwd_bowtie2_index_files,
        rev_converted_reference_fasta = rev_converted_reference_fasta,
        rev_bowtie2_index_files = rev_bowtie2_index_files,
        output_base_name = output_base_sample_name + ".paired_end",
        monitoring_script = monitoring_script
    }
  }
  if (!paired_end_run) {
    # map read 1 as single-end
    call Align as MapReadSingleEndRead1  {
      input:
        paired_end_run = false,
        directional = true,
        fastq_input_read_a = TrimSingleRead1.trimmed_fastqs[0],
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        fwd_converted_reference_fasta = fwd_converted_reference_fasta,
        fwd_bowtie2_index_files = fwd_bowtie2_index_files,
        rev_converted_reference_fasta = rev_converted_reference_fasta,
        rev_bowtie2_index_files = rev_bowtie2_index_files,
        output_base_name = output_base_sample_name + ".R1.single_end",
        monitoring_script = monitoring_script
    }

    # map read 2 as single-end
    call Align as MapReadSingleEndRead2  {
      input:
        paired_end_run = false,
        directional = false,
        fastq_input_read_a = TrimSingleRead2.trimmed_fastqs[0],
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        fwd_converted_reference_fasta = fwd_converted_reference_fasta,
        fwd_bowtie2_index_files = fwd_bowtie2_index_files,
        rev_converted_reference_fasta = rev_converted_reference_fasta,
        rev_bowtie2_index_files = rev_bowtie2_index_files,
        output_base_name = output_base_sample_name + ".R2.single_end",
        monitoring_script = monitoring_script
    }
  }

  # either the paired end aligned bam or...
  # the single end aligned read 1 bam and the single end aligned read 2 bam
  Array[Pair[File, String]] alignment_outputs = (
    if paired_end_run
    then [(select_first([MapReadsPairedEnd.mapped_bam_output]), "")]
    else [
      (select_first([MapReadSingleEndRead1.mapped_bam_output]), ".R1"),
      (select_first([MapReadSingleEndRead2.mapped_bam_output]), ".R2")
    ]
  )

  scatter (alignment_output in alignment_outputs) {
    # sort the bam in coordinate order to filter
    call Sort as SortAlignmentOutputBam {
      input:
        bam_input = alignment_output.left,
        output_base_name = output_base_sample_name + alignment_output.right,
        monitoring_script = monitoring_script
    }

    # remove duplicates from bam
    call FilterDuplicates {
      input:
        bam_input = SortAlignmentOutputBam.bam_sort_output,
        remove_duplicates = remove_duplicates,
        output_base_name = output_base_sample_name + alignment_output.right,
        monitoring_script = monitoring_script
    }

    # get a methylation report for the filtered duplicates bam
    call GetMethylationReport as GetMethylationReportForFilterDuplicates {
      input:
        bam_input = FilterDuplicates.bam_remove_dup_output,
        output_base_name = output_base_sample_name + alignment_output.right,
        monitoring_script = monitoring_script
    }

    if (defined(min_map_quality)) {
      # filter bam by map quality
      call FilterMapQuality {
      input:
        bam_input = FilterDuplicates.bam_remove_dup_output,
        min_map_quality = select_first([min_map_quality]),
        output_base_name = output_base_sample_name + alignment_output.right,
        monitoring_script = monitoring_script
      }

      # get methylation report for filtered bam
      call GetMethylationReport as GetMethylationReportForAboveMinMapQReads {
        input:
          bam_input = FilterMapQuality.bam_filter_above_min_mapq_output,
          output_base_name = output_base_sample_name + alignment_output.right + ".above_min_map_quality",
          monitoring_script = monitoring_script
      }

      # get methylation report for filtered bam
      call GetMethylationReport as GetMethylationReportForBelowMinMapQReads {
        input:
          bam_input = FilterMapQuality.bam_filter_below_min_mapq_output,
          output_base_name = output_base_sample_name + alignment_output.right + ".below_min_map_quality",
          monitoring_script = monitoring_script
      }
    }

    # if not filtering by map quality: the filtered duplicatets bam
    # else: the filtered map quality bam
    File filtered_bam = select_first([FilterMapQuality.bam_filter_above_min_mapq_output, FilterDuplicates.bam_remove_dup_output])
  }

  # if mapped in single end (2 alignment outputs)
  if (!paired_end_run) {
    File aligned_and_filtered_read1_bam = filtered_bam[0]
    File aligned_and_filtered_read2_bam = filtered_bam[1]

    if (extract_and_attach_barcodes_in_single_end_run) {
      # add barcodes from tagged unmapped bam to aligned bam
      call AttachBarcodes {
        input:
          mapped_bam_input = aligned_and_filtered_read1_bam,
          tagged_unmapped_bam_input = select_first([ExtractCellBarcodes.tagged_unmapped_bam_output]),
          reference_fasta = reference_fasta,
          reference_dictionary = reference_dictionary,
          cut_length = cut_length,
          output_base_name = output_base_sample_name + ".R1",
          monitoring_script = monitoring_script
      }
    }

    # merge read 1 and read 2 alligned, sorted and filtered bams
    call MergeBams {
      input:
        bam_input_a = select_first([AttachBarcodes.tagged_mapped_bam_output, aligned_and_filtered_read1_bam]),
        bam_input_b = aligned_and_filtered_read2_bam,
        output_base_name = output_base_sample_name,
        monitoring_script = monitoring_script
    }
  }

  # input bam is either the flitered and merged read 1 SE and read 2 SE bam or the filtered PE bam
  call AddReadGroup {
    input:
      bam_input = select_first([MergeBams.merged_bam_output, filtered_bam[0]]),
      read_group_library_name = read_group_library_name,
      read_group_platform_name = read_group_platform_name,
      read_group_platform_unit_name = read_group_platform_unit_name,
      read_group_platform_sample_name = output_base_sample_name,
      output_base_name = output_base_sample_name,
      monitoring_script = monitoring_script
  }

  # sort again in coordinate order after adding read group
  call Sort as SortFilteredBamWithReadGroup {
    input:
      bam_input = AddReadGroup.bam_with_read_group_output,
      output_base_name = output_base_sample_name + ".aligned.filtered",
      monitoring_script = monitoring_script
  }

  # index the outputted bams
  call IndexBam {
    input:
      bam_input = SortFilteredBamWithReadGroup.bam_sort_output,
      output_base_name = output_base_sample_name + ".aligned.filtered.sorted",
      monitoring_script = monitoring_script
  }

  # get methylated VCF
  call MethylationTypeCaller as GetMethylationSiteVCF {
    input:
      bam_input = SortFilteredBamWithReadGroup.bam_sort_output,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      reference_dictionary = reference_dictionary,
      output_base_name = output_base_sample_name,
      monitoring_script = monitoring_script
  }

  # convert VCF to ALL
  call VCFtoALLC {
    input:
      methylation_vcf_output_name = GetMethylationSiteVCF.methylation_vcf
  }

  # get number of sites that have a coverage greater than 1
  call ComputeCoverageDepth {
    input:
      bam = SortFilteredBamWithReadGroup.bam_sort_output,
      reference_fasta = reference_fasta,
      monitoring_script = monitoring_script
  }

  # output the bam, metrics and reports
  # select all will select an array consisting of 2 single end reports or an array of 1 paired end reports
  output {
    File aligned_and_filtered_bam = SortFilteredBamWithReadGroup.bam_sort_output
    File aligned_and_filtered_bam_index = IndexBam.bam_index_output

    File methylation_site_vcf = GetMethylationSiteVCF.methylation_vcf
    File methylation_site_vcf_index = GetMethylationSiteVCF.methylation_vcf_index

    File methylation_site_allc = VCFtoALLC.methylation_allc

    Int coverage_depth = ComputeCoverageDepth.total_depth_count

    ### TODO
    #
    # We know all the outputs below here either have 1 output or 2, depending on
    # whether or not the pipeline was run in paired-end mode.
    #
    # Ideally we could have 3 `File?` outputs with meaningfully distinct names for
    # each output instead of an array. Limitations in draft-2 WDL's handling of
    # optional types prevents us from doing this within a single workflow. We could
    # work around this by pulling pieces of the pipeline into a sub-workflow which
    # we explicitly call three times, but that'll only work in Terra if we can
    # make the sub-workflow public.

    Array[File] mapping_reports = select_all([
      MapReadSingleEndRead1.mapping_report_output,
      MapReadSingleEndRead2.mapping_report_output,
      MapReadsPairedEnd.mapping_report_output
    ])

    Array[File] duplicate_metrics = FilterDuplicates.metric_remove_dup_output

    Array[File] methylation_mbias_reports_filtered_duplicates = GetMethylationReportForFilterDuplicates.methylation_mbias_report_output
    Array[File] methylation_splitting_reports_filtered_duplicates = GetMethylationReportForFilterDuplicates.methylation_splitting_report_output
    Array[File] methylation_CpG_context_reports_filtered_duplicates = GetMethylationReportForFilterDuplicates.methylation_CpG_context_report_output
    Array[File] methylation_non_CpG_context_reports_filtered_duplicates = GetMethylationReportForFilterDuplicates.methylation_non_CpG_context_report_output

    Array[File] methylation_mbias_reports_filtered_above_min_mapq = select_all(
      GetMethylationReportForAboveMinMapQReads.methylation_mbias_report_output
    )
    Array[File] methylation_splitting_reports_filtered_above_min_mapq = select_all(
      GetMethylationReportForAboveMinMapQReads.methylation_splitting_report_output
    )
    Array[File] methylation_CpG_context_reports_filtered_above_min_mapq = select_all(
      GetMethylationReportForAboveMinMapQReads.methylation_CpG_context_report_output
    )
    Array[File] methylation_non_CpG_context_reports_filtered_above_min_mapq = select_all(
      GetMethylationReportForAboveMinMapQReads.methylation_non_CpG_context_report_output
    )

    Array[File] methylation_mbias_reports_min_filtered_below_mapq = select_all(
      GetMethylationReportForBelowMinMapQReads.methylation_mbias_report_output
    )
    Array[File] methylation_splitting_reports_filtered_below_min_mapq = select_all(
      GetMethylationReportForBelowMinMapQReads.methylation_splitting_report_output
    )
    Array[File] methylation_CpG_context_reports_filtered_below_min_mapq = select_all(
      GetMethylationReportForBelowMinMapQReads.methylation_CpG_context_report_output
    )
    Array[File] methylation_non_CpG_context_reports_filtered_below_min_mapq = select_all(
      GetMethylationReportForBelowMinMapQReads.methylation_non_CpG_context_report_output
    )
  }
}

task Trim {
    input {
      Int min_length
      # read 1 when trimming adapters
      File fastq_input_read_a
      # read 2 when trimming adapters
      File? fastq_input_read_b
      String output_base_name
      String? read1_adapter_seq
      String? read2_adapter_seq
      Int? quality_cutoff
      Int? cut_length
      File? monitoring_script
    }

      # input file size
      Float input_size = size(fastq_input_read_a, "GB") + (if defined(fastq_input_read_b) then size(fastq_input_read_b, "GB") else 0)

      # output names for trimmed reads
      String fastq_trimmed_adapter_output_name_a = output_base_name
                                                     + (if defined(fastq_input_read_b) then ".R1.trimmed_adapters.fastq.gz" else ".trimmed_single.fastq.gz")
      String fastq_trimmed_adapter_output_name_b = (if defined(fastq_input_read_b) then output_base_name + ".R2.trimmed_adapters.fastq.gz" else "")

  # using cutadapt to trim off sequence adapters in paired end mode and c/t adaptase and cell/well barcode in sinlge end mode
  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # fastq's, "-f", -A for paired adapters read 2"
    cutadapt \
      -f fastq \
      --minimum-length ~{min_length} \
      --output ~{fastq_trimmed_adapter_output_name_a} \
      ~{true="--paired-output"  false="" defined(fastq_input_read_b)} ~{fastq_trimmed_adapter_output_name_b} \
      ~{true="--quality-cutoff" false="" defined(quality_cutoff)} ~{quality_cutoff} \
      ~{true="--adapter" false="" defined(read1_adapter_seq)} ~{read1_adapter_seq} \
      ~{true="-A" false="" defined(read2_adapter_seq)} ~{read2_adapter_seq} \
      ~{true="--cut" false="" defined(cut_length)} ~{cut_length} \
      ~{true="--cut -" false="" defined(cut_length)}~{cut_length} \
      ~{fastq_input_read_a} ~{fastq_input_read_b}
  >>>

  # use docker image for given tool cutadapat
  runtime {
    docker: "quay.io/broadinstitute/cutadapt:1.18"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2 * input file size
    disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    Array[File] trimmed_fastqs = (if defined(fastq_input_read_b) then [fastq_trimmed_adapter_output_name_a, fastq_trimmed_adapter_output_name_b] else [fastq_trimmed_adapter_output_name_a])
    File monitoring_log = "monitoring.log"
  }
}

task CreateUnmappedBam {
    input {
      # read 1 when with adapter sequences trimmed off
      File fastq_input
      String output_base_name
      File? monitoring_script
    }

      # input file size
      Float input_size = size(fastq_input, "GB")

      # output name for unaligned bam
      String unmapped_bam_output_name = output_base_name + ".unmapped.bam"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # create an unmapped bam
    java -Xmx3000m -jar /usr/picard/picard.jar FastqToSam \
      FASTQ=~{fastq_input} \
      SAMPLE_NAME=~{output_base_name} \
      OUTPUT=~{unmapped_bam_output_name}
  >>>

  # use docker image for given tool cutadapat
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2.25 * input file size
    disks: "local-disk " + ceil(2.25 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3500 MiB"
  }

  output {
    File unmapped_bam_output = unmapped_bam_output_name
    File monitoring_log = "monitoring.log"
  }
}

task ExtractCellBarcodes {
    input {
      # read 1 when with adapter sequences trimmed off
      File fastq_input
      File unmapped_bam_input
      File? barcode_white_list
      Int barcode_start_pos
      Int barcode_length
      String output_base_name
      File? monitoring_script
    }

      # input file size
      Float input_size = size(fastq_input, "GB") + (if defined(barcode_white_list) then size(barcode_white_list, "GB") else 0)

      # output names for read1 bam with extracted barcodes
      String tagged_unmapped_bam_output_name = output_base_name + ".tagged_unmapped.bam"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # extract barcode and tag them to the unmapped bam
    AttachBarcodes \
      --r1 ~{fastq_input} \
      --u2 ~{unmapped_bam_input} \
      --cell-barcode-start-position ~{barcode_start_pos} \
      --cell-barcode-length ~{barcode_length} \
      ~{true="--whitelist" false="" defined(barcode_white_list)} ~{barcode_white_list} \
      --output-bamfile ~{tagged_unmapped_bam_output_name}
  >>>

  # use docker image for given tool sctools
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.4"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2.5 * input file size
    disks: "local-disk " + ceil(2.5 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    File tagged_unmapped_bam_output = tagged_unmapped_bam_output_name
    File monitoring_log = "monitoring.log"
  }
}

# map using bowtie2
task Align {
    input {
      Boolean paired_end_run
      Boolean directional
      File fastq_input_read_a
      File? fastq_input_read_b
      File reference_fasta
      File reference_fasta_index
      File fwd_converted_reference_fasta
      Array[File] fwd_bowtie2_index_files
      File rev_converted_reference_fasta
      Array[File] rev_bowtie2_index_files
      String output_base_name
      File? monitoring_script
    }

      # output name for mapped reads
      # output filenames are generated by Bismark
      # if/else used for paired end vs sinlge end alignment
      String bismark_mapped_bam_output_name = basename(fastq_input_read_a, ".fastq.gz") +
                                        (if paired_end_run then "_bismark_bt2_pe.bam" else "_bismark_bt2.bam")
      String bismark_mapping_report_output_name = basename(fastq_input_read_a, ".fastq.gz") +
                                            (if paired_end_run then "_bismark_bt2_PE_report.txt" else "_bismark_bt2_SE_report.txt")
      String mapped_bam_output_name = output_base_name + ".aligned.bam"
      String mapping_report_output_name = output_base_name + ".alignment_report.txt"

      # input file sizes
      Float reference_size = size(reference_fasta, "GB") +
                                size(fwd_converted_reference_fasta, "GB") +
                                size(rev_converted_reference_fasta, "GB")
      Float bam_size = size(fastq_input_read_a, "GB") +
                               (if paired_end_run then size(fastq_input_read_b, "GB") else 0)

  # align with bowtie2
  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # create dir struture required by bismark
    # move input files into corresponding sub dir
    declare -a REF_DIR=$(mktemp -d genome_refXXXXXX)
    mkdir $REF_DIR/Bisulfite_Genome
    mkdir $REF_DIR/Bisulfite_Genome/CT_conversion
    mkdir $REF_DIR/Bisulfite_Genome/GA_conversion

    mv ~{reference_fasta} $REF_DIR
    mv ~{reference_fasta_index} $REF_DIR
    mv ~{fwd_converted_reference_fasta} $REF_DIR/Bisulfite_Genome/CT_conversion
    mv ~{sep=' ' fwd_bowtie2_index_files} $REF_DIR/Bisulfite_Genome/CT_conversion
    mv ~{rev_converted_reference_fasta} $REF_DIR/Bisulfite_Genome/GA_conversion
    mv ~{sep=' ' rev_bowtie2_index_files} $REF_DIR/Bisulfite_Genome/GA_conversion

    # use bismark with bowtie2 to align the read
    # icpc to not attach read ID descriptions to read ID
    /Bismark/bismark $REF_DIR/ \
      --bowtie2 \
      --icpc \
      ~{true="-1" false="" paired_end_run} ~{fastq_input_read_a} \
      ~{true="-2" false="" paired_end_run} ~{fastq_input_read_b} \
      ~{true="-X 2000" false="" paired_end_run} \
      ~{true="--pbat" false="" directional}

    # rename defualt bismark output name to match basename
    mv ~{bismark_mapped_bam_output_name} ~{mapped_bam_output_name}
    mv ~{bismark_mapping_report_output_name} ~{mapping_report_output_name}
  >>>

  runtime {
    docker: "quay.io/broadinstitute/bismark:0.21.0"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 3 * input file size
    disks: "local-disk " + ceil(3 * (if (bam_size + reference_size) < 1 then 1 else (bam_size + reference_size))) + " HDD"
    cpu: 2
    memory: "15 GB"
  }

  output {
    File mapped_bam_output = mapped_bam_output_name
    File mapping_report_output = mapping_report_output_name
    File monitoring_log = "monitoring.log"
  }
}

task AttachBarcodes {
    input {
      File mapped_bam_input
      File tagged_unmapped_bam_input
      File reference_fasta
      File reference_dictionary
      Int cut_length
      String output_base_name
      File? monitoring_script
    }

      # input file size
      Float input_size = size(mapped_bam_input, "GB") +
        size(tagged_unmapped_bam_input, "GB") +
        size(reference_fasta, "GB") +
        size(reference_dictionary, "GB")

      # output names for bam with extracted barcodes
      String tagged_mapped_bam_output_name = output_base_name + ".tagged_mapped.bam"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # create an unmapped bam
    java -Xmx3000m -jar /usr/picard/picard.jar MergeBamAlignment \
      SORT_ORDER="unsorted" \
      ADD_MATE_CIGAR=true \
      R1_TRIM=~{cut_length} R2_TRIM=~{cut_length} \
      IS_BISULFITE_SEQUENCE=true \
      UNMAP_CONTAMINANT_READS=false \
      UNMAPPED_BAM=~{tagged_unmapped_bam_input} \
      ALIGNED_BAM=~{mapped_bam_input} \
      REFERENCE_SEQUENCE=~{reference_fasta} \
      OUTPUT=~{tagged_mapped_bam_output_name}
  >>>

  # use docker image for given tool cutadapat
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2 * input file size
    disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3500 MiB"
  }

  output {
    File tagged_mapped_bam_output = tagged_mapped_bam_output_name
    File monitoring_log = "monitoring.log"
  }
}

# merge read 1 and read 2 bams
task MergeBams {
    input {
      File bam_input_a
      File bam_input_b
      String output_base_name
      File? monitoring_script
    }

      # input file size
      Float input_size = size(bam_input_a, "GB") + size(bam_input_b, "GB")

      # output name for merged bam
      String merged_bam_output_name = output_base_name + ".merged.bam"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # merge 2 bams with samtools
    samtools merge \
      ~{merged_bam_output_name} \
      ~{bam_input_a} ~{bam_input_b}
  >>>

  runtime {
    docker: "quay.io/broadinstitute/samtools:1.9"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2 * input file size
    disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    File merged_bam_output = merged_bam_output_name
    File monitoring_log = "monitoring.log"
  }
}

# sort the mapped reads
task Sort {
    input {
      File bam_input
      String output_base_name
      File? monitoring_script
    }

      # input file size
      Float input_size = size(bam_input, "GB")

      # output name for sorted bam
      String bam_sort_output_name = output_base_name + ".sorted.bam"

  # sort with samtools
  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    java -Xmx3000m -jar /usr/picard/picard.jar SortSam \
      INPUT=~{bam_input} \
      SORT_ORDER=coordinate \
      MAX_RECORDS_IN_RAM=300000 \
      OUTPUT=~{bam_sort_output_name}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 3.25 * input file size
    disks: "local-disk " + ceil(3.25 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3500 MiB"
  }

  output {
    File bam_sort_output = bam_sort_output_name
    File monitoring_log = "monitoring.log"
  }
}

# filter bam by removing duplicates
task FilterDuplicates {
    input {
      File bam_input
      Boolean remove_duplicates
      String output_base_name
      File? monitoring_script
    }

      # input file size
      Float input_size = size(bam_input, "GB")

      # output name for filtered reads
      String bam_remove_dup_output_name = output_base_name + ".filtered.duplicates.bam"
      String metric_remove_dup_output_name = output_base_name + ".filtered.duplicate_metrics"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    java -Xmx3000m -jar /usr/picard/picard.jar MarkDuplicates \
      INPUT=~{bam_input} \
      OUTPUT=~{bam_remove_dup_output_name} \
      METRICS_FILE=~{metric_remove_dup_output_name} \
      REMOVE_DUPLICATES=~{remove_duplicates}
  >>>

  runtime {
     docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
     # if the input size is less than 1 GB adjust to min input size of 1 GB
     # disks should be set to 2 * input file size
     disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
     cpu: 1
     memory: "3500 MiB"
  }

  output {
    File bam_remove_dup_output = bam_remove_dup_output_name
    File metric_remove_dup_output = metric_remove_dup_output_name
    File monitoring_log = "monitoring.log"
  }
}

# filter bam by mapping quality
task FilterMapQuality {
    input {
      File bam_input
      Int min_map_quality
      String output_base_name
      File? monitoring_script
    }

      # input file size
      Float input_size = size(bam_input, "GB")

      # output name for filtered read
      String bam_filter_above_min_mapq_output_name = output_base_name + ".filtered.above_min_map_quality.bam"
      String bam_filter_below_min_mapq_output_name = output_base_name + ".filtered.below_min_map_quality.bam"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # filter for a map quality
    # -b output is bam, -h include header, -q reads with mapping quality >=
    samtools view \
      -bhq~{min_map_quality} \
      -U ~{bam_filter_below_min_mapq_output_name} \
      ~{bam_input} > ~{bam_filter_above_min_mapq_output_name}
  >>>

  runtime {
    docker: "quay.io/broadinstitute/samtools:1.9"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2 * input file size
    disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    File bam_filter_above_min_mapq_output = bam_filter_above_min_mapq_output_name
    File bam_filter_below_min_mapq_output = bam_filter_below_min_mapq_output_name
    File monitoring_log = "monitoring.log"
  }
}

# produce methylation report using bisamrk's methyation extractor
task GetMethylationReport {
    input {
      File bam_input
      String output_base_name
      File? monitoring_script
    }

  # input file size
  Float input_size = size(bam_input, "GB")

  # output name for reports
  String methylation_report_basename = basename(bam_input, ".bam")

  String methylation_mbias_report_bismark_name = methylation_report_basename + ".M-bias.txt"
  String methylation_mbias_report_output_name = output_base_name + ".Mbias_report.txt"

  String methylation_splitting_report_bismark_name = methylation_report_basename + "_splitting_report.txt"
  String methylation_splitting_report_output_name = output_base_name + ".splitting_report.txt"

  String methylation_CpG_context_report_bismark_name = "CpG_context_" + methylation_report_basename + ".txt"
  String methylation_CpG_context_report_output_name = output_base_name + ".CpG_context_report.txt"

  String methylation_non_CpG_context_report_bismark_name = "Non_CpG_context_" + methylation_report_basename + ".txt"
  String methylation_non_CpG_context_report_output_name = output_base_name + ".non_CpG_context_report.txt"



  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # get methylation report using Bismark
    /Bismark/bismark_methylation_extractor \
      --comprehensive \
      --merge_non_CpG \
      --report \
      ~{bam_input}

    # rename outputs
    mv ~{methylation_mbias_report_bismark_name} ~{methylation_mbias_report_output_name}
    mv ~{methylation_splitting_report_bismark_name} ~{methylation_splitting_report_output_name}
    mv ~{methylation_CpG_context_report_bismark_name} ~{methylation_CpG_context_report_output_name}
    mv ~{methylation_non_CpG_context_report_bismark_name} ~{methylation_non_CpG_context_report_output_name}
  >>>

  runtime {
    docker: "quay.io/broadinstitute/bismark:0.21.0"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    disks: "local-disk " + ceil(4.5 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    File methylation_mbias_report_output = methylation_mbias_report_output_name
    File methylation_splitting_report_output = methylation_splitting_report_output_name
    File methylation_CpG_context_report_output = methylation_CpG_context_report_output_name
    File methylation_non_CpG_context_report_output = methylation_non_CpG_context_report_output_name
    File monitoring_log = "monitoring.log"
  }
}

# add read group name for support of downstream tools including GATK
task AddReadGroup {
    input {
      File bam_input
      String read_group_library_name
      String read_group_platform_name
      String read_group_platform_unit_name
      String read_group_platform_sample_name
      String output_base_name
      File? monitoring_script
    }

      # input file size
      Float input_size = size(bam_input, "GB")

      # output name for bam with a read group
      String added_read_group_output_bam_name = output_base_name + ".with_added_read_group.bam"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
       chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    gatk --java-options "-Xms2500m -Xmx3000m" \
    AddOrReplaceReadGroups \
      --INPUT ~{bam_input} \
      --RGLB ~{read_group_library_name} \
      --RGPL ~{read_group_platform_name} \
      --RGPU ~{read_group_platform_unit_name} \
      --RGSM ~{read_group_platform_sample_name} \
      --OUTPUT ~{added_read_group_output_bam_name}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2 * input file size
    disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3500 MiB"
  }

  output {
    File bam_with_read_group_output = added_read_group_output_bam_name
    File monitoring_log = "monitoring.log"
  }
}

# create a VCF contain locus methylation information
task MethylationTypeCaller {
    input {
      File bam_input
      File reference_fasta
      File reference_fasta_index
      File reference_dictionary
      String output_base_name
      File? monitoring_script
    }

  # input file size
  Float input_size = size(bam_input, "GB") +
                      size(reference_fasta, "GB") +
                      size(reference_fasta_index, "GB") +
                      size(reference_dictionary, "GB")

  # output name for VCF and its index
  String methylation_vcf_output_name = output_base_name + ".vcf"
  String methylation_vcf_index_output_name = methylation_vcf_output_name + ".idx"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    gatk --java-options "-Xms2500m -Xmx3000m" \
    MethylationTypeCaller \
      --input ~{bam_input} \
      --reference ~{reference_fasta} \
      --output ~{methylation_vcf_output_name} \
      --create-output-variant-index 
  >>>

  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:4.3.0.0"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    disks: "local-disk " + ceil(4.5 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3500 MiB"
  }

  output {
    File methylation_vcf = methylation_vcf_output_name
    File methylation_vcf_index = methylation_vcf_index_output_name
    File monitoring_log = "monitoring.log"
  }
}

# create a ALLC from VCF 
task VCFtoALLC {
    input {
      File methylation_vcf_output_name
      Int disk_size_gib = if size(methylation_vcf_output_name, "GiB") < 1 then 5 else ceil(9 * size(methylation_vcf_output_name, "GiB"))
      Float mem_size_gib = 3.5
    }

  # input file size

  # output name for VCF and its index
  String methylation_allc_output_name = sub(basename(methylation_vcf_output_name), ".vcf$", ".allc")

  command <<<
    set -euo pipefail

    python3 /tools/convert-vcf-to-allc.py -i ~{methylation_vcf_output_name} -o ~{methylation_allc_output_name}
  >>>

  runtime {
    docker: "quay.io/humancellatlas/vcftoallc:v0.0.1"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    disks: "local-disk ~{disk_size_gib} HDD"
    cpu: 1
    memory: "~{mem_size_gib} GiB"
  }

  output {
    File methylation_allc = methylation_allc_output_name
  }
}


# get the number of sites the have coverage of 1 or more
task ComputeCoverageDepth {
    input {
      File bam
      File reference_fasta
      File? monitoring_script
    }

      Float input_size = size(bam, "GB") + size(reference_fasta, "GB")

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # get samtools output and pipe into word count
    samtools depth \
      --reference ~{reference_fasta} \
      ~{bam} \
      | wc -l
  >>>

  runtime {
    docker: "quay.io/broadinstitute/samtools:1.9"
    disks: "local-disk " + ceil(1.25 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    Int total_depth_count = read_int(stdout())
    File monitoring_log = "monitoring.log"
  }
}

# create bam index
task IndexBam {
    input {
      File bam_input
      String output_base_name
      File? monitoring_script
    }

      # input file size
      Float input_size = size(bam_input, "GB")

      # output name for indexed bam
      String bam_index_output_name = output_base_name + ".bam.bai"

  command <<<
    set -euo pipefail

    # if the WDL/task contains a monitoring script as input
    if [ ! -z "~{monitoring_script}" ]; then
      chmod a+x ~{monitoring_script}
      ~{monitoring_script} > monitoring.log &
    else
      echo "No monitoring script given as input" > monitoring.log &
    fi

    # index bam with samtools
    samtools index -b ~{bam_input} ~{bam_index_output_name}
  >>>

  runtime {
    docker: "quay.io/broadinstitute/samtools:1.9"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2 * input file size
    disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    File bam_index_output = bam_index_output_name
    File monitoring_log = "monitoring.log"
  }
}
