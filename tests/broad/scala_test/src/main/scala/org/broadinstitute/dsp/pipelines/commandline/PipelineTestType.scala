package org.broadinstitute.dsp.pipelines.commandline

import enumeratum.{Enum, EnumEntry}

import scala.collection.immutable

sealed abstract class PipelineTestType(val workflowName: String,
                                       val pipelineName: String,
                                       val workflowDir: String)
    extends EnumEntry

object PipelineTestType extends Enum[PipelineTestType] {
  override val values: immutable.IndexedSeq[PipelineTestType] = findValues

  // Broad Pipelines
  case object Arrays
      extends PipelineTestType(
        "TestArrays",
        "Arrays",
        "/broad/arrays/single_sample/"
      )
  case object BroadInternalRNAWithUMIs
      extends PipelineTestType(
        "TestBroadInternalRNAWithUMIs",
        "BroadInternalRNAWithUMIs",
        "/broad/internal/rna_seq/"
      )
  case object BroadInternalUltimaGenomics
      extends PipelineTestType(
        "TestBroadInternalUltimaGenomics",
        "BroadInternalUltimaGenomics",
        "/broad/internal/dna_seq/germline/single_sample/UltimaGenomics"
      )
  case object CheckFingerprint
      extends PipelineTestType(
        "TestCheckFingerprint",
        "CheckFingerprint",
        "/broad/qc"
      )
  case object CramToUnmappedBams
      extends PipelineTestType(
        "TestCramToUnmappedBams",
        "CramToUnmappedBams",
        "/broad/reprocessing/cram_to_unmapped_bams/"
      )
  case object ExternalExomeReprocessing
      extends PipelineTestType(
        "TestExternalExomeReprocessing",
        "ExternalExomeReprocessing",
        "/broad/reprocessing/external/exome/"
      )
  case object ExomeGermlineSingleSample
      extends PipelineTestType(
        "TestExomeGermlineSingleSample",
        "ExomeGermlineSingleSample",
        "/broad/dna_seq/germline/single_sample/exome/"
      )
  case object ExomeReprocessing
      extends PipelineTestType(
        "TestExomeReprocessing",
        "ExomeReprocessing",
        "/broad/reprocessing/exome/"
      )
  case object ExternalWholeGenomeReprocessing
      extends PipelineTestType(
        "TestExternalWholeGenomeReprocessing",
        "ExternalWholeGenomeReprocessing",
        "/broad/reprocessing/external/wgs/"
      )
  case object GDCWholeGenomeSomaticSingleSample
      extends PipelineTestType(
        "TestGDCWholeGenomeSomaticSingleSample",
        "GDCWholeGenomeSomaticSingleSample",
        "/broad/dna_seq/somatic/single_sample/wgs/gdc_genome/"
      )
  case object IlluminaGenotypingArray
      extends PipelineTestType(
        "TestIlluminaGenotypingArray",
        "IlluminaGenotypingArray",
        "/broad/genotyping/illumina/"
      )
  case object Imputation
      extends PipelineTestType(
        "TestImputation",
        "Imputation",
        "/broad/arrays/imputation"
      )
  case object JointGenotyping
      extends PipelineTestType(
        "TestJointGenotyping",
        "JointGenotyping",
        "/broad/dna_seq/germline/joint_genotyping/"
      )
  case object MultiSampleArrays
      extends PipelineTestType(
        "TestMultiSampleArrays",
        "MultiSampleArrays",
        "/broad/arrays/multi_sample/"
      )
  case object ReblockGVCF
      extends PipelineTestType(
        "TestReblockGVCF",
        "ReblockGVCF",
        "/broad/dna_seq/germline/joint_genotyping/reblocking/"
      )
  case object RNAWithUMIsPipeline
      extends PipelineTestType(
        "TestRNAWithUMIsPipeline",
        "RNAWithUMIsPipeline",
        "/broad/rna_seq/"
      )
  case object UltimaGenomicsJointGenotyping
      extends PipelineTestType(
        "TestUltimaGenomicsJointGenotyping",
        "UltimaGenomicsJointGenotyping",
        "/broad/dna_seq/germline/joint_genotyping/UltimaGenomics/"
      )
  case object UltimaGenomicsWholeGenomeGermline
      extends PipelineTestType(
        "TestUltimaGenomicsWholeGenomeGermline",
        "UltimaGenomicsWholeGenomeGermline",
        "/broad/dna_seq/germline/single_sample/ugwgs/"
      )
  case object ValidateChip
      extends PipelineTestType(
        "TestValidateChip",
        "ValidateChip",
        "/broad/arrays/validate_chip/"
      )
  case object VariantCalling
      extends PipelineTestType(
        "TestVariantCalling",
        "VariantCalling",
        "/broad/dna_seq/germline/variant_calling/"
      )
  case object WholeGenomeGermlineSingleSample
      extends PipelineTestType(
        "TestWholeGenomeGermlineSingleSample",
        "WholeGenomeGermlineSingleSample",
        "/broad/dna_seq/germline/single_sample/wgs/"
      )
  case object WholeGenomeReprocessing
      extends PipelineTestType(
        "TestWholeGenomeReprocessing",
        "WholeGenomeReprocessing",
        "/broad/reprocessing/wgs"
      )

  // Skylab Pipelines
  case object MultiSampleSmartSeq2
      extends PipelineTestType(
        "TestMultiSampleSmartSeq2",
        "MultiSampleSmartSeq2",
        "/skylab/smartseq2_multisample"
      )
  case object MultiSampleSmartSeq2SingleNucleus
      extends PipelineTestType(
        "TestMultiSampleSmartSeq2SingleNucleus",
        "MultiSampleSmartSeq2SingleNucleus",
        "/skylab/smartseq2_single_nucleus_multisample"
      )
  case object Optimus
      extends PipelineTestType(
        "TestOptimus",
        "Optimus",
        "/skylab/optimus/"
      )
  case object scATAC
      extends PipelineTestType(
        "TestscATAC",
        "scATAC",
        "/skylab/scATAC/"
      )
  case object SmartSeq2SingleSample
      extends PipelineTestType(
        "TestSmartSeq2SingleSample",
        "SmartSeq2SingleSample",
        "/skylab/smartseq2_single_sample"
      )

  // CEMBA Pipelines
  //case object CEMBA
  //    extends PipelineTestType(
  //      "TestCEMBA",
  //      "/cemba/cemba_methylcseq/"
  //    )
}
