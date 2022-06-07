package org.broadinstitute.dsp.pipelines.commandline

import enumeratum._

import scala.collection.immutable

sealed trait WorkflowTestType extends EnumEntry

object WorkflowTestType extends Enum[WorkflowTestType] {

  override val values: immutable.IndexedSeq[WorkflowTestType] = findValues

  case object AllOfUs extends WorkflowTestType
  case object AnnotationFiltration extends WorkflowTestType
  case object Arrays extends WorkflowTestType
  case object BroadInternalRNAWithUMIs extends WorkflowTestType
  case object BroadInternalUltimaGenomics extends WorkflowTestType
  case object CheckFingerprint extends WorkflowTestType
  case object CramToUnmappedBams extends WorkflowTestType
  case object CloudWorkflow extends WorkflowTestType
  case object Dummy extends WorkflowTestType
  case object ExternalReprocessing extends WorkflowTestType
  case object GDCWholeGenomeSomaticSingleSample extends WorkflowTestType
  case object GenotypeConcordance extends WorkflowTestType
  case object GermlineSingleSample extends WorkflowTestType
  case object IlluminaGenotypingArray extends WorkflowTestType
  case object Imputation extends WorkflowTestType
  case object JointGenotyping extends WorkflowTestType
  case object Reprocessing extends WorkflowTestType
  case object ReblockGvcf extends WorkflowTestType
  case object RNAWithUMIs extends WorkflowTestType
  case object SomaticSingleSample extends WorkflowTestType
  case object UltimaGenomicsWholeGenomeGermline extends WorkflowTestType
  case object UltimaGenomicsJointGenotyping extends WorkflowTestType
  case object ValidateChip extends WorkflowTestType
  case object VariantCalling extends WorkflowTestType
}
