package org.broadinstitute.dsp.pipelines.commandline

import org.broadinstitute.dsp.pipelines.config._

case class Config(
    test: WorkflowTestType = WorkflowTestType.Dummy,
    annotationFiltrationConfig: AnnotationFiltrationConfig =
      AnnotationFiltrationConfig(),
    arraysConfig: ArraysConfig = ArraysConfig(),
    checkFingerprintConfig: CheckFingerprintConfig = CheckFingerprintConfig(),
    cramToUnmappedBamsConfig: CramToUnmappedBamsConfig =
      CramToUnmappedBamsConfig(),
    genotypeConcordanceConfig: GenotypeConcordanceConfig =
      GenotypeConcordanceConfig(),
    gdcWholeGenomeSomaticSingleSampleConfig: GDCWholeGenomeSomaticSingleSampleConfig =
      GDCWholeGenomeSomaticSingleSampleConfig(),
    germlineCloudConfig: GermlineCloudWorkflowConfig =
      GermlineCloudWorkflowConfig(),
    illuminaGenotypingArrayConfig: IlluminaGenotypingArrayConfig =
      IlluminaGenotypingArrayConfig(),
    imputationConfig: ImputationConfig = ImputationConfig(),
    rnaWithUmisConfig: RNAWithUmisConfig = RNAWithUmisConfig(),
    somaticCloudWorkflowConfig: SomaticCloudWorkflowConfig =
      SomaticCloudWorkflowConfig(),
    validateChipConfig: ValidateChipConfig = ValidateChipConfig()
)
