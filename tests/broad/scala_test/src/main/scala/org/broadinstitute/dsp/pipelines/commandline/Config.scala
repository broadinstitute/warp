package org.broadinstitute.dsp.pipelines.commandline

import org.broadinstitute.dsp.pipelines.config._

case class Config(
    test: WorkflowTestType = WorkflowTestType.Dummy,
    annotationFiltrationConfig: AnnotationFiltrationConfig =
      AnnotationFiltrationConfig(),
    arraysConfig: ArraysConfig = ArraysConfig(),
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
    somaticCloudWorkflowConfig: SomaticCloudWorkflowConfig =
      SomaticCloudWorkflowConfig(),
    validateChipConfig: ValidateChipConfig = ValidateChipConfig()
)
