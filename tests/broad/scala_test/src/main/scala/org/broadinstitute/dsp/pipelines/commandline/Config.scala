package org.broadinstitute.dsp.pipelines.commandline

import org.broadinstitute.dsp.pipelines.config._

case class Config(
    test: WorkflowTestType = WorkflowTestType.Dummy,
    annotationFiltrationConfig: AnnotationFiltrationConfig =
      AnnotationFiltrationConfig(),
    arraysConfig: ArraysConfig = ArraysConfig(),
    broadInternalRNAWithUMIsConfig: BroadInternalRNAWithUMIsConfig =
      BroadInternalRNAWithUMIsConfig(),
    broadInternalUltimaGenomicsConfig: BroadInternalUltimaGenomicsConfig =
      BroadInternalUltimaGenomicsConfig(),
    checkFingerprintConfig: CheckFingerprintConfig = CheckFingerprintConfig(),
    cloudWorkflowConfig: CloudWorkflowConfig = CloudWorkflowConfig(),
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
    rnaWithUMIsConfig: RNAWithUMIsConfig = RNAWithUMIsConfig(),
    somaticCloudWorkflowConfig: SomaticCloudWorkflowConfig =
      SomaticCloudWorkflowConfig(),
    ultimaGenomicsWholeGenomeGermlineConfig: UltimaGenomicsWholeGenomeGermlineConfig =
      UltimaGenomicsWholeGenomeGermlineConfig(),
    ultimaGenomicsJointGenotypingConfig: UltimaGenomicsJointGenotypingConfig =
      UltimaGenomicsJointGenotypingConfig(),
    validateChipConfig: ValidateChipConfig = ValidateChipConfig()
)
