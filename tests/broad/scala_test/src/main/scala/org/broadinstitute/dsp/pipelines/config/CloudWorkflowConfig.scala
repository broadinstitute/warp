package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.{
  CromwellEnvironment,
  PipelineTestType,
  WorkflowTestCategory
}

case class CloudWorkflowConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing,
    env: CromwellEnvironment = CromwellEnvironment.Dev,
    pipeline: PipelineTestType = PipelineTestType.ExomeGermlineSingleSample,
    truthBranch: String = "develop",
    updateTruth: Boolean = false,
    useTimestamp: Option[String] = None,
    useCallCaching: Boolean = true
) extends BaseConfig
