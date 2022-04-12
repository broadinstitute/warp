package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.PapiVersion.PAPIv2
import org.broadinstitute.dsp.pipelines.commandline.{
  CromwellEnvironment,
  PapiVersion,
  PipelineTestType,
  WorkflowTestCategory
}

case class CloudWorkflowConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing,
    env: CromwellEnvironment = CromwellEnvironment.Dev,
    pipeline: PipelineTestType = PipelineTestType.ExomeGermlineSingleSample,
    truthBranch: String = "develop",
    papiVersion: PapiVersion = PAPIv2,
    updateTruth: Boolean = false,
    useTimestamp: Option[String] = None,
    useCallCaching: Boolean = true
) extends BaseConfig
