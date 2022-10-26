package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.util.DataType

import org.broadinstitute.dsp.pipelines.commandline.{
  CromwellEnvironment,
  WorkflowTestCategory
}

case class SomaticCloudWorkflowConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing,
    truthBranch: String = "master",
    dataType: DataType = DataType.Targeted,
    useTimestamp: Option[String] = None,
    useCallCaching: Boolean = true,
    updateTruth: Boolean = false,
    env: CromwellEnvironment = CromwellEnvironment.Dev
) extends BaseConfig
