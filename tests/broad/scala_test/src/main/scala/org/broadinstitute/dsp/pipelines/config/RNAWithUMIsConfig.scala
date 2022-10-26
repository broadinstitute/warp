package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.{
  CromwellEnvironment,
  WorkflowTestCategory
}

case class RNAWithUMIsConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing,
    truthBranch: String = "master",
    useTimestamp: Option[String] = None,
    useCallCaching: Boolean = true,
    updateTruth: Boolean = false,
    env: CromwellEnvironment = CromwellEnvironment.Dev
) extends BaseConfig
