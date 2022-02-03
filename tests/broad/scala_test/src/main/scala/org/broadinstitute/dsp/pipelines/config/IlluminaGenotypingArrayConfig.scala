package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.PapiVersion.PAPIv2
import org.broadinstitute.dsp.pipelines.commandline.{
  CromwellEnvironment,
  PapiVersion,
  WorkflowTestCategory
}

case class IlluminaGenotypingArrayConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing,
    env: CromwellEnvironment = CromwellEnvironment.Dev,
    useCallCaching: Boolean = true,
    updateTruth: Boolean = false,
    truthBranch: String = "master",
    papiVersion: PapiVersion = PAPIv2,
    useTimestamp: Option[String] = None
) extends BaseConfig
