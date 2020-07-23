package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.PapiVersion.PAPI
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
    truthBranch: String = "develop",
    papiVersion: PapiVersion = PAPI,
    useTimestamp: Option[String] = None
) extends BaseConfig
