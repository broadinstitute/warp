package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.util.DataType
import org.broadinstitute.dsp.pipelines.commandline.PapiVersion.PAPI
import org.broadinstitute.dsp.pipelines.commandline.{
  PapiVersion,
  CromwellEnvironment,
  WorkflowTestCategory
}

case class GermlineCloudWorkflowConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing,
    truthBranch: String = "develop",
    dataType: DataType = DataType.WGS,
    useTimestamp: Option[String] = None,
    useCallCaching: Boolean = true,
    updateTruth: Boolean = false,
    papiVersion: PapiVersion = PAPI,
    env: CromwellEnvironment = CromwellEnvironment.Dev
) extends BaseConfig
