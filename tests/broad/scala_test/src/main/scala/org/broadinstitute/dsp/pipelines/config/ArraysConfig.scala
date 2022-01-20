package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.PapiVersion.PAPIv2
import org.broadinstitute.dsp.pipelines.commandline.{
  CromwellEnvironment,
  PapiVersion,
  WorkflowTestCategory
}
import org.broadinstitute.dsp.pipelines.util.ArrayType

case class ArraysConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing,
    arrayType: ArrayType = ArrayType.Single,
    env: CromwellEnvironment = CromwellEnvironment.Dev,
    useCallCaching: Boolean = true,
    updateTruth: Boolean = false,
    truthBranch: String = "master",
    papiVersion: PapiVersion = PAPIv2,
    useTimestamp: Option[String] = None
) extends BaseConfig
