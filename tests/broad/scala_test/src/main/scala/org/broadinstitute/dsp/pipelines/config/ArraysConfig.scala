package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.PapiVersion.PAPI
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
    truthBranch: String = "develop",
    papiVersion: PapiVersion = PAPI,
    useTimestamp: Option[String] = None
) extends BaseConfig
