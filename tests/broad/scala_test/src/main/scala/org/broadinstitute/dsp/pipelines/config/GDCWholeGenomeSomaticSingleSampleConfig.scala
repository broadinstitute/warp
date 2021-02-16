package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.util.DataType
import org.broadinstitute.dsp.pipelines.commandline.PapiVersion.PAPIv2
import org.broadinstitute.dsp.pipelines.commandline.{
  CromwellEnvironment,
  PapiVersion,
  WorkflowTestCategory
}

case class GDCWholeGenomeSomaticSingleSampleConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing,
    truthBranch: String = "develop",
    dataType: DataType = DataType.WGS,
    useTimestamp: Option[String] = None,
    useCallCaching: Boolean = true,
    updateTruth: Boolean = false,
    papiVersion: PapiVersion = PAPIv2,
    env: CromwellEnvironment = CromwellEnvironment.Dev
) extends BaseConfig
