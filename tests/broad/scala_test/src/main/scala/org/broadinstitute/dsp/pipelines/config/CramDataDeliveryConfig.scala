package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.util.DataType
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment

case class CramDataDeliveryConfig(
    dataType: DataType = DataType.WGS,
    requester: String = "test.firec@gmail.com",
    leaveWorkspace: Boolean = false,
    env: CromwellEnvironment = CromwellEnvironment.Dev
)
