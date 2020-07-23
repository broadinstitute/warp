package org.broadinstitute.dsp.pipelines.config
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment

case class GenotypeConcordanceConfig(
    nistVersion: String = "3.3.2",
    env: CromwellEnvironment = CromwellEnvironment.Dev,
    useLatestNA12878: Boolean = false
)
