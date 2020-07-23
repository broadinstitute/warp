package org.broadinstitute.dsp.pipelines.commandline

case class ClioConfig(
    clioUrl: String,
    clioPort: Int,
    clioUseHttps: Boolean,
    serviceAccountJson: Option[String]
)
