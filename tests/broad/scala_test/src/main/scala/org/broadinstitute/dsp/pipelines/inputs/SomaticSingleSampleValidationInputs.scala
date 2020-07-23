package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}
import io.circe.{Decoder, Encoder}

object SomaticSingleSampleValidationInputs
    extends EncodableInputs[SomaticSingleSampleValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifySomaticSingleSample")
  override implicit val decoder: Decoder[SomaticSingleSampleValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[SomaticSingleSampleValidationInputs] =
    deriveConfiguredEncoder
}

case class SomaticSingleSampleValidationInputs(
    testMetrics: Seq[URI],
    truthMetrics: Seq[URI],
    testCram: URI,
    testCrai: URI,
    truthCram: URI,
    truthCrai: URI,
)
