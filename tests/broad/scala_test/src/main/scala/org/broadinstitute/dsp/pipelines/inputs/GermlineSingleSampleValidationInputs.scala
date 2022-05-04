package org.broadinstitute.dsp.pipelines.inputs
import java.net.URI

import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}
import io.circe.{Decoder, Encoder}

object GermlineSingleSampleValidationInputs
    extends EncodableInputs[GermlineSingleSampleValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifyGermlineSingleSample")
  override implicit val decoder: Decoder[GermlineSingleSampleValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[GermlineSingleSampleValidationInputs] =
    deriveConfiguredEncoder
}

case class GermlineSingleSampleValidationInputs(
    testMetrics: Seq[URI],
    truthMetrics: Seq[URI],
    testCram: URI,
    testCrai: URI,
    truthCram: URI,
    truthCrai: URI,
    testGvcf: URI,
    testGvcfIndex: URI,
    truthGvcf: URI,
    truthGvcfIndex: URI
)
