package org.broadinstitute.dsp.pipelines.inputs
import java.net.URI

import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}
import io.circe.{Decoder, Encoder}

object BroadInteralUGWrapperValidationInputs
    extends EncodableInputs[BroadInteralUGWrapperValidationInputs] {
  override def workflowNames: Seq[String] =
    Seq("VerifyBroadInternalUGWrapper")
  override implicit val decoder
    : Decoder[BroadInternalUGWrapperValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder
    : Encoder[BroadInternalUGWrapperValidationInputs] =
    deriveConfiguredEncoder
}

case class BroadInternalUGWrapperValidationInputs(
    testMetrics: Seq[URI],
    truthMetrics: Seq[URI],
    testCram: URI,
    truthCram: URI,
    testCrai: URI,
    truthCrai: URI,
    testVcf: URI,
    truthVcf: URI,
    testFilteredVcf: URI,
    truthFilteredVcf: URI,
    testGvcf: URI,
    truthGvcf: URI
)
