package org.broadinstitute.dsp.pipelines.inputs
import java.net.URI

import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}
import io.circe.{Decoder, Encoder}

object BroadInternalUltimaGenomicsValidationInputs
    extends EncodableInputs[BroadInternalUltimaGenomicsValidationInputs] {
  override def workflowNames: Seq[String] =
    Seq("VerifyBroadInternalUltimaGenomics")
  override implicit val decoder
    : Decoder[BroadInternalUltimaGenomicsValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder
    : Encoder[BroadInternalUltimaGenomicsValidationInputs] =
    deriveConfiguredEncoder
}

case class BroadInternalUltimaGenomicsValidationInputs(
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
