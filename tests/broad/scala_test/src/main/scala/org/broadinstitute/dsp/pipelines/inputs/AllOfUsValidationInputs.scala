package org.broadinstitute.dsp.pipelines.inputs
import java.net.URI

import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}
import io.circe.{Decoder, Encoder}

object AllOfUsValidationInputs
    extends EncodableInputs[AllOfUsValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifyAllOfUs")
  override implicit val decoder: Decoder[AllOfUsValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[AllOfUsValidationInputs] =
    deriveConfiguredEncoder
}

case class AllOfUsValidationInputs(
    testMetrics: Seq[URI],
    truthMetrics: Seq[URI],
    testCram: URI,
    testCrai: URI,
    truthCram: URI,
    truthCrai: URI,
    testGvcf: URI,
    truthGvcf: URI,
    testFiltrationReport: URI,
    truthFiltrationReport: URI,
    testSignificantVariantsVcf: URI,
    truthVariants: URI
)
