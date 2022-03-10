package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}
import io.circe.{Decoder, Encoder}

object VariantCallingValidationInputs
    extends EncodableInputs[VariantCallingValidationInputs] {

  override def workflowNames: Seq[String] = Seq("VerifyGvcf")

  override implicit val decoder: Decoder[VariantCallingValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[VariantCallingValidationInputs] =
    deriveConfiguredEncoder
}

case class VariantCallingValidationInputs(
    testGvcf: URI,
    truthGvcf: URI,
)
