package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.{Decoder, Encoder}
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}

object OptimusValidationInputs
    extends EncodableInputs[OptimusValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifyOptimus")
  override implicit val decoder: Decoder[OptimusValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[OptimusValidationInputs] =
    deriveConfiguredEncoder
}

case class OptimusValidationInputs(
    testBam: URI,
    truthBam: URI,
)
