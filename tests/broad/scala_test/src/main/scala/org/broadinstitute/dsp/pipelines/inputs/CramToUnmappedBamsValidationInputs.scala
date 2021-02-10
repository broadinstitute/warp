package org.broadinstitute.dsp.pipelines.inputs

import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}
import io.circe.{Decoder, Encoder}

object CramToUnmappedBamsValidationInputs
    extends EncodableInputs[CramToUnmappedBamsValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifyGermlineSingleSample")
  override implicit val decoder: Decoder[CramToUnmappedBamsValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[CramToUnmappedBamsValidationInputs] =
    deriveConfiguredEncoder
}
case class CramToUnmappedBamsValidationInputs(
    )
