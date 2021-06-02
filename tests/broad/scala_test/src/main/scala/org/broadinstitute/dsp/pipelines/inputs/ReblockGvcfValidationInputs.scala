package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.{Decoder, Encoder}
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}

object ReblockGvcfValidationInputs
    extends EncodableInputs[ReblockGvcfValidationInputs] {

  override def workflowNames: Seq[String] = Seq("VerifyGvcf")

  override implicit val decoder: Decoder[ReblockGvcfValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[ReblockGvcfValidationInputs] =
    deriveConfiguredEncoder

}

case class ReblockGvcfValidationInputs(
    testGvcf: URI,
    truthGvcf: URI,
)
