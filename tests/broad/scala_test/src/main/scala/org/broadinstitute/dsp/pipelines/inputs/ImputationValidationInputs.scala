package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.{Decoder, Encoder}
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}

object ImputationValidationInputs
    extends EncodableInputs[ImputationValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifyImputation")
  override implicit val decoder: Decoder[ImputationValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[ImputationValidationInputs] =
    deriveConfiguredEncoder
}
// TODO add in optional single_sample_vcfs if necessary
case class ImputationValidationInputs(
    truth_metrics: Seq[URI],
    test_metrics: Seq[URI],
    truth_vcf: URI,
    test_vcf: URI
)
