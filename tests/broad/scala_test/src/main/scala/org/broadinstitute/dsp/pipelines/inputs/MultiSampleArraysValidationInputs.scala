package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.{Decoder, Encoder}
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}

object MultiSampleArraysValidationInputs
    extends EncodableInputs[MultiSampleArraysValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifyMultiSampleArrays")
  override implicit val decoder: Decoder[MultiSampleArraysValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[MultiSampleArraysValidationInputs] =
    deriveConfiguredEncoder
}
case class MultiSampleArraysValidationInputs(
    truth_vcf: URI,
    test_vcf: URI
)
