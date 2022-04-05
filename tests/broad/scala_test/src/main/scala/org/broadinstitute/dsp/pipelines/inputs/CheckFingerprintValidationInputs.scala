package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.{Decoder, Encoder}
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}

object CheckFingerprintValidationInputs
    extends EncodableInputs[CheckFingerprintValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifyCheckFingerprint")
  override implicit val decoder: Decoder[CheckFingerprintValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[CheckFingerprintValidationInputs] =
    deriveConfiguredEncoder
}

case class CheckFingerprintValidationInputs(
    test_metrics: Seq[URI],
    truth_metrics: Seq[URI],
    truth_fingerprint_vcf: URI,
    test_fingerprint_vcf: URI
)
