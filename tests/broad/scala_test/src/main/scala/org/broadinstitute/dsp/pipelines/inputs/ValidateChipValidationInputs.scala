package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.{Decoder, Encoder}
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}

object ValidateChipValidationInputs
    extends EncodableInputs[ValidateChipValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifyValidateChip")
  override implicit val decoder: Decoder[ValidateChipValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[ValidateChipValidationInputs] =
    deriveConfiguredEncoder
}

case class ValidateChipValidationInputs(
    testGtc: URI,
    truthGtc: URI,
    beadPoolManifestFile: URI,
    testVcf: URI,
    truthVcf: URI,
    testGenotypeConcordanceVcf: URI,
    truthGenotypeConcordanceVcf: URI,
    testIndelGenotypeConcordanceVcf: URI,
    truthIndelGenotypeConcordanceVcf: URI,
    testMetrics: Seq[URI],
    truthMetrics: Seq[URI]
)
