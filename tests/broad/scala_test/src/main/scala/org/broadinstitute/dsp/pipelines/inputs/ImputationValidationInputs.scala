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

case class ImputationValidationInputs(
    haplotype_database: String,
    split_output_to_single_sample: Boolean,
    input_single_sample_vcfs: Option[String],
    input_single_sample_vcfs_indices: Option[String],
    truth_metrics: Seq[URI],
    test_metrics: Seq[URI],
    truth_vcf: URI,
    test_vcf: URI
)
