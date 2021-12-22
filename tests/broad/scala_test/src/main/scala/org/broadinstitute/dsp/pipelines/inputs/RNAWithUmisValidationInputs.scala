package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.{Decoder, Encoder}
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}

object RNAWithUmisValidationInputs
    extends EncodableInputs[RNAWithUmisValidationInputs] {
  override def workflowNames: Seq[String] = Seq("VerifyRNAWithUMIs")
  override implicit val decoder: Decoder[RNAWithUmisValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[RNAWithUmisValidationInputs] =
    deriveConfiguredEncoder
}

case class RNAWithUmisValidationInputs(
    truth_metrics: Seq[URI],
    test_metrics: Seq[URI],
    truth_text_metrics: Seq[URI],
    test_text_metrics: Seq[URI],
    test_output_bam: URI,
    truth_output_bam: URI,
    test_transcriptome_bam: URI,
    truth_transcriptome_bam: URI,
    test_gene_tpm: URI,
    truth_gene_tpm: URI,
    test_gene_counts: URI,
    truth_gene_counts: URI,
    test_exon_counts: URI,
    truth_exon_counts: URI
)
