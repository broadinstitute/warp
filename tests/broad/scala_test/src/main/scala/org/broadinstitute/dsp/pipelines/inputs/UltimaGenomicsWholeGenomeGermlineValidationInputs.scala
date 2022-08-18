package org.broadinstitute.dsp.pipelines.inputs
import java.net.URI

import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}
import io.circe.{Decoder, Encoder}

object UltimaGenomicsWholeGenomeGermlineValidationInputs
    extends EncodableInputs[UltimaGenomicsWholeGenomeGermlineValidationInputs] {
  override def workflowNames: Seq[String] =
    Seq("VerifyUltimaGenomicsWholeGenomeGermline")
  override implicit val decoder
    : Decoder[UltimaGenomicsWholeGenomeGermlineValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder
    : Encoder[UltimaGenomicsWholeGenomeGermlineValidationInputs] =
    deriveConfiguredEncoder
}

case class UltimaGenomicsWholeGenomeGermlineValidationInputs(
    testMetrics: Seq[URI],
    truthMetrics: Seq[URI],
    testCram: URI,
    truthCram: URI,
    testCrai: URI,
    truthCrai: URI,
    testVcf: URI,
    testVcfIndex: URI,
    truthVcf: URI,
    truthVcfIndex: URI,
    testFilteredVcf: URI,
    testFilteredVcfIndex: URI,
    truthFilteredVcf: URI,
    truthFilteredVcfIndex: URI,
    testGvcf: URI,
    testGvcfIndex: URI,
    truthGvcf: URI,
    truthGvcfIndex: URI
)
