package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.{Decoder, Encoder}
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}

object UltimaGenomicsJointGenotypingValidationInputs
    extends EncodableInputs[UltimaGenomicsJointGenotypingValidationInputs] {

  override def workflowNames: Seq[String] = Seq("VerifyJointGenotyping")

  override implicit val decoder
    : Decoder[UltimaGenomicsJointGenotypingValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder
    : Encoder[UltimaGenomicsJointGenotypingValidationInputs] =
    deriveConfiguredEncoder

}

case class UltimaGenomicsJointGenotypingValidationInputs(
    testMetrics: Seq[URI],
    truthMetrics: Seq[URI],
    testFingerprint: URI,
    truthFingerprint: URI,
    testVcfs: Seq[URI],
    testVcfIndexes: Seq[URI],
    truthVcfs: Seq[URI],
    truthVcfIndexes: Seq[URI],
    testIntervals: Seq[URI],
    truthIntervals: Seq[URI]
)
