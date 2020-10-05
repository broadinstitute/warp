package org.broadinstitute.dsp.pipelines.inputs

import java.net.URI

import io.circe.{Decoder, Encoder}
import io.circe.generic.extras.semiauto.{
  deriveConfiguredDecoder,
  deriveConfiguredEncoder
}

object JointGenotypingValidationInputs
    extends EncodableInputs[JointGenotypingValidationInputs] {

  override def workflowNames: Seq[String] = Seq("VerifyJointGenotyping")

  override implicit val decoder: Decoder[JointGenotypingValidationInputs] =
    deriveConfiguredDecoder
  override implicit val encoder: Encoder[JointGenotypingValidationInputs] =
    deriveConfiguredEncoder

}

case class JointGenotypingValidationInputs(
    testMetrics: Seq[URI],
    truthMetrics: Seq[URI],
    testFingerprint: URI,
    truthFingerprint: URI,
    testVcfs: Seq[URI],
    testVcfIndexes: Seq[URI],
    truthVcfs: Seq[URI],
    truthVcfIndexes: Seq[URI],
    testIntervals: Seq[URI],
    truthIntervals: Seq[URI],
    isExome: Boolean
)
