package org.broadinstitute.dsp.pipelines.inputs

import io.circe.Decoder
import io.circe.generic.extras.semiauto.deriveConfiguredDecoder

object ArraysDataDeliveryInputs
    extends CromwellWorkflowInputs[ArraysDataDeliveryInputs] {
  override def workflowNames = Seq("ArraysDataDelivery", "MultiSampleArrays")
  override implicit val decoder: Decoder[ArraysDataDeliveryInputs] =
    deriveConfiguredDecoder
}

case class ArraysDataDeliveryInputs(
    authDomain: Option[String],
    billingProject: String,
    chipwellbarcodeVersionTsv: String,
    clioUrl: String,
    environment: String,
    preemptibleTries: Int,
    requester: String,
    sampleSet: Option[String],
    samplesMetadata: Option[String],
    serviceAccountJson: String,
    workspaceName: String,
    overrideDeliveryCheck: Boolean,
    deliver_multi_sample_vcf: Boolean
) extends FireCloudWdlInputs
