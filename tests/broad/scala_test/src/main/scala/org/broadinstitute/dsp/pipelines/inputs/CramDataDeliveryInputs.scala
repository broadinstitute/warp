package org.broadinstitute.dsp.pipelines.inputs

import io.circe.generic.extras.semiauto.deriveConfiguredDecoder
import io.circe.Decoder

object CramDataDeliveryInputs
    extends CromwellWorkflowInputs[CramDataDeliveryInputs] {
  override def workflowNames =
    Seq("CramDataDelivery")
  override implicit val decoder: Decoder[CramDataDeliveryInputs] =
    deriveConfiguredDecoder
}

case class CramDataDeliveryInputs(
    sampleSet: Option[String],
    workspaceName: String,
    billingProject: String,
    serviceAccountJson: String,
    requester: String,
    bucket: Option[String],
    clioUrl: String,
    newWorkspace: Boolean,
    clinicalDataDelivery: Boolean,
    skipCramDeliveryCheck: Boolean,
    overrideCramDelivered: Boolean,
    feeForService: Boolean,
    fireCloudHost: Option[String],
    rawlsHost: Option[String],
    libraryCurator: Option[String],
    methodsDev: Option[String]
) extends FireCloudWdlInputs {
  val projectName: String = sampleSet.getOrElse(workspaceName)
}
