package org.broadinstitute.dsp.pipelines.inputs

import io.circe.generic.extras.semiauto.deriveConfiguredDecoder
import io.circe.Decoder

object DummyInputs extends CromwellWorkflowInputs[DummyInputs] {
  override def workflowNames: Seq[String] = Seq("DummyWorkflow")
  override implicit val decoder: Decoder[DummyInputs] = deriveConfiguredDecoder
}

case class DummyInputs(
    message: String,
    dummyInt: Int,
    dummyString: String,
    dummyBoolean: Boolean,
    dummyOption: Option[String]
)
