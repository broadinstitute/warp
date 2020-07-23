package org.broadinstitute.dsp.pipelines.inputs
import io.circe.{Encoder, Json, JsonObject}

/**
  * WDL 1.0 specifies the WorkflowNames must be at
  * @tparam A
  */
trait EncodableInputs[A] extends CromwellWorkflowInputs[A] {

  implicit val encoder: Encoder[A]

  private def augmentedEncoder(workflowName: String) = encoder.mapJson(
    json =>
      json.mapObject(
        jsonObj =>
          jsonObj.toList.foldLeft(JsonObject.empty)((json, tup) => {
            val (field, value) = tup
            json.add(s"$workflowName.$field", value)
          })
    )
  )

  def marshall(
      inputs: A,
      workflowName: String = workflowNames.head
  ): Json =
    augmentedEncoder(workflowName)(inputs)
}
