package org.broadinstitute.dsp.pipelines.inputs

import io.circe.{Decoder, Json}
import io.circe.parser.parse
import org.broadinstitute.clio.util.json.ModelAutoDerivation

trait CromwellWorkflowInputs[A] extends ModelAutoDerivation {

  def workflowNames: Seq[String]
  implicit val decoder: Decoder[A]

  def stripCommentsAndWorkflowNames(inputs: List[String]): String = {
    inputs
      .filter(!_.trim.startsWith("\"##"))
      .map(l => workflowNames.fold(l)((acc, cur) => acc.replace(s"$cur.", "")))
      .mkString("\n")
  }

  def unmarshall(parsedJson: Json): Decoder.Result[A] = {
    parsedJson.as[A]
  }

  def apply(text: List[String]): A = {
    val rawJson = stripCommentsAndWorkflowNames(text)
    val json = parse(rawJson)
      .flatMap(unmarshall)
    json match {
      case Left(e) =>
        throw new RuntimeException(
          "There was an error parsing the inputs to the workflow",
          e
        )
      case Right(inputs) => inputs
    }
  }
}
