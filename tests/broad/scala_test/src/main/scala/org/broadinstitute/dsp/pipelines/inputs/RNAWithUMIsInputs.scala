package org.broadinstitute.dsp.pipelines.inputs

import cats.syntax.either._
import io.circe.parser._

class RNAWithUMIsInputs(inputs: String) {
  import org.broadinstitute.clio.JsonUtils.JsonOps

  private val parsed = parse(inputs).valueOr { e =>
    throw new RuntimeException("Could not get JSON object from inputs", e)
  }

  def getBasename(workflowName: String): String =
    parsed.unsafeGet[String](s"$workflowName.output_basename")
}
