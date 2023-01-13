package org.broadinstitute.dsp.pipelines.inputs

import cats.syntax.either._
import io.circe.parser._

class CheckFingerprintInputs(inputs: String) {
  import org.broadinstitute.clio.JsonUtils.JsonOps

  private val parsed = parse(inputs).valueOr { e =>
    throw new RuntimeException("Could not get JSON object from inputs", e)
  }

  def getBasename(workflowName: String): String =
    parsed.unsafeGet[String](s"$workflowName.output_basename")

  def getSampleAlias(workflowName: String): String =
    parsed.unsafeGet[String](s"$workflowName.sample_alias")
}
