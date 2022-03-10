package org.broadinstitute.dsp.pipelines.inputs

import better.files.File
import cats.syntax.either._
import io.circe.parser._

class VariantCallingInputs(inputs: String) {

  import org.broadinstitute.clio.JsonUtils.JsonOps

  private val parsed = parse(inputs).valueOr { e =>
    throw new RuntimeException("Could not get JSON object from inputs", e)
  }

  def getGvcfBasename(workflowName: String): String =
    File(
      parsed.unsafeGet[String](
        s"$workflowName.final_vcf_base_name"
      )
    ).name
}
