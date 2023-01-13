package org.broadinstitute.dsp.pipelines.inputs

import cats.syntax.either._
import io.circe.Json
import io.circe.parser._

class GermlineSingleSampleInputs(inputs: String) {
  import org.broadinstitute.clio.JsonUtils.JsonOps

  private val parsed = parse(inputs).valueOr { e =>
    throw new RuntimeException("Could not get JSON object from inputs", e)
  }

  private def sampleAndUnmappedBams(workflowName: String) =
    parsed.unsafeGet[Json](
      s"$workflowName.sample_and_unmapped_bams"
    )

  def getBaseFileName(workflowName: String): String =
    sampleAndUnmappedBams(workflowName).unsafeGet[String]("base_file_name")
}
