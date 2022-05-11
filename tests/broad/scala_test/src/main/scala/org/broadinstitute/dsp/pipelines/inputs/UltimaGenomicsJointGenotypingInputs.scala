package org.broadinstitute.dsp.pipelines.inputs

import cats.syntax.either._
import io.circe.parser._

class UltimaGenomicsJointGenotypingInputs(inputs: String) {

  import org.broadinstitute.clio.JsonUtils.JsonOps

  private val parsed = parse(inputs).valueOr { e =>
    throw new RuntimeException("Could not get JSON object from inputs", e)
  }

  def getCallsetName(workflowName: String): String =
    parsed.unsafeGet[String](
      s"$workflowName.callset_name"
    )

  def getFinalGvcfBaseName(workflowName: String): String =
    parsed.unsafeGet[String](
      s"$workflowName.final_gvcf_base_name"
    )
}
