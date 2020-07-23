package org.broadinstitute.dsp.pipelines.inputs

import cats.syntax.either._
import io.circe.parser.parse

class ValidateChipInputs(inputs: String) {
  import org.broadinstitute.clio.JsonUtils.JsonOps

  private val parsed = parse(inputs).valueOr { e =>
    throw new RuntimeException("Could not get JSON object from inputs", e)
  }

  def chipWellBarcode(workflowName: String) =
    parsed.unsafeGet[String](s"$workflowName.chip_well_barcode")

  def getBeadPoolManifestFile(workflowName: String): String =
    parsed.unsafeGet[String](s"$workflowName.bead_pool_manifest_file")
}
