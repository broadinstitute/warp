package org.broadinstitute.dsp.pipelines.inputs

import cats.syntax.either._
import io.circe.parser._

class SingleSampleArraysInputs(inputs: String) {
  import org.broadinstitute.clio.JsonUtils.JsonOps

  private val parsed = parse(inputs).valueOr { e =>
    throw new RuntimeException("Could not get JSON object from inputs", e)
  }

  def getChipwellBarcode(workflowName: String): String =
    parsed.unsafeGet[String](s"$workflowName.chip_well_barcode")

  def getArraysMetadataPath(workflowName: String): Option[String] =
    parsed.unsafeGet[Option[String]](s"$workflowName.arrays_metadata_path")

  def getBeadPoolManifestFile(workflowName: String): Option[String] =
    parsed.unsafeGet[Option[String]](s"$workflowName.bead_pool_manifest_file")

  def getBeadPoolManifestFilename(workflowName: String): Option[String] =
    parsed.unsafeGet[Option[String]](
      s"$workflowName.bead_pool_manifest_filename")
}
