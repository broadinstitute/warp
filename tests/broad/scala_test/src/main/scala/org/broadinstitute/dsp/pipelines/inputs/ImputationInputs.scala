package org.broadinstitute.dsp.pipelines.inputs

import cats.syntax.either._
import io.circe.parser._

class ImputationInputs(inputs: String) {
  import org.broadinstitute.clio.JsonUtils.JsonOps

  private val parsed = parse(inputs).valueOr { e =>
    throw new RuntimeException("Could not get JSON object from inputs", e)
  }

  def getBasename(workflowName: String): String =
    parsed.unsafeGet[String](s"$workflowName.output_callset_name")

  def getHaplotypeDatabase(workflowName: String): String =
    parsed.unsafeGet[String](s"$workflowName.haplotype_database")

  def getSplitOutputToSingleSample(workflowName: String): String =
    parsed.unsafeGet[String](s"$workflowName.split_output_to_single_sample")

  def getSingleSampleVcfs(workflowName: String): Array =
    parsed.unsafeGet[String](s"$workflowName.single_sample_vcfs")

  def getSingleSampleVcfsIndices(workflowName: String): Array =
    parsed.unsafeGet[String](s"$workflowName.single_sample_vcf_indices")
}
