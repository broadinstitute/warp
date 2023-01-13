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

  def getOutputCallsetName(workflowName: String): String =
    parsed.unsafeGet[String](s"$workflowName.output_callset_name")

  def getSplitOutputToSingleSample(workflowName: String): Boolean =
    parsed.unsafeGet[Boolean](s"$workflowName.split_output_to_single_sample")

  def getSingleSampleVcfs(workflowName: String): Option[Seq[String]] =
    parsed.unsafeGet[Option[Seq[String]]](s"$workflowName.single_sample_vcfs")

  def getSingleSampleVcfsIndices(workflowName: String): Option[Seq[String]] =
    parsed.unsafeGet[Option[Seq[String]]](
      s"$workflowName.single_sample_vcf_indices")

  def getMultiSampleVcf(workflowName: String): Option[String] =
    parsed.unsafeGet[Option[String]](s"$workflowName.multi_sample_vcf")

  def getMultiSampleVcfIndices(workflowName: String): Option[String] =
    parsed.unsafeGet[Option[String]](s"$workflowName.multi_sample_vcf_index")
}
