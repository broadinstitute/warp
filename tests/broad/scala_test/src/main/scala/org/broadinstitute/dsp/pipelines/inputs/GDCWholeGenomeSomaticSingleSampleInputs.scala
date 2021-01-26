package org.broadinstitute.dsp.pipelines.inputs


class GDCWholeGenomeSomaticSingleSampleInputs(inputs: String) {

  private val parsed = parse(inputs).valueOr { e =>
    throw new RuntimeException("Could not get JSON object from inputs", e)
  }

  def getBaseFileName(workflowName: String): String =
    parsed.unsafeGet[String](
      s"$workflowName.base_file_name")
}
