package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.util.DataType.Exome
import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRunParameters,
  WorkflowTest
}
import org.broadinstitute.dsp.pipelines.commandline.WorkflowTestCategory
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  GermlineSingleSampleInputs,
  GermlineSingleSampleValidationInputs
}

class GermlineSingleSampleTester(testerConfig: GermlineCloudWorkflowConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends GermlineCloudWorkflowTester(testerConfig) {

  override val workflowName: String = s"${dataTypePrefix}GermlineSingleSample"

  override lazy val workflowDir: File =
    CromwellWorkflowTester.PipelineRoot / "broad" / "dna_seq" / "germline" / "single_sample" / dataTypeString

  override protected val validationWorkflowName: String =
    "VerifyGermlineSingleSample"

  protected lazy val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/germline_single_sample/$dataTypeString/$testTypeString/$timestamp/"
    )
  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/germline_single_sample/$dataTypeString/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override def generateRunParameters: Seq[WorkflowRunParameters] = {

    if (testerConfig.category.equals(WorkflowTestCategory.Load) &&
        testerConfig.dataType.equals(Exome)) {
      generateExomeLoadRunParameters(envString, resultsPrefix, timestamp)
    } else {
      logger.info(s"workflowInputRoot: $workflowInputRoot")
      inputFileNames.map { fileName =>
        val projectSampleAlias = fileName.replace(".json", "")
        val resultsPath =
          resultsPrefix.resolve(s"$projectSampleAlias/")
        val truthPath = truthPrefix.resolve(s"$projectSampleAlias/")

        WorkflowRunParameters(
          id = s"${envString}_$projectSampleAlias",
          workflowInputs = getInputContents(fileName),
          resultsCloudPath = resultsPath,
          truthCloudPath = truthPath
        )
      }
    }
  }

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val germlineSingleSampleInputs = new GermlineSingleSampleInputs(
      workflowTest.runParameters.workflowInputs
    )
    val outputBaseName =
      germlineSingleSampleInputs.getBaseFileName(workflowName)
    val gvcfBaseName =
      germlineSingleSampleInputs.getFinalGvcfBaseName(workflowName)
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath
    val metricsFileNames = ioUtil
      .listGoogleObjects(truthCloudPath)
      .filter(_.getPath.endsWith("metrics"))
      .map(uriToFilename)
    val validationInputs = GermlineSingleSampleValidationInputs(
      testMetrics = metricsFileNames.map(resultsCloudPath.resolve),
      truthMetrics = metricsFileNames.map(truthCloudPath.resolve),
      testCram = resultsCloudPath.resolve(s"$outputBaseName.cram"),
      testCrai = resultsCloudPath.resolve(s"$outputBaseName.cram.crai"),
      truthCram = truthCloudPath.resolve(s"$outputBaseName.cram"),
      truthCrai = truthCloudPath.resolve(s"$outputBaseName.cram.crai"),
      testGvcf = resultsCloudPath.resolve(s"$gvcfBaseName.g.vcf.gz"),
      truthGvcf = truthCloudPath.resolve(s"$gvcfBaseName.g.vcf.gz")
    )
    GermlineSingleSampleValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }

  private def generateExomeLoadRunParameters(
      envString: String,
      resultsPrefix: URI,
      timestamp: String
  ): Seq[WorkflowRunParameters] = {
    case class LoadSample(barcodeLane: String, sampleAlias: String)
    val loadInputTemplate =
      (workflowDir / "test_inputs" / "load_template.json").contentAsString

    val loadSamples = (workflowDir / "test_inputs" / "load_samples.tsv").lines
      .map(line => {
        val lineArray = line.split("\t")
        LoadSample(lineArray(0), lineArray(1))
      })

    loadSamples.map { s =>
      val id = s"${s.sampleAlias}.${s.barcodeLane}"
      val resultsPath = resultsPrefix.resolve(s"$timestamp/$id")
      val inputContents = loadInputTemplate
        .replace("{SAMPLE_ALIAS}", s.sampleAlias)
        .replace("{BARCODE_LANE}", s.barcodeLane)
      WorkflowRunParameters(
        id = s"${envString}_$id",
        workflowInputs = inputContents,
        resultsCloudPath = resultsPath,
        truthCloudPath = URI.create("notused")
      )
    }.toSeq
  }
}
