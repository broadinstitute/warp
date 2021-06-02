package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  SomaticSingleSampleInputs,
  SomaticSingleSampleValidationInputs
}

class SomaticSingleSampleTester(testerConfig: SomaticCloudWorkflowConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  val dataTypePrefix: String = dataTypePrefix(testerConfig.dataType)
  val dataTypeString: String = testerConfig.dataType.entryName.toLowerCase

  override val workflowName: String = s"${dataTypePrefix}SomaticSingleSample"

  val workflowDir
    : File = CromwellWorkflowTester.WarpRoot / "beta_pipelines" / "broad" / "somatic" / "single_sample" / dataTypeString

  override protected val validationWorkflowName: String =
    "VerifySomaticSingleSample"

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/somatic_single_sample/$dataTypeString/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/somatic_single_sample/$dataTypeString/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val somaticSingleSampleInputs = new SomaticSingleSampleInputs(
      workflowTest.runParameters.workflowInputs
    )
    val outputBaseName = somaticSingleSampleInputs.getBaseFileName(workflowName)
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath
    val metricsFileNames = ioUtil
      .listGoogleObjects(truthCloudPath)
      .filter(_.getPath.endsWith("metrics"))
      .map(uriToFilename)
    val validationInputs = SomaticSingleSampleValidationInputs(
      testMetrics = metricsFileNames.map(resultsCloudPath.resolve),
      truthMetrics = metricsFileNames.map(truthCloudPath.resolve),
      testCram = resultsCloudPath.resolve(s"$outputBaseName.cram"),
      testCrai = resultsCloudPath.resolve(s"$outputBaseName.cram.crai"),
      truthCram = truthCloudPath.resolve(s"$outputBaseName.cram"),
      truthCrai = truthCloudPath.resolve(s"$outputBaseName.cram.crai")
    )
    SomaticSingleSampleValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
