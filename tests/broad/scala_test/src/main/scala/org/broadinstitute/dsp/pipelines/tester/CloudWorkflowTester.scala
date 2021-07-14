package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.config._

class CloudWorkflowTester(testerConfig: CloudWorkflowConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override protected def workflowName: String =
    testerConfig.pipeline.workflowName

  override protected def workflowDir: File =
    File(
      CromwellWorkflowTester.PipelineRoot + testerConfig.pipeline.workflowDir)

  protected val pipeline: String = testerConfig.pipeline.getClass.getName

  override protected lazy val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/$pipeline/$testTypeString/$timestamp/"
    )

  override protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/$pipeline/truth/$testTypeString/${testerConfig.truthBranch}/"
    )

  //override protected def buildValidationWdlInputs(
  //    workflowTest: WorkflowTest
  //): String = {
  //  val somaticSingleSampleInputs = new SomaticSingleSampleInputs(
  //    workflowTest.runParameters.workflowInputs
  //  )
  //  val outputBaseName = somaticSingleSampleInputs.getBaseFileName(workflowName)
  //  val resultsCloudPath =
  //    workflowTest.runParameters.resultsCloudPath
  //  val truthCloudPath = workflowTest.runParameters.truthCloudPath
  //  val metricsFileNames = ioUtil
  //    .listGoogleObjects(truthCloudPath)
  //    .filter(_.getPath.endsWith("metrics"))
  //    .map(uriToFilename)
  //  val validationInputs = SomaticSingleSampleValidationInputs(
  //    testMetrics = metricsFileNames.map(resultsCloudPath.resolve),
  //    truthMetrics = metricsFileNames.map(truthCloudPath.resolve),
  //    testCram = resultsCloudPath.resolve(s"$outputBaseName.cram"),
  //    testCrai = resultsCloudPath.resolve(s"$outputBaseName.cram.crai"),
  //    truthCram = truthCloudPath.resolve(s"$outputBaseName.cram"),
  //    truthCrai = truthCloudPath.resolve(s"$outputBaseName.cram.crai")
  //  )
  //  SomaticSingleSampleValidationInputs
  //    .marshall(validationInputs)
  //    .printWith(implicitly)
  //}
}
