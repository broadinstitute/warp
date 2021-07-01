package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  OptimusInputs,
  OptimusValidationInputs
}

class OptimusTester(testerConfig: OptimusConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override val workflowName: String = s"Optimus"

  val workflowDir
    : File = CromwellWorkflowTester.WarpRoot / "pipelines" / "skylab" / "optimus"

  override protected val validationWorkflowName: String =
    "VerifyOptimus"

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/optimus/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/optimus/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val optimusInputs =
      new OptimusInputs(
        workflowTest.runParameters.workflowInputs
      )
    val outputBaseName =
      optimusInputs.getBaseFileName(workflowName)
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath
    val validationInputs = OptimusValidationInputs(
      testBam = resultsCloudPath.resolve(s"$outputBaseName.bam"),
      truthBam = truthCloudPath.resolve(s"$outputBaseName.bam")
    )
    OptimusValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
