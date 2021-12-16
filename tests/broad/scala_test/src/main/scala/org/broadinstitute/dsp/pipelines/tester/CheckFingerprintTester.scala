package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  CheckFingerprintInputs,
  CheckFingerprintValidationInputs
}

class CheckFingerprintTester(testerConfig: CheckFingerprintConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override val workflowName: String = "CheckFingerprint"

  val workflowDir
    : File = CromwellWorkflowTester.PipelineRoot / "broad" / "infrastructure"

  protected val vaultTokenPath: String =
    s"gs://broad-dsp-gotc-arrays-$envString-tokens/arrayswdl.token"

  override protected val validationWorkflowName: String =
    "VerifyCheckFingerprint"

  // Validation uses the same options as the arrays workflow
  override protected lazy val validationWdlOptions: String = {
    readTestOptions(releaseDir, env)
  }

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/infrastructure/check_fingerprint/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/infrastructure/check_fingerprint/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override def getInputContents(fileName: String): String = {
    super
      .getInputContents(fileName)
      .lines
      .map(
        _.replace("{ENV}", envString)
          .replace("{VAULT_TOKEN_PATH}", vaultTokenPath)
      )
      .mkString
  }

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val checkFingerprintInputs = new CheckFingerprintInputs(
      workflowTest.runParameters.workflowInputs
    )
    val fileNameSafeSampleAlias =
      checkFingerprintInputs.getSampleAlias(workflowName).replaceAll(" ", "_")

    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath

    val metricsFileNames = ioUtil
      .listGoogleObjects(truthCloudPath)
      .filter(_.getPath.endsWith("metrics"))
      .map(uriToFilename)

    val validationInputs = CheckFingerprintValidationInputs(
      test_metrics = metricsFileNames.map(resultsCloudPath.resolve),
      truth_metrics = metricsFileNames.map(truthCloudPath.resolve),
      test_fingerprint_vcf = resultsCloudPath.resolve(
        s"$fileNameSafeSampleAlias.reference.fingerprint.vcf.gz"),
      truth_fingerprint_vcf = truthCloudPath.resolve(
        s"$fileNameSafeSampleAlias.reference.fingerprint.vcf.gz"),
    )
    CheckFingerprintValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
