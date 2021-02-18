package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import io.circe.parser.parse
import io.circe.syntax._
import io.circe.Json
import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRunParameters,
  WorkflowTest
}
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.ReprocessingInputs
import org.broadinstitute.dsp.pipelines.tester.CromwellWorkflowTester.WarpGitHash

import scala.concurrent.Future

class ExternalReprocessingTester(testerConfig: GermlineCloudWorkflowConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ReprocessingTester(testerConfig) {

  override val workflowName: String = s"External${dataTypePrefix}Reprocessing"

  override lazy val workflowDir: File =
    CromwellWorkflowTester.PipelineRoot / "broad" / "reprocessing" / "external" / dataTypeString

  override protected lazy val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/external_reprocessing/$dataTypeString/$testTypeString/$timestamp/"
    )
  override protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/external_reprocessing/$dataTypeString/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  protected val destinationCloudPath: String =
    s"gs://broad-gotc-$envString/external_reprocessing/$dataTypeString/$testTypeString/$timestamp/copied/"

  protected val googleAccountVaultPath: String =
    s"secret/dsde/gotc/$envString/picard/picard-account.pem"

  protected val vaultTokenPath: String =
    s"gs://broad-dsp-gotc-$envString-tokens/picardsa.token"

  override def generateRunParameters: Seq[WorkflowRunParameters] = {
    super.generateRunParameters.map(
      rp =>
        rp.copy(
          workflowInputs = rp.workflowInputs
            .replace("{TRUTH_BRANCH}", testerConfig.truthBranch)
            .replace("{DESTINATION_CLOUD_PATH}", destinationCloudPath)
            .replace("{GOOGLE_ACCOUNT_VAULT_PATH}", googleAccountVaultPath)
            .replace("{VAULT_TOKEN_PATH}", vaultTokenPath)
      )
    )
  }

  // Note - we are explicitly setting the google_project here so that when running in a non-dev environment,
  // The workflow can still access the test data AND can then read from the vault
  override def readTestOptions(
      releaseDir: File,
      environment: CromwellEnvironment
  ): String = {
    val defaultOptions = Array(
      "read_from_cache" -> testerConfig.useCallCaching.asJson,
      "backend" -> testerConfig.papiVersion.entryName.asJson,
      "monitoring_script" -> "gs://broad-gotc-test-storage/cromwell_monitoring_script.sh".asJson,
      "google_project" -> "broad-exomes-dev1".asJson
    )

    val optionsJson = defaultOptions ++ environment.environmentOptions

    parse(
      (releaseDir / workflowName / s"${workflowName}_${WarpGitHash}.options.json").contentAsString)
      .fold(
        e => throw new RuntimeException("Could not create options json", e),
        _.deepMerge(Json.obj(optionsJson: _*)).noSpaces
      )
  }

  override protected def testerValidation(
      finishedRun: WorkflowTest): Future[Unit] = {
    val uriCloudPath = URI.create(destinationCloudPath)

    val reprocessingInputs = new ReprocessingInputs(
      finishedRun.runParameters.workflowInputs
    )
    val outputBaseName = reprocessingInputs.getBaseFileName(workflowName)

    val suffixes =
      Seq(".cram",
          ".cram.crai",
          ".g.vcf.gz",
          ".g.vcf.gz.tbi",
          ".preBqsr.selfSM")
    val files =
      suffixes.map(suffix => uriCloudPath.resolve(s"$outputBaseName$suffix"))

    val allFilesCopied = files.map(ioUtil.isGoogleObject).reduce(_ && _)

    if (allFilesCopied) {
      Future.successful(())
    } else {
      val msg =
        s"Workflow ${finishedRun.workflow.id} did not succeed in copying all files to the final directory. Please make sure ${files
          .mkString(", ")} are being exist in the output directory"
      Future.failed(new CromwellWorkflowTester.TestFailedException(msg))
    }
  }
}
