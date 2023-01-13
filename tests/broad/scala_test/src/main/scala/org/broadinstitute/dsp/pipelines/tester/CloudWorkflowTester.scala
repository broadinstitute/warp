package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer

import better.files.File
import io.circe.Json
import io.circe.parser.parse
import io.circe.syntax._
import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRunParameters,
  WorkflowTest
}
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment
import org.broadinstitute.dsp.pipelines.config.CloudWorkflowConfig
import org.broadinstitute.dsp.pipelines.tester.CromwellWorkflowTester.WarpGitHash

import scala.concurrent.Future
import scala.util.matching.Regex

class CloudWorkflowTester(testerConfig: CloudWorkflowConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends CromwellWorkflowTester {

  // Name of the workflow that we are running
  // i.e. TestExomeGermlineSingleSample
  override def workflowName: String =
    testerConfig.pipeline.workflowName

  // Name of the pipeline that is being tested
  // i.e. ExomeGermlineSingleSample
  protected val pipeline: String =
    testerConfig.pipeline.pipelineName

  // Directory in WARP where the main workflow lives
  // i.e. pipelines/broad/dna_seq/germline/single_sample/exome/
  protected def workflowDir: File =
    File(
      CromwellWorkflowTester.PipelineRoot + testerConfig.pipeline.workflowDir)

  // Bundle everything up into a single WDL
  // Grab the test workflow that we are running from verification/test-wdls
  protected lazy val releaseDir: File =
    CromwellWorkflowTester.runReleaseWorkflow(
      CromwellWorkflowTester.TestsRoot / s"$workflowName.wdl",
      env
    )

  // Bucket to copy test results to
  protected lazy val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/$pipeline/$testTypeString/$timestamp/"
    )

  // Bucket to copy results to if we are updating truth
  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/$pipeline/truth/$testTypeString/${testerConfig.truthBranch}/"
    )

  // Location for the test inputs in WARP
  // i.e. /Plumbing or /Scientific
  protected def workflowInputRoot: File =
    workflowDir / "test_inputs" / testerConfig.category.entryName

  // All of our plumbing or scientific test inputs
  protected lazy val inputFileNames: Seq[String] =
    workflowInputRoot.list
      .filter(_.name.endsWith(".json"))
      .toSeq
      .map(_.name.toString)

  // plumbing or scientific
  protected val testTypeString: String =
    testerConfig.category.entryName.toLowerCase

  // Environment settings
  protected def env: CromwellEnvironment = testerConfig.env

  // Updating truth
  protected lazy val updateTruth: Boolean = testerConfig.updateTruth

  // Timestamp for results bucket
  protected val timestamp: String =
    testerConfig.useTimestamp.getOrElse(CromwellWorkflowTester.getTimestamp)

  // Get string contents of test wdl
  override lazy val wdlContents: String =
    readWdlFromReleaseDir(releaseDir, workflowName)

  // ExternalWholeGenomeReprocessing destination
  protected val destinationCloudPathWGS: String =
    s"gs://broad-gotc-$envString/external_reprocessing/WGS/$testTypeString/$timestamp/copied/"

  // ExternalExomeReprocessing destination
  protected val destinationCloudPathExome: String =
    s"gs://broad-gotc-$envString/external_reprocessing/Exome/$testTypeString/$timestamp/copied/"

  // Final destination for external reprocessing copy step
  val destinationCloudPath =
    if (pipeline.contains("Exome")) destinationCloudPathExome
    else destinationCloudPathWGS

  // Values below are needed so we can auth as picard service account to copy output files
  protected val vaultTokenPath: String =
    s"gs://broad-dsp-gotc-$envString-tokens/picardsa.token"
  protected val googleAccountVaultPath: String =
    s"secret/dsde/gotc/$envString/picard/picard-account.pem"
  // Specific to arrays workflows to allow Arrays authenticate with Vault
  protected val vaultTokenPathArrays: String =
    s"gs://broad-dsp-gotc-arrays-$envString-tokens/arrayswdl.token"

  // Always run in broad-exomes-dev1 google project
  protected lazy val googleProject: String = {
    s"broad-exomes-dev1"
  }

  protected def testerValidation(finishedRun: WorkflowTest): Future[Unit] = {
    val _ = finishedRun
    Future.successful(())
  }

  protected def uriToFilename(uri: URI): String = {
    val path = uri.getPath
    path.substring(path.lastIndexOf('/') + 1)
  }

  /**
    * Generate the run parameters for each testing sample
    */
  def generateRunParameters: Seq[WorkflowRunParameters] = {
    inputFileNames.map { fileName =>
      val inputsName = fileName.replace(".json", "")
      val resultsPath = resultsPrefix.resolve(s"$inputsName/")
      val truthPath = truthPrefix.resolve(s"$inputsName/")

      logger.info(s"Generating WDL inputs for -> $fileName")

      WorkflowRunParameters(
        id = s"${envString}_$inputsName",
        workflowInputs = getInputContents(fileName, resultsPath, truthPath),
        resultsCloudPath = resultsPath,
        truthCloudPath = truthPath
      )
    }
  }

  /**
    * Format the inputs for the test WDL
    */
  def getInputContents(fileName: String,
                       resultsPath: URI,
                       truthPath: URI): String = {
    var defaultInputs = Array(
      workflowName + ".truth_path" -> truthPath.asJson,
      workflowName + ".results_path" -> resultsPath.asJson,
      workflowName + ".update_truth" -> updateTruth.asJson,
      workflowName + ".vault_token_path" -> vaultTokenPath.asJson,
      workflowName + ".google_account_vault_path" -> googleAccountVaultPath.asJson,
    )

    val tokenPipelines =
      List("Arrays",
           "BroadInternalRNAWithUMIs",
           "BroadInternalUltimaGenomics",
           "CheckFingerprint")
    val externalPipelines =
      List("ExternalWholeGenomeReprocessing", "ExternalExomeReprocessing")

    // Only add the vault token for the pipelines above
    // These allow these pipelines to auth with vault/gcp
    if (tokenPipelines.contains(pipeline)) {
      defaultInputs = defaultInputs :+ workflowName + ".vault_token_path_arrays" -> vaultTokenPathArrays.asJson
      defaultInputs = defaultInputs :+ workflowName + ".environment" -> envString.asJson
    }

    // Only add the destinationCloudPath for the external pipelines
    if (externalPipelines.contains(pipeline)) {
      defaultInputs = defaultInputs :+ workflowName + ".destination_cloud_path" -> destinationCloudPath.asJson
    }

    /**
      * If we have nested inputs in our test inputs we need to push them down a level
      * e.g.
      * ExomeGermlineSingleSample.AggregatedBamQC.CollectReadgroupBamQualityMetrics.collect_gc_bias_metrics ->
      * TestExomeGermlineSingleSample.ExomeGermlineSingleSample.AggregatedBamQC.CollectReadgroupBamQualityMetrics.collect_gc_bias_metrics
      */
    val pattern = new Regex(s"($workflowName).([A-Z]\\w+).")

    /** Find any instance of the pipeline followed by . and replace with wrapper workflow
      * e.g.
      * Arrays. -> TestArrays.
      *
      * This handles the case where the wrapper workflow is a substring of a nested input (CheckFingerprint CheckFingerprintTask)
      */
    var inputsString = (workflowInputRoot / fileName).contentAsString
      .replace(s"$pipeline.", s"$workflowName.")

    // Replace all of the injected {} value in the test inputs file
    inputsString = inputsString
      .replaceAll("\\{TRUTH_BRANCH}", testerConfig.truthBranch)
      .replaceAll("\\DESTINATION_CLOUD_PATH", destinationCloudPath)
      .replaceAll("\\VAULT_TOKEN_PATH", vaultTokenPath)
      .replaceAll("\\GOOGLE_ACCOUNT_VAULT_PATH", googleAccountVaultPath)

    // If wrapper workflow name is follow by [A-Z] then we know its a nested input
    inputsString = pattern.replaceAllIn(
      inputsString,
      m => s"$workflowName.$pipeline." + m.group(2) + ".")

    parse(inputsString).fold(
      e => throw new RuntimeException("Could not create inputs json", e),
      _.deepMerge(Json.obj(defaultInputs: _*)).noSpaces
    )
  }

  /**
    * Collect all the samples to run
    */
  override def runTest: Future[Unit] = {
    logger.info(
      s"Running the $workflowName workflow using ${testerConfig.category} data"
    )
    val samples = generateRunParameters
    runFullTest(samples)
  }

  /**
    * Run the full tests. This is the normal mode of operation.
    */
  private def runFullTest(samples: Seq[WorkflowRunParameters]): Future[Unit] = {
    for {
      submittedSamples <- submitBatchWorkflows(
        samples,
        Some(readTestOptions(releaseDir, env)),
        dependenciesZipFromReleaseDir(releaseDir, workflowName)
      )
      finishedSamples <- awaitBatchCromwellWorkflowCompletion(submittedSamples)
    } yield ()
  }

  /**
    * Generate options.json for test sample
    */
  def readTestOptions(releaseDir: File,
                      environment: CromwellEnvironment): String = {
    val defaultOptions = Array(
      "read_from_cache" -> testerConfig.useCallCaching.asJson,
      "backend" -> testerConfig.papiVersion.entryName.asJson,
      "monitoring_script" -> "gs://broad-gotc-test-storage/cromwell_monitoring_script.sh".asJson,
      "google_project" -> googleProject.asJson
    )

    val optionsJson = defaultOptions ++ environment.environmentOptions

    parse(
      (releaseDir / workflowName / s"${workflowName}_${WarpGitHash}.options.json").contentAsString)
      .fold(
        e => throw new RuntimeException("Could not create options json", e),
        _.deepMerge(Json.obj(optionsJson: _*)).noSpaces
      )
  }

}
