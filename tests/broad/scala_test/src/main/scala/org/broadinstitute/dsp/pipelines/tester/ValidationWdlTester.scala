package org.broadinstitute.dsp.pipelines.tester

import java.net.URI
import java.util.UUID

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import akka.stream.scaladsl.{Sink, Source}
import better.files.File
import cromwell.api.model.{
  SubmittedWorkflow,
  Succeeded,
  WorkflowId,
  WorkflowSingleSubmission
}
import io.circe.Json
import io.circe.parser.parse
import io.circe.syntax._
import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRunParameters,
  WorkflowTest
}
import org.broadinstitute.dsp.pipelines.commandline.{
  CromwellEnvironment,
  WorkflowTestCategory
}
import org.broadinstitute.dsp.pipelines.config.BaseConfig
import org.broadinstitute.dsp.pipelines.tester.CromwellWorkflowTester.WarpGitHash

import scala.collection.immutable.Iterable
import scala.concurrent.Future

abstract class ValidationWdlTester(testerConfig: BaseConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends CromwellWorkflowTester {

  protected def validationWorkflowName: String

  protected def buildValidationWdlInputs(
      cloudWorkflowTest: WorkflowTest): String

  protected val resultsPrefix: URI
  protected val truthPrefix: URI
  protected def workflowDir: File
  protected def workflowInputRoot: File =
    workflowDir / "test_inputs" / testerConfig.category.entryName
  protected val testTypeString: String =
    testerConfig.category.entryName.toLowerCase

  protected val timestamp: String =
    testerConfig.useTimestamp.getOrElse(CromwellWorkflowTester.getTimestamp)

  override lazy val wdlContents: String =
    readWdlFromReleaseDir(releaseDir, workflowName)

  protected def env: CromwellEnvironment = testerConfig.env

  protected lazy val releaseDir: File =
    CromwellWorkflowTester.runReleaseWorkflow(localWdlPath, env)

  protected def localValidationWdlPath: File =
    CromwellWorkflowTester.WarpRoot / "verification" / s"$validationWorkflowName.wdl"

  lazy val validationWdlContents: String =
    readWdlFromReleaseDir(validationReleaseDir, validationWorkflowName)

  protected lazy val validationReleaseDir: File =
    CromwellWorkflowTester.runReleaseWorkflow(localValidationWdlPath, env)

  protected def localWdlPath: File = workflowDir / s"$workflowName.wdl"

  protected lazy val inputFileNames: Seq[String] =
    workflowInputRoot.list.toSeq.map(_.name.toString)

  /**
    * If we're not updating the truth data, just validate the runs.
    * Else, use the provided run data as new truth data
    *
    * @param finishedRuns Completed workflow test runs
    * @return A Future of the work
    */
  def validateRunsOrUpdateTruth(
      finishedRuns: Seq[WorkflowTest],
      updateTruth: Boolean,
      testCategory: WorkflowTestCategory
  ): Future[Unit] = {
    if (updateTruth) {
      logger.info("Updating truth data instead of running validation")
      updateTruthData(finishedRuns)
    } else if (testCategory == WorkflowTestCategory.Load) {
      logger.info("Skipping validation for load test")
      Future.successful(())
    } else {
      logger.info("Validating runs against known truth data")
      for {
        _ <- Future.sequence(finishedRuns.map(testerValidation))
        submittedValidationRuns <- submitBatchValidationWorkflows(
          finishedRuns,
          validationWdlContents,
          Option(validationWdlOptions),
          buildValidationWdlInputs,
          zippedImports = dependenciesZipFromReleaseDir(validationReleaseDir,
                                                        validationWorkflowName)
        )
        _ <- awaitBatchCromwellWorkflowCompletion(submittedValidationRuns)
      } yield ()
    }
  }

  protected def testerValidation(finishedRun: WorkflowTest): Future[Unit] = {
    val _ = finishedRun
    Future.successful(())
  }

  protected lazy val googleProject: String = {
    s"broad-gotc-${env.picardEnv}"
  }

  protected lazy val validationWdlOptions: String = Json
    .obj(
      Seq("read_from_cache" -> true.asJson, "write_to_cache" -> true.asJson)
        ++ parse(
          readTestOptions(
            releaseDir,
            env
          )
        ).toOption
          .flatMap(_.asObject)
          .flatMap(_("google_project"))
          .map(project => "google_project" -> project)
          .toSeq
        ++ env.environmentOptions: _*
    )
    .noSpaces

  /**
    * Update the truth data by deleting the old truth data and putting the new run data in its place
    *
    * @param tests the test runs to use as new truth
    * @return A future of the operation
    */
  def updateTruthData(tests: Seq[WorkflowTest]): Future[Unit] = {
    val parallelism = 8
    Source(Iterable(tests: _*))
      .flatMapConcat { test =>
        logger.info(
          s"Replacing truth data at ${test.runParameters.truthCloudPath} with data at ${test.runParameters.resultsCloudPath}"
        )
        // Delete the current truth data
        stagingIoUtil
          .deleteCloudObjects(
            Iterable(
              stagingIoUtil
                .listGoogleObjects(test.runParameters.truthCloudPath): _*
            ).flatMap(
              uri =>
                if (uri.toString.endsWith("logs/"))
                  Iterable(stagingIoUtil.listGoogleObjects(uri): _*)
                else Iterator(uri)
            )
          )
          // Get all the new files
          .flatMapMerge(
            parallelism, { _ =>
              Source(
                Iterable(
                  stagingIoUtil
                    .listGoogleObjects(test.runParameters.resultsCloudPath): _*
                ).flatMap(
                  uri =>
                    if (uri.toString.endsWith("logs/"))
                      Iterable(stagingIoUtil.listGoogleObjects(uri): _*)
                    else Iterator(uri)
                )
              )
              // Shove all the new files into the truth location
                .flatMapMerge(
                  parallelism, { result =>
                    Source.single(
                      stagingIoUtil.copyGoogleObject(
                        result,
                        test.runParameters.truthCloudPath
                          .resolve(uriToFilename(result))
                      )
                    )
                  }
                )
            }
          )
      }
      .runWith(Sink.ignore)
      .map(_ => ())
  }

  protected def uriToFilename(uri: URI): String = {
    val path = uri.getPath
    path.substring(path.lastIndexOf('/') + 1)
  }

  def generateRunParameters: Seq[WorkflowRunParameters] = {
    logger.info(s"workflowInputRoot: $workflowInputRoot")
    workflowInputRoot.list.toSeq.map(_.name.toString).map { fileName =>
      val inputsName = fileName.replace(".json", "")
      val resultsPath =
        resultsPrefix.resolve(s"$inputsName/")
      val truthPath = truthPrefix.resolve(s"$inputsName/")

      WorkflowRunParameters(
        id = s"${envString}_$inputsName",
        workflowInputs = getInputContents(fileName),
        resultsCloudPath = resultsPath,
        truthCloudPath = truthPath
      )
    }
  }

  def getInputContents(fileName: String): String =
    (workflowInputRoot / fileName).contentAsString

  override def runTest: Future[Unit] = {
    logger.info(
      s"Running the $workflowName workflow using ${testerConfig.category} data"
    )
    val samples = generateRunParameters

    if (testerConfig.useTimestamp.isDefined) {
      usePreviousRun(samples)
    } else {
      runFullTest(samples)
    }
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
      _ <- copyBatchCromwellWorkflowResults(finishedSamples)
      _ <- validateRunsOrUpdateTruth(
        finishedSamples,
        testerConfig.updateTruth,
        testerConfig.category
      )

    } yield ()
  }

  /**
    * Run only the validation portion of the tests. This is triggered by providing a "validation" timestamp input
    * as a command line argument.
    */
  private def usePreviousRun(
      samples: Seq[WorkflowRunParameters]): Future[Unit] = {
    testerConfig.useTimestamp.foreach { timestamp =>
      logger.info(
        s"Only running validation workflows with timestamp $timestamp"
      )
    }
    val validationRuns = samples.map(
      run =>
        WorkflowTest(
          runParameters = run,
          workflow = SubmittedWorkflow(
            WorkflowId(UUID.randomUUID()),
            env.cromwellUrl,
            WorkflowSingleSubmission(None,
                                     None,
                                     None,
                                     None,
                                     None,
                                     None,
                                     None,
                                     None,
                                     None)
          ),
          workflowStatus = Succeeded
      )
    )
    validateRunsOrUpdateTruth(
      validationRuns,
      testerConfig.updateTruth,
      testerConfig.category
    )
  }

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
