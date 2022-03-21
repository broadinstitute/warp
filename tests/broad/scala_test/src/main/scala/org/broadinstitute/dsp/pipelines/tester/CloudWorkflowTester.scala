package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import akka.stream.scaladsl.{Sink, Source}

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

import scala.collection.immutable.Iterable
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

  // Directory in WARP where the test workflow lives
  // i.e. pipelines/broad/dna_seq/germline/single_sample/exome/TestExomeGermlineSingleSample.wdl
  protected def workflowDir: File =
    File(
      CromwellWorkflowTester.PipelineRoot + testerConfig.pipeline.workflowDir)

  // Bundle everything up into a single WDL
  protected lazy val releaseDir: File =
    CromwellWorkflowTester.runReleaseWorkflow(
      workflowDir / s"$workflowName.wdl",
      env
    )

  // Store the results
  protected lazy val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/$pipeline/$testTypeString/$timestamp/"
    )

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
    workflowInputRoot.list.toSeq.map(_.name.toString)

  protected val testTypeString: String =
    testerConfig.category.entryName.toLowerCase

  protected val timestamp: String =
    testerConfig.useTimestamp.getOrElse(CromwellWorkflowTester.getTimestamp)

  override lazy val wdlContents: String =
    readWdlFromReleaseDir(releaseDir, workflowName)

  protected def env: CromwellEnvironment = testerConfig.env

  protected lazy val updateTruth: Boolean = testerConfig.updateTruth

  protected lazy val useTimestamp: Option[String] = testerConfig.useTimestamp

  /**
    * If we're not updating the truth data, just validate the runs.
    * Else, use the provided run data as new truth data
    *
    * @param finishedRuns Completed workflow test runs
    * @return A Future of the work
    */
  // TODO: rename and refactor
  /*
  def validateRunsOrUpdateTruth(
      finishedRuns: Seq[WorkflowTest],
      updateTruth: Boolean,
      testCategory: WorkflowTestCategory
  ): Future[Unit] = {
    if (updateTruth) {
      logger.info("Updating truth data instead of running validation")
      updateTruthData(finishedRuns)
    }
  }
   */

  protected def testerValidation(finishedRun: WorkflowTest): Future[Unit] = {
    val _ = finishedRun
    Future.successful(())
  }

  protected lazy val googleProject: String = {
    s"broad-gotc-${env.picardEnv}"
  }

  /*
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

   */

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
      //val metricsFileNames = ioUtil
      //.listGoogleObjects(truthPath)
      //.filter(_.getPath.endsWith("metrics"))
      //.map(uriToFilename)
      val stupidInputs = getInputContents(fileName, resultsPath, truthPath)
      logger.info(stupidInputs)

      WorkflowRunParameters(
        id = s"${envString}_$inputsName",
        workflowInputs = getInputContents(fileName, resultsPath, truthPath),
        resultsCloudPath = resultsPath,
        truthCloudPath = truthPath
      )
    }
  }

  def getInputContents(fileName: String,
                       resultsPath: URI,
                       truthPath: URI): String = {
    val defaultInputs = Array(
      workflowName + ".truth_path" -> truthPath.asJson,
      workflowName + ".results_path" -> resultsPath.asJson,
      workflowName + ".update_truth" -> updateTruth.asJson,
      workflowName + ".use_timestamp" -> useTimestamp.asJson,
      workflowName + ".timestamp" -> timestamp.asJson,
      workflowName + ".cromwell_env" -> envString.asJson
    )

    /**
      * If we have nested inputs in our test inputs we need to push them down a level
      * e.g.
      * ExomeGermlineSingleSample.AggregatedBamQC.CollectReadgroupBamQualityMetrics.collect_gc_bias_metrics ->
      * TestExomeGermlineSingleSample.ExomeGermlineSingleSample.AggregatedBamQC.CollectReadgroupBamQualityMetrics.collect_gc_bias_metrics
      */
    val pattern = new Regex("(TestExomeGermlineSingleSample).([A-Z]\\w+).")

    var inputsString = (workflowInputRoot / fileName).contentAsString
      .replace(pipeline, workflowName)

    inputsString = pattern.replaceAllIn(
      inputsString,
      m => s"$workflowName.$pipeline." + m.group(2) + ".")

    parse(inputsString).fold(
      e => throw new RuntimeException("Could not create inputs json", e),
      _.deepMerge(Json.obj(defaultInputs: _*)).noSpaces
    )
  }

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
      _ <- copyBatchCromwellWorkflowResults(finishedSamples)
    } yield ()
  }

// TODO: either fix this so it will work, or disable it. Question: Do we need / want this functionality

  /**
    * Run only the validation portion of the tests. This is triggered by providing a "validation" timestamp input
    * as a command line argument.
    */
  /*
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
