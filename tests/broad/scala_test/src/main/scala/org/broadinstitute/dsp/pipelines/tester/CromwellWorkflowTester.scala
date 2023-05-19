package org.broadinstitute.dsp.pipelines.tester

import java.net.URI
import java.nio.ByteBuffer
import java.text.SimpleDateFormat
import java.util.Date

import akka.actor.{ActorSystem, Scheduler}
import akka.http.scaladsl.model.headers.OAuth2BearerToken
import akka.http.scaladsl.unmarshalling.Unmarshaller.UnsupportedContentTypeException
import akka.pattern.after
import akka.stream.ActorMaterializer
import akka.stream.scaladsl.{Sink, Source}
import better.files.File
import com.typesafe.scalalogging.StrictLogging
import cromwell.api.CromwellClient
import cromwell.api.model._
import de.heikoseeberger.akkahttpcirce.FailFastCirceSupport
import io.circe.parser._
import org.apache.http.client.fluent.Request
import org.broadinstitute.clio.client.util.IoUtil
import org.broadinstitute.clio.util.auth.ClioCredentials
import org.broadinstitute.clio.util.json.ModelAutoDerivation
import org.broadinstitute.dsp.pipelines.util.DataType
import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRun,
  WorkflowRunParameters,
  WorkflowTest,
  WorkflowValidation
}
import org.broadinstitute.dsp.pipelines.commandline.{
  Config,
  CromwellEnvironment,
  WorkflowTestType
}
import org.broadinstitute.dsp.pipelines.tester.CromwellWorkflowTester.WarpGitHash
import org.broadinstitute.dsp.pipelines.util.{CromwellUtils, VaultUtil}

import scala.collection.immutable.Iterable
import scala.collection.JavaConverters._
import scala.concurrent.duration._
import scala.concurrent.{ExecutionContext, Future}
import scala.sys.process.Process
import scala.util.{Failure, Success, Try}

object CromwellWorkflowTester {
  class TestFailedException(msg: String) extends RuntimeException(msg)

  lazy val WarpRoot: File =
    File(Process(Seq("git", "rev-parse", "--show-toplevel")).!!.trim)

  lazy val WarpGitHash: String =
    Process(Seq("git", "rev-parse", "--short", "HEAD")).!!.trim

  lazy val GotcRoot: File = WarpRoot / "genomes_in_the_cloud"

  lazy val PipelineRoot: File = WarpRoot / "pipelines"

  lazy val TestsRoot: File = WarpRoot / "verification" / "test-wdls"

  def apply(config: Config)(
      implicit system: ActorSystem,
      mat: ActorMaterializer
  ): CromwellWorkflowTester = {
    import WorkflowTestType._

    config.test match {
      case AllOfUs =>
        new AllOfUsTester(config.germlineCloudConfig)
      case AnnotationFiltration =>
        new AnnotationFiltrationTester(config.annotationFiltrationConfig)
      case BroadInternalRNAWithUMIs =>
        new BroadInternalRNAWithUMIsTester(
          config.broadInternalRNAWithUMIsConfig)
      case BroadInternalUltimaGenomics =>
        new BroadInternalUltimaGenomicsTester(
          config.broadInternalUltimaGenomicsConfig)
      case CheckFingerprint =>
        new CheckFingerprintTester(config.checkFingerprintConfig)
      case CramToUnmappedBams =>
        new CramToUnmappedBamsTester(config.cramToUnmappedBamsConfig)
      case CloudWorkflow =>
        new CloudWorkflowTester(config.cloudWorkflowConfig)
      case Dummy => new DummyTester()
      case ExternalReprocessing =>
        new ExternalReprocessingTester(config.germlineCloudConfig)
      case GDCWholeGenomeSomaticSingleSample =>
        new GDCWholeGenomeSomaticSingleSampleTester(
          config.gdcWholeGenomeSomaticSingleSampleConfig)
      case GenotypeConcordance =>
        new GenotypeConcordanceTester(config.genotypeConcordanceConfig)
      case GermlineSingleSample =>
        new GermlineSingleSampleTester(config.germlineCloudConfig)
      case ReblockGvcf =>
        new ReblockGvcfTester(config.germlineCloudConfig)
      case Reprocessing =>
        new ReprocessingTester(config.germlineCloudConfig)
      case JointGenotyping =>
        new JointGenotypingTester(config.germlineCloudConfig)
      case ValidateChip =>
        new ValidateChipTester(config.validateChipConfig)
      case Arrays =>
        new ArraysTester(config.arraysConfig)
      case IlluminaGenotypingArray =>
        new IlluminaGenotypingArrayTester(config.illuminaGenotypingArrayConfig)
      case Imputation =>
        new ImputationTester(config.imputationConfig)
      case RNAWithUMIs =>
        new RNAWithUMIsTester(config.rnaWithUMIsConfig)
      case SomaticSingleSample =>
        new SomaticSingleSampleTester(config.somaticCloudWorkflowConfig)
      case UltimaGenomicsJointGenotyping =>
        new UltimaGenomicsJointGenotypingTester(
          config.ultimaGenomicsJointGenotypingConfig)
      case UltimaGenomicsWholeGenomeGermline =>
        new UltimaGenomicsWholeGenomeGermlineTester(
          config.ultimaGenomicsWholeGenomeGermlineConfig)
      case VariantCalling =>
        new VariantCallingTester(config.germlineCloudConfig)
    }
  }

  private lazy val pipelinesReleaseScript: File =
    WarpRoot / "scripts" / "build_pipeline_release.sh"

  def runReleaseWorkflow(wdlPath: File, env: CromwellEnvironment): File = {
    val tmp = File.newTemporaryDirectory().deleteOnExit()
    Process(
      pipelinesReleaseScript.pathAsString,
      Seq(
        "-w",
        wdlPath.pathAsString,
        "-o",
        tmp.pathAsString,
        "-v",
        WarpGitHash,
        "-e",
        env.picardEnv
      )
    ).!!
    tmp
  }

  def getTimestamp: String =
    new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss").format(new Date())
}

abstract class CromwellWorkflowTester(
    implicit am: ActorMaterializer,
    as: ActorSystem,
) extends StrictLogging
    with FailFastCirceSupport
    with ModelAutoDerivation
    with CromwellUtils {

  implicit val ec: ExecutionContext = as.dispatcher
  implicit val scheduler: Scheduler = as.scheduler

  protected def env: CromwellEnvironment

  protected lazy val envString: String = env.picardEnv

  def workflowName: String

  def wdlContents: String

  def runTest: Future[Unit]

  protected def cromwellClient: () => CromwellClient = () => getCromwellClient

  protected lazy val vaultToken: String = sys.env
    .get("VAULT_TOKEN_PATH")
    .map(File(_))
    .orElse(Some(File.home / ".vault-token").filter(_.exists))
    .fold(throw new RuntimeException("Could not get vault token"))(
      _.contentAsString.trim)

  protected lazy val serviceAccountJson: File = {
    val json = VaultUtil.getPicardServiceAccount(env, vaultToken)
    File.newTemporaryFile().deleteOnExit().write(json.spaces2)
  }

  protected lazy val googleCredentials: ClioCredentials = new ClioCredentials(
    Some(serviceAccountJson)
  )

  protected lazy val ioUtil: IoUtil = IoUtil(googleCredentials)

  // Used for truth data copying
  protected lazy val stagingIoUtil: IoUtil = env match {
    case CromwellEnvironment.Staging => ioUtil
    case _ =>
      val json =
        VaultUtil.getPicardServiceAccount(CromwellEnvironment.Staging,
                                          vaultToken)
      val jsonFile = File.newTemporaryFile().deleteOnExit().write(json.spaces2)
      val stagingGoogleCredentials: ClioCredentials = new ClioCredentials(
        Some(jsonFile))
      IoUtil(stagingGoogleCredentials)
  }

  /**
    * Start a single workflow in Cromwell using the contents of the files provided
    *
    * @param wdlPath        Path to the WDL
    * @param wdlInputsPath  Path to the inputs for this submission
    * @param wdlOptionsPath Optional WDL options for this submission
    * @return A Future containing the submitted workflow
    */
  protected def runWorkflowFromFiles(
      wdlPath: File,
      wdlInputsPath: File,
      wdlOptionsPath: File
  ): Future[SubmittedWorkflow] = {
    val wdl = wdlPath.contentAsString
    val wdlInputs = wdlInputsPath.contentAsString
    val wdlOptions = wdlOptionsPath.contentAsString
    runWorkflow(wdl, wdlInputs, Option(wdlOptions))
  }

  /**
    * Start a single workflow in Cromwell
    *
    * @param wdl        The wdl to run
    * @param wdlInputs  Inputs to the wdl
    * @param wdlOptions Optional wdl options
    * @return A Future containing the submitted workflow
    */
  protected def runWorkflow(
      wdl: String,
      wdlInputs: String,
      wdlOptions: Option[String],
      wdlLabels: Option[List[Label]] = None,
      zippedImports: Option[File] = None,
      workflowName: String = this.workflowName
  ): Future[SubmittedWorkflow] = {
    val wss =
      prepareSingleSubmission(wdl,
                              wdlInputs,
                              wdlOptions,
                              wdlLabels,
                              zippedImports)
    submitWorkflow(wss, workflowName)
  }

  /**
    * Prepare a single submission, but don't submit it to Cromwell
    *
    * @param wdl           The WDL
    * @param wdlInputs     The inputs to the WDL
    * @param wdlOptions    Optional options to the WDL
    * @param wdlLabels a list of WDL labels to label a job with
    * @param zippedImports Optional imports for subworkflows
    * @return A WorkflowSingleSubmission ready to be sent to Cromwell
    */
  protected def prepareSingleSubmission(
      wdl: String,
      wdlInputs: String,
      wdlOptions: Option[String] = None,
      wdlLabels: Option[List[Label]] = None,
      zippedImports: Option[File] = None
  ): WorkflowSingleSubmission = {
    WorkflowSingleSubmission(
      workflowSource = Option(wdl),
      workflowType = Option("wdl"),
      workflowRoot = None,
      workflowUrl = None,
      workflowTypeVersion = None,
      inputsJson = Option(wdlInputs),
      options = wdlOptions,
      labels = wdlLabels,
      zippedImports = zippedImports
    )
  }

  /**
    * Submit a single workflow to Cromwell. Adds a shutdown hook to abort the workflow
    * if the program is stopped
    *
    * @param wss The WorkflowSingleSubmission to submit
    * @return A Future containing the submitted workflow
    */
  protected def submitWorkflow(
      wss: WorkflowSingleSubmission,
      workflowName: String
  ): Future[SubmittedWorkflow] = {
    logger.info(s"Submitting $workflowName job to Cromwell")
    cromwellClient().submit(wss).value.unsafeToFuture().map {
      case Left(value) =>
        throw new RuntimeException(s"an error occurred: " + value)
      case Right(submitted) => {
        addAbortTerminationHook(submitted)
        submitted
      }
    }
  }

  protected def retry[A](
      op: => Future[A],
      endWhen: A => Boolean,
      d: FiniteDuration
  ): Future[A] = op.flatMap { a =>
    if (endWhen(a)) Future(a) else after(d, scheduler)(retry(op, endWhen, d))
  }

  /**
    * Wait for a cromwell workflow to complete
    *
    * @param submission The submission to wait for
    * @return A successful or failed Future which represents a successful or failed workflow
    */
  protected def awaitCromwellWorkflowCompletion(
      submission: SubmittedWorkflow
  ): Future[Unit] = {
    awaitBatchCromwellWorkflowCompletion(
      collection.immutable.Seq(
        WorkflowTest.fromSubmittedWorkflow(submission, Submitted))
    ).map(_ => ())
  }

  /**
    * Submit a batch of workflows
    *
    * @param workflowRuns Runs to submit
    * @param workflowOptions Optional options for the workflow
    * @param zippedImports Optional zipped dependencies for the workflow
    * @return A Future of submitted workflow tests.
    */
  protected def submitBatchWorkflows(
      workflowRuns: Seq[WorkflowRunParameters],
      workflowOptions: Option[String],
      zippedImports: Option[File] = None
  ): Future[Seq[WorkflowTest]] = {
    Source(Iterable(workflowRuns: _*))
      .flatMapMerge(
        8, { run =>
          Source
            .fromFuture(
              runWorkflow(
                wdlContents,
                run.workflowInputs,
                workflowOptions,
                Some(List(Label("test-run-id", run.id))),
                zippedImports
              )
            )
            .map { submittedWorkflow =>
              WorkflowTest(
                runParameters = run,
                workflow = submittedWorkflow,
                workflowStatus = Submitted
              )
            }
        }
      )
      .runWith(Sink.collection)
  }

  protected def submitBatchValidationWorkflows(
      workflowTests: Seq[WorkflowTest],
      validationWorkflowContents: String,
      workflowOptions: Option[String],
      inputBuilder: WorkflowTest => String,
      zippedImports: Option[File] = None
  ): Future[Seq[WorkflowValidation]] = {

    Source(Iterable(workflowTests: _*))
      .flatMapMerge(
        8, { test =>
          Source
            .fromFuture(
              runWorkflow(
                validationWorkflowContents,
                inputBuilder(test),
                workflowOptions,
                Some(
                  List(Label("validation-run-id", test.getRunParameters.id))),
                zippedImports,
                this.workflowName + " Validation"
              )
            )
            .map { submittedWorkflow =>
              WorkflowValidation(
                workflowTest = test,
                workflow = submittedWorkflow,
                workflowStatus = Submitted
              )
            }
        }
      )
      .runWith(Sink.collection)
  }

  /**
    * Copy results from a batch of cromwell workflows to a cloud destination
    *
    * @param finishedWorkflows finished workflows that need to have their results copied
    * @return Future denoting the work to be done
    */
  protected def copyBatchCromwellWorkflowResults(
      finishedWorkflows: Seq[WorkflowTest]
  ): Future[Seq[Unit]] = {
    Source(Iterable(finishedWorkflows: _*))
      .flatMapMerge(
        8, { finishedWorkflow =>
          val resultsCloudPath = finishedWorkflow.runParameters.resultsCloudPath
          val id = finishedWorkflow.runParameters.id

          val fut = for {
            cromwellMetadata <- cromwellClient()
              .metadata(
                finishedWorkflow.workflow.id,
                Option(Map("includeKey" -> List("backendLogs"),
                           "expandSubWorkflows" -> List("true")))
              )
              .value
              .unsafeToFuture()
            cromwellOutputs <- cromwellClient()
              .outputs(finishedWorkflow.workflow.id)
              .value
              .unsafeToFuture()
          } yield {
            logger.info(s"Copying logs for $id to ${resultsCloudPath}logs/")
            cromwellMetadata match {
              case Left(value) =>
                throw new RuntimeException(s"Cromwell error: $value")
              case Right(value) =>
                value.getLogs.foreach { log =>
                  logger.debug(s"Copying $log to ${resultsCloudPath}logs/")
                  Try(
                    ioUtil.copyGoogleObject(log,
                                            resultsCloudPath.resolve("logs/")))
                    .recover {
                      case _ =>
                        logger.info(
                          s"$log failed to be copied. Probably a call-cached task. Continuing..."
                        )
                    }
                }
            }
            logger.info(s"Copying outputs for $id to $resultsCloudPath")
            cromwellOutputs match {
              case Left(value) =>
                throw new RuntimeException(s"Cromwell error: $value")
              case Right(value) =>
                value.getAllCloudOutputs.foreach { output =>
                  val dest = resultsCloudPath.resolve(File(output.getPath).name)
                  logger.debug(s"Copying $output to $dest")
                  ioUtil.copyGoogleObject(output,
                                          resultsCloudPath.resolve(dest))
                }
            }
          }
          Source.fromFuture(fut)
        }
      )
      .runWith(Sink.collection)
  }

  /**
    * Wait for multiple Cromwell workflows to complete
    *
    * @param workflowRuns A Seq of workflow runs
    * @return A future containing the terminal status of the workflows
    */
  protected def awaitBatchCromwellWorkflowCompletion[A <: WorkflowRun[A]](
      workflowRuns: Seq[A],
  ): Future[Seq[A]] = {
    logger.info(s"Awaiting completion of ${workflowRuns.size} workflows")
    for (workflowRun <- workflowRuns) {
      val workflowId = workflowRun.workflow.id
      logger.info(
        s"A timing diagram for workflow ${workflowRun.getRunParameters.id} will be available at " +
          s"${cromwellClient().statusEndpoint(workflowId).toString().replace("status", "timing")}"
      )
    }
    def op: Future[Seq[A]] = {
      Source(Iterable(workflowRuns: _*))
        .flatMapMerge(
          8, { workflowRun =>
            if (workflowRun.workflowStatus.isInstanceOf[TerminalStatus]) {
              Source.single(workflowRun)
            } else {
              Source.fromFuture(
                cromwellClient()
                  .status(workflowRun.workflow.id)
                  .value
                  .unsafeToFuture()
                  .map {
                    case Left(value) =>
                      throw new RuntimeException(s"an error occurred: " + value)
                    case Right(status) => workflowRun.withWorkflowStatus(status)
                  }
                  // Sometimes cromwell barfs, so don't abort everything.
                  .recoverWith {
                    case ex: UnsupportedContentTypeException =>
                      logger.warn(
                        "Encountered an exception while decoding response from cromwell. Trying again.",
                        ex
                      )
                      Future.successful(workflowRun)
                  }
              )
            }
          }
        )
        .runWith(Sink.collection)
    }
    def endWhen(workflowRuns: Seq[A]): Boolean = {
      workflowRuns.groupBy(_.workflowStatus).foreach { group =>
        val (status, seq) = group
        logger.info(
          s"Workflows with status of $status: ${seq.map(_.workflow.id).mkString("[", ", ", "]")}"
        )
      }
      workflowRuns.forall(_.workflowStatus.isInstanceOf[TerminalStatus])
    }

    // Need to wait 5 seconds for cromwell to register the job
    after(5.seconds, scheduler)(retry[Seq[A]](op, endWhen, 30.seconds))
      .flatMap {
        case list
            if list.forall(sample => sample.workflowStatus.equals(Succeeded)) =>
          for (sample <- list) {
            logger.info(
              s"Workflow ${sample.getRunParameters.id} (${sample.workflow.id}) completed successfully"
            )

          }
          Future.successful(list)
        case containsFailures =>
          logger.error("Not all workflows succeeded")
          for (sample <- containsFailures) {
            logger.error(
              s"Workflow ${sample.getRunParameters.id} (${sample.workflow.id}) ended with status: ${sample.workflowStatus}"
            )
          }
          val failingIds =
            containsFailures.filter(!_.workflowStatus.equals(Succeeded))
          val msg =
            s"Workflows with ids ${failingIds.map(_.workflow.id).mkString("[", ",", "]")} failed"
          val metadataFutures: Seq[Future[WorkflowMetadata]] = failingIds.map({
            failed =>
              val endpoint =
                cromwellClient().metadataEndpoint(failed.workflow.id)
              logger.error(s"Look for failures at $endpoint")
              for {
                metadataResponse <- cromwellClient()
                  .metadata(
                    failed.workflow.id,
                    args = Option(Map("expandSubWorkflows" -> List("true")))
                  )
                  .value
                  .unsafeToFuture()
                metadata <- Future.fromTry(
                  metadataResponse match {
                    case Left(value) =>
                      Failure(new RuntimeException(s"Cromwell error: $value"))
                    case Right(value) => Try.apply(value)
                  }
                )
              } yield {
                metadata
              }
          })
          Future
            .sequence(metadataFutures)
            .map(
              metadataFuture =>
                metadataFuture.foreach(metadata => {
                  trawlMetadataPrintingFailures(metadata)
                })
            )
            .flatMap(
              _ =>
                Future.failed(
                  new CromwellWorkflowTester.TestFailedException(msg))
            )
      }
  }

  /**
    * On system shutdown, terminate the workflow in Cromwell if it is still
    * running to save on costs for larger tests
    *
    * @param submitted Workflow to abort on system shutdown
    */
  protected def addAbortTerminationHook(submitted: SubmittedWorkflow): Unit = {
    // We need to use synchronous http requests here because the
    // actor system is shutting down and won't wait for futures.
    val _ = as.registerOnTermination {
      logger.info("Searching for any running workflows in Cromwell to abort")
      val content =
        Request
          .Get(cromwellClient().statusEndpoint(submitted.id).toString())
          .execute()
          .returnContent()
      for {
        json <- parse(content.asString)
        cromwellStatus <- json.as[CromwellStatus]
      } yield {
        val workflowStatus = WorkflowStatus(cromwellStatus.status)
        if (!workflowStatus.isInstanceOf[TerminalStatus]) {
          logger.info(s"Aborting cromwell workflow with id: ${submitted.id}")
          val _ = Request
            .Post(cromwellClient().abortEndpoint(submitted.id).toString())
            .execute()
        }
      }
    }
  }

  /**
    * Convenience method to make testing easier
    *
    * @return A CromwellClient
    */
  protected def getCromwellClient: CromwellClient = {

    /** Scopes needed to get user info from a Google account to prove "valid user" identity. */
    val scopes = Seq(
      "https://www.googleapis.com/auth/userinfo.profile",
      "https://www.googleapis.com/auth/userinfo.email",
      "https://www.googleapis.com/auth/devstorage.read_write"
    )
    val googleCreds = googleCredentials.userInfo().createScoped(scopes.asJava)
    googleCreds.refresh()
    new CromwellClient(
      env.cromwellUrl,
      "1",
      Option(
        OAuth2BearerToken(
          googleCreds.getAccessToken.getTokenValue
        )
      )
    )
  }

  /**
    * Given the metadata for the Cromwell workflow run by the tester, find all task failures and print them to
    * console so we don't have to go looking for them manually
    *
    * @param metadata The metadata object to parse.
    */
  protected def trawlMetadataPrintingFailures(
      metadata: WorkflowMetadata): Unit = {
    parse(metadata.value) match {
      case Left(parsingFailure) => throw parsingFailure
      case Right(json)          =>
        /* Calls is an object of lists of objects of lists of call data, e.g.:
         * calls: {
         *   LookAtAllTheseScatters: [
         *     {
         *       "thisScattersALot": true,
         *       "lotsOfObjects": true
         *     },
         *     { ... },
         *     {
         *       "lastObject": true
         *     }
         *   ],
         *   SomethingWithoutAnyScatters: [
         *     {
         *       "onlyOneObject": true,
         *       "butStillAList": true
         *     }
         *   ]
         * }
         */
        json
        // This finds even subworkflow calls too - it returns all key/value pairs with 'calls' as the key.
          .findAllByKey("calls")
          // Everything in this json framework is a top-level Json object or an Option so we have to de-Json and
          // de-Option it to iterate. Here, we start with a List[Json], turn each Json into an Option[JsonObject],
          // then turn each Option[JsonObject] into a JsonObject. We end up with a List[JsonObject].
          .map(
            callNames =>
              callNames.asObject.getOrElse(
                throw new RuntimeException("Could not get calls as JsonObject")
            )
          )
          .foreach(
            // And we have to do this same thing WITH EACH LAYER OF NESTING.
            callsObject =>
              callsObject.keys.foreach(callName => {
                // This is now the lists of calls. This results in a List[JsonObject].
                val calls = callsObject(callName)
                  .flatMap(_.asArray)
                  .getOrElse(
                    throw new RuntimeException(
                      s"Could not get calls for $callName")
                  )
                  .flatMap(_.asObject)
                // We are finally iterating over the actual innermost call objects.
                calls.foreach(callObject => {
                  val executionStatus =
                    callObject("executionStatus").flatMap(_.asString)
                  val failed = executionStatus.map(_.equals("Failed"))
                  failed match {
                    case Some(true) => {
                      // Accessing the stderr object also involves de-Jsoning.
                      val stderrLocation = callObject("stderr")
                        .flatMap(_.asString)
                      stderrLocation match {
                        case Some(location) => {
                          val stderrUri = new URI(location)
                          // Just in case this is an absurd like 30 GB cloud file, we don't want to get the entire thing. So we fork.
                          val stderrContents =
                            if (ioUtil.googleObjectExists(stderrUri)) {
                              // For cloud files, we want to get only the last megabyte of file.
                              val stderrBlob = ioUtil.requireBlob(stderrUri)
                              val size = stderrBlob.getSize
                              try {
                                val reader = stderrBlob.reader()
                                val bytes = math.min(size, 1048576).toInt
                                val position = size - bytes
                                reader.seek(position)
                                val buffer = ByteBuffer.allocate(bytes)
                                reader.read(buffer)
                                new String(buffer.array(), "UTF-8")
                              } catch {
                                case e: Throwable => throw e
                              }
                            } else if (ioUtil.isGoogleObject(stderrUri)) {
                              // If it's a cloud object but doesn't exist, then further debugging is needed
                              val jobId =
                                callObject("jobId").flatMap(_.asString).orNull
                              s"gcloud stderr file should exist but doesn't - check worker logs with gcloud alpha genomics operations describe $jobId"
                            } else {
                              // For on-prem files, dump the entire thing.
                              ioUtil.readFile(stderrUri)
                            }

                          logger.error(
                            callName + " failed with the following stderr:")
                          logger.error(stderrContents)
                          logger.error("stderr location: " + location)
                        }
                        case None => {
                          // If there is no stderr, try getting the failures block(s).
                          logger.error(s"Could not get stderr for $callName.")
                          val failures = callObject("failures")
                            .flatMap(_.asArray)
                            .getOrElse(
                              throw new RuntimeException(
                                s"Could not get stderr or failures for $callName!"
                              )
                            )
                            .flatMap(_.asObject)
                          // Accessing the failures object also also involves de-Jsoning.git status
                          failures.foreach(failureObject => {
                            val message = failureObject("message")
                              .flatMap(_.asString)
                            message match {
                              case Some(messageText) => {
                                logger.error(
                                  callName + " failed with the following failure message:"
                                )
                                logger.error(messageText)
                              }
                              case None => {
                                logger.error(
                                  s"Could not get failure message(s) for $callName!"
                                )
                              }
                            }
                          })
                        }
                      }
                    }
                    case _ =>
                      logger.info(
                        s"executionStatus for $callName is $executionStatus")
                  }
                })
              })
          )
    }
  }

  protected def assertEqual[A](actual: A,
                               expected: A,
                               description: String): Try[Unit] =
    if (actual == expected) {
      Success(())
    } else {
      val err =
        s"""Unexpected value for $description:
           |Expected: $expected
           |Actual: $actual
         """.stripMargin
      Failure(new CromwellWorkflowTester.TestFailedException(err))
    }

  /**
    * Read the main workflow WDL from the provided release dir
    *
    * @param releaseDir The release dir with WDLs and dependencies
    * @return The WDL contents
    */
  protected def readWdlFromReleaseDir(releaseDir: File,
                                      workflowName: String): String =
    (releaseDir / workflowName / s"${workflowName}_${WarpGitHash}.wdl").contentAsString

  /**
    * Get the dependencies zip from the provided release dir
    *
    * @param releaseDir The release dir with WDLs and dependencies
    * @param workflowName The name of the workflow
    * @return A File of the dependencies ZIP
    */
  protected def dependenciesZipFromReleaseDir(
      releaseDir: File,
      workflowName: String): Option[File] = {
    val zipFile = releaseDir / workflowName / s"${workflowName}_${WarpGitHash}.zip"
    if (zipFile.exists) Option(zipFile) else None
  }

  /**
    * Get the data type string prefix from the DataType enum
    *
    * @param dataType DataType to get string name of
    * @return The string name of the DataType
    */
  protected def dataTypePrefix(dataType: DataType): String = dataType match {
    case DataType.WGS      => "WholeGenome"
    case DataType.Exome    => "Exome"
    case DataType.RNA      => "RNA"
    case DataType.Targeted => "Targeted"
  }
}
