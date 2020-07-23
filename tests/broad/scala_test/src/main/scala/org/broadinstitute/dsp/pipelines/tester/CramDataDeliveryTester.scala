package org.broadinstitute.dsp.pipelines.tester

import java.net.URI
import java.time.ZoneId
import java.util.UUID

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import akka.stream.scaladsl.Sink
import better.files.File
import org.broadinstitute.clio.client.commands.DeliverCram
import org.broadinstitute.clio.client.dispatch.DeliverExecutor
import org.broadinstitute.clio.transfer.model.cram.{CramKey, CramMetadata}
import org.broadinstitute.clio.util.model.Location
import org.broadinstitute.clio.util.model.{DataType => ClioDataType}
import org.broadinstitute.dsp.pipelines.util.DataType
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment
import org.broadinstitute.dsp.pipelines.config.CramDataDeliveryConfig
import org.broadinstitute.dsp.pipelines.inputs.{
  CramDataDeliveryInputs,
  FireCloudWdlInputs
}

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

class CramDataDeliveryTester(testerConfig: CramDataDeliveryConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends FireCloudDeliveryTester[CramKey, CramMetadata] {

  override val env: CromwellEnvironment = testerConfig.env

  val workflowName = s"CramDataDelivery"

  private lazy val dataTypeString = testerConfig.dataType.entryName.toLowerCase

  override val inputsDir
    : File = CromwellWorkflowTester.DsdePipelinesRoot / "delivery" / "cram" / "test_inputs" / dataTypeString

  protected lazy val timestamp: String = CromwellWorkflowTester.getTimestamp

  protected lazy val inputsCloudPath: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/cram_data_delivery/$dataTypeString/$timestamp/"
    )

  override val workflowTimezoneId: ZoneId = ZoneId.systemDefault()

  lazy val workflowDir: File =
    CromwellWorkflowTester.DsdePipelinesRoot / "delivery" / "cram"

  private lazy val fakeReleaseDir =
    CromwellWorkflowTester.runReleaseWorkflow(
      workflowDir / s"$workflowName.wdl",
      env
    )

  override val wdlContents: String = readWdlFromReleaseDir(fakeReleaseDir)

  val inputsTsv: File = inputsDir / "expected" / "deliverable_samples.tsv"

  override lazy val differentiator: String = UUID.randomUUID().toString

  val wdlInputContents: List[String] =
    (inputsDir / "cramDataDeliveryInputs.json").lines
      .map(_.replace("{DIFF}", differentiator))
      .map(_.replace("{PEM_JSON}", serviceAccountJson.pathAsString))
      .map(_.replace("{REQUESTER}", testerConfig.requester))
      .map(_.replace("{INPUTS_BUCKET}", inputsCloudPath.toString))
      .toList

  val wdlRedeliveryInputContents: List[String] =
    (inputsDir / "cramDataRedeliveryInputs.json").lines
      .map(_.replace("{DIFF}", differentiator))
      .map(_.replace("{PEM_JSON}", serviceAccountJson.pathAsString))
      .map(_.replace("{REQUESTER}", testerConfig.requester))
      .map(_.replace("{INPUTS_BUCKET}", inputsCloudPath.toString))
      .toList

  override val wdlInputsDecoded: FireCloudWdlInputs =
    CramDataDeliveryInputs(wdlInputContents)

  override lazy val preparedMetadata: CramMetadata = CramMetadata(
    workspaceName = Option(wdlInputsDecoded.workspaceName)
  )

  private val clioDataType: ClioDataType = testerConfig.dataType match {
    case DataType.WGS   => ClioDataType.WGS
    case DataType.Exome => ClioDataType.Exome
    case DataType.RNA =>
      throw new UnsupportedOperationException("RNA data type not supported")
    case DataType.Targeted =>
      throw new UnsupportedOperationException(
        "Targeted data type not supported")
  }

  def uploadFireCloudMetadata(bucket: String): Future[Unit] = {
    def writeCloudFile(resource: String) = {
      val data = (inputsDir / "expected" / resource).lines
        .map(transformExpected)
        .map(_.replace("{FC_BUCKET}", bucket))
        .mkString("\n")
      stagingIoUtil.writeGoogleObjectData(data,
                                          inputsCloudPath.resolve(resource))
    }

    writeCloudFile("samples.tsv")
    writeCloudFile("participants.tsv")
    writeCloudFile("sample_set_entity.tsv")
    writeCloudFile("sample_set_membership.tsv")
    writeCloudFile("libAttrs.json")
    writeCloudFile("tags.json")
    writeCloudFile("deliverable_samples.tsv")
    stagingIoUtil.writeGoogleObjectData(
      serviceAccountJson.contentAsString,
      inputsCloudPath.resolve("picard-account.json")
    )

    Future.successful(())
  }

  /**
    * Run the test. This method is the entrypoint for all tests
    * and must be implemented
    *
    * @return
    */
  override def runTest: Future[Unit] = {
    val workflowOptions =
      (fakeReleaseDir / workflowName / s"$workflowName.options.json").contentAsString
    val test = for {
      createdWorkspace <- fireCloudClient.createWorkspace
      _ = logger.info(
        s"Created Firecloud workspace: ${createdWorkspace.namespace}/${createdWorkspace.name}"
      )
      _ <- uploadFireCloudMetadata(createdWorkspace.bucketName)
      submittedWorkflow <- runWorkflow(
        wdlContents,
        wdlInputContents
          .map(_.replace("{FC_BUCKET}", createdWorkspace.bucketName))
          .mkString,
        Option(workflowOptions),
        None,
        Option(dependenciesZipFromReleaseDir(fakeReleaseDir))
      )
      _ <- awaitCromwellWorkflowCompletion(submittedWorkflow)
      _ = logger.info("Verifying workflow results")
      _ <- verifyWorkflow(createdWorkspace.bucketName)
      // we redeliver the workspace here to make sure that cram data deliveries behave the same when delivering to an existing workspace
      _ = logger.info("Redelivering workspace")
      resubmittedWorkflow <- runWorkflow(
        wdlContents,
        wdlRedeliveryInputContents
          .map(_.replace("{FC_BUCKET}", createdWorkspace.bucketName))
          .mkString,
        Option(workflowOptions),
        None,
        Option(dependenciesZipFromReleaseDir(fakeReleaseDir))
      )
      _ <- awaitCromwellWorkflowCompletion(resubmittedWorkflow)
      _ = logger.info("Verifying workflow results")
      _ <- verifyWorkflow(createdWorkspace.bucketName)
    } yield ()
    test.transformWith {
      case Failure(ex) =>
        cleanup().transform(_ => throw ex)
      case Success(_) => cleanup()
    }
  }

  /**
    * Clean up after the test. The crams + metrics need to be moved back to their
    * original location and the FireCloud workspace needs to be deleted
    *
    * @return Success or failure of cleanup
    */
  private[tester] def cleanup(): Future[Unit] = {
    logger.info("Deleting picard account json from staging path")
    stagingIoUtil.deleteGoogleObject(
      inputsCloudPath.resolve("picard-account.json"))

    logger.info("Moving things back in Clio and deleting FireCloud workspace")

    val keys = readKeys()

    val moved: Seq[Future[Try[Unit]]] = keys.map { key =>
      val bucket = testerConfig.dataType match {
        case DataType.WGS   => "broad-gotc-dev-storage"
        case DataType.Exome => "broad-exomes-dev-storage"
        case DataType.RNA =>
          throw new UnsupportedOperationException("RNA data type not supported")
        case DataType.Targeted =>
          throw new UnsupportedOperationException(
            "Targeted data type not supported")
      }
      val devStorage =
        s"gs://$bucket/pipeline/${key.project}/${key.sampleAlias}/v${key.version}/"
      val cram = URI.create(s"$devStorage/${key.sampleAlias}.cram")
      logger.info(s"Moving $key back to $devStorage")
      if (!ioUtil.googleObjectExists(cram)) {
        new DeliverExecutor(
          DeliverCram(key, "", "", URI.create(devStorage), force = true)
        ).execute(clioWebClient, ioUtil)
          .runWith(Sink.ignore)
          .transformWith {
            case Failure(ex) =>
              Future.successful(
                Failure(
                  new RuntimeException(
                    s"Failed to move $key back to $devStorage",
                    ex)
                )
              )
            case Success(_) =>
              Future.successful(Success(()))
          }
      } else {
        logger.warn(s"$key has already been moved back. Skipping.")
        Future.successful(Success(()))
      }
    }

    if (!testerConfig.leaveWorkspace) {
      cleanupFireCloudWorkspace(keys, moved)
    } else {
      Future.sequence(moved).map(_ => ())
    }
  }

  /** Get Clio keys for the crams specified in the delivery input TSV. */
  override def readKeys(): Seq[CramKey] = {
    inputsTsv.lineIterator.toList
      .map(_.split("\t"))
      .map(a => (a(0), a(1), a(2), a(3)))
      .map { t =>
        val (project, sampleAlias, version, _) = t
        CramKey(
          project = project,
          dataType = clioDataType,
          sampleAlias = sampleAlias,
          version = version.toInt,
          location = Location.GCP
        )
      }
  }
}
