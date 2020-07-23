package org.broadinstitute.dsp.pipelines.tester

import java.net.URI
import java.time.{Instant, ZoneId, ZoneOffset}

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import akka.stream.scaladsl.Sink
import better.files.File
import io.circe.parser._
import org.broadinstitute.clio.client.commands.MoveArrays
import org.broadinstitute.clio.client.dispatch.Executor
import org.broadinstitute.clio.client.dispatch.MoveExecutor
import org.broadinstitute.clio.transfer.model.ArraysIndex
import org.broadinstitute.clio.transfer.model.arrays.{ArraysKey, ArraysMetadata}
import org.broadinstitute.clio.util.model.Location
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment
import org.broadinstitute.dsp.pipelines.config.ArraysDataDeliveryConfig
import org.broadinstitute.dsp.pipelines.file.TsvParser
import org.broadinstitute.dsp.pipelines.firecloud.model.autogen.{
  WorkspaceACL,
  WorkspaceCatalog
}
import org.broadinstitute.dsp.pipelines.inputs.ArraysDataDeliveryInputs

import scala.collection.immutable
import scala.concurrent.Future
import scala.util.{Failure, Success, Try}

class ArraysDataDeliveryTester(testerConfig: ArraysDataDeliveryConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends FireCloudDeliveryTester[ArraysKey, ArraysMetadata] {

  import CromwellWorkflowTester.DsdePipelinesRoot

  override val env = CromwellEnvironment.Dev

  private val user = System.getProperty("user.name", "greenteam")
  private val requester = "test.firec@gmail.com"
  private val timestamp = Instant.now.toString.replaceAll("[:.]", "-")

  override val workflowName: String = "ArraysDataDelivery"
  override val differentiator = s"test-arrays-$user-$timestamp"

  private val workspaceName = s"$differentiator-workspace"
  private val picardServiceAccountName = "picard-service-account-pem.json"

  private val storage =
    s"broad-gotc-$envString-storage/test/arrays/delivery/$timestamp"
  private val browserPrefix = "https://console.cloud.google.com/storage/browser"
  private val browser =
    s"$browserPrefix/$storage/?project=broad-gotc-$envString-storage"
  private val gsStorage = s"gs://$storage"
  private val gsDir = URI.create(s"$gsStorage/")
  private val libraryCuratorEmail = "GROUP_DSP_library_curator@firecloud.org"
  private val methodsDevEmail = "GROUP_DSP_library_curator@firecloud.org"

  override val inputsDir
    : File = DsdePipelinesRoot / "delivery" / "arrays" / "test_inputs"

  override val workflowTimezoneId: ZoneId = ZoneOffset.UTC

  logger.info(s"Find uploaded test files under here in browser:\n$browser")

  protected lazy val fakeReleaseDir: File =
    CromwellWorkflowTester.runReleaseWorkflow(
      DsdePipelinesRoot / "delivery" / "arrays" / s"$workflowName.wdl",
      env
    )

  override lazy val wdlContents: String = readWdlFromReleaseDir(fakeReleaseDir)

  private val wdlOptions =
    (fakeReleaseDir / workflowName / s"ArraysDataDelivery.options.json").contentAsString

  private val cloudPicardServiceAccountJson =
    gsDir.resolve(picardServiceAccountName)
  private val cloudSamplesMetadataFile =
    gsDir.resolve(s"$differentiator.metadata")

  private val chipwellBarcodeVersionTsv = inputsDir / "ArraysDataDeliveryChipwellBarcodeVersion.tsv"

  private val cloudChipwellBarcodeVersionTsv =
    gsDir.resolve("ArraysDataDeliveryChipwellBarcodeVersion.tsv")

  private val wdlInputs: List[String] =
    (inputsDir / "ArraysDataDeliveryWf.inputs.json").lines
      .map(_.replace("{REQUESTER}", requester))
      .map(_.replace("{SAMPLE_SET}", s"$differentiator-sampleset"))
      .map(_.replace("{WORKSPACE_NAME}", workspaceName))
      .map(_.replace("{INPUT_TSV}", cloudChipwellBarcodeVersionTsv.toString))
      .map(_.replace("{SERVICE_ACCOUNT_JSON}",
                     cloudPicardServiceAccountJson.toString))
      .map(_.replace("{ENV}", envString.toUpperCase))
      .map(_.replace("{SAMPLES_METADATA}", cloudSamplesMetadataFile.toString))
      .map(
        _.replace(
          "{COMPUTE_SERVICE_ACCOUNT_JSON}",
          cloudPicardServiceAccountJson.toString
        )
      )
      .toList

  override val wdlInputsDecoded = ArraysDataDeliveryInputs(wdlInputs)

  override lazy val preparedMetadata: ArraysMetadata = ArraysMetadata(
    workspaceName = Option(wdlInputsDecoded.workspaceName)
  )

  override def runTest: Future[Unit] = {
    val inputs = parse(wdlInputs.mkString).fold(throw _, _.spaces2)
    val options = parse(wdlOptions).fold(throw _, _.spaces2)

    val keys = readKeys()

    // Copy inputs to the cloud.
    ioUtil.writeGoogleObjectData(
      serviceAccountJson.contentAsString,
      cloudPicardServiceAccountJson
    )
    ioUtil.writeGoogleObjectData(
      chipwellBarcodeVersionTsv.contentAsString,
      cloudChipwellBarcodeVersionTsv
    )
    ioUtil.writeGoogleObjectData(
      differentiator,
      cloudSamplesMetadataFile
    )

    val initialMetadata = getInitialMetadata(keys)

    val test = for {
      _ <- initialMetadata
      submittedWorkflow <- runWorkflow(
        wdlContents,
        inputs,
        Some(options),
        None,
        Some(dependenciesZipFromReleaseDir(fakeReleaseDir))
      )
      id = submittedWorkflow.id
      cloudExecutionDir = URI.create(
        s"gs://broad-gotc-$envString-cromwell-execution/$workflowName/$id/"
      )

      _ = logger.info(s"Cromwell executes here: $cloudExecutionDir")
      _ = logger.info(
        s"Browse it here: $browserPrefix/${cloudExecutionDir.toString.stripPrefix("gs://")}"
      )

      _ <- awaitCromwellWorkflowCompletion(submittedWorkflow)
      cromwellOutputs <- cromwellClient().outputs(id).value.unsafeToFuture()
      workspaceBucket <- Future.fromTry(
        cromwellOutputs match {
          case Left(value) =>
            Failure(new RuntimeException(s"Cromwell error: $value"))
          case Right(value) =>
            value.getStringFromOutputs(s"$workflowName.workspace_bucket")
        }
      )

      _ = logger.info("Verifying workflow results.")

      expectedSamplesTsv = TsvParser.parseFromInput(
        inputsDir / "expected" / "samples.tsv",
        transformExpected(_).replace("{FC_BUCKET}", workspaceBucket)
      )
      _ <- Future.sequence(
        keys.map { key =>
          expectedSamplesTsv
            .find { row =>
              row("chip_well_barcode") == key.chipwellBarcode.name &&
              row("version") == key.version.toString
            }
            .fold(
              Future.failed[Unit](
                new IllegalStateException(
                  s"No expected sample row for input $key")
              )
            )(checkSingleSampleFileDelivery(key, _))
        }
      )

      expectedSampleSetTsv = TsvParser.parseFromInput(
        inputsDir / "expected" / "sample_set_entity.tsv",
        transformExpected(_).replace("{FC_BUCKET}", workspaceBucket)
      )
      _ <- checkMultiSampleFileDelivery(expectedSampleSetTsv.head)

      // Check that the requester has read access to the workspace
      _ = logger.info(s"Verifying access control lists.")
      wacl <- fireCloudClient.getWorkspaceAcl
      _ <- checkAccessControl(requester, "READER", wacl)
      _ <- checkAccessControl(methodsDevEmail, "READER", wacl)
      _ <- checkAccessControl(libraryCuratorEmail, "READER", wacl)
      _ <- checkCataloger
      _ <- verifyWorkflow(workspaceBucket)
    } yield ()

    lazy val runCleanup = initialMetadata.flatMap(cleanup)

    test.transformWith {
      case Failure(ex) => runCleanup.transform(_ => throw ex)
      case Success(_)  => runCleanup
    }
  }

  private def checkCataloger: Future[Unit] = {
    fireCloudClient.getWorkspaceCatalogers.flatMap { catalogers =>
      Future.fromTry(
        assertEqual(
          catalogers.contains(
            WorkspaceCatalog(libraryCuratorEmail, catalog = true)),
          true,
          s"$libraryCuratorEmail should be a cataloger"
        )
      )
    }
  }

  private def checkAccessControl(
      email: String,
      accessLevel: String,
      workspaceAcl: WorkspaceACL
  ): Future[Unit] = {
    Future.fromTry(
      assertEqual(
        workspaceAcl.acl(email).accessLevel,
        accessLevel,
        s"The $email should have $accessLevel access"
      )
    )
  }

  private def getInitialMetadata(
      keys: Seq[ArraysKey]
  ): Future[Map[ArraysKey, ArraysMetadata]] = {
    Future
      .sequence {
        keys.map { k =>
          clioWebClient
            .getMetadataForKey(ArraysIndex)(k, includeDeleted = false)
            .map(k -> _)
            .runWith(Sink.head)
        }
      }
      .map(_.toMap)
  }

  private def checkSingleSampleFileDelivery(
      key: ArraysKey,
      expectedSamplesRow: Map[String, String]
  ): Future[Unit] = {
    logger.info(s"Verifying that files for $key were actually delivered")

    def checkDelivered(tsvHeader: String, clioUri: URI): Try[Unit] = {
      if (!ioUtil.googleObjectExists(clioUri)) {
        Failure(
          new CromwellWorkflowTester.TestFailedException(
            s"Dangling path found in clio after delivery: $clioUri"
          )
        )
      } else {
        val expectedPath = expectedSamplesRow(tsvHeader)
        if (clioUri != URI.create(expectedPath)) {
          Failure(
            new CromwellWorkflowTester.TestFailedException(
              s"Cloud file delivered to wrong path: Expected '$expectedPath' but found '$clioUri'"
            )
          )
        } else {
          Success(())
        }
      }
    }

    clioWebClient
      .getMetadataForKey(ArraysIndex)(key, includeDeleted = false)
      .mapConcat { metadata =>
        immutable.Iterable.concat(
          metadata.vcfPath.map("vcf" -> _),
          metadata.vcfIndexPath.map("vcf_index" -> _),
          metadata.gtcPath.map("gtc" -> _),
          metadata.grnIdatPath.map("green_idat" -> _),
          metadata.redIdatPath.map("red_idat" -> _)
        )
      }
      .mapAsyncUnordered(5) {
        case (header, path) => Future(checkDelivered(header, path))
      }
      .runFold(Seq.empty[String]) { (errs, check) =>
        check.fold(err => err.getMessage +: errs, _ => errs)
      }
      .flatMap { maybeErrs =>
        if (maybeErrs.isEmpty) {
          Future.unit
        } else {
          Future.failed(
            new CromwellWorkflowTester.TestFailedException(
              maybeErrs.mkString("\n"))
          )
        }
      }
  }

  private def checkMultiSampleFileDelivery(
      expectedSampleSetRow: Map[String, String]
  ): Future[Unit] = {
    logger.info(
      "Verifying that the multi-sample arrays files were actually delivered.")
    val vcfPath = expectedSampleSetRow("vcf")
    val vcf = URI.create(vcfPath)
    val multiSampleFolder = vcf.resolve(".")

    val allExist = immutable
      .Iterable(
        vcf,
        URI.create(s"$vcfPath.tbi"),
        URI.create(expectedSampleSetRow("sample_metadata")),
        multiSampleFolder.resolve(
          "PsychChip_v1-1_15073391_A1.1.3.extended.csv"),
        multiSampleFolder.resolve("PsychChip_v1-1_15073391_A1.csv"),
        multiSampleFolder.resolve(
          "PsychChip_v1-1_15073391_A1.1.3.bad_assays.csv"),
        multiSampleFolder.resolve("PsychChip_v1-1_15073391_A1_ClusterFile.egt")
      )
      .filterNot(ioUtil.googleObjectExists)

    if (allExist.isEmpty) {
      Future.unit
    } else {
      Future.failed(
        new RuntimeException(
          s"Delivered files not found! ${allExist.map(_.toString).mkString("\n")}"
        )
      )
    }
  }

  /** Get Clio keys for the crams specified in the delivery input TSV. */
  override def readKeys(): Seq[ArraysKey] = {
    val lines = chipwellBarcodeVersionTsv.lines.toList
    lines
      .map(_.split("\t"))
      .map(a => (a(0), a(1)))
      .map { t =>
        val (chipwellBarcode, version) = t
        ArraysKey(
          chipwellBarcode = Symbol(chipwellBarcode),
          version = version.toInt,
          location = Location.GCP
        )
      }
  }

  /**
    * Clean up after the test. The crams need to be moved back to their
    * original location and the FireCloud workspace needs to be deleted
    *
    * @return Success or failure of cleanup
    */
  private[tester] def cleanup(
      initialMetadata: Map[ArraysKey, ArraysMetadata]
  ): Future[Unit] = {
    import Executor.SourceMonadOps

    logger.info("Deleting files that were staged in the cloud.")

    ioUtil.deleteGoogleObject(cloudPicardServiceAccountJson)
    ioUtil.deleteGoogleObject(cloudChipwellBarcodeVersionTsv)

    logger.info("Moving things back in Clio and deleting FireCloud workspace")

    val moved: Iterable[Future[Try[Unit]]] = initialMetadata.map {
      case (key, metadata) =>
        val originalLocation = metadata.vcfPath
          .orElse(metadata.gtcPath)
          .fold(
            throw new IllegalStateException(
              s"Just ran a pointless test; no VCF or GTC for $key"
            )
          ) { uri =>
            val str = uri.toString
            URI.create(str.splitAt(str.lastIndexOf('/') + 1)._1)
          }

        val resetStream = for {
          _ <- new MoveExecutor(MoveArrays(key, originalLocation))
            .execute(clioWebClient, ioUtil)
          idatPatch = ArraysMetadata(
            redIdatPath = metadata.redIdatPath,
            grnIdatPath = metadata.grnIdatPath,
            workspaceName = Some("")
          )
          _ = logger.info(
            s"Resetting idat paths and blanking workspace name for $key")
          _ <- clioWebClient.upsert(ArraysIndex)(key, idatPatch, force = true)
        } yield ()

        resetStream.runWith(Sink.ignore).transformWith {
          case Failure(ex) =>
            Future.successful(
              Failure(
                new RuntimeException(
                  s"""Failed to move $key back to its original location.
                     |All the files should be moved back to the same location as the params_path""".stripMargin,
                  ex
                )
              )
            )
          case Success(_) =>
            Future.successful(Success(()))
        }
    }
    if (!testerConfig.leaveWorkspace) {
      cleanupFireCloudWorkspace(initialMetadata.keys, moved)
    } else {
      Future.sequence(moved).map(_ => ())
    }
  }
}
