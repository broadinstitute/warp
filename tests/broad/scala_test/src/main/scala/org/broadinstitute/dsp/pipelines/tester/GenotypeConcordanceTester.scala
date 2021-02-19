package org.broadinstitute.dsp.pipelines.tester

import java.net.{InetAddress, URI}

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import akka.stream.alpakka.ftp.FtpSettings
import akka.stream.alpakka.ftp.scaladsl.Ftp
import akka.stream.scaladsl.Sink
import better.files.{File, Resource}
import cromwell.api.model.SubmittedWorkflow
import org.broadinstitute.clio.client.webclient.ClioWebClient
import org.broadinstitute.clio.transfer.model.GvcfIndex
import org.broadinstitute.clio.util.json.ModelAutoDerivation
import org.broadinstitute.dsp.pipelines.util.DataType.WGS
import org.broadinstitute.dsp.pipelines.commandline.WorkflowTestCategory
import org.broadinstitute.dsp.pipelines.config.{
  GenotypeConcordanceConfig,
  GermlineCloudWorkflowConfig
}
import org.broadinstitute.dsp.pipelines.file.TsvParser

import scala.concurrent.Future
import scala.util.Failure

class GenotypeConcordanceTester(testerConfig: GenotypeConcordanceConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends CromwellWorkflowTester
    with ModelAutoDerivation {

  override val env = testerConfig.env

  override val workflowName: String = "GenotypeConcordanceWorkflow"
  private val wdlName = "GenotypeConcordanceWf"

  val nistFtpSettings: FtpSettings = FtpSettings(
    host = InetAddress.getByName("ftp-trace.ncbi.nlm.nih.gov"),
    passiveMode = true
  )

  val geneConcordDir
    : File = CromwellWorkflowTester.GotcRoot / "genotype_concordance"

  override val wdlContents: String =
    (geneConcordDir / s"$wdlName.wdl").contentAsString

  val truthUri: URI = URI.create(
    s"gs://broad-gotc-test-storage/genotype_concordance/truth/v${testerConfig.nistVersion}/"
  )

  val wdlOptionsContents: String =
    (geneConcordDir / s"$wdlName.options.json").contentAsString

  override def runTest: Future[Unit] = {
    val (truthVcf, truthVcfIndex): (String, String) = findTruthFiles()

    for {
      _ <- assertLatestNistRelease(truthVcf, truthVcfIndex)
      gvcfPath <- getOrGenerateGvcf()
      submittedWorkflow <- runWorkflow(
        wdlContents,
        buildWorkflowInputs(gvcfPath, truthVcf, truthVcfIndex),
        Some(wdlOptionsContents)
      )
      _ <- awaitCromwellWorkflowCompletion(submittedWorkflow)
      _ <- validateResults(submittedWorkflow)
    } yield ()
  }

  private def buildWorkflowInputs(
      testGvcf: String,
      truthVcf: String,
      truthVcfIndex: String
  ): String = {
    Resource
      .getAsString("genotypeConcordance/genotypeConcordanceInputs.json")
      .lines
      .map(
        _.replace("{GVCF}", testGvcf)
          .replace("{GVCF_INDEX}", s"$testGvcf.tbi")
          .replace("{TRUTH_VCF}", truthVcf)
          .replace("{TRUTH_VCF_INDEX}", truthVcfIndex)
      )
      .mkString
  }

  private def validateResults(workflow: SubmittedWorkflow): Future[Unit] = {
    for {
      outputs <- cromwellClient().outputs(workflow.id).value.unsafeToFuture()
      outputPath <- Future.fromTry(
        outputs match {
          case Left(value) =>
            Failure(new RuntimeException(s"Cromwell error: $value"))
          case Right(value) =>
            value.getStringFromOutputs(s"$workflowName.summary_metrics")
        }
      )
    } yield {
      val summaryMetricsPath = URI.create(outputPath)
      logger.info(s"Looking at summary metrics at $summaryMetricsPath")
      val summaryMetrics =
        ioUtil.readGoogleObjectData(summaryMetricsPath).split("\n")
      val parsed = TsvParser.parseExcludeComments(summaryMetrics, "#")
      val passed = parsed.forall { line =>
        line("GENOTYPE_CONCORDANCE").toDouble > 0.95 &&
        line("NON_REF_GENOTYPE_CONCORDANCE").toDouble > 0.95
      }

      if (!passed) {
        throw new CromwellWorkflowTester.TestFailedException(
          "Concordance did not pass quality threshold."
        )
      }
    }
  }

  private def findTruthFiles(): (String, String) = {
    val truthFiles = ioUtil
      .listGoogleObjects(truthUri)
      .map(_.toString)

    val maybeFiles = for {
      vcf <- truthFiles.find(_.endsWith("vcf.gz"))
      idx <- truthFiles.find(_.endsWith("vcf.gz.tbi"))
    } yield {
      (vcf, idx)
    }

    maybeFiles.getOrElse {
      throw new IllegalStateException(
        s"Could not find truth VCF / index in $truthUri; found: $truthFiles"
      )
    }
  }

  private def assertLatestNistRelease(
      vcfPath: String,
      vcfIndexPath: String
  ): Future[Unit] = {
    val gotcTruthFiles = List(vcfPath, vcfIndexPath)
      .map(_.split("/").last)

    val latestNistFiles = Ftp
      .ls("giab/ftp/release/NA12878_HG001/latest/GRCh38/", nistFtpSettings)
      .map(_.name)
      .runWith(Sink.seq)

    latestNistFiles.map { nistFiles =>
      val upToDate = gotcTruthFiles.forall(nistFiles.contains)

      if (!upToDate) {
        throw new CromwellWorkflowTester.TestFailedException(
          s"The latest files in $truthUri do not match the latest NIST release: $latestNistFiles"
        )
      } else {
        logger.info(
          "Samples in gcloud are up to date with the latest NIST release")
      }
    }
  }

  private def getOrGenerateGvcf(): Future[String] = {
    if (testerConfig.useLatestNA12878) {
      queryClioForLatestGvcf()
    } else {
      generateFreshGvcf()
    }
  }

  private def queryClioForLatestGvcf(): Future[String] = {
    import io.circe.literal._

    logger.info("Retrieving latest nightly GVCF from Clio")
    val clioClient = ClioWebClient(
      googleCredentials,
      clioHost = s"clio.gotc-$envString.broadinstitute.org",
      clioPort = 443,
      useHttps = true
    )

    val gvcfField = "gvcf_path"

    val query =
      json"""{
            "_source": $gvcfField,
            "query": {
              "bool": {
                "must": [
                  {"term": {"project.exact": "G96830"}},
                  {"term": {"sample_alias.exact": "NA12878"}}
                ]
              }
            },
            "sort": [ {"version": "desc"} ],
            "size": 1
          }"""

    clioClient
      .query(GvcfIndex)(query, raw = true)
      .runWith(Sink.headOption)
      .map {
        _.fold(
          throw new IllegalStateException(
            s"No latest GVCF registered in Clio in $envString"
          )
        ) { json =>
          import org.broadinstitute.clio.JsonUtils.JsonOps
          json.unsafeGet[String](gvcfField)
        }
      }
  }

  private def generateFreshGvcf(): Future[String] = {
    logger.info("Generating fresh germline single-sample GVCF")

    val singleSampleTest = new GermlineSingleSampleTester(
      GermlineCloudWorkflowConfig(
        dataType = WGS,
        category = WorkflowTestCategory.Scientific,
        env = testerConfig.env
      )
    )

    val germlineWdlName = singleSampleTest.workflowName
    val tmpReleaseDir = CromwellWorkflowTester.runReleaseWorkflow(
      singleSampleTest.workflowDir / s"$germlineWdlName.wdl",
      env
    )

    for {
      submission <- runWorkflow(
        (tmpReleaseDir / germlineWdlName / s"$germlineWdlName.wdl").contentAsString,
        singleSampleTest.getInputContents("G96830.NA12878.json"),
        Some(
          singleSampleTest.readTestOptions(tmpReleaseDir, env)
        ),
        None,
        dependenciesZipFromReleaseDir(tmpReleaseDir, workflowName),
        singleSampleTest.workflowName
      )
      _ <- awaitCromwellWorkflowCompletion(submission)
      outputs <- cromwellClient().outputs(submission.id).value.unsafeToFuture()
    } yield {
      val expectedOutput =
        s"${singleSampleTest.workflowName}.output_vcf"

      outputs match {
        case Left(value) =>
          throw new RuntimeException(s"Cromwell error: $value")
        case Right(value) =>
          value
            .getStringFromOutputs(expectedOutput)
            .getOrElse(
              throw new IllegalStateException(
                s"Fresh single-sample run with ID '${submission.id.id}' has no output '$expectedOutput'"
              )
            )
      }
    }
  }

}
