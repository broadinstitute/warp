package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import cromwell.api.model.SubmittedWorkflow
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment
import org.broadinstitute.dsp.pipelines.config.AnnotationFiltrationConfig
import org.broadinstitute.dsp.pipelines.file.TsvParser
import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRunParameters,
  WorkflowTest
}
import org.broadinstitute.dsp.pipelines.tester.CromwellWorkflowTester.WarpGitHash

import scala.concurrent.Future
import scala.sys.process._
import scala.util.{Failure, Success, Try}

class AnnotationFiltrationTester(testConfig: AnnotationFiltrationConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends CromwellWorkflowTester {

  override val env: CromwellEnvironment = testConfig.env

  override val workflowName = "AnnotationFiltration"

  val workflowDir
    : File = CromwellWorkflowTester.PipelineRoot / "broad" / "annotation_filtration"

  protected lazy val releaseDir: File =
    CromwellWorkflowTester.runReleaseWorkflow(
      workflowDir / s"$workflowName.wdl",
      env
    )

  override lazy val wdlContents: String =
    (releaseDir / workflowName / s"${workflowName}_${WarpGitHash}.wdl").contentAsString

  val testInputs: String = (workflowDir / "test_inputs" / {
    testConfig.category.entryName
  } / "hg38.json").contentAsString

  override def runTest: Future[Unit] = {
    val workflowOptions =
      (releaseDir / workflowName / s"${workflowName}_${WarpGitHash}.options.json").contentAsString
    val workflowRuns = Seq(
      WorkflowRunParameters(
        id = s"${envString}_${testConfig.category.entryName}",
        workflowInputs = testInputs,
        resultsCloudPath = URI.create(
          s"gs://broad-gotc-test-results/$envString/annotation_filtration/${testConfig.category.entryName.toLowerCase}/${CromwellWorkflowTester.getTimestamp}/"
        ),
        truthCloudPath = URI.create("")
      )
    )
    for {
      submittedWorkflows <- submitBatchWorkflows(
        workflowRuns,
        Some(workflowOptions),
        dependenciesZipFromReleaseDir(releaseDir, workflowName)
      )
      finishedSamples <- awaitBatchCromwellWorkflowCompletion(
        submittedWorkflows)
      _ <- copyBatchCromwellWorkflowResults(finishedSamples)
      outputReports <- Future.sequence(
        submittedWorkflows
          .map(x => readTsvFromCromwellOutput(x.workflow, "filtration_report"))
      )
      _ <- Future.fromTry(checkReport(outputReports.flatten))
      costTuples <- getCosts(submittedWorkflows)
    } yield {
      costTuples.foreach { tup =>
        val (workflowTest, cost) = tup
        logger.info(s"${workflowTest.runParameters.id} costs: $cost")
      }
    }
  }

  private def getCosts(
      workflows: Seq[WorkflowTest]
  ): Future[Seq[(WorkflowTest, String)]] = {
    Future.sequence(workflows.map { workflow =>
      val workflowId = workflow.workflow.id
      val tempFile = File.newTemporaryFile().deleteOnExit()
      for {
        metadata <- cromwellClient()
          .metadata(workflowId,
                    args = Option(Map("expandSubWorkflows" -> List("true"))))
          .value
          .unsafeToFuture()
        value <- metadata match {
          case Left(value) =>
            Future.failed(new RuntimeException(s"Cromwell error: $value"))
          case Right(value) => Future.successful(value.value)
        }
      } yield {
        tempFile.overwrite(value)
        val cost =
          s"${CromwellWorkflowTester.WarpRoot}/scripts/calculate_cost.py --only_total --metadata $tempFile"
            .!!
        (workflow, cost)
      }
    })
  }

  /**
    * Read tsv from a cromwell output
    *
    * @param submittedWorkflow Workflow to get output from
    * @param outputName        Name of the output to get
    * @return The contents of the cromwell output as Tsv output
    */
  protected def readTsvFromCromwellOutput(
      submittedWorkflow: SubmittedWorkflow,
      outputName: String
  ): Future[Seq[Map[String, String]]] =
    readFileFromCromwellOutput(submittedWorkflow, outputName).map(
      s => TsvParser.parse(s.split("\n"))
    )

  /**
    * Read json from a cromwell output
    *
    * @param submittedWorkflow Workflow to get output from
    * @param outputName        Name of the output to get
    * @return The contents of the cromwell output as a string
    */
  protected def readFileFromCromwellOutput(
      submittedWorkflow: SubmittedWorkflow,
      outputName: String
  ): Future[String] =
    for {
      outputs <- cromwellClient()
        .outputs(submittedWorkflow.id)
        .value
        .unsafeToFuture()
      outputPath <- Future.fromTry(
        outputs match {
          case Left(value) =>
            Failure(new RuntimeException(s"Cromwell error: $value"))
          case Right(value) =>
            value.getStringFromOutputs(s"$workflowName.$outputName")
        }
      )
    } yield ioUtil.readFile(URI.create(outputPath))

  private def checkReport(tsvLines: Seq[Map[String, String]]): Try[Unit] = {
    val errors = tsvLines.foldLeft(Seq.empty[String]) { (errs, tsvLine) =>
      val vcfPath = tsvLine("vcf")
      val markedAsSignificant = tsvLine("has_significant_variants").toBoolean
      if (vcfPath.contains("positive") && !markedAsSignificant) {
        s"VCF $vcfPath incorrectly marked as having no significant variants" +: errs
      } else if (vcfPath.contains("negative") && markedAsSignificant) {
        s"VCF $vcfPath incorrectly marked as having significant variants" +: errs
      } else {
        errs
      }
    }

    if (errors.isEmpty) {
      Success(())
    } else {
      Failure(
        new CromwellWorkflowTester.TestFailedException(
          s"Some VCFs were incorrectly classified:\n${errors.mkString("\n")}"
        )
      )
    }
  }
}
