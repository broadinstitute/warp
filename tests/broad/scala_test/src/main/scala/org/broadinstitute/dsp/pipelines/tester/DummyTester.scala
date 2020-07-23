package org.broadinstitute.dsp.pipelines.tester

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.Resource
import cromwell.api.model.SubmittedWorkflow
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment
import org.broadinstitute.dsp.pipelines.inputs.DummyInputs

import scala.concurrent.Future
import scala.util.Failure

class DummyTester(
    implicit am: ActorMaterializer,
    as: ActorSystem,
) extends CromwellWorkflowTester {

  override val env = CromwellEnvironment.Dev

  override val workflowName: String = "DummyWorkflow"

  override val wdlContents: String =
    Resource.getAsString("dummy/dummyWorkflow.wdl")

  private val wdlInputs: List[String] =
    Resource.getAsString("dummy/dummyInputs.json").lines.toList

  private val wdlInputsDecoded: DummyInputs = DummyInputs(wdlInputs)

  override def runTest: Future[Unit] = {
    val wdlOptions = None

    logger.info(s"Submitting dummy workflow to ${env.cromwellUrl}")
    for {
      submittedWorkflow <- runWorkflow(wdlContents,
                                       wdlInputs.mkString,
                                       wdlOptions)
      _ = logger.info(s"Submitted workflow with id: ${submittedWorkflow.id}")
      completedWorkflowStatus <- awaitCromwellWorkflowCompletion(
        submittedWorkflow)
      _ = logger.info(
        s"Workflow with id ${submittedWorkflow.id} completed with status $completedWorkflowStatus"
      )
      _ = logger.info(s"Validating workflow with id: ${submittedWorkflow.id}")
      validation <- validateWorkflowResults(submittedWorkflow)
    } yield {
      validation
    }
  }

  def validateWorkflowResults(
      submittedWorkflow: SubmittedWorkflow): Future[Unit] = {
    for {
      outputs <- cromwellClient()
        .outputs(submittedWorkflow.id)
        .value
        .unsafeToFuture()
      echoed <- Future.fromTry(
        outputs match {
          case Left(value) =>
            Failure(new RuntimeException(s"Cromwell error: $value"))
          case Right(value) =>
            value.getStringFromOutputs("DummyWorkflow.echoed")
        }
      )
    } yield {
      if (!wdlInputsDecoded.dummyBoolean || wdlInputsDecoded.dummyOption.isEmpty) {
        throw new RuntimeException("The dummy inputs didn't decode correctly")
      }
      if (echoed equals "This is a dummy workflow!") {
        ()
      } else {
        throw new RuntimeException(
          "The dummy workflow didn't echo the right thing")
      }
    }
  }
}
