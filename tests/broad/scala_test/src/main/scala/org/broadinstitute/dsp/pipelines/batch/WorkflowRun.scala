package org.broadinstitute.dsp.pipelines.batch

import java.net.URI

import cromwell.api.model.{SubmittedWorkflow, WorkflowStatus}

/**
  * A trait to represent the status of workflows.
  * This allows for the batch running of workflows.
  * @tparam A The type parameter for denoting which type of workflow run this is.
  *           The run could be the 'test' portion, or the 'validation' portion of a workflow test.
  *           The single sample and arrays workflows follow this pattern, where the workflow is run,
  *           an then another workflow is used to validate the outputs of the first.
  */
sealed trait WorkflowRun[A <: WorkflowRun[_]] {
  val workflow: SubmittedWorkflow
  val workflowStatus: WorkflowStatus

  def getRunParameters: WorkflowRunParameters
  def withWorkflow(workflow: SubmittedWorkflow): A
  def withWorkflowStatus(workflowStatus: WorkflowStatus): A
}

object WorkflowTest {

  def fromSubmittedWorkflow(
      submittedWorkflow: SubmittedWorkflow,
      status: WorkflowStatus
  ): WorkflowTest = {
    WorkflowTest(
      WorkflowRunParameters(
        submittedWorkflow.id.id.toString,
        submittedWorkflow.workflow.inputsJson.getOrElse(""),
        URI.create(""),
        URI.create("")
      ),
      submittedWorkflow,
      status
    )
  }
}

/**
  * A representation of the 'test' portion of a workflow run.
  * This is the workflow being tested, like the single sample or arrays workflow.
  * @param runParameters The initial parameters for this run. Contains information like ID and results cloud path.
  * @param workflow The SubmittedWorkflow that is sent to Cromwell
  * @param workflowStatus The status of the workflow. This is updated by copying the case class with a new status.
  */
case class WorkflowTest(
    runParameters: WorkflowRunParameters,
    workflow: SubmittedWorkflow,
    workflowStatus: WorkflowStatus
) extends WorkflowRun[WorkflowTest] {
  override def getRunParameters: WorkflowRunParameters = runParameters
  override def withWorkflow(workflow: SubmittedWorkflow): WorkflowTest =
    copy(workflow = workflow)
  override def withWorkflowStatus(
      workflowStatus: WorkflowStatus
  ): WorkflowTest = copy(workflowStatus = workflowStatus)
}

/**
  * A representation of the 'validation' portion of a workflow run.
  * This is the validation workflow that is testing the original workflow.
  * @param workflowTest The WorkflowTest that this WorkflowValidation is validating.
  * @param workflow The validation workflow that was submitted to Cromwell
  * @param workflowStatus The status of the validation workflow. Like the WorkflowTest, this is updated by
  *                       copying the case class with a new status.
  */
case class WorkflowValidation(
    workflowTest: WorkflowTest,
    workflow: SubmittedWorkflow,
    workflowStatus: WorkflowStatus
) extends WorkflowRun[WorkflowValidation] {
  override def getRunParameters: WorkflowRunParameters =
    workflowTest.getRunParameters
  override def withWorkflow(workflow: SubmittedWorkflow): WorkflowValidation =
    copy(workflow = workflow)
  override def withWorkflowStatus(
      workflowStatus: WorkflowStatus
  ): WorkflowValidation = copy(workflowStatus = workflowStatus)
}
