package org.broadinstitute.dsp.pipelines.batch

import java.net.URI

/**
  * A case class to represent a run of one sample through cromwell
  *
  * @param id ID to use in logging
  * @param workflowInputs The inputs to the workflow
  * @param resultsCloudPath The cloud location where the results should be moved to after the workflow finishes
  * @param truthCloudPath The cloud location where the truth data exists
  */
case class WorkflowRunParameters(
    id: String,
    workflowInputs: String,
    resultsCloudPath: URI,
    truthCloudPath: URI
) {

  override def toString: String = s"Workflow with ID: $id"
}
