package org.broadinstitute.dsp.pipelines.util

import better.files.Resource
import cromwell.api.model.WorkflowMetadata
import org.scalatest.{FlatSpec, Matchers}

class CromwellUtilsTest extends FlatSpec with Matchers with CromwellUtils {

  behavior of "CromwellUtilsTest.FromWorkflowMetadata implicit class"

  it should "properly retrieve all logs from metadata" in {
    val metadata =
      WorkflowMetadata(Resource.getAsString("cromwellMetadata.json"))
    metadata.getLogs.size should be(579)

  }

}
