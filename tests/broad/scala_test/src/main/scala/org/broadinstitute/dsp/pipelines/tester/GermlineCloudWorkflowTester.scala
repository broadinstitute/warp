package org.broadinstitute.dsp.pipelines.tester

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import org.broadinstitute.dsp.pipelines.batch.WorkflowRunParameters
import org.broadinstitute.dsp.pipelines.util.DataType
import org.broadinstitute.dsp.pipelines.config.GermlineCloudWorkflowConfig

abstract class GermlineCloudWorkflowTester(
    testerConfig: GermlineCloudWorkflowConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  /**
    * @return The WorkflowRunParameters to be used in this test
    */
  protected def generateRunParameters: Seq[WorkflowRunParameters]

  protected val dataTypePrefix: String = dataTypePrefix(testerConfig.dataType)

  protected val dataTypeString: String =
    testerConfig.dataType.entryName.toLowerCase

  protected def dataType: DataType = testerConfig.dataType
}
