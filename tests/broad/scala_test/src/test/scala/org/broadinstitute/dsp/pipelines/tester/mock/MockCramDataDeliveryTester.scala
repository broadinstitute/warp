package org.broadinstitute.dsp.pipelines.tester.mock

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import cromwell.api.CromwellClient
import cromwell.api.model.SubmittedWorkflow
import org.broadinstitute.clio.client.webclient.ClioWebClient
import org.broadinstitute.dsp.pipelines.config.CramDataDeliveryConfig
import org.broadinstitute.dsp.pipelines.firecloud.FireCloudClient
import org.broadinstitute.dsp.pipelines.tester.CramDataDeliveryTester

import scala.concurrent.Future

class MockCramDataDeliveryTester(
  config: CramDataDeliveryConfig,
  mockCromwellClient: CromwellClient,
  mockFireCloudClient: FireCloudClient,
  mockClioWebClient: ClioWebClient,
  mockDiff: String
)(implicit am: ActorMaterializer, as: ActorSystem)
    extends CramDataDeliveryTester(config) {

  override lazy val differentiator: String = mockDiff

  override lazy val clioWebClient: ClioWebClient = mockClioWebClient

  override lazy val fireCloudClient: FireCloudClient = mockFireCloudClient

  override def cleanup(): Future[Unit] = Future.unit

  override def addAbortTerminationHook(submitted: SubmittedWorkflow): Unit = ()

  override def getCromwellClient: CromwellClient = mockCromwellClient

  override def awaitCromwellWorkflowCompletion(
    submission: SubmittedWorkflow
  ): Future[Unit] = Future.unit
}
