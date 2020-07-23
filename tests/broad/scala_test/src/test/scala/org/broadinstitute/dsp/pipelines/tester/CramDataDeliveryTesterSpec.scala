package org.broadinstitute.dsp.pipelines.tester

import java.time.ZoneId
import java.util.UUID

import better.files.Resource
import cats.data.EitherT
import cromwell.api.CromwellClient
import cromwell.api.model._
import io.circe.Decoder
import io.circe.parser._
import org.broadinstitute.clio.client.webclient.MockClioWebClient
import org.broadinstitute.dsp.pipelines.TestKitSuite
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment
import org.broadinstitute.dsp.pipelines.config.CramDataDeliveryConfig
import org.broadinstitute.dsp.pipelines.file.TsvParser
import org.broadinstitute.dsp.pipelines.firecloud.FireCloudClient
import org.broadinstitute.dsp.pipelines.firecloud.model.autogen.{Entity, Workspace}
import org.broadinstitute.dsp.pipelines.firecloud.model.{
  SampleEntity,
  SampleSetEntity,
  WorkspaceResponse
}
import org.broadinstitute.dsp.pipelines.tester.mock.MockCramDataDeliveryTester
import org.broadinstitute.dsp.pipelines.util.DateFormatters.{
  FireCloudWorkspaceDateFormatters => DateFormatters
}
import org.broadinstitute.dsp.pipelines.util.FireCloudDecoders
import org.scalamock.scalatest.AsyncMockFactory
import spray.json.{JsObject, JsString}

import scala.concurrent.{ExecutionContext, Future}
import scala.util.Success

class CramDataDeliveryTesterSpec
    extends TestKitSuite("workflow-test")
    with AsyncMockFactory {

  import FireCloudDecoders._

  val submittedWorkflow = SubmittedWorkflow(
    WorkflowId(UUID.randomUUID()),
    CromwellEnvironment.Dev.cromwellUrl,
    WorkflowSingleSubmission(None, None, None, None, None, None, None, None, None)
  )

  def transform(line: String, diff: String): String = {
    line
      .replace("{DIFF}", diff)
      .replace("{FC_BUCKET}", s"fc-$diff")
      .replace("{WORKSPACE_DATE}", DateFormatters.dateForTags(ZoneId.systemDefault()))
      .replace("{DATE}", DateFormatters.dateForTsv(ZoneId.systemDefault()))
  }

  behavior of "cramDataDeliveryTester"

  it should "check workspace properties correctly" in {
    val tester = getMockTester
    val workspaceResponse = parse(
      Resource
        .getAsString("cramDataDelivery/workspaceResponse.json")
        .lines
        .map(transform(_, tester.differentiator))
        .mkString
    ).flatMap(_.as[WorkspaceResponse])
      .getOrElse(sys.error("Error decoding json workspace response"))

    tester.checkWorkspaceProperties(workspaceResponse.workspace) should be(Success(()))
  }

  it should "check samples correctly" in {
    val tester = getMockTester
    val tsvSamples =
      TsvParser.parseFromInput(
        (tester.inputsDir / "expected" / "samples.tsv"),
        transform(_, tester.differentiator)
      )
    val jsonSamples =
      parse(
        Resource
          .getAsString("cramDataDelivery/samples.json")
          .lines
          .map(transform(_, tester.differentiator))
          .mkString
      ).flatMap(_.as[Seq[SampleEntity]])
        .getOrElse(sys.error("Error decoding json samples"))

    tester.checkEntities(jsonSamples, tsvSamples, "entity:sample_id") should
      be(Success(()))
  }

  it should "check participants correctly" in {
    val tester = getMockTester
    val tsvParticipants =
      TsvParser.parseFromInput(
        (tester.inputsDir / "expected" / "participants.tsv"),
        transform(_, tester.differentiator)
      )
    val jsonParticipants =
      parse(
        Resource
          .getAsString("cramDataDelivery/participants.json")
          .lines
          .map(transform(_, tester.differentiator))
          .mkString
      ).flatMap(_.as[Seq[Entity]])
        .getOrElse(sys.error("Error decoding json participants"))

    tester.checkEntities(jsonParticipants, tsvParticipants, "entity:participant_id") should
      be(Success(()))
  }

  it should "check samples sets correctly" in {
    val tester = getMockTester
    val sampleSetEntities = TsvParser.parseFromInput(
      (tester.inputsDir / "expected" / "sample_set_entity.tsv"),
      transform(_, tester.differentiator)
    )
    val sampleSetMembership = TsvParser.parseFromInput(
      (tester.inputsDir / "expected" / "sample_set_membership.tsv"),
      transform(_, tester.differentiator)
    )
    val jsonSampleSets = parse(
      Resource
        .getAsString("cramDataDelivery/sampleSets.json")
        .lines
        .map(transform(_, tester.differentiator))
        .mkString
    ).flatMap(_.as[Seq[SampleSetEntity]])
      .getOrElse(sys.error("Error decoding json participants"))

    tester.checkSampleSets(jsonSampleSets, sampleSetEntities, sampleSetMembership) should
      be(Success(()))
  }

  it should "kick off a workflow and verify it" in {
    val tester = getMockTester
    tester.runTest.map(_ should be(()))
  }

  // Need to create a new tester for each test.
  // See https://github.com/paulbutcher/ScalaMock/issues/117
  def getMockTester: CramDataDeliveryTester = {
    val cramDataDeliveryConfig = CramDataDeliveryConfig()

    val uuid = UUID.randomUUID().toString

    def setupFireCloudResponse[A: Decoder](fileName: String): Future[A] =
      Future.fromTry(
        parse(
          Resource
            .getAsString(s"cramDataDelivery/$fileName")
            .lines
            .map(transform(_, uuid))
            .mkString
        ).flatMap(_.as[A]).toTry
      )

    val mockFireCloudClient = stub[FireCloudClient]

    (mockFireCloudClient.createWorkspace _)
      .when()
      .returns(setupFireCloudResponse[Workspace]("workspace.json"))
    (mockFireCloudClient.getWorkspace _)
      .when()
      .returns(setupFireCloudResponse[WorkspaceResponse]("workspaceResponse.json"))
    (mockFireCloudClient.getWorkspaceSampleSets _)
      .when()
      .returns(setupFireCloudResponse[Seq[SampleSetEntity]]("sampleSets.json"))
    (mockFireCloudClient.getWorkspaceSamples _)
      .when()
      .returns(setupFireCloudResponse[Seq[SampleEntity]]("samples.json"))
    (mockFireCloudClient.getWorkspaceParticipants _)
      .when()
      .returns(setupFireCloudResponse[Seq[Entity]]("participants.json"))

    val mockCromwellClient: CromwellClient = stub[CromwellClient]

    (mockCromwellClient
      .submit(_: WorkflowSingleSubmission)(_: ExecutionContext))
      .when(*, *)
      .returns(EitherT.rightT(submittedWorkflow))
    (mockCromwellClient
      .status(_: WorkflowId)(_: ExecutionContext))
      .when(*, *)
      .returning(EitherT.rightT(Succeeded))
    (mockCromwellClient
      .abort(_: WorkflowId)(_: ExecutionContext))
      .when(*, *)
      .returns(EitherT.rightT(Aborting))
    (mockCromwellClient
      .outputs(_: WorkflowId)(_: ExecutionContext))
      .when(*, *)
      .returns(
        EitherT.rightT(
          WorkflowOutputs(
            submittedWorkflow.id.toString(),
            JsObject(
              Map(
                "WholeGenomeCramDataDelivery.workspace_bucket" ->
                  JsString(s"fc-$uuid")
              )
            )
          )
        )
      )

    val mockClioWebClient = new MockClioWebClient("cramDataDelivery/clio")

    new MockCramDataDeliveryTester(
      cramDataDeliveryConfig,
      mockCromwellClient,
      mockFireCloudClient,
      mockClioWebClient,
      uuid
    )
  }

}
