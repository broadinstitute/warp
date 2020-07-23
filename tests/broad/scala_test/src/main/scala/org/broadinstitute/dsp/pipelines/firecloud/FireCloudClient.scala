package org.broadinstitute.dsp.pipelines.firecloud

import java.net.URL

import akka.actor.ActorSystem
import akka.http.scaladsl.Http
import akka.http.scaladsl.coding.{Deflate, Gzip, NoCoding}
import akka.http.scaladsl.model._
import akka.http.scaladsl.model.headers.HttpEncodings
import akka.http.scaladsl.unmarshalling.{Unmarshal, Unmarshaller}
import akka.stream.ActorMaterializer
import cromwell.api.CromwellClient.UnsuccessfulRequestException
import de.heikoseeberger.akkahttpcirce.FailFastCirceSupport
import io.circe.Json
import io.circe.syntax._
import org.broadinstitute.clio.client.webclient.GoogleCredentialsGenerator
import org.broadinstitute.clio.util.auth.ClioCredentials
import org.broadinstitute.dsp.pipelines.firecloud.model._
import org.broadinstitute.dsp.pipelines.firecloud.model.autogen._
import org.broadinstitute.dsp.pipelines.util.FireCloudDecoders

import scala.concurrent.Future
import scala.util.{Failure, Try}

class FireCloudClient(
    url: URL,
    workspaceNamespace: String,
    workspaceName: String,
    credentials: ClioCredentials
)(
    implicit as: ActorSystem,
    am: ActorMaterializer,
) {

  import FailFastCirceSupport._
  import FireCloudDecoders._
  import as.dispatcher

  lazy val credentialsGenerator = GoogleCredentialsGenerator(credentials)

  def getWorkspace: Future[WorkspaceResponse] =
    httpGetRequest[WorkspaceResponse](existingWorkspaceEndpoint)

  def createWorkspace: Future[Workspace] =
    httpPostRequest[Workspace](workspaceEndpoint)

  def getWorkspaceAcl: Future[WorkspaceACL] =
    httpGetRequest[WorkspaceACL](workspaceAclEndpoint)

  def getWorkspaceCatalogers: Future[Seq[WorkspaceCatalog]] =
    httpGetRequest[Seq[WorkspaceCatalog]](workspaceCatalogEndpoint)

  def getWorkspaceTags: Future[Seq[String]] =
    httpGetRequest[Seq[String]](workspaceTagsEndpoint)

  def deleteWorkspace: Future[Unit] =
    httpDeleteRequest(existingWorkspaceEndpoint)

  def getWorkspaceParticipants: Future[Seq[Entity]] =
    httpGetRequest[Seq[Entity]](workspaceEntityWithTypeEndpoint("participant"))

  def getWorkspaceSamples: Future[Seq[SampleEntity]] =
    httpGetRequest[Seq[SampleEntity]](workspaceEntityWithTypeEndpoint("sample"))

  def getWorkspaceSampleSets: Future[Seq[SampleSetEntity]] =
    httpGetRequest[Seq[SampleSetEntity]](
      workspaceEntityWithTypeEndpoint("sample_set"))

  def workspaceEntityWithTypeEndpoint(entityType: String): String =
    s"$existingWorkspaceEndpoint/entities/$entityType"

  def workspaceAclEndpoint = s"$existingWorkspaceEndpoint/acl"

  def workspaceCatalogEndpoint = s"$existingWorkspaceEndpoint/catalog"

  def existingWorkspaceEndpoint =
    s"$url/api/workspaces/$workspaceNamespace/$workspaceName"

  def workspaceEndpoint =
    s"$url/api/workspaces"

  def workspaceTagsEndpoint = s"$existingWorkspaceEndpoint/tags"

  private def httpGetRequest[A](uri: String)(
      implicit um: Unmarshaller[HttpResponse, A]
  ): Future[A] = {
    makeRequest[A](
      HttpRequest(uri = uri)
        .addCredentials(credentialsGenerator.generateCredentials())
    )
  }

  private def httpPostRequest[A](uri: String)(
      implicit um: Unmarshaller[HttpResponse, A]
  ): Future[A] =
    makeRequest[A](
      HttpRequest(
        method = HttpMethods.POST,
        uri = uri,
        entity = HttpEntity(
          ContentTypes.`application/json`,
          Json
            .obj(
              "namespace" -> workspaceNamespace.asJson,
              "name" -> workspaceName.asJson,
              "attributes" -> Json.obj(),
              "authorizationDomain" -> Json.arr()
            )
            .spaces2
        )
      ).addCredentials(credentialsGenerator.generateCredentials())
    )

  private def httpDeleteRequest(uri: String): Future[Unit] = {
    val request =
      HttpRequest(
        method = HttpMethods.DELETE,
        uri = uri
      ).addCredentials(credentialsGenerator.generateCredentials())
    executeRequest(request).flatMap { response =>
      if (!response.status.isSuccess()) {
        Future.failed(
          UnsuccessfulRequestException(s"Failed to DELETE $uri", response))
      } else {
        Future.unit
      }
    }
  }

  private def makeRequest[A](request: HttpRequest)(
      implicit um: Unmarshaller[HttpResponse, A]
  ): Future[A] =
    for {
      response <- executeRequest(request)
      decoded <- Future.fromTry(decodeResponse(response))
      unmarshalled <- Unmarshal(decoded).to[A].recover {
        case ex =>
          val msg =
            s"""Failed to parse response as JSON. Response was:
             |$decoded""".stripMargin
          throw new RuntimeException(msg, ex)
      }
    } yield unmarshalled

  private def executeRequest(
      httpRequest: HttpRequest
  ): Future[HttpResponse] =
    Http().singleRequest(httpRequest)

  private val decoders =
    Map(
      HttpEncodings.gzip -> Gzip,
      HttpEncodings.deflate -> Deflate,
      HttpEncodings.identity -> NoCoding
    )

  private def decodeResponse(response: HttpResponse): Try[HttpResponse] = {
    decoders.get(response.encoding) map { decoder =>
      Try(decoder.decodeMessage(response))
    } getOrElse Failure(
      UnsuccessfulRequestException(s"No decoder for ${response.encoding}",
                                   response)
    )
  }
}
