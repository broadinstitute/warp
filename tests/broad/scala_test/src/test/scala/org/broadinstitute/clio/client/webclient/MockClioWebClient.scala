package org.broadinstitute.clio.client.webclient

import akka.NotUsed
import akka.http.scaladsl.model.headers.OAuth2BearerToken
import akka.http.scaladsl.model.{HttpRequest, HttpResponse}
import akka.stream.scaladsl.{Flow, Source}
import better.files.File
import io.circe.Json
import io.circe.parser._
import org.broadinstitute.clio.client.ClioClientConfig
import org.broadinstitute.clio.transfer.model.{ClioIndex, CramIndex, GvcfIndex, UbamIndex}

// FIXME: We shouldn't need this; replace it with scalamock.
class MockClioWebClient(dataResourceDir: String)
    extends ClioWebClient(
      Flow[HttpRequest].map(_ => HttpResponse()),
      ClioClientConfig.responseTimeout,
      ClioClientConfig.maxRequestRetries,
      () => OAuth2BearerToken("fake")
    ) {

  override def simpleQuery[I](clioIndex: ClioWebClient.QueryAux[I])(
    input: I,
    includeDeleted: Boolean
  ): Source[Json, NotUsed] = {
    clioIndex match {
      case CramIndex =>
        parse(File(pathMapper(clioIndex)(input)).contentAsString)
          .fold(Source.failed, Source.single)
      case _ => sys.error("Unexpected clio query")
    }
  }

  private def pathMapper[CI <: ClioIndex](
    clioIndex: CI
  ): clioIndex.QueryInputType => String = { qi =>
    clioIndex match {
      case CramIndex =>
        val queryInput = qi.asInstanceOf[CramIndex.QueryInputType]
        s"$dataResourceDir/${queryInput.project.get}_${queryInput.sampleAlias.get}_${queryInput.version.get}.json"
      case UbamIndex =>
        val queryInput = qi.asInstanceOf[UbamIndex.QueryInputType]
        s"$dataResourceDir/${queryInput.flowcellBarcode.get}_${queryInput.lane.get}_${queryInput.libraryName.get}.json"
      case GvcfIndex =>
        val queryInput = qi.asInstanceOf[GvcfIndex.QueryInputType]
        s"$dataResourceDir/${queryInput.project.get}_${queryInput.sampleAlias.get}_${queryInput.version.get}.json"
    }

  }
}
