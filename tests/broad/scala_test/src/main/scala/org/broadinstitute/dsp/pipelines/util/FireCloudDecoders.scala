package org.broadinstitute.dsp.pipelines.util

import java.time.OffsetDateTime
import java.util.UUID

import io.circe.{Decoder, HCursor}
import org.broadinstitute.dsp.pipelines.firecloud.model._
import org.broadinstitute.dsp.pipelines.firecloud.model.autogen._

object FireCloudDecoders extends CirceDecoders {

  //TODO Could probably combine sample and sample set decoders, but focusing on getting
  // a working test suite now
  implicit val sampleEntityDecoder: Decoder[SampleEntity] =
    new Decoder[SampleEntity] {

      def apply(c: HCursor): Decoder.Result[SampleEntity] =
        for {
          name <- c.downField("name").as[String]
          entityType <- c.downField("entityType").as[String]
          participant <- c
            .downField("attributes")
            .downField("participant")
            .as[EntityReference]
          otherFields <- c
            .downField("attributes")
            .downField("participant")
            .delete
            .as[Map[String, String]]
        } yield {
          SampleEntity(name, entityType, participant, otherFields)
        }
    }

  implicit val sampleSetEntityDecoder: Decoder[SampleSetEntity] =
    new Decoder[SampleSetEntity] {
      override def apply(c: HCursor): Decoder.Result[SampleSetEntity] = {
        for {
          name <- c.downField("name").as[String]
          entityType <- c.downField("entityType").as[String]
          samples <- c
            .downField("attributes")
            .downField("samples")
            .as[AttributeValue[Seq[EntityReference]]]
          otherFields <- c
            .downField("attributes")
            .downField("samples")
            .delete
            .as[Map[String, String]]
        } yield {
          SampleSetEntity(name, entityType, samples.items.toSet, otherFields)
        }
      }
    }

  implicit val workspaceDecoder: Decoder[Workspace] = new Decoder[Workspace] {

    def apply(c: HCursor): Decoder.Result[Workspace] =
      for {
        workspaceId <- c.downField("workspaceId").as[UUID]
        name <- c.downField("name").as[String]
        isLocked <- c.downField("isLocked").as[Boolean]
        lastModified <- c.downField("lastModified").as[OffsetDateTime]
        createdBy <- c.downField("createdBy").as[String]
        bucketName <- c.downField("bucketName").as[String]
        namespace <- c.downField("namespace").as[String]
        authorizationDomain <- c
          .downField("authorizationDomain")
          .as[Seq[ManagedGroupRef]]
        createdDate <- c.downField("createdDate").as[OffsetDateTime]
        attributesJson <- c
          .downField("attributes")
          .focus
          .toRight(sys.error("Error getting json from workspace attributes"))
      } yield {
        val attributes = attributesJson.asObject
          .map { obj =>
            obj.toList.toMap
          }
          .getOrElse(sys.error("Could not turn attributes list into Map"))
        Workspace(
          namespace = namespace,
          name = name,
          authorizationDomain = authorizationDomain,
          workspaceId = workspaceId,
          bucketName = bucketName,
          createdDate = createdDate,
          lastModified = lastModified,
          createdBy = createdBy,
          isLocked = isLocked,
          attributes = attributes
        )
      }
  }
}
