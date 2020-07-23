package org.broadinstitute.dsp.pipelines.util

import java.net.URI
import java.time.OffsetDateTime
import java.util.UUID

import io.circe.generic.AutoDerivation
import io.circe.{Decoder, Encoder}

trait CirceDecoders extends AutoDerivation {

  implicit val encodeUri: Encoder[URI] =
    Encoder.encodeString.contramap(_.toString)

  implicit val decodeUri: Decoder[URI] =
    Decoder.decodeString.map(URI.create)

  implicit val encodeUUID: Encoder[UUID] =
    Encoder.encodeString.contramap(_.toString)

  implicit val decodeUUID: Decoder[UUID] =
    Decoder.decodeString.map(UUID.fromString)

  implicit val encodeOffsetDateTime: Encoder[OffsetDateTime] =
    Encoder.encodeString.contramap(_.toString)

  implicit val decodeOffsetDateTime: Decoder[OffsetDateTime] =
    Decoder.decodeString.map(OffsetDateTime.parse)
}
