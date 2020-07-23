package org.broadinstitute.dsp.pipelines.firecloud.model

/**
  * These files are used in conjunction with the edited swagger-codegened files.
  * They make circe parsing and decoding much easier.
  */
case class SampleEntity(
    name: String,
    entityType: String,
    participant: EntityReference,
    attributes: Map[String, String]
) extends FireCloudEntity
