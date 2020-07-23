package org.broadinstitute.dsp.pipelines.firecloud.model

/**
  * These files are used in conjunction with the edited swagger-codegened files.
  * They make circe parsing and decoding much easier.
  */
trait FireCloudEntity {
  val name: String
  val entityType: String
  val attributes: Map[String, String]
}
