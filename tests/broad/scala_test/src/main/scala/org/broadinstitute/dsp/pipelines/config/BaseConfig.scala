package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.{
  PapiVersion,
  CromwellEnvironment,
  WorkflowTestCategory
}

/**
  * This exists so that GermlineCloudWorkflowTester can accept either a GermlineCloudWorkflowConfig or an
  * AllOfUsConfig.
  */
trait BaseConfig {
  def category: WorkflowTestCategory
  def truthBranch: String
  def useTimestamp: Option[String]
  def useCallCaching: Boolean
  def updateTruth: Boolean
  def papiVersion: PapiVersion
  def env: CromwellEnvironment
}
