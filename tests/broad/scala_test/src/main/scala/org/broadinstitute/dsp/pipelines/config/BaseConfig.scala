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
  def env: CromwellEnvironment
  def truthBranch: String
  def updateTruth: Boolean
  def useTimestamp: Option[String]
  def useCallCaching: Boolean
  def papiVersion: PapiVersion
}
