package org.broadinstitute.dsp.pipelines.inputs

trait FireCloudWdlInputs {

  def billingProject: String
  def requester: String
  def workspaceName: String
}
