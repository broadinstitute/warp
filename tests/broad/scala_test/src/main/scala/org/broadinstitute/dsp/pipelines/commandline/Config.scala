package org.broadinstitute.dsp.pipelines.commandline

import org.broadinstitute.dsp.pipelines.config._

case class Config(
    test: WorkflowTestType = WorkflowTestType.Dummy,
    germlineCloudConfig: GermlineCloudWorkflowConfig = GermlineCloudWorkflowConfig(),
    cloudWorkflowConfig: CloudWorkflowConfig = CloudWorkflowConfig()
)
