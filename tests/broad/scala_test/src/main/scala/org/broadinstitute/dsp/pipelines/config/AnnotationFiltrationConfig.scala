package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.{
  CromwellEnvironment,
  WorkflowTestCategory
}

case class AnnotationFiltrationConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing,
    env: CromwellEnvironment = CromwellEnvironment.Dev
)
