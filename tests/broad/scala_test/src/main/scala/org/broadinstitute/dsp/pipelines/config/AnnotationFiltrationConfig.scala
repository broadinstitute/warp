package org.broadinstitute.dsp.pipelines.config

import org.broadinstitute.dsp.pipelines.commandline.WorkflowTestCategory

case class AnnotationFiltrationConfig(
    category: WorkflowTestCategory = WorkflowTestCategory.Plumbing
)
