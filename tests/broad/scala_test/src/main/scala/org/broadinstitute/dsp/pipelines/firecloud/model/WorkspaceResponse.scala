package org.broadinstitute.dsp.pipelines.firecloud.model

import org.broadinstitute.dsp.pipelines.firecloud.model.autogen.Workspace

/**
  * These files are used in conjunction with the edited swagger-codegened files.
  * They make circe parsing and decoding much easier.
  */
case class WorkspaceResponse(
    catalog: Boolean,
    workspaceSubmissionStats: WorkspaceSubmissionStats,
    accessLevel: String,
    owners: Seq[String],
    canShare: Boolean,
    canCompute: Boolean,
    workspace: Workspace
)
