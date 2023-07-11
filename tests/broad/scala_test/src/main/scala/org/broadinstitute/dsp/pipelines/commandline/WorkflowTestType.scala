package org.broadinstitute.dsp.pipelines.commandline

import enumeratum._

import scala.collection.immutable

sealed trait WorkflowTestType extends EnumEntry

object WorkflowTestType extends Enum[WorkflowTestType] {

  override val values: immutable.IndexedSeq[WorkflowTestType] = findValues
  case object CloudWorkflow extends WorkflowTestType
  case object Dummy extends WorkflowTestType
  case object GermlineSingleSample extends WorkflowTestType
}
