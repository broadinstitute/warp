package org.broadinstitute.dsp.pipelines.commandline
import enumeratum.{Enum, EnumEntry}

import scala.collection.immutable

sealed trait WorkflowTestCategory extends EnumEntry

object WorkflowTestCategory extends Enum[WorkflowTestCategory] {
  override val values: immutable.IndexedSeq[WorkflowTestCategory] = findValues

  case object Plumbing extends WorkflowTestCategory
  case object Scientific extends WorkflowTestCategory
  case object Load extends WorkflowTestCategory
}
