package org.broadinstitute.dsp.pipelines.util

import enumeratum.{Enum, EnumEntry}

import scala.collection.immutable

sealed trait ArrayType extends EnumEntry {
  def pathName: String
  def resourceName: String
  def workflowBaseName: String
}

object ArrayType extends Enum[ArrayType] {
  override val values: immutable.IndexedSeq[ArrayType] = findValues

  case object Single extends ArrayType {
    override def pathName: String = "single_sample"
    override def resourceName: String = "singleSample"
    override def workflowBaseName: String = "Arrays"
  }
  case object Multi extends ArrayType {
    override def pathName: String = "multi_sample"
    override def resourceName: String = "multiSample"
    override def workflowBaseName: String = "MultiSampleArrays"
  }
}
