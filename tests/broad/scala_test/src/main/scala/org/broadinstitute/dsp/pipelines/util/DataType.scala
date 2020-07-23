package org.broadinstitute.dsp.pipelines.util

import enumeratum.{Enum, EnumEntry}

import scala.collection.immutable

sealed trait DataType extends EnumEntry

object DataType extends Enum[DataType] {
  override val values: immutable.IndexedSeq[DataType] = findValues

  case object WGS extends DataType
  case object Exome extends DataType
  case object RNA extends DataType
  case object Targeted extends DataType
}
