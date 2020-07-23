package org.broadinstitute.dsp.pipelines.commandline
import enumeratum.{Enum, EnumEntry}

import scala.collection.immutable

sealed trait PapiVersion extends EnumEntry

object PapiVersion extends Enum[PapiVersion] {

  override val values: immutable.IndexedSeq[PapiVersion] = findValues

  case object PAPI extends PapiVersion
  case object PAPIv2 extends PapiVersion
}
