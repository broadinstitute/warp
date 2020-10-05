package org.broadinstitute.dsp.pipelines.commandline

import java.net.URL

import enumeratum.{Enum, EnumEntry}
import io.circe.Json
import io.circe.syntax._

import scala.collection.immutable

sealed trait CromwellEnvironment extends EnumEntry {
  def cromwellUrl: URL
  def environmentOptions: Seq[(String, Json)] = Seq.empty
  def picardEnv: String
}

object CromwellEnvironment extends Enum[CromwellEnvironment] {
  override val values: immutable.IndexedSeq[CromwellEnvironment] = findValues

  def optionsString: String =
    s"[${CromwellEnvironment.lowerCaseNamesToValuesMap.keys.mkString("|")}]"

  case object Dev extends CromwellEnvironment {
    override val cromwellUrl = new URL(
      "https://cromwell.gotc-dev.broadinstitute.org")
    override val picardEnv: String = "dev"
  }

  case object Staging extends CromwellEnvironment {
    override val cromwellUrl = new URL(
      "https://cromwell.gotc-staging.broadinstitute.org")
    override val picardEnv: String = "staging"
  }

  case object Prod extends CromwellEnvironment {
    override val cromwellUrl = new URL(
      "https://cromwell.gotc-prod.broadinstitute.org")
    override val picardEnv: String = "prod"
  }

  case object Test extends CromwellEnvironment {
    override val cromwellUrl = new URL(
      "https://cromwell-test.gotc-dev.broadinstitute.org"
    )
    override val picardEnv: String = "dev"
    override val environmentOptions: Seq[(String, Json)] = Seq(
      "jes_gcs_root" -> "gs://broad-gotc-dev-execution1".asJson,
      "google_project" -> "broad-exomes-dev1".asJson,
    )
  }

  case object Pharma5 extends CromwellEnvironment {
    override val cromwellUrl = new URL(
      "https://cromwell-pharma5.gotc-prod.broadinstitute.org"
    )
    override val environmentOptions: Seq[(String, Json)] = Seq(
      "jes_gcs_root" -> "gs://broad-pharma5-execution1/".asJson,
      "google_project" -> "broad-pharma5-compute1".asJson,
    )
    override val picardEnv: String = "prod"
  }
  case object JGDev extends CromwellEnvironment {
    override val cromwellUrl = new URL(
      "https://cromwell-jg.gotc-dev.broadinstitute.org")
    override val picardEnv: String = "dev"
    override val environmentOptions: Seq[(String, Json)] = Seq(
      "jes_gcs_root" -> "gs://broad-gotc-dev-execution1".asJson,
      "google_project" -> "broad-exomes-dev1".asJson,
    )
  }

  case object JGProd extends CromwellEnvironment {
    override val cromwellUrl = new URL(
      "https://cromwell-jg.gotc-prod.broadinstitute.org")
    override val picardEnv: String = "prod"
    override val environmentOptions: Seq[(String, Json)] = Seq(
      "jes_gcs_root" -> "gs://broad-gotc-prod-execution1".asJson,
      "google_project" -> "broad-exomes-prod1".asJson,
    )
  }
}
