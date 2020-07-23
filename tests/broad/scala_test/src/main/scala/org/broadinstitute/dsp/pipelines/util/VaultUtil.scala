package org.broadinstitute.dsp.pipelines.util

import java.util.{Base64, Map => JavaMap}

import com.bettercloud.vault.{Vault, VaultConfig}
import io.circe.Json
import io.circe.syntax._
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment

import scala.collection.JavaConverters._

object VaultUtil {

  private val config: VaultConfig = new VaultConfig()
    .address("https://clotho.broadinstitute.org:8200")

  def getDatabaseCredentials(
      path: String,
      token: String
  ): (String, String) = {
    val data = getVaultData(path, token)
    (data.get("username"), data.get("password"))
  }

  def getDatabaseKeystores(
      path: String,
      token: String
  ): (Array[Byte], Array[Byte]) = {
    val data = getVaultData(path, token)
    val keystore = Base64.getDecoder.decode(data.get("keystore").trim())
    val truststore = Base64.getDecoder.decode(data.get("truststore").trim())
    (keystore, truststore)
  }

  def getVaultData(
      path: String,
      token: String
  ): JavaMap[String, String] = {
    val vault = new Vault(config.token(token).build())
    vault.logical().read(path).getData
  }

  def getPicardServiceAccount(environment: CromwellEnvironment,
                              token: String): Json = {
    getVaultData(
      s"secret/dsde/gotc/${environment.picardEnv}/picard/picard-account.pem",
      token
    ).asScala.asJson
  }
}
