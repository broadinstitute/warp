import com.typesafe.sbt.SbtGit.git
import Def.Initialize

name := "workflow_tests"

scalaVersion := "2.12.6"

resolvers ++= Seq(
  "Broad Artifactory Releases" at "https://broadinstitute.jfrog.io/broadinstitute/libs-release/",
  "Broad Artifactory Snapshots" at "https://broadinstitute.jfrog.io/broadinstitute/libs-snapshot/"
)

// For local clio docker compose
val elasticSearchVersion = "5.4.0_6"

val akkaVersion = "2.5.25"
val akkaHttpCirceVersion = "1.29.1"
val akkaHttpVersion = "10.1.10"
val alpakkaFtpVersion = "0.19"
val betterFilesVersion = "3.6.0"
val circeVersion = "0.12.3"
val circeGenericVersion = "0.12.2"
val clioVersion = "96c5dfb93e29e85d89bc2e3cfb5dba6a177c8255"
val commonsIoVersion = "2.5"
val cromwellVersion = "47"
val enumeratumVersion = "1.5.13"
val googleCloudNio = "0.52.0-alpha"
val httpClientVersion = "4.5.5"
val logbackClassicVersion = "1.2.3"
val mysqlJdbcDriverVersion = "5.1.46"
val oracleOjdbcDriverVersion = "12.2.0.1"
val picardPrivateVersion = "95f781b6847e80fb4dddd35320fc4235cd545b15"
val scoptVersion = "3.7.0"
val testContainersScalaVersion = "0.8.0"
val scalaFmtVersion = "1.4.0"
val scalaLoggingVersion = "3.9.0"
val vaultVersion = "3.1.0"
val ztZipVersion = "1.13"

// Test dependencies
val scalaTestVersion = "3.0.5"
val scalamockVersion = "3.6.0"

// Git Versioning
enablePlugins(GitVersioning)
lazy val gitShaVersion: Initialize[Option[String]] = Def.setting {
  val hash = git.gitHeadCommit.value match {
    case Some(sha) => s"${sha.take(7)}"
    case None      => "UNKNOWN"
  }
  Option(hash)
}
git.formattedShaVersion := gitShaVersion.value

// ScalaFmt
inThisBuild(
  Seq(
    scalafmtVersion := scalaFmtVersion,
    scalafmtOnCompile := true,
    ignoreErrors in scalafmt := false
  )
)

// SbtBuildInfo
val workflowBuildInfoKeys = Def.setting {
  Seq[BuildInfoKey](
    name,
    version,
    BuildInfoKey.constant(("elasticsearchVersion", elasticSearchVersion)),
    BuildInfoKey.constant(("clioVersion", clioVersion)),
    BuildInfoKey.constant(
      ("confDir",
       ((resourceDirectory in Compile).value / "clio").getAbsolutePath)
    ),
    BuildInfoKey.constant("cromwellVersion", cromwellVersion)
  )
}

lazy val `workflow_tests` = (project in file("."))
  .enablePlugins(BuildInfoPlugin)
  .settings(
    fork := true,
    buildInfoKeys := workflowBuildInfoKeys.value,
    buildInfoPackage := "org.broadinstitute.dsp.pipelines.util",
    buildInfoObject := "WorkflowTestBuildInfo",
    test in assembly := {},
    assemblyMergeStrategy in assembly := {
      case PathList("logback.xml") => MergeStrategy.first
      case x =>
        val oldStrategy = (assemblyMergeStrategy in assembly).value
        oldStrategy(x)
    },
    assemblyJarName in assembly := s"${name.value}-${version.value}.jar",
    libraryDependencies ++= Seq(
      "com.google.cloud" % "google-cloud-nio" % googleCloudNio,
      "ch.qos.logback" % "logback-classic" % logbackClassicVersion,
      "com.bettercloud" % "vault-java-driver" % vaultVersion,
      "com.beachape" %% "enumeratum" % enumeratumVersion,
      "com.dimafeng" %% "testcontainers-scala" % testContainersScalaVersion
        exclude ("org.slf4j", "jcl-over-slf4j"),
      "com.github.pathikrit" %% "better-files" % betterFilesVersion,
      "com.github.scopt" %% "scopt" % scoptVersion,
      "com.lightbend.akka" %% "akka-stream-alpakka-ftp" % alpakkaFtpVersion,
      "mysql" % "mysql-connector-java" % mysqlJdbcDriverVersion,
      "com.typesafe.akka" %% "akka-stream" % akkaVersion,
      "com.typesafe.scala-logging" %% "scala-logging" % scalaLoggingVersion,
      "de.heikoseeberger" %% "akka-http-circe" % akkaHttpCirceVersion,
      "io.circe" %% "circe-core" % circeVersion,
      "io.circe" %% "circe-generic" % circeVersion,
      "io.circe" %% "circe-generic-extras" % circeGenericVersion,
      "io.circe" %% "circe-literal" % circeVersion,
      "io.circe" %% "circe-parser" % circeVersion,
      "org.apache.httpcomponents" % "httpclient" % httpClientVersion,
      "org.apache.httpcomponents" % "fluent-hc" % httpClientVersion,
      "org.broadinstitute" %% "clio-client" % clioVersion,
      "org.broadinstitute" %% "cromwell-api-client" % cromwellVersion
        exclude ("com.typesafe.akka", "akka-http-spray-json"),
      "com.typesafe.akka" %% "akka-http-spray-json" % akkaHttpVersion,
      "org.zeroturnaround" % "zt-zip" % ztZipVersion
    ),
    libraryDependencies ++= Seq(
      "org.scalatest" %% "scalatest" % scalaTestVersion,
      "com.typesafe.akka" %% "akka-http-testkit" % akkaHttpVersion,
      "org.scalamock" %% "scalamock-scalatest-support" % scalamockVersion
    ).map(_ % Test),
    dependencyOverrides ++= Set(
      "com.typesafe.akka" %% "akka-actor" % akkaVersion,
      "com.typesafe.akka" %% "akka-stream" % akkaVersion,
      "commons-io" % "commons-io" % commonsIoVersion,
      "com.github.pathikrit" %% "better-files" % betterFilesVersion
    ),
    scalacOptions in Compile ++= Seq(
      "-deprecation",
      "-encoding",
      "UTF-8",
      "-explaintypes",
      "-feature",
      "-target:jvm-1.8",
      "-unchecked",
      "-Xfatal-warnings",
      "-Xfuture",
      "-Xlint",
      "-Xmax-classfile-name",
      "200",
      "-Yno-adapted-args",
      "-Ywarn-dead-code",
      "-Ywarn-extra-implicit",
      "-Ywarn-inaccessible",
      "-Ywarn-infer-any",
      "-Ywarn-nullary-override",
      "-Ywarn-nullary-unit",
      "-Ywarn-numeric-widen",
      "-Ywarn-unused",
      "-Ywarn-unused-import",
      "-Ywarn-value-discard"
    ),
    scalacOptions in (Compile, console) := (scalacOptions in Compile).value filterNot Set(
      "-Xfatal-warnings",
      "-Xlint",
      "-Ywarn-unused",
      "-Ywarn-unused-import"
    )
  )
