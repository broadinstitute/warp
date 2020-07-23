val assemblyVersion = "0.14.7"
val sbtScalafmtVersion = "1.15"
val sbtBuildInfoVersion = "0.7.0"
val sbtGitVersion = "0.9.3"

addSbtPlugin("com.eed3si9n" % "sbt-assembly" % assemblyVersion)
addSbtPlugin("com.lucidchart" % "sbt-scalafmt" % sbtScalafmtVersion)
addSbtPlugin("com.eed3si9n" % "sbt-buildinfo" % sbtBuildInfoVersion)
addSbtPlugin("com.typesafe.sbt" % "sbt-git" % sbtGitVersion)

// Silence a build warning.
val Slf4jVersion = "1.7.25"

libraryDependencies += "org.slf4j" % "slf4j-nop" % Slf4jVersion
