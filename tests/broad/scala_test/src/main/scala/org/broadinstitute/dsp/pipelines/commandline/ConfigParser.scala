package org.broadinstitute.dsp.pipelines.commandline

import java.net.URL

import better.files.File
import enumeratum.{Enum, EnumEntry}
import org.broadinstitute.dsp.pipelines.util.DataType
import org.broadinstitute.dsp.pipelines.commandline.WorkflowTestType._
import org.broadinstitute.dsp.pipelines.config.GermlineCloudWorkflowConfig
import org.broadinstitute.dsp.pipelines.util.WorkflowTestBuildInfo
import scopt.Read.reads
import scopt.{OptionDef, OptionParser, Read, RenderingMode}

object ConfigParser {

  implicit val fileRead: Read[File] = reads { s =>
    val file = File(s)
    if (file.notExists) {
      throw new IllegalArgumentException(s"'$file' does not exist")
    }
    file
  }

  implicit val urlRead: Read[URL] = reads {
    new URL(_)
  }

  implicit def enumRead[E <: EnumEntry: Enum]: Read[E] = reads {
    implicitly[Enum[E]].withNameInsensitive
  }
}

class ConfigParser
    extends OptionParser[Config](
      s"java [options] -jar ${WorkflowTestBuildInfo.name}-${WorkflowTestBuildInfo.version}.jar"
    ) {

  import ConfigParser._

  override def showUsageOnError: Boolean = true

  override def renderingMode: RenderingMode = RenderingMode.OneColumn

  def germlineCloudPipelineCommandLineConfig(
      workflowTestType: WorkflowTestType,
      checkConfigBlock: Config => Either[String, Unit],
      defaultGermlineCloudWorkflowConfig: Option[GermlineCloudWorkflowConfig] =
        None
  ): OptionDef[Unit, Config] = {
    cmd(workflowTestType.entryName)
      .text(s"Test the ${workflowTestType.entryName} workflow")
      .action { (_, config) =>
        defaultGermlineCloudWorkflowConfig
          .fold(config)(
            germlineConfig => config.copy(germlineCloudConfig = germlineConfig)
          )
          .copy(test = workflowTestType)
      }
      .children(
        opt[WorkflowTestCategory]('t', "test")
          .text("The type of test to run")
          .required()
          .action { (test, config) =>
            config.copy(
              germlineCloudConfig =
                config.germlineCloudConfig.copy(category = test)
            )
          },
        opt[DataType]('d', "data-type")
          .text(
            s"The data type to test ${DataType.values.mkString("[", ",", "]")}")
          .required()
          .action { (dataType, config) =>
            config.copy(
              germlineCloudConfig = config.germlineCloudConfig.copy(
                dataType = dataType
              )
            )
          },
        opt[String]('b', "branch")
          .text("The branch of truth data to test against (Defaults to master)")
          .optional()
          .action { (branch, config) =>
            config.copy(
              germlineCloudConfig =
                config.germlineCloudConfig.copy(truthBranch = branch)
            )
          },
        opt[String]("use-timestamp")
          .text(
            "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
          )
          .optional()
          .action { (timestamp, config) =>
            config.copy(
              germlineCloudConfig = config.germlineCloudConfig
                .copy(useTimestamp = Option(timestamp))
            )
          },
        opt[Unit]('u', "uncached")
          .text("Disable call-caching for the main workflow run")
          .optional()
          .action { (_, config) =>
            config.copy(
              germlineCloudConfig =
                config.germlineCloudConfig.copy(useCallCaching = false)
            )
          },
        opt[CromwellEnvironment]('e', "env")
          .text(
            s"The environment that this should run in ${CromwellEnvironment.optionsString}"
          )
          .required()
          .action { (env, config) =>
            config.copy(
              germlineCloudConfig = config.germlineCloudConfig.copy(env = env)
            )
          },
        opt[Unit]("update-truth")
          .text(
            "Update the truth data with the results of this run."
          )
          .optional()
          .action { (_, config) =>
            config.copy(
              germlineCloudConfig =
                config.germlineCloudConfig.copy(updateTruth = true)
            )
          },
        checkConfig(checkConfigBlock)
      )
  }

  // NOTE: All the `note("")`s are injecting newlines. Without them, scopt bunches everything
  // up and it's very hard to read.
  help("help").text("Show this output")

  note("")
  cmd(Dummy.entryName)
    .text("Run a dummy test to check the plumbing of the test harness")
    .action((_, config) => config.copy(test = Dummy))


  note("")
  cmd(CloudWorkflow.entryName)
    .text("Test a cloud workflow")
    .action(
      (_, config) =>
        config.copy(
          test = CloudWorkflow
      )
    )
    .children(
      opt[PipelineTestType]('p', "pipeline")
        .text("The pipeline to test")
        .required()
        .action { (pipeline, config) =>
          config.copy(
            cloudWorkflowConfig =
              config.cloudWorkflowConfig.copy(pipeline = pipeline)
          )
        },
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            cloudWorkflowConfig =
              config.cloudWorkflowConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to develop)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            cloudWorkflowConfig =
              config.cloudWorkflowConfig.copy(truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            cloudWorkflowConfig = config.cloudWorkflowConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            cloudWorkflowConfig =
              config.cloudWorkflowConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            cloudWorkflowConfig = config.cloudWorkflowConfig
              .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            cloudWorkflowConfig =
              config.cloudWorkflowConfig.copy(useCallCaching = false)
          )
        }
    )
}
