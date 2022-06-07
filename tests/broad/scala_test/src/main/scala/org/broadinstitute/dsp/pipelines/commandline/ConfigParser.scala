package org.broadinstitute.dsp.pipelines.commandline

import java.net.URL

import better.files.File
import enumeratum.{Enum, EnumEntry}
import org.broadinstitute.dsp.pipelines.util.DataType
import org.broadinstitute.dsp.pipelines.commandline.PapiVersion.PAPIv2
import org.broadinstitute.dsp.pipelines.commandline.WorkflowTestType._
import org.broadinstitute.dsp.pipelines.config.GermlineCloudWorkflowConfig
import org.broadinstitute.dsp.pipelines.util.{ArrayType, WorkflowTestBuildInfo}
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
        opt[PapiVersion]("papi-version")
          .text("The version of Pipelines API to use")
          .optional()
          .action { (papiVersion, config) =>
            config.copy(
              germlineCloudConfig =
                config.germlineCloudConfig.copy(papiVersion = papiVersion)
            )
          },
        checkConfig(checkConfigBlock)
      )
  }

  // NOTE: All the `note("")`s are injecting newlines. Without them, scopt bunches everything
  // up and it's very hard to read.
  help("help").text("Show this output")

  note("")
  germlineCloudPipelineCommandLineConfig(
    AllOfUs, { config =>
      (config.test, config.germlineCloudConfig.category) match {
        case (AllOfUs, WorkflowTestCategory.Load) =>
          failure(
            "The AllOfUs test is not configured to run load."
          )
        case _ => success
      }
    },
    Some(GermlineCloudWorkflowConfig(papiVersion = PAPIv2))
  )

  note("")
  cmd(AnnotationFiltration.entryName)
    .text("Test the WGS clinical annotation/filtration workflow")
    .action { (_, config) =>
      config.copy(test = AnnotationFiltration)
    }
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            annotationFiltrationConfig =
              config.annotationFiltrationConfig.copy(category = test)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            annotationFiltrationConfig =
              config.annotationFiltrationConfig.copy(env = env)
          )
        }
    )

  note("")
  cmd(Arrays.entryName)
    .text("Test the Arrays (Single- or Multi-sample) workflow")
    .action(
      (_, config) =>
        config.copy(
          test = Arrays
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .optional()
        .action { (test, config) =>
          config.copy(
            arraysConfig = config.arraysConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            arraysConfig = config.arraysConfig.copy(truthBranch = branch)
          )
        },
      opt[ArrayType]('a', "array-type")
        .text(
          s"The Array type of test to run ${ArrayType.values.mkString("[", ", ", "]")}"
        )
        .optional()
        .action { (arrayType, config) =>
          config.copy(
            arraysConfig = config.arraysConfig.copy(arrayType = arrayType)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            arraysConfig = config.arraysConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            arraysConfig = config.arraysConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            arraysConfig = config.arraysConfig
              .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            arraysConfig = config.arraysConfig.copy(useCallCaching = false)
          )
        },
      opt[PapiVersion]("papi-version")
        .text("The version of Pipelines API to use")
        .optional()
        .action { (papiVersion, config) =>
          config.copy(
            arraysConfig = config.arraysConfig.copy(papiVersion = papiVersion)
          )
        }
    )

  note("")
  cmd(IlluminaGenotypingArray.entryName)
    .text("Test the IlluminaGenotypingArray workflow")
    .action(
      (_, config) =>
        config.copy(
          test = IlluminaGenotypingArray
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            illuminaGenotypingArrayConfig =
              config.illuminaGenotypingArrayConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            illuminaGenotypingArrayConfig =
              config.illuminaGenotypingArrayConfig.copy(truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            illuminaGenotypingArrayConfig =
              config.illuminaGenotypingArrayConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            illuminaGenotypingArrayConfig =
              config.illuminaGenotypingArrayConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            illuminaGenotypingArrayConfig = config.illuminaGenotypingArrayConfig
              .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            illuminaGenotypingArrayConfig =
              config.illuminaGenotypingArrayConfig.copy(useCallCaching = false)
          )
        },
      opt[PapiVersion]("papi-version")
        .text("The version of Pipelines API to use")
        .optional()
        .action { (papiVersion, config) =>
          config.copy(
            illuminaGenotypingArrayConfig = config.illuminaGenotypingArrayConfig
              .copy(papiVersion = papiVersion)
          )
        }
    )

  note("")
  cmd(Dummy.entryName)
    .text("Run a dummy test to check the plumbing of the test harness")
    .action((_, config) => config.copy(test = Dummy))

  note("")
  germlineCloudPipelineCommandLineConfig(
    ExternalReprocessing, { config =>
      (config.test, config.germlineCloudConfig.category) match {
        case (ExternalReprocessing, WorkflowTestCategory.Load) =>
          failure(
            "The external reprocessing test is not configured to run load."
          )
        case (ExternalReprocessing, WorkflowTestCategory.Scientific) =>
          failure(
            "The external reprocessing test is not configured to run scientific."
          )
        case _ => success
      }
    }
  )

  note("")
  cmd(GenotypeConcordance.entryName)
    .text("Run genotype concordance on a GVCF from NA12878")
    .action(
      (_, config) => config.copy(test = GenotypeConcordance)
    )
    .children(
      opt[String]('n', "nist-version")
        .text("The NIST release version to run against")
        .optional()
        .action(
          (version, config) =>
            config.copy(
              genotypeConcordanceConfig =
                config.genotypeConcordanceConfig.copy(nistVersion = version)
          )
        ),
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            genotypeConcordanceConfig =
              config.genotypeConcordanceConfig.copy(env = env)
          )
        },
      opt[Unit]('l', "use-latest")
        .text(
          "Use the latest GVCF generated in the environment, instead of generating a fresh one"
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            genotypeConcordanceConfig =
              config.genotypeConcordanceConfig.copy(useLatestNA12878 = true)
          )
        }
    )

  note("")
  germlineCloudPipelineCommandLineConfig(
    GermlineSingleSample, { config =>
      val c = config.germlineCloudConfig
      (config.test, c.category, c.useTimestamp, c.updateTruth) match {
        case (GermlineSingleSample, WorkflowTestCategory.Load, Some(_), _) =>
          failure(
            "Cannot use previous run for Load test, since it does not run validation."
          )
        case (GermlineSingleSample, WorkflowTestCategory.Load, _, true) =>
          failure(
            "Cannot update truth data for Load test, since it does not run validation."
          )
        case _ => success
      }
    }
  )

  note("")
  germlineCloudPipelineCommandLineConfig(
    JointGenotyping, { config =>
      // Tuple is needed because Plumbing is the default category
      (config.test, config.germlineCloudConfig.category) match {
        case (JointGenotyping, WorkflowTestCategory.Load) =>
          failure("The joint genotyping test is not configured to run load")
        case _ =>
          success
      }
    }
  )

  note("")
  germlineCloudPipelineCommandLineConfig(
    Reprocessing, { config =>
      (config.test, config.germlineCloudConfig.category) match {
        case (Reprocessing, WorkflowTestCategory.Load) =>
          failure(
            "The reprocessing test is not configured to run load."
          )
        case (Reprocessing, WorkflowTestCategory.Scientific) =>
          failure(
            "The reprocessing test is not configured to run scientific."
          )
        case _ => success
      }
    }
  )
  note("")
  germlineCloudPipelineCommandLineConfig(
    ReblockGvcf, { config =>
      (config.test, config.germlineCloudConfig.category) match {
        case (ReblockGvcf, WorkflowTestCategory.Load) =>
          failure(
            "The ReblockGvcf test is not configured to run load."
          )
        case _ => success
      }
    }
  )
  note("")
  cmd(ValidateChip.entryName)
    .text(s"Test the ${ValidateChip.entryName} workflow")
    .action { (_, config) =>
      config.copy(test = ValidateChip)
    }
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            validateChipConfig = config.validateChipConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            validateChipConfig =
              config.validateChipConfig.copy(truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            validateChipConfig = config.validateChipConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            validateChipConfig =
              config.validateChipConfig.copy(updateTruth = true)
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            validateChipConfig =
              config.validateChipConfig.copy(useCallCaching = false)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            validateChipConfig = config.validateChipConfig
              .copy(useTimestamp = Option(timestamp))
          )
        }
    )
  note("")
  cmd(SomaticSingleSample.entryName)
    .text("Test the SomaticSingleSample workflow")
    .action(
      (_, config) =>
        config.copy(
          test = SomaticSingleSample
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            somaticCloudWorkflowConfig =
              config.somaticCloudWorkflowConfig.copy(category = test)
          )
        },
      opt[DataType]('d', "data-type")
        .text(
          s"The data type to test ${DataType.values.mkString("[", ",", "]")}")
        .required()
        .action { (dataType, config) =>
          config.copy(
            somaticCloudWorkflowConfig = config.somaticCloudWorkflowConfig.copy(
              dataType = dataType
            )
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            somaticCloudWorkflowConfig =
              config.somaticCloudWorkflowConfig.copy(truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            somaticCloudWorkflowConfig =
              config.somaticCloudWorkflowConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            somaticCloudWorkflowConfig =
              config.somaticCloudWorkflowConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            somaticCloudWorkflowConfig = config.somaticCloudWorkflowConfig
              .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            somaticCloudWorkflowConfig =
              config.somaticCloudWorkflowConfig.copy(useCallCaching = false)
          )
        },
      opt[PapiVersion]("papi-version")
        .text("The version of Pipelines API to use")
        .optional()
        .action { (papiVersion, config) =>
          config.copy(
            somaticCloudWorkflowConfig =
              config.somaticCloudWorkflowConfig.copy(papiVersion = papiVersion)
          )
        }
    )

  note("")
  cmd(CramToUnmappedBams.entryName)
    .text("Test the CramToUnmappedBams workflow")
    .action(
      (_, config) =>
        config.copy(
          test = CramToUnmappedBams
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            cramToUnmappedBamsConfig =
              config.cramToUnmappedBamsConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            cramToUnmappedBamsConfig =
              config.cramToUnmappedBamsConfig.copy(truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            cramToUnmappedBamsConfig =
              config.cramToUnmappedBamsConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            cramToUnmappedBamsConfig =
              config.cramToUnmappedBamsConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            cramToUnmappedBamsConfig = config.cramToUnmappedBamsConfig
              .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            cramToUnmappedBamsConfig =
              config.cramToUnmappedBamsConfig.copy(useCallCaching = false)
          )
        },
      opt[PapiVersion]("papi-version")
        .text("The version of Pipelines API to use")
        .optional()
        .action { (papiVersion, config) =>
          config.copy(
            cramToUnmappedBamsConfig =
              config.cramToUnmappedBamsConfig.copy(papiVersion = papiVersion)
          )
        }
    )

  note("")
  cmd(GDCWholeGenomeSomaticSingleSample.entryName)
    .text("Test the GDCWholeGenomeSomaticSingleSample workflow")
    .action(
      (_, config) =>
        config.copy(
          test = GDCWholeGenomeSomaticSingleSample
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            gdcWholeGenomeSomaticSingleSampleConfig =
              config.gdcWholeGenomeSomaticSingleSampleConfig.copy(
                category = test)
          )
        },
      opt[DataType]('d', "data-type")
        .text(
          s"The data type to test ${DataType.values.mkString("[", ",", "]")}")
        .required()
        .action { (dataType, config) =>
          config.copy(
            gdcWholeGenomeSomaticSingleSampleConfig =
              config.gdcWholeGenomeSomaticSingleSampleConfig.copy(
                dataType = dataType
              )
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            gdcWholeGenomeSomaticSingleSampleConfig =
              config.gdcWholeGenomeSomaticSingleSampleConfig.copy(
                truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            gdcWholeGenomeSomaticSingleSampleConfig =
              config.gdcWholeGenomeSomaticSingleSampleConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            gdcWholeGenomeSomaticSingleSampleConfig =
              config.gdcWholeGenomeSomaticSingleSampleConfig.copy(
                updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            gdcWholeGenomeSomaticSingleSampleConfig =
              config.gdcWholeGenomeSomaticSingleSampleConfig
                .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            gdcWholeGenomeSomaticSingleSampleConfig =
              config.gdcWholeGenomeSomaticSingleSampleConfig.copy(
                useCallCaching = false)
          )
        },
      opt[PapiVersion]("papi-version")
        .text("The version of Pipelines API to use")
        .optional()
        .action { (papiVersion, config) =>
          config.copy(
            gdcWholeGenomeSomaticSingleSampleConfig =
              config.gdcWholeGenomeSomaticSingleSampleConfig.copy(
                papiVersion = papiVersion)
          )
        }
    )
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

  note("")
  germlineCloudPipelineCommandLineConfig(
    VariantCalling, { config =>
      (config.test, config.germlineCloudConfig.category) match {
        case (VariantCalling, WorkflowTestCategory.Scientific) =>
          failure(
            "The VariantCalling test is not configured to run scientific. This case is covered by the GermlineSingleSample tests"
          )
        case (VariantCalling, WorkflowTestCategory.Load) =>
          failure(
            "The VariantCalling test is not configured to run load."
          )
        case _ => success
      }
    },
    Some(GermlineCloudWorkflowConfig(papiVersion = PAPIv2))
  )

  note("")
  cmd(Imputation.entryName)
    .text("Test the Imputation workflow")
    .action(
      (_, config) =>
        config.copy(
          test = Imputation
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .optional()
        .action { (test, config) =>
          config.copy(
            imputationConfig = config.imputationConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            imputationConfig =
              config.imputationConfig.copy(truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            imputationConfig = config.imputationConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            imputationConfig = config.imputationConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            imputationConfig = config.imputationConfig
              .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            imputationConfig =
              config.imputationConfig.copy(useCallCaching = false)
          )
        },
      opt[PapiVersion]("papi-version")
        .text("The version of Pipelines API to use")
        .optional()
        .action { (papiVersion, config) =>
          config.copy(
            imputationConfig =
              config.imputationConfig.copy(papiVersion = papiVersion)
          )
        }
    )

  note("")
  cmd(RNAWithUMIs.entryName)
    .text("Test the RNA with UMIs workflow")
    .action(
      (_, config) =>
        config.copy(
          test = RNAWithUMIs
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .optional()
        .action { (test, config) =>
          config.copy(
            rnaWithUMIsConfig = config.rnaWithUMIsConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            rnaWithUMIsConfig =
              config.rnaWithUMIsConfig.copy(truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            rnaWithUMIsConfig = config.rnaWithUMIsConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            rnaWithUMIsConfig =
              config.rnaWithUMIsConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            rnaWithUMIsConfig = config.rnaWithUMIsConfig
              .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            rnaWithUMIsConfig =
              config.rnaWithUMIsConfig.copy(useCallCaching = false)
          )
        }
    )

  note("")
  cmd(CheckFingerprint.entryName)
    .text("Test the CheckFingerprint workflow")
    .action(
      (_, config) =>
        config.copy(
          test = CheckFingerprint
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .optional()
        .action { (test, config) =>
          config.copy(
            checkFingerprintConfig =
              config.checkFingerprintConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            checkFingerprintConfig =
              config.checkFingerprintConfig.copy(truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            checkFingerprintConfig =
              config.checkFingerprintConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            checkFingerprintConfig =
              config.checkFingerprintConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            checkFingerprintConfig = config.checkFingerprintConfig
              .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            checkFingerprintConfig =
              config.checkFingerprintConfig.copy(useCallCaching = false)
          )
        }
    )

  note("")
  cmd(BroadInternalRNAWithUMIs.entryName)
    .text("Test the Broad Internal RNA with UMIs workflow")
    .action(
      (_, config) =>
        config.copy(
          test = BroadInternalRNAWithUMIs
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .optional()
        .action { (test, config) =>
          config.copy(
            broadInternalRNAWithUMIsConfig =
              config.broadInternalRNAWithUMIsConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            broadInternalRNAWithUMIsConfig =
              config.broadInternalRNAWithUMIsConfig.copy(truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            broadInternalRNAWithUMIsConfig =
              config.broadInternalRNAWithUMIsConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            broadInternalRNAWithUMIsConfig =
              config.broadInternalRNAWithUMIsConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            broadInternalRNAWithUMIsConfig =
              config.broadInternalRNAWithUMIsConfig
                .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            broadInternalRNAWithUMIsConfig =
              config.broadInternalRNAWithUMIsConfig.copy(useCallCaching = false)
          )
        }
    )

  note("")
  cmd(UltimaGenomicsWholeGenomeGermline.entryName)
    .text(s"Test the ${UltimaGenomicsWholeGenomeGermline.entryName} workflow")
    .action(
      (_, config) =>
        config.copy(
          test = UltimaGenomicsWholeGenomeGermline
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            ultimaGenomicsWholeGenomeGermlineConfig =
              config.ultimaGenomicsWholeGenomeGermlineConfig.copy(
                category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            ultimaGenomicsWholeGenomeGermlineConfig =
              config.ultimaGenomicsWholeGenomeGermlineConfig.copy(
                truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            ultimaGenomicsWholeGenomeGermlineConfig =
              config.ultimaGenomicsWholeGenomeGermlineConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            ultimaGenomicsWholeGenomeGermlineConfig =
              config.ultimaGenomicsWholeGenomeGermlineConfig.copy(
                updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            ultimaGenomicsWholeGenomeGermlineConfig =
              config.ultimaGenomicsWholeGenomeGermlineConfig
                .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            ultimaGenomicsWholeGenomeGermlineConfig =
              config.ultimaGenomicsWholeGenomeGermlineConfig.copy(
                useCallCaching = false)
          )
        }
    )

  note("")
  cmd(BroadInternalUltimaGenomics.entryName)
    .text(s"Test the ${BroadInternalUltimaGenomics.entryName} workflow")
    .action(
      (_, config) =>
        config.copy(
          test = BroadInternalUltimaGenomics
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            broadInternalUltimaGenomicsConfig =
              config.broadInternalUltimaGenomicsConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            broadInternalUltimaGenomicsConfig =
              config.broadInternalUltimaGenomicsConfig.copy(
                truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            broadInternalUltimaGenomicsConfig =
              config.broadInternalUltimaGenomicsConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            broadInternalUltimaGenomicsConfig =
              config.broadInternalUltimaGenomicsConfig.copy(updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            ultimaGenomicsWholeGenomeGermlineConfig =
              config.ultimaGenomicsWholeGenomeGermlineConfig
                .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            broadInternalUltimaGenomicsConfig =
              config.broadInternalUltimaGenomicsConfig.copy(
                useCallCaching = false)
          )
        }
    )

  note("")
  cmd(UltimaGenomicsJointGenotyping.entryName)
    .text(s"Test the ${UltimaGenomicsJointGenotyping.entryName} workflow")
    .action(
      (_, config) =>
        config.copy(
          test = UltimaGenomicsJointGenotyping
      )
    )
    .children(
      opt[WorkflowTestCategory]('t', "test")
        .text("The type of test to run")
        .required()
        .action { (test, config) =>
          config.copy(
            ultimaGenomicsJointGenotypingConfig =
              config.ultimaGenomicsJointGenotypingConfig.copy(category = test)
          )
        },
      opt[String]('b', "branch")
        .text("The branch of truth data to test against (Defaults to master)")
        .optional()
        .action { (branch, config) =>
          config.copy(
            ultimaGenomicsJointGenotypingConfig =
              config.ultimaGenomicsJointGenotypingConfig.copy(
                truthBranch = branch)
          )
        },
      opt[CromwellEnvironment]('e', "env")
        .text(
          s"The environment that this should run in ${CromwellEnvironment.optionsString}"
        )
        .required()
        .action { (env, config) =>
          config.copy(
            ultimaGenomicsJointGenotypingConfig =
              config.ultimaGenomicsJointGenotypingConfig.copy(env = env)
          )
        },
      opt[Unit]("update-truth")
        .text(
          "Update the truth data with the results of this run."
        )
        .optional()
        .action { (_, config) =>
          config.copy(
            ultimaGenomicsJointGenotypingConfig =
              config.ultimaGenomicsJointGenotypingConfig.copy(
                updateTruth = true)
          )
        },
      opt[String]("use-timestamp")
        .text(
          "Do not run the workflows. Instead, just use a previous runs timestamp (yyyy-MM-dd-HH-mm-ss)"
        )
        .optional()
        .action { (timestamp, config) =>
          config.copy(
            ultimaGenomicsJointGenotypingConfig =
              config.ultimaGenomicsJointGenotypingConfig
                .copy(useTimestamp = Option(timestamp))
          )
        },
      opt[Unit]('u', "uncached")
        .text("Disable call-caching for the main workflow run")
        .optional()
        .action { (_, config) =>
          config.copy(
            ultimaGenomicsJointGenotypingConfig =
              config.ultimaGenomicsJointGenotypingConfig.copy(
                useCallCaching = false)
          )
        }
    )
}
