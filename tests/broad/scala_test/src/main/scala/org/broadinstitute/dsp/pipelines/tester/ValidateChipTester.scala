package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config.ValidateChipConfig
import org.broadinstitute.dsp.pipelines.inputs.{
  ValidateChipInputs,
  ValidateChipValidationInputs
}

class ValidateChipTester(testerConfig: ValidateChipConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override val workflowName: String = "ValidateChip"

  val workflowDir
    : File = CromwellWorkflowTester.PipelineRoot / "broad" / "arrays" / "validate_chip"

  protected val resultsPrefix: URI = URI.create(
    s"gs://broad-gotc-test-results/$envString/validate_chip/$testTypeString/$timestamp/"
  )
  protected val truthPrefix: URI = URI.create(
    s"gs://broad-gotc-test-storage/validate_chip/$testTypeString/truth/${testerConfig.truthBranch}/"
  )

  override protected val validationWorkflowName: String = "VerifyValidateChip"

  override protected lazy val googleProject: String = {
    if (env.picardEnv.equals("dev")) {
      s"broad-gotc-${env.picardEnv}"
    } else {
      s"broad-arrays-${env.picardEnv}"
    }
  }

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val validateChipInputs = new ValidateChipInputs(
      workflowTest.runParameters.workflowInputs
    )
    val outputBaseName = validateChipInputs.chipWellBarcode(workflowName)
    val resultsCloudPath = workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath
    val metricsFileNames = ioUtil
      .listGoogleObjects(resultsCloudPath)
      .filter(_.getPath.endsWith("metrics"))
      .map(uriToFilename)
    val validationInputs = ValidateChipValidationInputs(
      testGtc = resultsCloudPath.resolve(s"$outputBaseName.gtc"),
      truthGtc = truthCloudPath.resolve(s"$outputBaseName.gtc"),
      beadPoolManifestFile =
        new URI(validateChipInputs.getBeadPoolManifestFile("ValidateChip")),
      testVcf = resultsCloudPath.resolve(s"$outputBaseName.vcf.gz"),
      truthVcf = truthCloudPath.resolve(s"$outputBaseName.vcf.gz"),
      testGenotypeConcordanceVcf = resultsCloudPath.resolve(
        s"$outputBaseName.gc.genotype_concordance.vcf.gz"),
      truthGenotypeConcordanceVcf = truthCloudPath.resolve(
        s"$outputBaseName.gc.genotype_concordance.vcf.gz"),
      testIndelGenotypeConcordanceVcf = resultsCloudPath.resolve(
        s"$outputBaseName.indels.gc.genotype_concordance.vcf.gz"
      ),
      truthIndelGenotypeConcordanceVcf = truthCloudPath.resolve(
        s"$outputBaseName.indels.gc.genotype_concordance.vcf.gz"),
      testMetrics = metricsFileNames.map(resultsCloudPath.resolve),
      truthMetrics = metricsFileNames.map(truthCloudPath.resolve)
    )
    ValidateChipValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
