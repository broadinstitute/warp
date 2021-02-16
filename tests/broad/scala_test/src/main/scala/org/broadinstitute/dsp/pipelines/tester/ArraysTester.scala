package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  MultiSampleArraysInputs,
  MultiSampleArraysValidationInputs,
  SingleSampleArraysInputs,
  SingleSampleArraysValidationInputs
}
import org.broadinstitute.dsp.pipelines.util.ArrayType

class ArraysTester(testerConfig: ArraysConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override val workflowName: String = testerConfig.arrayType.workflowBaseName

  val workflowDir
    : File = CromwellWorkflowTester.PipelineRoot / "broad" / "arrays" / testerConfig.arrayType.pathName

  protected val vaultTokenPath: String =
    s"gs://broad-dsp-gotc-arrays-$envString-tokens/arrayswdl.token"

  override protected val validationWorkflowName: String = s"Verify$workflowName"

  override protected lazy val googleProject: String = {
    if (env.picardEnv.equals("dev")) {
      s"broad-gotc-${env.picardEnv}"
    } else {
      s"broad-arrays-${env.picardEnv}"
    }
  }

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/arrays/${testerConfig.arrayType.pathName}/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/arrays/${testerConfig.arrayType.pathName}/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override def getInputContents(fileName: String): String = {
    super
      .getInputContents(fileName)
      .lines
      .map(
        _.replace("{ENV}", envString)
          .replace("{VAULT_TOKEN_PATH}", vaultTokenPath)
      )
      .mkString
  }

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest): String = {
    testerConfig.arrayType match {
      case ArrayType.Single =>
        buildSingleSampleValidationWdlInputs(workflowTest)
      case ArrayType.Multi =>
        buildMultiSampleArraysValidationWdlInputs(workflowTest)
    }
  }

  private def buildSingleSampleValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val singleSampleArraysInputs = new SingleSampleArraysInputs(
      workflowTest.runParameters.workflowInputs
    )
    val outputBaseName = singleSampleArraysInputs.getChipwellBarcode("Arrays")
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath

    val metricsFileNames = ioUtil
      .listGoogleObjects(truthCloudPath)
      .filter(_.getPath.endsWith("metrics"))
      .map(uriToFilename)
    val validationInputs = SingleSampleArraysValidationInputs(
      test_metrics = metricsFileNames.map(resultsCloudPath.resolve),
      truth_metrics = metricsFileNames.map(truthCloudPath.resolve),
      bead_pool_manifest_file =
        new URI(singleSampleArraysInputs.getBeadPoolManifestFile("Arrays")),
      test_gtc = resultsCloudPath.resolve(s"$outputBaseName.gtc"),
      truth_gtc = truthCloudPath.resolve(s"$outputBaseName.gtc"),
      test_vcf = resultsCloudPath.resolve(s"$outputBaseName.vcf.gz"),
      truth_vcf = truthCloudPath.resolve(s"$outputBaseName.vcf.gz"),
      test_fp_vcf =
        resultsCloudPath.resolve(s"$outputBaseName.fingerprint.vcf.gz"),
      truth_fp_vcf =
        truthCloudPath.resolve(s"$outputBaseName.fingerprint.vcf.gz"),
      truth_green_idat_md5 =
        truthCloudPath.resolve(s"${outputBaseName}_Grn.idat.md5sum"),
      test_green_idat_md5 =
        resultsCloudPath.resolve(s"${outputBaseName}_Grn.idat.md5sum"),
      truth_red_idat_md5 =
        truthCloudPath.resolve(s"${outputBaseName}_Red.idat.md5sum"),
      test_red_idat_md5 =
        resultsCloudPath.resolve(s"${outputBaseName}_Red.idat.md5sum")
    )
    SingleSampleArraysValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }

  private def buildMultiSampleArraysValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val callsetName = new MultiSampleArraysInputs(
      workflowTest.runParameters.workflowInputs
    ).getCallsetName("MultiSampleArrays")

    val validationInputs = MultiSampleArraysValidationInputs(
      test_vcf = workflowTest.runParameters.resultsCloudPath
        .resolve(s"$callsetName.vcf.gz"),
      truth_vcf = workflowTest.runParameters.truthCloudPath
        .resolve(s"$callsetName.vcf.gz")
    )

    MultiSampleArraysValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
