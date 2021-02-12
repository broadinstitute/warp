package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  SingleSampleArraysInputs,
  SingleSampleArraysValidationInputs
}

class IlluminaGenotypingArrayTester(testerConfig: IlluminaGenotypingArrayConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override val workflowName: String = "IlluminaGenotypingArray"

  val workflowDir
    : File = CromwellWorkflowTester.PipelineRoot / "broad" / "genotyping" / "illumina"

  override protected val validationWorkflowName: String = "VerifyArrays"

  override protected lazy val googleProject: String = {
    if (env.picardEnv.equals("dev")) {
      s"broad-gotc-${env.picardEnv}"
    } else {
      s"broad-arrays-${env.picardEnv}"
    }
  }

  // Validation uses the same options as the arrays workflow
  override protected lazy val validationWdlOptions: String = {
    readTestOptions(releaseDir, env)
  }

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/genotyping/illumina/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/genotyping/illumina/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val singleSampleArraysInputs = new SingleSampleArraysInputs(
      workflowTest.runParameters.workflowInputs
    )
    val outputBaseName =
      singleSampleArraysInputs.getChipwellBarcode(workflowName)
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
        new URI(singleSampleArraysInputs.getBeadPoolManifestFile(workflowName)),
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
}
