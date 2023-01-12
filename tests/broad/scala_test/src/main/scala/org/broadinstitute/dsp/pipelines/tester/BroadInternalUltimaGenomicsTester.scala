package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  UltimaGenomicsWholeGenomeGermlineInputs,
  UltimaGenomicsWholeGenomeGermlineValidationInputs
}

class BroadInternalUltimaGenomicsTester(
    testerConfig: BroadInternalUltimaGenomicsConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override val workflowName: String = "BroadInternalUltimaGenomics"

  val workflowDir
    : File = CromwellWorkflowTester.PipelineRoot / "broad" / "internal" / "dna_seq" / "germline" / "single_sample" / "UltimaGenomics"

  protected val vaultTokenPath: String =
    s"gs://broad-dsp-gotc-arrays-$envString-tokens/arrayswdl.token"

  override protected val validationWorkflowName: String =
    "VerifyUltimaGenomicsWholeGenomeGermline"

  override protected lazy val googleProject: String = {
    if (env.picardEnv.equals("dev")) {
      s"broad-gotc-${env.picardEnv}"
    } else {
      s"broad-arrays-${env.picardEnv}"
    }
  }

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/broad_internal/dna_seq/germline/single_sample/UltimaGenomics/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/broad_internal/dna_seq/germline/single_sample/UltimaGenomics/$testTypeString/truth/${testerConfig.truthBranch}/"
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
      workflowTest: WorkflowTest
  ): String = {
    val ultimaGenomicsWholeGenomeGermlineInputs =
      new UltimaGenomicsWholeGenomeGermlineInputs(
        workflowTest.runParameters.workflowInputs
      )
    val outputBaseName =
      ultimaGenomicsWholeGenomeGermlineInputs.getBaseFileName(workflowName)
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath
    val metricsFileNames = ioUtil
      .listGoogleObjects(truthCloudPath)
      .filter(_.getPath.endsWith("metrics"))
      .map(uriToFilename)
    val validationInputs = UltimaGenomicsWholeGenomeGermlineValidationInputs(
      testMetrics = metricsFileNames.map(resultsCloudPath.resolve),
      truthMetrics = metricsFileNames.map(truthCloudPath.resolve),
      testCram = resultsCloudPath.resolve(s"$outputBaseName.cram"),
      truthCram = truthCloudPath.resolve(s"$outputBaseName.cram"),
      testCrai = resultsCloudPath.resolve(s"$outputBaseName.cram.crai"),
      truthCrai = truthCloudPath.resolve(s"$outputBaseName.cram.crai"),
      testVcf = resultsCloudPath.resolve(s"$outputBaseName.vcf.gz"),
      truthVcf = truthCloudPath.resolve(s"$outputBaseName.vcf.gz"),
      testFilteredVcf =
        resultsCloudPath.resolve(s"$outputBaseName.filtered.vcf.gz"),
      truthFilteredVcf =
        truthCloudPath.resolve(s"$outputBaseName.filtered.vcf.gz"),
      testGvcf =
        resultsCloudPath.resolve(s"$outputBaseName.annotated.rb.g.vcf.gz"),
      truthGvcf =
        truthCloudPath.resolve(s"$outputBaseName.annotated.rb.g.vcf.gz")
    )
    UltimaGenomicsWholeGenomeGermlineValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
