package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File

import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRunParameters,
  WorkflowTest
}
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  UltimaGenomicsWholeGenomeGermlineInputs,
  UltimaGenomicsWholeGenomeGermlineValidationInputs
}

class UltimaGenomicsWholeGenomeGermlineTester(
    testerConfig: UltimaGenomicsWholeGenomeGermlineConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override val workflowName: String = "UltimaGenomicsWholeGenomeGermline"

  override lazy val workflowDir: File =
    CromwellWorkflowTester.PipelineRoot / "broad" / "dna_seq" / "germline" / "single_sample" / "ugwgs"

  override protected val validationWorkflowName: String =
    "VerifyUltimaGenomicsWholeGenomeGermline"

  protected lazy val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/ugwgs/$testTypeString/$timestamp/"
    )
  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/ugwgs/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override def generateRunParameters: Seq[WorkflowRunParameters] = {

    logger.info(s"workflowInputRoot: $workflowInputRoot")
    inputFileNames.map { fileName =>
      val projectSampleAlias = fileName.replace(".json", "")
      val resultsPath =
        resultsPrefix.resolve(s"$projectSampleAlias/")
      val truthPath = truthPrefix.resolve(s"$projectSampleAlias/")
      logger.info("Running: " + fileName)
      logger.info("Results will be written to: " + resultsPath)

      WorkflowRunParameters(
        id = s"${envString}_$projectSampleAlias",
        workflowInputs = getInputContents(fileName),
        resultsCloudPath = resultsPath,
        truthCloudPath = truthPath
      )
    }
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
