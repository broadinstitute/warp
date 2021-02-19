package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.commandline.CromwellEnvironment
import org.broadinstitute.dsp.pipelines.config.GermlineCloudWorkflowConfig
import org.broadinstitute.dsp.pipelines.inputs.{
  AllOfUsInputs,
  AllOfUsValidationInputs
}

class AllOfUsTester(testerConfig: GermlineCloudWorkflowConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends GermlineSingleSampleTester(testerConfig) {
  override protected def env: CromwellEnvironment = testerConfig.env

  override val workflowName: String = "AllOfUsWholeGenome"

  override protected val validationWorkflowName: String = "VerifyAllOfUs"

  override lazy val workflowDir: File =
    CromwellWorkflowTester.WarpRoot / "projects" / "all_of_us"

  override lazy val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/all_of_us/$testTypeString/$timestamp/"
    )

  override lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/all_of_us/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val allOfUsInputs = new AllOfUsInputs(
      workflowTest.runParameters.workflowInputs
    )
    val outputBaseName = allOfUsInputs.getBaseFileName(workflowName)
    val gvcfBaseName = allOfUsInputs.getFinalGvcfBaseName(workflowName)

    val resultsCloudPath = workflowTest.runParameters.resultsCloudPath
    val sampleTruthCloudPath = workflowTest.runParameters.truthCloudPath
    val allOfUsTruthCloudPath = new URI(
      s"gs://broad-gotc-test-storage/all_of_us/$testTypeString/truth/"
    )

    val metricsFileNames = ioUtil
      .listGoogleObjects(sampleTruthCloudPath)
      .filter(_.getPath.endsWith("metrics"))
      .map(uriToFilename)

    val validationInputs = AllOfUsValidationInputs(
      testMetrics = metricsFileNames.map(resultsCloudPath.resolve),
      truthMetrics = metricsFileNames.map(sampleTruthCloudPath.resolve),
      testCram = resultsCloudPath.resolve(s"$outputBaseName.cram"),
      testCrai = resultsCloudPath.resolve(s"$outputBaseName.cram.crai"),
      truthCram = sampleTruthCloudPath.resolve(s"$outputBaseName.cram"),
      truthCrai = sampleTruthCloudPath.resolve(s"$outputBaseName.cram.crai"),
      testGvcf = resultsCloudPath.resolve(s"$gvcfBaseName.g.vcf.gz"),
      truthGvcf = sampleTruthCloudPath.resolve(s"$gvcfBaseName.g.vcf.gz"),
      testFiltrationReport = resultsCloudPath.resolve("filtration_report.tsv"),
      truthFiltrationReport =
        allOfUsTruthCloudPath.resolve("filtration_report.tsv"),
      testSignificantVariantsVcf = resultsCloudPath.resolve(
        s"$gvcfBaseName.filtered.clinsig-variants.vcf.gz"),
      truthVariants = allOfUsTruthCloudPath.resolve(s"mutations_to_find.csv")
    )
    AllOfUsValidationInputs.marshall(validationInputs).printWith(implicitly)
  }
}
