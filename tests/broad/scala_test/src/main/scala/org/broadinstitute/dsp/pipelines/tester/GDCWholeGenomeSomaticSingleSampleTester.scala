package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  GDCWholeGenomeSomaticSingleSampleInputs,
  GDCWholeGenomeSomaticSingleSampleValidationInputs
}

class GDCWholeGenomeSomaticSingleSampleTester(
    testerConfig: GDCWholeGenomeSomaticSingleSampleConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  val dataTypePrefix: String = dataTypePrefix(testerConfig.dataType)
  val dataTypeString: String = testerConfig.dataType.entryName.toLowerCase

  override val workflowName: String = s"GDC${dataTypePrefix}SomaticSingleSample"

  val workflowDir
    : File = CromwellWorkflowTester.WarpRoot / "pipelines" / "broad" / "dna_seq" / "somatic" / "single_sample" / dataTypeString / "gdc_genome"

  override protected val validationWorkflowName: String =
    "VerifyGDCSomaticSingleSample"

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/somatic_single_sample/gdc/$dataTypeString/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/somatic_single_sample/gdc/$dataTypeString/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val gdcWholeGenomeSomaticSingleSampleInputs =
      new GDCWholeGenomeSomaticSingleSampleInputs(
        workflowTest.runParameters.workflowInputs
      )
    val outputBaseName =
      gdcWholeGenomeSomaticSingleSampleInputs.getBaseFileName(workflowName)
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath
    val metricsFileNames = ioUtil
      .listGoogleObjects(truthCloudPath)
      .filter(_.getPath.endsWith("metrics"))
      .map(uriToFilename)
    val validationInputs = GDCWholeGenomeSomaticSingleSampleValidationInputs(
      testMetrics = metricsFileNames.map(resultsCloudPath.resolve),
      truthMetrics = metricsFileNames.map(truthCloudPath.resolve),
      testBam = resultsCloudPath.resolve(s"$outputBaseName.aln.mrkdp.bam"),
      testBai = resultsCloudPath.resolve(s"$outputBaseName.aln.mrkdp.bai"),
      truthBam = truthCloudPath.resolve(s"$outputBaseName.aln.mrkdp.bam"),
      truthBai = truthCloudPath.resolve(s"$outputBaseName.aln.mrkdp.bai")
    )
    GDCWholeGenomeSomaticSingleSampleValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
