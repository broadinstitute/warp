package org.broadinstitute.dsp.pipelines.tester
import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRunParameters,
  WorkflowTest
}
import org.broadinstitute.dsp.pipelines.config.UltimaGenomicsJointGenotypingConfig
import org.broadinstitute.dsp.pipelines.inputs.{
  UltimaGenomicsJointGenotypingInputs,
  UltimaGenomicsJointGenotypingValidationInputs
}

class UltimaGenomicsJointGenotypingTester(
    testerConfig: UltimaGenomicsJointGenotypingConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override def workflowDir: File =
    CromwellWorkflowTester.PipelineRoot / "broad" / "dna_seq" / "germline" / "joint_genotyping" / "UltimaGenomics"

  override def workflowName: String = "UltimaGenomicsJointGenotyping"

  override protected def workflowInputRoot: File =
    workflowDir / "test_inputs" / testerConfig.category.entryName

  protected val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/joint_genotyping/UltimaGenomics/$testTypeString/$timestamp/"
    )
  protected val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/joint_genotyping/UltimaGenomics/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override def generateRunParameters: Seq[WorkflowRunParameters] = {
    inputFileNames.map { fileName =>
      val inputsName = fileName.replace(".json", "")
      val resultsPath =
        resultsPrefix.resolve(s"$inputsName/")
      val truthPath = truthPrefix.resolve(s"$inputsName/")

      WorkflowRunParameters(
        id = s"${envString}_$inputsName",
        workflowInputs = getInputContents(fileName),
        resultsCloudPath = resultsPath,
        truthCloudPath = truthPath
      )
    }
  }

  override protected val validationWorkflowName: String =
    "VerifyUltimaGenomicsJointGenotyping"

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath

    def getMatchingFiles(suffix: String): Seq[String] = {
      ioUtil
        .listGoogleObjects(truthCloudPath)
        .filter(_.getPath.endsWith(suffix))
        .map(uriToFilename)
    }
    val callsetName =
      new UltimaGenomicsJointGenotypingInputs(
        workflowTest.runParameters.workflowInputs)
        .getCallsetName(workflowName)

    val metricsFileNames = Seq(
      s"$callsetName.variant_calling_detail_metrics",
      s"$callsetName.variant_calling_summary_metrics"
    )

    val interval_lists = getMatchingFiles("interval_list")
    val vcfs = getMatchingFiles("vcf.gz")
    val indexes = getMatchingFiles("vcf.gz.tbi")

    val validationInputs = UltimaGenomicsJointGenotypingValidationInputs(
      testMetrics = metricsFileNames.map(resultsCloudPath.resolve),
      truthMetrics = metricsFileNames.map(truthCloudPath.resolve),
      testFingerprint =
        resultsCloudPath.resolve(s"$callsetName.fingerprintcheck"),
      truthFingerprint =
        truthCloudPath.resolve(s"$callsetName.fingerprintcheck"),
      testVcfs = vcfs.map(resultsCloudPath.resolve),
      truthVcfs = vcfs.map(truthCloudPath.resolve),
      testVcfIndexes = indexes.map(resultsCloudPath.resolve),
      truthVcfIndexes = indexes.map(truthCloudPath.resolve),
      testIntervals = interval_lists.map(resultsCloudPath.resolve),
      truthIntervals = interval_lists.map(truthCloudPath.resolve)
    )
    UltimaGenomicsJointGenotypingValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
