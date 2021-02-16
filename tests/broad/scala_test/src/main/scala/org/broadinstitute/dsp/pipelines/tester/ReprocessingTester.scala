package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import io.circe.{Json, JsonObject}
import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRunParameters,
  WorkflowTest
}
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  GermlineSingleSampleValidationInputs,
  ReprocessingInputs
}

class ReprocessingTester(testerConfig: GermlineCloudWorkflowConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends GermlineSingleSampleTester(testerConfig) {

  override val workflowName: String = s"${dataTypePrefix}Reprocessing"

  override lazy val workflowDir: File =
    CromwellWorkflowTester.PipelineRoot / "broad" / "reprocessing" / dataTypeString

  override protected val validationWorkflowName: String = "VerifyReprocessing"

  override protected lazy val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/reprocessing/$dataTypeString/$testTypeString/$timestamp/"
    )
  override protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/reprocessing/$dataTypeString/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override def generateRunParameters: Seq[WorkflowRunParameters] = {
    super.generateRunParameters.map(
      rp =>
        rp.copy(
          workflowInputs = rp.workflowInputs.replace("{TRUTH_BRANCH}",
                                                     testerConfig.truthBranch)
      )
    )
  }

  override protected def buildValidationWdlInputs(
      germlineSingleSampleTest: WorkflowTest
  ): String = {
    val reprocessingInputs = new ReprocessingInputs(
      germlineSingleSampleTest.runParameters.workflowInputs
    )
    val outputBaseName = reprocessingInputs.getBaseFileName(workflowName)
    val gvcfBaseName = reprocessingInputs.getFinalGvcfBaseName(workflowName)
    val resultsCloudPath =
      germlineSingleSampleTest.runParameters.resultsCloudPath
    val truthCloudPath = germlineSingleSampleTest.runParameters.truthCloudPath
    val metricsFileNames = ioUtil
      .listGoogleObjects(truthCloudPath)
      .filter(_.getPath.endsWith("metrics"))
      .map(uriToFilename)

    val validationInputs = GermlineSingleSampleValidationInputs.marshall(
      GermlineSingleSampleValidationInputs(
        testMetrics = metricsFileNames.map(resultsCloudPath.resolve),
        truthMetrics = metricsFileNames.map(truthCloudPath.resolve),
        testCram = resultsCloudPath.resolve(s"$outputBaseName.cram"),
        testCrai = resultsCloudPath.resolve(s"$outputBaseName.cram.crai"),
        truthCram = truthCloudPath.resolve(s"$outputBaseName.cram"),
        truthCrai = truthCloudPath.resolve(s"$outputBaseName.cram.crai"),
        testGvcf = resultsCloudPath.resolve(s"$gvcfBaseName.g.vcf.gz"),
        truthGvcf = truthCloudPath.resolve(s"$gvcfBaseName.g.vcf.gz")
      ),
      validationWorkflowName
    )

    val revertedBams = ioUtil
      .listGoogleObjects(resultsCloudPath)
      .filter(_.toString.endsWith(".bam"))
      .sorted
    val truthBams = revertedBams
      .map(uri => File(uri.getPath).name)
      .map { fileName =>
        getCloudUnmappedBam(outputBaseName, fileName)
      }
      .sorted

    val added = validationInputs.asObject.fold(JsonObject.empty)(
      _.add(
        s"$validationWorkflowName.bam_pairs",
        Json.arr(
          revertedBams
            .zip(truthBams)
            .map(
              pair =>
                Json.obj(
                  "test_bam" -> Json.fromString(pair._1.toString),
                  "truth_bam" -> Json.fromString(pair._2.toString)
              )
            ): _*
        )
      )
    )
    Json.fromJsonObject(added).toString()
  }

  private def getCloudUnmappedBam(projectSampleAlias: String,
                                  fileName: String): URI =
    URI
      .create(
        s"gs://broad-gotc-test-storage/germline_single_sample/$dataTypeString/$testTypeString/bams/$projectSampleAlias/"
      )
      .resolve(fileName)
}
