package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import io.circe.{Json, JsonObject}
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{CramToUnmappedBamsInputs, CramToUnmappedBamsValidationInputs, GermlineSingleSampleValidationInputs, ReprocessingInputs}

class CramToUnmappedBamsTester(testerConfig: CramToUnmappedBamsConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  val dataTypePrefix: String = dataTypePrefix(testerConfig.dataType)
  val dataTypeString: String = testerConfig.dataType.entryName.toLowerCase

  override val workflowName: String = s"CramToUnmappedBams"

  val workflowDir
    : File = CromwellWorkflowTester.DsdePipelinesRoot / "pipelines" / "broad" / "reprocessing" / "cram_to_unmapped_bams"

  lazy val localValidationWdlPath: File =
    CromwellWorkflowTester.DsdePipelinesRoot / "verification" / "VerifyCramToUnmappedBams.wdl"

  private val verifyWorkflowName = "VerifyCramToUnmappedBams"

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/cram_to_unmapped_bams/$dataTypeString/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/cram_to_unmapped_bams/$dataTypeString/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val cramToUnmappedBamsInputs = new CramToUnmappedBamsInputs(
      workflowTest.runParameters.workflowInputs
    )
    val outputBaseName = cramToUnmappedBamsInputs.getBaseFileName(workflowName)
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val validationInputs = CramToUnmappedBamsValidationInputs(
    )
    CramToUnmappedBamsValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)

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
        s"$verifyWorkflowName.bam_pairs",
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
        s"gs://broad-gotc-test-storage/cram_to_unmapped_bams/$dataTypeString/$testTypeString/truth/${testerConfig.truthBranch}/$projectSampleAlias/"
      )
      .resolve(fileName)

}
