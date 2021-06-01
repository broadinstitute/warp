package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.{
  WorkflowRunParameters,
  WorkflowTest
}
import org.broadinstitute.dsp.pipelines.config.GermlineCloudWorkflowConfig
import org.broadinstitute.dsp.pipelines.inputs.{
  ReblockGvcfInputs,
  ReblockGvcfValidationInputs
}

class ReblockGvcfTester(testerConfig: GermlineCloudWorkflowConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends GermlineCloudWorkflowTester(testerConfig) {

  override def workflowDir: File =
    CromwellWorkflowTester.PipelineRoot / "broad" / "dna_seq" / "germline" / "joint_genotyping" / "reblocking"

  override def workflowName: String = "ReblockGVCF"

  override protected val validationWorkflowName: String = s"VerifyGvcf"

  override protected def workflowInputRoot: File =
    workflowDir / "test_inputs" / dataTypeString / testerConfig.category.entryName

  protected val resultsPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-results/$envString/reblock_gvcf/$dataTypeString/$testTypeString/$timestamp/"
    )
  protected val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/reblock_gvcf/$dataTypeString/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override def generateRunParameters: Seq[WorkflowRunParameters] = {
    inputFileNames.map { fileName =>
      val inputsName = fileName.replace(".json", "")
      val resultsPath =
        resultsPrefix.resolve(s"$inputsName/")
      val truthPath = truthPrefix.resolve(s"$inputsName/")
      val inputs =
        getInputContents(fileName).replace("{TRUTH_BRANCH}",
                                           testerConfig.truthBranch)

      WorkflowRunParameters(
        id = s"${envString}_$inputsName",
        workflowInputs = inputs,
        resultsCloudPath = resultsPath,
        truthCloudPath = truthPath
      )
    }
  }

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath

    val gvcfBasename =
      new ReblockGvcfInputs(workflowTest.runParameters.workflowInputs)
        .getGvcfBasename(workflowName)

    val validationInputs = ReblockGvcfValidationInputs(
      testGvcf = resultsCloudPath.resolve(s"$gvcfBasename.reblocked.g.vcf.gz"),
      truthGvcf = truthCloudPath.resolve(s"$gvcfBasename.reblocked.g.vcf.gz"),
    )
    ReblockGvcfValidationInputs.marshall(validationInputs).printWith(implicitly)
  }
}
