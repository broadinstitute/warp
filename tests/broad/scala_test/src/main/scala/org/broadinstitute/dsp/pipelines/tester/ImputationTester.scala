package org.broadinstitute.dsp.pipelines.tester

import java.net.URI

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import better.files.File
import org.broadinstitute.dsp.pipelines.batch.WorkflowTest
import org.broadinstitute.dsp.pipelines.config._
import org.broadinstitute.dsp.pipelines.inputs.{
  ImputationInputs,
  ImputationValidationInputs
}

class ImputationTester(testerConfig: ImputationConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {

  override val workflowName: String = "Imputation"

  val workflowDir
    : File = CromwellWorkflowTester.PipelineRoot / "broad" / "arrays" / "imputation"

  override protected val validationWorkflowName: String =
    "VerifyImputation"

  // Validation uses the same options as the arrays workflow
  override protected lazy val validationWdlOptions: String = {
    readTestOptions(releaseDir, env)
  }

  protected lazy val resultsPrefix: URI = {
    URI.create(
      s"gs://broad-gotc-test-results/$envString/imputation/$testTypeString/$timestamp/"
    )
  }

  protected lazy val truthPrefix: URI =
    URI.create(
      s"gs://broad-gotc-test-storage/imputation/$testTypeString/truth/${testerConfig.truthBranch}/"
    )

  override protected def buildValidationWdlInputs(
      workflowTest: WorkflowTest
  ): String = {
    val imputationInputs = new ImputationInputs(
      workflowTest.runParameters.workflowInputs
    )
    val outputBaseName =
      imputationInputs.getBasename(workflowName)
    val outputCallsetName =
      imputationInputs.getOutputCallsetName(workflowName)
    val splitOutputToSingleSample =
      imputationInputs.getSplitOutputToSingleSample(workflowName)
    val singleSampleVcfs =
      imputationInputs.getSingleSampleVcfs(workflowName)
    val singleSampleVcfsIndices =
      imputationInputs.getSingleSampleVcfsIndices(workflowName)
    val multiSampleVcf =
      imputationInputs.getMultiSampleVcf(workflowName)
    val multiSampleVcfIndices =
      imputationInputs.getMultiSampleVcfIndices(workflowName)
    val resultsCloudPath =
      workflowTest.runParameters.resultsCloudPath
    val truthCloudPath = workflowTest.runParameters.truthCloudPath

    val validationInputs = ImputationValidationInputs(
      split_output_to_single_sample = splitOutputToSingleSample,
      output_callset_name = outputCallsetName,
      input_single_sample_vcfs = singleSampleVcfs,
      input_single_sample_vcfs_indices = singleSampleVcfsIndices,
      input_multi_sample_vcf = multiSampleVcf,
      input_multi_sample_vcf_index = multiSampleVcfIndices,
      test_metrics = Array(
        resultsCloudPath.resolve(s"${outputBaseName}_chunk_info.tsv"),
        resultsCloudPath.resolve(s"${outputBaseName}_failed_chunks.tsv"),
        resultsCloudPath.resolve(s"n_failed_chunks.txt"),
        resultsCloudPath.resolve(
          s"${outputBaseName}_aggregated_imputation_metrics.tsv")
      ),
      truth_metrics = Array(
        truthCloudPath.resolve(s"${outputBaseName}_chunk_info.tsv"),
        truthCloudPath.resolve(s"${outputBaseName}_failed_chunks.tsv"),
        truthCloudPath.resolve(s"n_failed_chunks.txt"),
        truthCloudPath.resolve(
          s"${outputBaseName}_aggregated_imputation_metrics.tsv")
      ),
      test_vcf = resultsCloudPath.resolve(s"$outputBaseName.vcf.gz"),
      truth_vcf = truthCloudPath.resolve(s"$outputBaseName.vcf.gz"),
      test_vcf_index = resultsCloudPath.resolve(s"$outputBaseName.vcf.gz.tbi"),
      truth_vcf_index = truthCloudPath.resolve(s"$outputBaseName.vcf.gz.tbi")
    )
    ImputationValidationInputs
      .marshall(validationInputs)
      .printWith(implicitly)
  }
}
