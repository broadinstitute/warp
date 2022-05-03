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

class BroadInternalUGWrapperTester(
    testerConfig: BroadInternalUGWrapperConfig)(
    implicit am: ActorMaterializer,
    as: ActorSystem
) extends ValidationWdlTester(testerConfig) {
    
}