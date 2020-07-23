package org.broadinstitute.dsp.pipelines.tester

import org.broadinstitute.dsp.pipelines.TestKitSuite

class DummyTesterSpec extends TestKitSuite("workflow-test") {

  behavior of "DummyTester"

  it should "run a trivial workflow" in {
    new DummyTester().runTest.map(_ should be(()))
  }

}
