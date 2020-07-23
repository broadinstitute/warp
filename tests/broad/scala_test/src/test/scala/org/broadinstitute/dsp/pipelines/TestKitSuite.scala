package org.broadinstitute.dsp.pipelines

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import akka.testkit.TestKit
import org.scalatest.{AsyncFlatSpecLike, BeforeAndAfterAll, Matchers}

/**
  * A mix of Akka TestKit with ScalaTest mixed in to clean up the actor system.
  *
  * @param actorSystemName The name of the actor system.
  */
abstract class TestKitSuite(actorSystemName: String)
    extends TestKit(ActorSystem(actorSystemName))
    with AsyncFlatSpecLike
    with Matchers
    with BeforeAndAfterAll {

  /**
    * Most Akka-streams APIs require an implicit Materializer, so we provide
    * one here tied to the ActorSystem managed by the TestKit for convenience.
    */
  implicit val materializer: ActorMaterializer = ActorMaterializer()

  override protected def afterAll(): Unit = {
    shutdown()
  }
}
