package org.broadinstitute.dsp.pipelines

import akka.actor.ActorSystem
import akka.stream.ActorMaterializer
import com.typesafe.scalalogging.StrictLogging
import org.broadinstitute.dsp.pipelines.commandline.{Config, ConfigParser}
import org.broadinstitute.dsp.pipelines.tester.CromwellWorkflowTester

import scala.concurrent.Future
import scala.util.{Failure, Success}

object WorkflowTest extends StrictLogging {

  private val parser = new ConfigParser

  def main(args: Array[String]): Unit = {

    parser.parse(args, Config()).fold(sys.exit(1)) { config =>
      implicit val system: ActorSystem = ActorSystem("workflow-test")
      implicit val materializer: ActorMaterializer = ActorMaterializer()
      import system.dispatcher

      sys.addShutdownHook({ val _ = system.terminate() })

      Future(CromwellWorkflowTester(config)).flatMap(_.runTest).onComplete {
        case Success(_) =>
          logger.info("The tests were successful")
          sys.exit(0)
        case Failure(ex) =>
          logger.error("The tests were unsuccessful", ex)
          sys.exit(1)
      }
    }
  }
}
