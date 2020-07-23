package org.broadinstitute.dsp.pipelines.util

import java.net.URI

import cromwell.api.model.{WorkflowMetadata, WorkflowOutputs}
import io.circe.{Json, JsonObject}
import io.circe.parser.parse
import spray.json.{JsArray, JsObject, JsString, JsValue}

import scala.util.Try

trait CromwellUtils {

  implicit class FromWorkflowOutputs(workflowOutputs: WorkflowOutputs) {

    def getStringFromOutputs(output: String): Try[String] = {
      (for {
        parsed <- parse(workflowOutputs.outputs.compactPrint).toOption
        foo <- parsed.asObject
        bar <- foo(output)
        baz <- bar.asString
      } yield baz)
        .toRight(
          new IllegalArgumentException(
            s"$output was not found in the cromwell outputs ${workflowOutputs.outputs}"
          )
        )
        .toTry
    }

    def getAllCloudOutputs: Seq[URI] = {
      def getCloudOutput(jsValue: JsValue): Option[String] =
        jsValue match {
          case string: JsString if string.value.startsWith("gs://") =>
            Option(string.value)
          case _ => None
        }
      def getCloudOutputs(jsValue: JsValue): Iterable[Option[String]] =
        jsValue match {
          case string: JsString => Iterable(getCloudOutput(string))
          case array: JsArray   => array.elements.map(getCloudOutput)
          case obj: JsObject    => obj.fields.values.flatMap(getCloudOutputs)
          case _                => Iterable.empty
        }

      workflowOutputs.outputs.asJsObject.fields.values
        .flatMap(getCloudOutputs)
        .flatten
        .map(URI.create)
        .toSeq
    }
  }

  implicit class FromWorkflowMetadata(workflowMetadata: WorkflowMetadata) {

    def getLogs: Seq[URI] = {
      parse(workflowMetadata.value)
        .map(getAllLogs)
        .getOrElse(Seq.empty)
        .map(URI.create)
    }

    private def getAllLogs(json: Json): Seq[String] = {
      def handleObject(jsonObject: JsonObject) =
        Seq.empty ++
          jsonObject("log").flatMap(_.asString) ++
          jsonObject("stderr").flatMap(_.asString) ++
          jsonObject("stdout").flatMap(_.asString) ++
          jsonObject.values.flatMap(getAllLogs)

      def handleArray(jsonArray: Vector[Json]) = jsonArray.flatMap(getAllLogs)

      json.fold(
        jsonNull = Seq.empty,
        jsonBoolean = _ => Seq.empty,
        jsonNumber = _ => Seq.empty,
        jsonString = _ => Seq.empty,
        jsonObject = handleObject,
        jsonArray = handleArray
      )
    }
  }
}
