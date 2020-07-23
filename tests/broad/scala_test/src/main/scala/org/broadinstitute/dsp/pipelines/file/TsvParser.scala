package org.broadinstitute.dsp.pipelines.file

import better.files.File

object TsvParser {

  def parse(tsv: Seq[String]): Seq[Map[String, String]] = {
    val emptyLinesStripped = tsv.filter(_.length != 0)
    val headers: Seq[String] = emptyLinesStripped.head.split("\t").toSeq
    emptyLinesStripped.tail.map(line => headers.zip(line.split("\t")).toMap)
  }

  def parseExcludeComments(
      tsv: Seq[String],
      commentChar: String
  ): Seq[Map[String, String]] = {
    parse(tsv.filter(!_.startsWith(commentChar)))
  }

  def parseFromInput(
      resourcePath: File,
      transform: String => String = identity
  ): Seq[Map[String, String]] = {
    parse(resourcePath.lines.map(transform).toSeq)
  }
}
