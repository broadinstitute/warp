package org.broadinstitute.dsp.pipelines.util

import java.time._
import java.time.format.DateTimeFormatter

object DateFormatters {

  object FireCloudWorkspaceDateFormatters {

    def dateForTsv(zoneId: ZoneId): String =
      DateTimeFormatter
        .ofPattern("yyyy/MM/dd")
        .format(OffsetDateTime.now(zoneId))

    def dateForTags(zoneId: ZoneId): String =
      DateTimeFormatter
        .ofPattern("MMM d, yyyy")
        .format(LocalDateTime.now(zoneId))
  }
}
