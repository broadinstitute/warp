# 1.0.2
2026-11-24 (Date of Last Commit)

* Renamed the PeakCalling task to PeakCallingTask so it no longer shares a name with the PeakCalling workflow; nested task inputs set as `PeakCalling.PeakCalling.<input>` are now `PeakCalling.PeakCallingTask.<input>`
* Fixed WDL type errors reported by the Sprocket linter in shared tasks (for example, typed `cpu` runtime values as `Int` and unwrapped optional values with `select_first`); no functional impact on pipeline outputs

# 1.0.1
2025-08-11 (Date of Last Commit)

* Reorganized pipelines into the wdl directory

# 1.0.0

2025-03-24 (Date of Last Commit)

* The first major version release for the PeakCalling pipeline.