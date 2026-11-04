# aou_9.1.0

2026-11-04 (Date of Last Commit)

* Removed the predefined machine type workflow input, breaking compatibility with input configurations that set it, and changed the default GPU type to nvidia-tesla-t4

# aou_9.0.0

2026-09-30 (Date of Last Commit)

* Updated the WDL to the WDL 1.0 specification
* Added pipeline version tracking
* Added a predefined_machine_type input with a default of g2-standard-16
* Added a gpu_type input with a default of nvidia-l4