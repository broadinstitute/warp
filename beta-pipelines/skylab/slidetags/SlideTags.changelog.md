# 1.1.0
2025-02-25 (Date of Last Commit)

* Removed boolean variable is_slidetags; no longer needed with new updates
* Refactored the STAR alignment step (STARsoloFastq) in Optimus and removed tasks FastqProcessing and MergeSortBamFiles; we are no longer sharding. We are now running one instance of STAR

