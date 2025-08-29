# 1.0.4
2025-08-26 (Date of Last Commit)

* Update tasks to use maxRetries of 1 instead of 2. 1 retry is sufficient for transient errors and helps reduce costs.

# 1.0.3
2025-07-31 (Date of Last Commit)

* Add maxRetries to tasks to help with performance on Google Batch.

# 1.0.2
2025-04-07 (Date of Last Commit)

* Update Imputation Tasks to not request external IP addresses (use runtime attribute 'noAddress: true').

# 1.0.1
2025-04-01 (Date of Last Commit)

* Update Imputation Tasks to use dockers (ubuntu and tidyverse) that have been moved to GAR.

# 1.0.0
2025-02-24 (Date of Last Commit)

* Initial release of pipeline to calculate the number of samples, i.e. quota used by an imputation service that uses ImputationBeagle.wdl.
