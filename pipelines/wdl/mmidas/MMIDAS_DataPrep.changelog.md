# 1.0.1
2026-08-06 (Date of Last Commit)

* Fixed a misleading normalization sanity check in 01_data_prep.py. It printed per-cell sums of log1p values against an expectation of "~log(1e6+1) ~ 13.8", which is not what that sum equals, so a correct run reported values in the thousands and looked broken. It now checks the invariant that actually holds -- inverting log1p recovers CPM, so each non-empty cell sums back to 1e6 -- warns if that deviates by more than 0.1%, and reports the maximum single log-CPM value against its log1p(1e6) bound plus the all-zero cell count. No change to pipeline outputs.
* Corrected the selected_genes gene count in the WDL comments and parameter_meta from ~1,252 to 5,032, matching genes_SS_ALM-VISp.csv.
* The Docker image now installs the mmidas package from a pinned revision of the fork rather than a copy of a local directory, so any published image can be rebuilt from source. See "Docker image provenance" in dashboard.md.
* Updated the Docker image to us.gcr.io/broad-gotc-prod/mmidas:1.0.0-0.1.0-1786046379

# 1.0.0
2026-06-24 (Date of Last Commit)

* Initial release. Wraps 01_data_prep.py to load raw Mouse Smart-seq ALM/VISp count matrices from the Allen Brain Atlas, filter to neuronal cells, normalize to log-CPM, and write a single AnnData .h5ad file for downstream MMIDAS training.
