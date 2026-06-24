# MMIDAS_DataPrep Changelog

## 1.0.0
2026-06-24

- Initial release. Wraps `01_data_prep.py` to load raw Mouse Smart-seq ALM/VISp count matrices from the Allen Brain Atlas, filter to neuronal cells, normalize to log-CPM, and write a single AnnData `.h5ad` file for downstream MMIDAS training.
