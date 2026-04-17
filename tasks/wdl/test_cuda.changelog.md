# 1.0.0
2026-04-14 (Date of Last Commit)

* Initial release of GPU/CUDA smoke test WDL for verifying CUDA availability on Cromwell GPU nodes
* Tests PyTorch CUDA availability using the configured docker image, GPU type, and NVIDIA driver version
* Returns stdout log; exits non-zero if CUDA is not available
