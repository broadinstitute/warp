"""
WARP Test Harness - Module entry point.

Allows running the test harness as a Python module:
  python -m scripts.test_harness mt_coverage_merge
"""

import sys
from .cli import main

if __name__ == "__main__":
    sys.exit(main())
