#!/usr/bin/env python3
"""Build an on-disk mtDNA coverage database (v2) from per-sample coverage TSVs.

Contract (matches v1 `coverage_merging` + `add_coverage_annotations` outputs):

Inputs
------
- A TSV (the pipeline "processed TSV") that contains at least:
  - `s`: sample ID
  - `coverage`: path to per-sample base-level coverage metrics TSV

- Each per-sample coverage TSV is expected to be a tab-delimited table with
  at least columns:
  - chrom (string)
  - pos   (int; 1-based)
  - x     (coverage int)

The v1 pipeline reads these files via `hl.import_matrix_table` and renames
column `x` to `coverage`.

Outputs
-------
- `coverage.h5`
  - dataset `/coverage`: shape (n_samples, n_positions) dtype uint16/uint32
  - dataset `/sample_id`: shape (n_samples,) variable-length string
  - dataset `/pos`: shape (n_positions,) int32 (positions)
  - attrs: reference, chrom, position_start, position_end

- `coverage_summary.tsv`
  - columns: pos, mean, median, over_100, over_1000

Notes
-----
- Exactness: mean, median, and threshold fractions are computed exactly.
- Median definition: we match the usual median of a multiset. For even N, we
  use the average of the two middle values and then cast to int32 to match the
  Hail schema (`median: int32`).
  This is the most likely behavior of `hl.median` on integer arrays.
  If parity testing reveals a different rounding rule, we can adjust the
  even-N path.

This script intentionally avoids Hail/Spark.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import gcsfs  # type: ignore
except Exception:  # pragma: no cover
    gcsfs = None

try:
    import h5py  # type: ignore
except Exception as e:  # pragma: no cover
    raise RuntimeError(
        "Missing dependency 'h5py'. Install it in the runtime image/env used for v2."
    ) from e


@dataclass(frozen=True)
class CoverageDBConfig:
    reference: str = "GRCh38"
    chrom: str = "chrM"
    position_start: int = 1
    position_end: int = 16569

    @property
    def n_positions(self) -> int:
        return self.position_end - self.position_start + 1


logger = logging.getLogger("coverage_db")


def _read_processed_tsv(input_tsv: str, sample_col: str, coverage_col: str) -> List[Tuple[str, str]]:
    with open(input_tsv, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in (sample_col, coverage_col) if c not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"Input TSV is missing required columns {missing}. Found columns={reader.fieldnames}"
            )

        out: List[Tuple[str, str]] = []
        for row in reader:
            s = row[sample_col]
            cov_path = row[coverage_col]
            if s is None or s == "":
                raise ValueError("Encountered empty sample id in processed TSV")
            if cov_path is None or cov_path == "":
                # v1 filters out rows with missing coverage metrics earlier; keep strict here.
                raise ValueError(f"Missing coverage path for sample {s}")
            out.append((s, cov_path))

    if len(out) == 0:
        raise ValueError("No rows found in processed TSV")
    return out


def _infer_uint_dtype(max_value: int) -> np.dtype:
    if max_value <= np.iinfo(np.uint16).max:
        return np.uint16
    return np.uint32


def _parse_coverage_tsv(
    path: str,
    expected_chrom: Optional[str],
    position_start: int,
    position_end: int,
    *,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
    cov_col: str = "x",
) -> np.ndarray:
    """Parse one per-sample coverage TSV into a dense 1D array of length n_positions."""
    n_pos = position_end - position_start + 1
    out = np.zeros(n_pos, dtype=np.int64)

    # Support gs:// URIs (Terra often leaves these as remote paths).
    if path.startswith("gs://"):
        if gcsfs is None:
            raise RuntimeError(
                "Encountered a gs:// coverage path but dependency 'gcsfs' is not installed in the runtime image. "
                "Install gcsfs or localize the TSVs before running."
            )
        fs = gcsfs.GCSFileSystem()  # relies on application default credentials in Terra/Batch
        f = fs.open(path, "rt")
    else:
        f = open(path, "r", newline="")

    with f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Coverage TSV {path} appears to have no header")

        # Some AoU pipelines name coverage column 'coverage' instead of 'x'.
        # Allow either; prefer the explicitly requested cov_col if present.
        if cov_col in reader.fieldnames:
            cov_col_eff = cov_col
        elif cov_col == "x" and "coverage" in reader.fieldnames:
            cov_col_eff = "coverage"
        else:
            raise ValueError(
                f"Coverage TSV {path} missing columns {['x']} (or fallback 'coverage'). Found columns={reader.fieldnames}"
            )

        required = [chrom_col, pos_col, cov_col_eff]
        missing = [c for c in required if c not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"Coverage TSV {path} missing columns {missing}. Found columns={reader.fieldnames}"
            )

        for row in reader:
            chrom = row[chrom_col]
            if expected_chrom is not None and chrom != expected_chrom:
                # Be strict: v1 uses chrom from file and then maps to GRCh38 locus.
                raise ValueError(
                    f"Coverage TSV {path} has chrom={chrom}, expected {expected_chrom}"
                )
            pos = int(row[pos_col])
            if pos < position_start or pos > position_end:
                continue
            cov = int(row[cov_col_eff])
            out[pos - position_start] = cov

    return out


def _median_int32_from_int_array(values: np.ndarray) -> np.int32:
    """Compute the median of 1D integer array as an int32.

    This matches Hail's `hl.median` semantics for collections of numbers.

    In particular, for even N, Hail returns the *lower* of the two middle
    values (a discrete median), NOT the average. See Hail docs example:
    `hl.median([1, 3, 5, 6, 7, 9]) == 5`.
    """
    if values.ndim != 1:
        raise ValueError("values must be 1D")
    n = int(values.shape[0])
    if n == 0:
        raise ValueError("Cannot compute median of empty array")

    # Work on a copy because partition mutates.
    v = values.astype(np.int64, copy=True)

    # Hail uses a discrete median: element at floor((n-1)/2) in the sorted array.
    k = (n - 1) // 2
    v_part = np.partition(v, k)
    return np.int32(v_part[k])


def build_coverage_db(
    *,
    samples: Sequence[Tuple[str, str]],
    out_h5: str,
    out_summary_tsv: str,
    config: CoverageDBConfig,
    batch_size: int,
    position_block_size: int,
    compression: Optional[str],
    compression_level: Optional[int],
    dtype: Optional[str],
    expected_chrom: Optional[str],
    skip_summary: bool = False,
) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(out_h5)) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(out_summary_tsv)) or ".", exist_ok=True)

    n_samples = len(samples)
    n_pos = config.n_positions

    logger.info(
        "Starting build_coverage_db: n_samples=%d n_positions=%d out_h5=%s out_summary_tsv=%s skip_summary=%s",
        n_samples,
        n_pos,
        out_h5,
        out_summary_tsv,
        skip_summary,
    )

    # We determine dtype either from user request or by scanning maxima per-batch.
    user_dtype: Optional[np.dtype] = None
    if dtype is not None:
        if dtype.lower() in ("uint16", "u16"):
            user_dtype = np.uint16
        elif dtype.lower() in ("uint32", "u32"):
            user_dtype = np.uint32
        else:
            raise ValueError(f"Unsupported dtype: {dtype}")

    pos_vec = np.arange(config.position_start, config.position_end + 1, dtype=np.int32)

    # Streaming accumulators for summary stats (only if requested).
    sum_cov: Optional[np.ndarray]
    count_gt_100: Optional[np.ndarray]
    count_gt_1000: Optional[np.ndarray]
    if skip_summary:
        sum_cov = None
        count_gt_100 = None
        count_gt_1000 = None
    else:
        sum_cov = np.zeros(n_pos, dtype=np.float64)
        count_gt_100 = np.zeros(n_pos, dtype=np.int64)
        count_gt_1000 = np.zeros(n_pos, dtype=np.int64)

    # Create HDF5 file and datasets.
    with h5py.File(out_h5, "w") as h5:
        h5.attrs["reference"] = config.reference
        h5.attrs["chrom"] = config.chrom
        h5.attrs["position_start"] = config.position_start
        h5.attrs["position_end"] = config.position_end

        str_dt = h5py.string_dtype(encoding="utf-8")
        h5.create_dataset("/sample_id", shape=(n_samples,), dtype=str_dt)
        h5.create_dataset("/pos", data=pos_vec, dtype=np.int32)

        # We'll create /coverage with the chosen dtype after the first batch if needed.
        cov_ds = None

        max_seen = 0
        for start in range(0, n_samples, batch_size):
            end = min(n_samples, start + batch_size)
            batch = samples[start:end]

            logger.info(
                "Parsing coverage TSVs: samples %d..%d of %d",
                start,
                end,
                n_samples,
            )

            # Parse batch into an int64 matrix (batch_samples, n_pos)
            mat = np.zeros((len(batch), n_pos), dtype=np.int64)
            for i, (sid, cov_path) in enumerate(batch):
                mat[i, :] = _parse_coverage_tsv(
                    cov_path,
                    expected_chrom,
                    config.position_start,
                    config.position_end,
                )
                h5["/sample_id"][start + i] = sid

            batch_max = int(mat.max()) if mat.size else 0
            max_seen = max(max_seen, batch_max)

            if cov_ds is None:
                final_dtype = user_dtype if user_dtype is not None else _infer_uint_dtype(max_seen)
                logger.info("Creating HDF5 dataset /coverage dtype=%s", str(final_dtype))
                cov_ds = h5.create_dataset(
                    "/coverage",
                    shape=(n_samples, n_pos),
                    dtype=final_dtype,
                    chunks=(min(batch_size, n_samples), min(position_block_size, n_pos)),
                    compression=compression,
                    compression_opts=compression_level,
                )

            cov_ds[start:end, :] = mat.astype(cov_ds.dtype, copy=False)

            if (start == 0) or (end == n_samples) or (start // batch_size) % 10 == 0:
                logger.info("Wrote HDF5 rows %d..%d", start, end)

            # Update streaming summary stats.
            if not skip_summary:
                assert sum_cov is not None
                assert count_gt_100 is not None
                assert count_gt_1000 is not None
                sum_cov += mat.sum(axis=0)
                count_gt_100 += (mat > 100).sum(axis=0)
                count_gt_1000 += (mat > 1000).sum(axis=0)

        assert cov_ds is not None

        if not skip_summary:
            assert sum_cov is not None
            assert count_gt_100 is not None
            assert count_gt_1000 is not None

            logger.info("Computing per-position summary stats (mean/median/over_100/over_1000)")

            # Compute mean and threshold fractions.
            mean_cov = sum_cov / float(n_samples)
            over_100 = count_gt_100 / float(n_samples)
            over_1000 = count_gt_1000 / float(n_samples)

            # Compute exact medians blockwise across positions.
            medians = np.zeros(n_pos, dtype=np.int32)
            for p0 in range(0, n_pos, position_block_size):
                p1 = min(n_pos, p0 + position_block_size)
                if (p0 == 0) or (p1 == n_pos) or (p0 // position_block_size) % 50 == 0:
                    logger.info("Median progress: positions %d..%d of %d", p0, p1, n_pos)

                block = cov_ds[:, p0:p1]
                # Convert to int64 for exact selection.
                block = block.astype(np.int64, copy=False)
                for j in range(p1 - p0):
                    medians[p0 + j] = _median_int32_from_int_array(block[:, j])

            logger.info("Writing summary TSV: %s", out_summary_tsv)
            with open(out_summary_tsv, "w", newline="") as out:
                w = csv.writer(out, delimiter="\t")
                w.writerow(["pos", "mean", "median", "over_100", "over_1000"])
                for i, pos in enumerate(pos_vec.tolist()):
                    w.writerow([
                        pos,
                        float(mean_cov[i]),
                        int(medians[i]),
                        float(over_100[i]),
                        float(over_1000[i]),
                    ])
        else:
            logger.info("skip_summary=true; not computing/writing coverage_summary.tsv")

    logger.info("Done build_coverage_db")


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-tsv", required=True, help="Processed TSV containing sample and coverage paths")
    p.add_argument("--out-h5", default="coverage.h5", help="Output HDF5 path")
    p.add_argument("--out-summary-tsv", default="coverage_summary.tsv", help="Output summary TSV")

    p.add_argument(
        "--skip-summary",
        action="store_true",
        help="If set, skip computing/writing per-position summary metrics (coverage_summary.tsv).",
    )

    p.add_argument("--sample-col", default="s", help="Sample ID column in processed TSV")
    p.add_argument("--coverage-col", default="coverage", help="Coverage path column in processed TSV")

    p.add_argument("--position-start", type=int, default=1)
    p.add_argument("--position-end", type=int, default=16569)
    p.add_argument("--reference", default="GRCh38")
    p.add_argument("--chrom", default="chrM")

    p.add_argument("--expected-chrom", default=None, help="If set, validate coverage TSV chrom equals this")

    p.add_argument("--batch-size", type=int, default=2000, help="Samples per write batch")
    p.add_argument(
        "--position-block-size",
        type=int,
        default=256,
        help="Positions per median computation block",
    )

    p.add_argument(
        "--compression",
        default="lzf",
        help="HDF5 compression (e.g. lzf, gzip, or empty for none)",
    )
    p.add_argument(
        "--compression-level",
        type=int,
        default=None,
        help="Compression level (only for some codecs like gzip)",
    )
    p.add_argument(
        "--dtype",
        default=None,
        help="Force coverage dtype (uint16 or uint32). If omitted, inferred from data.",
    )

    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    )

    compression = args.compression
    if compression is not None and compression.strip() == "":
        compression = None

    cfg = CoverageDBConfig(
        reference=args.reference,
        chrom=args.chrom,
        position_start=args.position_start,
        position_end=args.position_end,
    )

    samples = _read_processed_tsv(args.input_tsv, args.sample_col, args.coverage_col)

    build_coverage_db(
        samples=samples,
        out_h5=args.out_h5,
        out_summary_tsv=args.out_summary_tsv,
        config=cfg,
        batch_size=args.batch_size,
        position_block_size=args.position_block_size,
        compression=compression,
        compression_level=args.compression_level,
        dtype=args.dtype,
        expected_chrom=args.expected_chrom,
        skip_summary=args.skip_summary,
    )

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
