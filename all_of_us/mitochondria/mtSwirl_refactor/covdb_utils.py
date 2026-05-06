"""Utilities for working with the v2 HDF5 mtDNA coverage DB (coverage.h5)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

try:
    import h5py  # type: ignore
except Exception as e:  # pragma: no cover
    raise RuntimeError(
        "Missing dependency 'h5py'. Install it in the runtime image/env used for v2."
    ) from e


@dataclass(frozen=True)
class CovDBIndex:
    sample_to_index: Dict[str, int]
    pos_to_index: Dict[int, int]


def open_covdb_index(h5_path: str) -> CovDBIndex:
    """Open `coverage.h5` and build lookup maps.

    Expected datasets:
    - /sample_id: array[str] length n_samples
    - /pos: int array length n_positions
    """
    with h5py.File(h5_path, "r") as h5:
        sample_ids = h5["sample_id"][...]
        pos = h5["pos"][...]

    # h5py may return bytes for vlen strings in some builds; normalize.
    sample_ids_norm: List[str] = []
    for x in sample_ids:
        if isinstance(x, bytes):
            sample_ids_norm.append(x.decode("utf-8"))
        else:
            sample_ids_norm.append(str(x))

    sample_to_index = {s: i for i, s in enumerate(sample_ids_norm)}
    pos_to_index = {int(p): i for i, p in enumerate(pos.tolist())}

    return CovDBIndex(sample_to_index=sample_to_index, pos_to_index=pos_to_index)


def read_covdb_block(
    *,
    h5_path: str,
    sample_indices: np.ndarray,
    pos_indices: np.ndarray,
) -> np.ndarray:
    """Read coverage for a set of samples and positions.

    Returns a dense array with shape (n_samples, n_pos).

    Notes:
    - h5py requires index arrays to be strictly increasing. We therefore
      sort/unique the requested indices and restore the original order
      after reading.
    - This uses numpy advanced indexing in-memory. For large arrays you want to
      keep pos_indices to a contiguous block and sample_indices sorted.
    """
    sample_indices = np.asarray(sample_indices)
    pos_indices = np.asarray(pos_indices)

    # h5py requires strictly increasing index arrays. Use unique (sorted) indices
    # and re-expand to the original order after reading.
    sample_unique, sample_inv = np.unique(sample_indices, return_inverse=True)
    pos_unique, pos_inv = np.unique(pos_indices, return_inverse=True)

    with h5py.File(h5_path, "r") as h5:
        dset = h5["coverage"]
        # Advanced indexing: dset[sample_indices, :][:, pos_indices]
        # h5py supports 1D index arrays for each axis.
        arr = dset[sample_unique, :][:, pos_unique]

    # Restore requested order (including any duplicate indices).
    arr = arr[sample_inv, :][:, pos_inv]
    return np.asarray(arr)
