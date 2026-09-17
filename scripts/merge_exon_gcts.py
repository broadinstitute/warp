#!/usr/bin/env python3
"""Merge exon GCT files using the exon IDs shared by every input file."""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import TextIO


class MergeError(ValueError):
    """Report an invalid or incompatible GCT file."""


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def read_header(handle: TextIO, path: Path) -> tuple[int, tuple[str, ...]]:
    if handle.readline().rstrip("\r\n") != "#1.2":
        raise MergeError(f"{path}: expected GCT version #1.2")

    dimensions = handle.readline().rstrip("\r\n").split("\t")
    if len(dimensions) != 2:
        raise MergeError(f"{path}: invalid GCT dimensions line")
    try:
        row_count, sample_count = (int(value) for value in dimensions)
    except ValueError as error:
        raise MergeError(f"{path}: GCT dimensions must be integers") from error
    if row_count < 0 or sample_count < 1:
        raise MergeError(f"{path}: invalid GCT dimensions")

    columns = handle.readline().rstrip("\r\n").split("\t")
    if len(columns) != sample_count + 2 or columns[:2] != ["Name", "Description"]:
        raise MergeError(f"{path}: invalid GCT column header")
    sample_names = tuple(columns[2:])
    if any(not sample_name for sample_name in sample_names):
        raise MergeError(f"{path}: empty sample name")
    return row_count, sample_names


def read_features(path: Path) -> tuple[int, tuple[str, ...], dict[str, str]]:
    feature_ids: list[str] = []
    descriptions: dict[str, str] = {}
    with open_text(path, "rt") as handle:
        row_count, sample_names = read_header(handle, path)
        for row_number in range(1, row_count + 1):
            fields = handle.readline().rstrip("\r\n").split("\t")
            if len(fields) != len(sample_names) + 2:
                raise MergeError(f"{path}: invalid field count at row {row_number}")
            feature_id, description = fields[:2]
            if not feature_id:
                raise MergeError(f"{path}: empty feature ID at row {row_number}")
            if feature_id in descriptions:
                raise MergeError(f"{path}: duplicate feature ID: {feature_id}")
            feature_ids.append(feature_id)
            descriptions[feature_id] = description
        if handle.readline():
            raise MergeError(f"{path}: contains rows beyond its declared row count")
    return row_count, tuple(feature_ids), descriptions


def read_values(path: Path, wanted_features: set[str]) -> dict[str, tuple[str, ...]]:
    values: dict[str, tuple[str, ...]] = {}
    with open_text(path, "rt") as handle:
        row_count, sample_names = read_header(handle, path)
        for row_number in range(1, row_count + 1):
            fields = handle.readline().rstrip("\r\n").split("\t")
            if len(fields) != len(sample_names) + 2:
                raise MergeError(f"{path}: invalid field count at row {row_number}")
            feature_id = fields[0]
            if feature_id in wanted_features:
                values[feature_id] = tuple(fields[2:])
        if handle.readline():
            raise MergeError(f"{path}: contains rows beyond its declared row count")
    if len(values) != len(wanted_features):
        missing = wanted_features.difference(values)
        raise MergeError(f"{path}: missing shared feature IDs: {sorted(missing)[:5]}")
    return values


def merge(paths: list[Path], output: Path) -> None:
    if not paths:
        raise MergeError("at least one input GCT is required")

    metadata = [read_features(path) for path in paths]
    sample_names = []
    for path, (_, names, _) in zip(paths, metadata):
        sample_names.extend(names)
    if len(sample_names) != len(set(sample_names)):
        raise MergeError("duplicate sample names across input GCT files")

    shared_features = set(metadata[0][1])
    for _, feature_ids, _ in metadata[1:]:
        shared_features.intersection_update(feature_ids)
    if not shared_features:
        raise MergeError("input GCT files have no shared feature IDs")

    ordered_features = [
        feature_id for feature_id in metadata[0][1] if feature_id in shared_features
    ]
    descriptions = metadata[0][2]
    for path, (_, _, current_descriptions) in zip(paths[1:], metadata[1:]):
        for feature_id in ordered_features:
            if descriptions[feature_id] != current_descriptions[feature_id]:
                raise MergeError(
                    f"{path}: description differs for feature ID {feature_id}"
                )

    values_by_file = [
        read_values(path, shared_features)
        for path in paths
    ]
    with open_text(output, "wt") as handle:
        handle.write("#1.2\n")
        handle.write(f"{len(ordered_features)}\t{len(sample_names)}\n")
        handle.write("\t".join(("Name", "Description", *sample_names)) + "\n")
        for feature_id in ordered_features:
            values = [
                value
                for file_values in values_by_file
                for value in file_values[feature_id]
            ]
            handle.write(
                "\t".join((feature_id, descriptions[feature_id], *values)) + "\n"
            )

    dropped_counts = [row_count - len(shared_features) for row_count, _, _ in metadata]
    print(
        f"Merged {len(paths)} files: {len(ordered_features)} shared exons; "
        f"dropped {sum(dropped_counts)} file-specific rows.",
        file=sys.stderr,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="Input .gct or .gct.gz files")
    parser.add_argument("-o", "--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        merge(args.inputs, args.output)
    except (MergeError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())