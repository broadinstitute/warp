#!/usr/bin/env python3
"""Merge RNA-SeQC tables without loading complete matrices into memory."""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
import tempfile
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterator, Sequence


class MergeError(ValueError):
    """Report an invalid or incompatible RNA-SeQC input."""


@dataclass(frozen=True)
class GctHeader:
    row_count: int
    sample_names: tuple[str, ...]


@dataclass(frozen=True)
class SampleRecord:
    sample_id: str
    tpm_gct: str
    count_gct: str
    exon_count_gct: str
    metrics_tsv: str
    insert_size_hist: str | None = None


@dataclass(frozen=True)
class SampleManifest:
    records: tuple[SampleRecord, ...]
    include_insert_sizes: bool


REQUIRED_MANIFEST_COLUMNS = (
    "sample_id",
    "tpm_gct",
    "count_gct",
    "exon_count_gct",
    "metrics_tsv",
)
INSERT_SIZE_COLUMN = "insert_size_hist"
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9._-]+$")
SAFE_OUTPUT_PREFIX = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
VALID_GCS_URI = re.compile(r"^gs://[^/\s]+/.+$")


def open_text(path: Path, mode: str) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def read_path_list(path: Path) -> list[Path]:
    paths = [
        Path(line.strip()) for line in path.read_text().splitlines() if line.strip()
    ]
    if not paths:
        raise MergeError(f"Input list is empty: {path}")
    return paths


def read_sample_ids(path: Path) -> tuple[str, ...]:
    sample_ids = tuple(path.read_text().splitlines())
    if not sample_ids:
        raise MergeError(f"Sample ID list is empty: {path}")
    for line_number, sample_id in enumerate(sample_ids, start=1):
        if not sample_id:
            raise MergeError(f"{path}:{line_number}: invalid sample ID: empty value")
        if "\t" in sample_id:
            raise MergeError(f"{path}:{line_number}: invalid sample ID: tab character")
        if sample_id != sample_id.strip():
            raise MergeError(
                f"{path}:{line_number}: invalid sample ID: surrounding whitespace"
            )
    return sample_ids


def parse_gct_header(handle: IO[str], path: Path) -> GctHeader:
    version = handle.readline().rstrip("\r\n")
    if version != "#1.2":
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
    if columns[:2] != ["Name", "Description"]:
        raise MergeError(f"{path}: expected Name and Description columns")
    if len(columns) != sample_count + 2:
        raise MergeError(f"{path}: GCT sample count does not match its header")
    sample_names = tuple(columns[2:])
    if any(not sample_name for sample_name in sample_names):
        raise MergeError(f"{path}: empty sample name")
    return GctHeader(row_count=row_count, sample_names=sample_names)


@contextmanager
def atomic_text_writer(path: Path) -> Iterator[IO[str]]:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".gz" if path.suffix == ".gz" else ".tmp",
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        with open_text(temporary_path, "wt") as handle:
            yield handle
        os.replace(temporary_path, path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def write_lines_atomic(path: Path, values: Sequence[str]) -> None:
    with atomic_text_writer(path) as handle:
        for value in values:
            handle.write(f"{value}\n")


def ensure_unique(values: Sequence[str], label: str) -> None:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    if duplicates:
        joined = ", ".join(sorted(duplicates))
        raise MergeError(f"Duplicate {label}: {joined}")


def read_sample_manifest(path: Path) -> SampleManifest:
    with path.open("rt", newline="") as handle:
        lines = [line.rstrip("\r\n") for line in handle]
    if not lines:
        raise MergeError(f"Sample manifest is empty: {path}")

    header = tuple(lines[0].split("\t"))
    valid_headers = (
        REQUIRED_MANIFEST_COLUMNS,
        (*REQUIRED_MANIFEST_COLUMNS, INSERT_SIZE_COLUMN),
    )
    if header not in valid_headers:
        raise MergeError(
            f"{path}: invalid header; expected five required columns and an optional "
            f"{INSERT_SIZE_COLUMN} column"
        )
    include_insert_sizes = len(header) == 6

    records: list[SampleRecord] = []
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(header):
            raise MergeError(
                f"{path}:{line_number}: expected {len(header)} tab-separated columns"
            )
        sample_id, *uris = fields
        if not SAFE_SAMPLE_ID.fullmatch(sample_id):
            raise MergeError(
                f"{path}:{line_number}: invalid sample ID; use only letters, numbers, "
                "periods, underscores, and hyphens"
            )
        for column, uri in zip(header[1:], uris):
            if uri != uri.strip() or not VALID_GCS_URI.fullmatch(uri):
                raise MergeError(
                    f"{path}:{line_number}: {column} must be a valid "
                    "gs://bucket/object URI without surrounding whitespace"
                )
        records.append(
            SampleRecord(
                sample_id=sample_id,
                tpm_gct=uris[0],
                count_gct=uris[1],
                exon_count_gct=uris[2],
                metrics_tsv=uris[3],
                insert_size_hist=uris[4] if include_insert_sizes else None,
            )
        )

    if not records:
        raise MergeError(f"Sample manifest has no data rows: {path}")
    ensure_unique(tuple(record.sample_id for record in records), "sample IDs")
    return SampleManifest(tuple(records), include_insert_sizes)


def validate_output_prefix(path: Path) -> None:
    values = path.read_text().splitlines()
    if len(values) != 1 or not SAFE_OUTPUT_PREFIX.fullmatch(values[0]):
        raise MergeError(
            "invalid output prefix; start with a letter or number and use only "
            "letters, numbers, periods, underscores, and hyphens"
        )


def validate_sample_manifest(
    path: Path, batch_size: int, prefix_file: Path | None = None
) -> None:
    if batch_size < 1:
        raise MergeError("Batch size must be at least 1")
    if prefix_file is not None:
        validate_output_prefix(prefix_file)
    manifest = read_sample_manifest(path)
    sample_count = len(manifest.records)
    batch_count = (sample_count + batch_size - 1) // batch_size
    write_lines_atomic(Path("sample_count.txt"), (str(sample_count),))
    write_lines_atomic(Path("batch_count.txt"), (str(batch_count),))
    write_lines_atomic(
        Path("include_insert_sizes.txt"),
        ("true" if manifest.include_insert_sizes else "false",),
    )


def local_stage_path(
    staging_directory: Path,
    record_number: int,
    sample_id: str,
    suffix: str,
    uri: str,
) -> str:
    filename = f"{record_number:06d}.{sample_id}.{suffix}"
    if uri.endswith(".gz"):
        filename += ".gz"
    return str(staging_directory / filename)


def prepare_manifest_batch(
    path: Path,
    batch_index: int,
    batch_size: int,
    staging_directory: Path,
    merge_exons: bool = True,
) -> None:
    if batch_index < 0:
        raise MergeError("Batch index must be nonnegative")
    if batch_size < 1:
        raise MergeError("Batch size must be at least 1")
    manifest = read_sample_manifest(path)
    start = batch_index * batch_size
    records = manifest.records[start : start + batch_size]
    if not records:
        raise MergeError(f"Batch {batch_index + 1} is empty")

    file_specs = [
        ("tpm_gct", "gene_tpm.gct", "local_tpm.list"),
        ("count_gct", "gene_reads.gct", "local_count.list"),
        ("metrics_tsv", "metrics.tsv", "local_metrics.list"),
    ]
    if merge_exons:
        file_specs.insert(2, ("exon_count_gct", "exon_reads.gct", "local_exon.list"))
    local_lists: dict[str, list[str]] = {name: [] for _, _, name in file_specs}
    if manifest.include_insert_sizes:
        local_lists["local_insert.list"] = []
    transfers: list[str] = []
    sample_ids: list[str] = []

    for offset, record in enumerate(records):
        record_number = start + offset + 1
        sample_ids.append(record.sample_id)
        for attribute, suffix, list_name in file_specs:
            uri = getattr(record, attribute)
            destination = local_stage_path(
                staging_directory, record_number, record.sample_id, suffix, uri
            )
            transfers.append(f"{uri}\t{destination}")
            local_lists[list_name].append(destination)
        if record.insert_size_hist is not None:
            destination = local_stage_path(
                staging_directory,
                record_number,
                record.sample_id,
                "fragmentSizes.txt",
                record.insert_size_hist,
            )
            transfers.append(f"{record.insert_size_hist}\t{destination}")
            local_lists["local_insert.list"].append(destination)

    write_lines_atomic(Path("transfers.tsv"), transfers)
    write_lines_atomic(Path("batch_sample_ids.list"), sample_ids)
    for list_name, destinations in local_lists.items():
        write_lines_atomic(Path(list_name), destinations)


def merge_gcts(
    input_paths: Sequence[Path],
    output_path: Path,
    sample_output_path: Path | None = None,
    expected_samples: Sequence[str] | None = None,
    sample_names_override: Sequence[str] | None = None,
) -> tuple[str, ...]:
    with ExitStack() as stack:
        handles = [stack.enter_context(open_text(path, "rt")) for path in input_paths]
        headers = [
            parse_gct_header(handle, path) for handle, path in zip(handles, input_paths)
        ]

        row_count = headers[0].row_count
        for header, path in zip(headers[1:], input_paths[1:]):
            if header.row_count != row_count:
                raise MergeError(f"{path}: GCT row count differs from the first input")

        input_sample_names = tuple(
            sample_name for header in headers for sample_name in header.sample_names
        )
        if sample_names_override is not None:
            if len(sample_names_override) != len(input_sample_names):
                raise MergeError(
                    "Assigned sample-name count does not match the GCT "
                    "sample-column count"
                )
            sample_names = tuple(sample_names_override)
        else:
            sample_names = input_sample_names
        ensure_unique(sample_names, "sample names")
        if expected_samples is not None and sample_names != tuple(expected_samples):
            raise MergeError("GCT sample names do not match the expected sample order")

        with atomic_text_writer(output_path) as output:
            output.write("#1.2\n")
            output.write(f"{row_count}\t{len(sample_names)}\n")
            output.write("\t".join(("Name", "Description", *sample_names)) + "\n")

            for row_index in range(1, row_count + 1):
                combined_values: list[str] = []
                reference_fields: tuple[str, str] | None = None
                for handle, path, header in zip(handles, input_paths, headers):
                    line = handle.readline()
                    if not line:
                        raise MergeError(f"{path}: ended before GCT row {row_index}")
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) != len(header.sample_names) + 2:
                        raise MergeError(
                            f"{path}: invalid field count at GCT row {row_index}"
                        )
                    identifiers = (fields[0], fields[1])
                    if reference_fields is None:
                        reference_fields = identifiers
                    elif identifiers != reference_fields:
                        raise MergeError(
                            f"{path}: feature differs at GCT row {row_index}; "
                            f"expected {reference_fields[0]}, found {identifiers[0]}"
                        )
                    combined_values.extend(fields[2:])

                assert reference_fields is not None
                output.write("\t".join((*reference_fields, *combined_values)) + "\n")

            for handle, path in zip(handles, input_paths):
                if handle.readline():
                    raise MergeError(
                        f"{path}: contains rows beyond its declared GCT row count"
                    )

    if sample_output_path is not None:
        write_lines_atomic(sample_output_path, sample_names)
    return sample_names


def read_individual_metrics(path: Path) -> tuple[str, tuple[str, ...], tuple[str, ...]]:
    with open_text(path, "rt") as handle:
        rows = [line.rstrip("\r\n").split("\t") for line in handle]
    if not rows or rows[0][0] != "Sample" or len(rows[0]) != 2 or not rows[0][1]:
        raise MergeError(f"{path}: expected a Sample row with one sample name")
    for row_number, row in enumerate(rows[1:], 2):
        if len(row) != 2 or not row[0]:
            raise MergeError(f"{path}: invalid metrics row {row_number}")
    metric_names = tuple(row[0] for row in rows[1:])
    ensure_unique(metric_names, f"metric names in {path}")
    metric_values = tuple(row[1] for row in rows[1:])
    return rows[0][1], metric_names, metric_values


def merge_individual_metrics(
    input_paths: Sequence[Path],
    expected_samples: Sequence[str] | None,
    output_path: Path,
    sample_names_override: Sequence[str] | None = None,
) -> None:
    parsed = [read_individual_metrics(path) for path in input_paths]
    input_sample_names = tuple(item[0] for item in parsed)
    if sample_names_override is not None:
        if len(sample_names_override) != len(input_sample_names):
            raise MergeError(
                "Assigned sample-name count does not match the metrics file count"
            )
        sample_names = tuple(sample_names_override)
    else:
        sample_names = input_sample_names
    ensure_unique(sample_names, "metrics sample names")
    if expected_samples is not None and sample_names != tuple(expected_samples):
        raise MergeError("Metrics sample names do not match the expected sample order")

    metric_names = parsed[0][1]
    for item, path in zip(parsed[1:], input_paths[1:]):
        if item[1] != metric_names:
            raise MergeError(
                f"{path}: metric names or order differ from the first input"
            )

    with atomic_text_writer(output_path) as output:
        output.write("\t".join(("sample_id", *metric_names)) + "\n")
        for sample_name, (_, _, metric_values) in zip(sample_names, parsed):
            output.write("\t".join((sample_name, *metric_values)) + "\n")


def merge_aggregated_metrics(
    input_paths: Sequence[Path],
    expected_samples: Sequence[str],
    output_path: Path,
) -> None:
    sample_names: list[str] = []
    expected_header: tuple[str, ...] | None = None
    with atomic_text_writer(output_path) as output:
        for path in input_paths:
            with open_text(path, "rt") as handle:
                header = tuple(handle.readline().rstrip("\r\n").split("\t"))
                if len(header) < 2 or header[0] != "sample_id":
                    raise MergeError(f"{path}: invalid aggregated metrics header")
                if expected_header is None:
                    expected_header = header
                    output.write("\t".join(header) + "\n")
                elif header != expected_header:
                    raise MergeError(
                        f"{path}: metrics header differs from the first input"
                    )

                for row_number, line in enumerate(handle, 2):
                    fields = line.rstrip("\r\n").split("\t")
                    if len(fields) != len(header) or not fields[0]:
                        raise MergeError(
                            f"{path}: invalid aggregated metrics row {row_number}"
                        )
                    sample_names.append(fields[0])
                    output.write(line if line.endswith("\n") else f"{line}\n")

        ensure_unique(sample_names, "metrics sample names")
        if tuple(sample_names) != tuple(expected_samples):
            raise MergeError(
                "Metrics sample names do not match the expected sample order"
            )


def read_individual_insert_sizes(path: Path) -> dict[int, int]:
    distribution: dict[int, int] = {}
    with open_text(path, "rt") as handle:
        for row_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2:
                raise MergeError(f"{path}: invalid insert-size row {row_number}")
            try:
                insert_size, count = (int(value) for value in fields)
            except ValueError as error:
                raise MergeError(
                    f"{path}: insert-size row {row_number} must contain integers"
                ) from error
            if insert_size < 0 or count < 0:
                raise MergeError(f"{path}: insert sizes and counts must be nonnegative")
            if insert_size in distribution:
                raise MergeError(f"{path}: duplicate insert-size bin {insert_size}")
            distribution[insert_size] = count
    if not distribution:
        raise MergeError(f"{path}: insert-size file is empty")
    return distribution


def merge_individual_insert_sizes(
    input_paths: Sequence[Path],
    expected_samples: Sequence[str],
    output_path: Path,
) -> None:
    if len(input_paths) != len(expected_samples):
        raise MergeError(
            "Insert-size file count does not match the expected sample count"
        )
    distributions = [read_individual_insert_sizes(path) for path in input_paths]
    bins = sorted(
        {insert_size for distribution in distributions for insert_size in distribution}
    )
    with atomic_text_writer(output_path) as output:
        output.write("\t".join(("insert_size", *expected_samples)) + "\n")
        for insert_size in bins:
            counts = (
                str(distribution.get(insert_size, 0))
                for distribution in distributions
            )
            output.write("\t".join((str(insert_size), *counts)) + "\n")


def read_aggregated_insert_row(
    handle: IO[str],
    path: Path,
    sample_count: int,
    row_number: int,
    previous_bin: int | None,
) -> tuple[int, tuple[str, ...]] | None:
    line = handle.readline()
    if not line:
        return None
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != sample_count + 1:
        raise MergeError(f"{path}: invalid insert-size row {row_number}")
    try:
        insert_size = int(fields[0])
        counts = tuple(int(value) for value in fields[1:])
    except ValueError as error:
        raise MergeError(
            f"{path}: insert-size row {row_number} must contain integers"
        ) from error
    if insert_size < 0 or any(count < 0 for count in counts):
        raise MergeError(f"{path}: insert sizes and counts must be nonnegative")
    if previous_bin is not None and insert_size <= previous_bin:
        raise MergeError(f"{path}: insert-size bins must be strictly increasing")
    return insert_size, tuple(str(count) for count in counts)


def merge_aggregated_insert_sizes(
    input_paths: Sequence[Path],
    expected_samples: Sequence[str],
    output_path: Path,
) -> None:
    with ExitStack() as stack:
        handles = [stack.enter_context(open_text(path, "rt")) for path in input_paths]
        sample_groups: list[tuple[str, ...]] = []
        for handle, path in zip(handles, input_paths):
            fields = handle.readline().rstrip("\r\n").split("\t")
            if len(fields) < 2 or fields[0] not in ("", "insert_size"):
                raise MergeError(f"{path}: invalid aggregated insert-size header")
            sample_groups.append(tuple(fields[1:]))

        sample_names = tuple(sample for group in sample_groups for sample in group)
        ensure_unique(sample_names, "insert-size sample names")
        if sample_names != tuple(expected_samples):
            raise MergeError(
                "Insert-size sample names do not match the expected sample order"
            )

        previous_bins: list[int | None] = [None] * len(handles)
        row_numbers = [2] * len(handles)
        current_rows = [
            read_aggregated_insert_row(handle, path, len(samples), 2, None)
            for handle, path, samples in zip(handles, input_paths, sample_groups)
        ]

        with atomic_text_writer(output_path) as output:
            output.write("\t".join(("insert_size", *sample_names)) + "\n")
            while any(row is not None for row in current_rows):
                insert_size = min(row[0] for row in current_rows if row is not None)
                combined_counts: list[str] = []
                for index, sample_group in enumerate(sample_groups):
                    row = current_rows[index]
                    if row is not None and row[0] == insert_size:
                        combined_counts.extend(row[1])
                        previous_bins[index] = row[0]
                        row_numbers[index] += 1
                        current_rows[index] = read_aggregated_insert_row(
                            handles[index],
                            input_paths[index],
                            len(sample_group),
                            row_numbers[index],
                            previous_bins[index],
                        )
                    else:
                        combined_counts.extend("0" for _ in sample_group)
                output.write("\t".join((str(insert_size), *combined_counts)) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_manifest_parser = subparsers.add_parser(
        "validate-manifest", help="Validate one combined sample manifest"
    )
    validate_manifest_parser.add_argument("--input", type=Path, required=True)
    validate_manifest_parser.add_argument("--batch-size", type=int, required=True)
    validate_manifest_parser.add_argument("--prefix-file", type=Path)

    prepare_batch_parser = subparsers.add_parser(
        "prepare-batch", help="Write transfer and local-file lists for one batch"
    )
    prepare_batch_parser.add_argument("--input", type=Path, required=True)
    prepare_batch_parser.add_argument("--batch-index", type=int, required=True)
    prepare_batch_parser.add_argument("--batch-size", type=int, required=True)
    prepare_batch_parser.add_argument("--staging-directory", type=Path, required=True)
    prepare_batch_parser.add_argument(
        "--skip-exons", action="store_true",
        help="Omit exon files from transfer and local-file lists; keep QC metrics",
    )

    gct_parser = subparsers.add_parser("gct", help="Merge GCT sample columns")
    gct_parser.add_argument("--input-list", type=Path, required=True)
    gct_parser.add_argument("--output", type=Path, required=True)
    gct_parser.add_argument("--sample-output", type=Path)
    gct_parser.add_argument("--expected-samples", type=Path)
    gct_parser.add_argument("--sample-names", type=Path)

    individual_metrics_parser = subparsers.add_parser(
        "metrics-individual", help="Merge individual metrics files"
    )
    individual_metrics_parser.add_argument("--input-list", type=Path, required=True)
    individual_metrics_names = individual_metrics_parser.add_mutually_exclusive_group(
        required=True
    )
    individual_metrics_names.add_argument("--expected-samples", type=Path)
    individual_metrics_names.add_argument("--sample-names", type=Path)
    individual_metrics_parser.add_argument("--output", type=Path, required=True)

    aggregated_metrics_parser = subparsers.add_parser(
        "metrics-aggregated", help="Merge batch-level metrics tables"
    )
    aggregated_metrics_parser.add_argument("--input-list", type=Path, required=True)
    aggregated_metrics_parser.add_argument(
        "--expected-samples", type=Path, required=True
    )
    aggregated_metrics_parser.add_argument("--output", type=Path, required=True)

    individual_insert_parser = subparsers.add_parser(
        "insert-sizes-individual", help="Merge individual insert-size distributions"
    )
    individual_insert_parser.add_argument("--input-list", type=Path, required=True)
    individual_insert_parser.add_argument(
        "--expected-samples", type=Path, required=True
    )
    individual_insert_parser.add_argument("--output", type=Path, required=True)

    aggregated_insert_parser = subparsers.add_parser(
        "insert-sizes-aggregated", help="Merge batch-level insert-size distributions"
    )
    aggregated_insert_parser.add_argument("--input-list", type=Path, required=True)
    aggregated_insert_parser.add_argument(
        "--expected-samples", type=Path, required=True
    )
    aggregated_insert_parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "validate-manifest":
            validate_sample_manifest(args.input, args.batch_size, args.prefix_file)
        elif args.command == "prepare-batch":
            prepare_manifest_batch(
                args.input,
                args.batch_index,
                args.batch_size,
                args.staging_directory,
                merge_exons=not args.skip_exons,
            )
        elif args.command == "gct":
            if args.expected_samples and args.sample_names:
                raise MergeError(
                    "Use either --expected-samples or --sample-names, not both"
                )
            merge_gcts(
                read_path_list(args.input_list),
                args.output,
                args.sample_output,
                read_sample_ids(args.expected_samples)
                if args.expected_samples
                else None,
                read_sample_ids(args.sample_names) if args.sample_names else None,
            )
        elif args.command == "metrics-individual":
            merge_individual_metrics(
                read_path_list(args.input_list),
                read_sample_ids(args.expected_samples)
                if args.expected_samples
                else None,
                args.output,
                read_sample_ids(args.sample_names) if args.sample_names else None,
            )
        elif args.command == "metrics-aggregated":
            merge_aggregated_metrics(
                read_path_list(args.input_list),
                read_sample_ids(args.expected_samples),
                args.output,
            )
        elif args.command == "insert-sizes-individual":
            merge_individual_insert_sizes(
                read_path_list(args.input_list),
                read_sample_ids(args.expected_samples),
                args.output,
            )
        elif args.command == "insert-sizes-aggregated":
            merge_aggregated_insert_sizes(
                read_path_list(args.input_list),
                read_sample_ids(args.expected_samples),
                args.output,
            )
    except (MergeError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())