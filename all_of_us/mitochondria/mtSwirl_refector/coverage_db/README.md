# coverage_db (v2)

This folder contains the v2 coverage database builder used by the `mt_coverage_merge` rewrite.

## What it does

`build_coverage_db.py` reads per-sample base-level coverage TSVs, stores them in a single HDF5 file, and computes exact per-position summary metrics:

- mean coverage
- median coverage
- fraction of samples with coverage > 100
- fraction of samples with coverage > 1000

These correspond to the row fields in the v1 coverage MT schema.

## Files

- `build_coverage_db.py`: main implementation
- `smoke_test_build_coverage_db.py`: small self-contained correctness smoke test

## Dependencies

- `numpy`
- `h5py`

The repo’s current dev environment may not have `h5py` installed.
The intended v2 runtime image / WDL task should install/include it.

## Docker

A minimal Dockerfile is provided at `generate_mtdna_call_mt/coverage_db/Dockerfile`.

Build from the **repo root** (so `generate_mtdna_call_mt/` is in the build context):

```bash
docker build \
	-f generate_mtdna_call_mt/coverage_db/Dockerfile \
	-t mt-coverage-db:dev \
	.
```

Quick sanity check (runs the smoke test inside the container):

```bash
docker run --rm mt-coverage-db:dev \
	python -m generate_mtdna_call_mt.coverage_db.smoke_test_build_coverage_db
```

For WDL integration, you’ll typically override the ENTRYPOINT and run:

```bash
python -m generate_mtdna_call_mt.coverage_db.build_coverage_db \
	--input-tsv processed_data.tsv \
	--output-h5 coverage.h5 \
	--output-summary coverage_summary.tsv
```

## Notes on exact median

The output schema expects `median` to be an `int32`.
For even N, the builder computes the average of the two middle values and casts to int.
If parity testing against Hail indicates a different rounding convention, adjust that conversion.
