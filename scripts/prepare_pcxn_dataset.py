#!/usr/bin/env python3
"""
Scale-safe ETL pipeline that turns one or more raw PCxN pairwise pathway
co-expression edge lists - one per tissue - into the artifacts needed for the
web deployment architecture. Each input file is processed fully independently
into its own subtree:

  <output-dir>/<tissue_id>/pcxn_master.parquet
                           - single canonical archive (ZSTD), bidirectional,
                             for bulk download / DuckDB / containerized querying.
  <output-dir>/<tissue_id>/pathways/<id>.json|parquet
                           - one small file per pathway holding its full edge
                             neighborhood, sorted by significance, for R2 /
                             Cloudflare Worker serving at a deterministic key
                             (R2 key convention: <tissue_id>/pathways/<id>.json).
  <output-dir>/<tissue_id>/pathway_index.json
                           - lightweight per-tissue directory (id, name, edge
                             counts, min/max correlation) for frontend
                             autocomplete and routing.
  <output-dir>/tissues_index.json
                           - top-level manifest listing every tissue (id,
                             source file, pathway/edge counts, and paths to
                             its own pcxn_master.parquet / pathway_index.json)
                             for a frontend tissue picker.

Each input filename is expected to look like ``<tissue>_pathway_estimates.tsv``
(or ``.tsv.gz``); the tissue id is derived by stripping that known suffix (see
``--tissue-suffix``) and slugifying what remains.

Design notes (why it's built this way)
---------------------------------------
* Input is streamed in bounded-size batches with PyArrow's CSV reader
  (``pyarrow.csv.open_csv`` over ``pyarrow.input_stream(..., compression="detect")``),
  so a 100+ GB(.gz) file is never materialized in memory. Only one batch
  (``--batch-size-mb``) is resident at a time.
* Symmetry (A->B implies B->A) is enforced by emitting both directed rows for
  every input row. Rows are then routed by hash(pathway_a) into one of
  ``--num-buckets`` on-disk partitions. Because both directions of an edge are
  generated *before* bucketing, and bucketing keys off each row's own
  ``pathway_a``, all edges for a given pathway always land in exactly one
  bucket - regardless of how large the full dataset is. This is what makes
  the pathway-grouping step (building per-pathway slices) possible without
  ever sorting or grouping the whole dataset in memory.
* Each bucket is sized (via ``--num-buckets``) to comfortably fit in RAM, so
  phase 2 (dedup, sort, slice + index generation, append to master) can use
  ordinary in-memory Polars operations per bucket while the *overall* pipeline
  stays bounded by (batch size + one bucket) rather than by total input size.
* Pathway identifiers stay UTF-8 strings end to end. Parquet's built-in
  dictionary encoding already gives categorical-equivalent compression on
  disk for repeated string columns, without the cross-batch categorical
  consistency issues (global string cache) that an explicit ``Categorical``
  dtype would introduce across independently-processed streaming batches.
* ``correlation`` is downcast to float32 as soon as a batch is read. ``p_value``
  and ``fdr`` are, by default, transformed to -log10(x) float32 (see
  ``--pvalue-mode``) before downcasting, since real PCxN p-values/FDRs routinely
  underflow float32 (~1.18e-38 floor) and would otherwise collapse to 0.0.

Requirements: polars>=0.20, pyarrow>=14 (see requirements.txt).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import re
import shutil
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import polars as pl
import pyarrow as pa
import pyarrow.csv as pa_csv
import pyarrow.parquet as pq

logger = logging.getLogger("prepare_pcxn_dataset")

# normalized (lowercase, alnum-only) header aliases we auto-detect against
COLUMN_ALIASES = {
    "pathway_a": ["pathwaya", "pathway1", "pathwayasource", "sourcepathway", "nodea", "source"],
    "pathway_b": ["pathwayb", "pathway2", "pathwaybtarget", "targetpathway", "nodeb", "target"],
    "correlation": ["correlation", "pathcor", "cor", "r", "weight", "pcxn"],
    "p_value": ["pvalue", "p"],
    "fdr": ["fdr", "padj", "qvalue", "adjp", "adjustedp"],
}


def _normalize_header(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def slugify(name: str) -> str:
    """Deterministic, filesystem/URL-safe id derived from an arbitrary raw name."""
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    slug = re.sub(r"_+", "_", slug).strip("_")
    return slug or "pathway"


def derive_tissue_id(path: Path, suffix: str) -> str:
    """Derive a tissue id from an input filename like <tissue>_pathway_estimates.tsv[.gz]."""
    stem = path.name
    for ext in (".tsv.gz", ".tsv"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break
    if suffix and stem.lower().endswith(suffix.lower()):
        stem = stem[: -len(suffix)]
    return slugify(stem)


def _short_hash(text: str) -> str:
    return hashlib.md5(text.encode("utf-8")).hexdigest()[:8]


def _fmt_float(v: Optional[float]) -> Optional[float]:
    """Render a float32 value with ~6 significant digits, no float64 noise tail."""
    if v is None:
        return None
    try:
        return float(f"{v:.6g}")
    except (ValueError, TypeError):
        return None


@dataclass
class ColumnMap:
    pathway_a: str
    pathway_b: str
    correlation: str
    p_value: str
    fdr: str
    extra: List[str] = field(default_factory=list)


def resolve_columns(header: List[str], args: argparse.Namespace) -> ColumnMap:
    normalized = {_normalize_header(h): h for h in header}
    overrides = {
        "pathway_a": args.col_pathway_a,
        "pathway_b": args.col_pathway_b,
        "correlation": args.col_correlation,
        "p_value": args.col_pvalue,
        "fdr": args.col_fdr,
    }
    # Raises ValueError (not SystemExit) on failure: this runs per-tissue inside
    # run_pipeline_for_tissue, and a single malformed/mismatched tissue file
    # should fail that tissue without aborting the rest of the batch.
    resolved: Dict[str, str] = {}
    for canon, override in overrides.items():
        if override:
            if override not in header:
                raise ValueError(f"--col-{canon.replace('_', '-')}={override!r} not found in input header: {header}")
            resolved[canon] = override
            continue
        candidates = [canon] + COLUMN_ALIASES[canon]
        found = None
        for cand in candidates:
            key = _normalize_header(cand)
            if key in normalized:
                found = normalized[key]
                break
        if found is None:
            raise ValueError(
                f"Could not auto-detect the '{canon}' column in header {header}. "
                f"Pass --col-{canon.replace('_', '-')} explicitly."
            )
        resolved[canon] = found

    extra = [c.strip() for c in (args.extra_cols or "").split(",") if c.strip()]
    for c in extra:
        if c not in header:
            raise ValueError(f"--extra-cols column {c!r} not found in input header: {header}")

    logger.info("Resolved columns: %s (extra=%s)", resolved, extra)
    return ColumnMap(**resolved, extra=extra)


def peek_header(input_path: Path) -> List[str]:
    opener = gzip.open if input_path.suffix == ".gz" else open
    with opener(input_path, "rt", encoding="utf-8", errors="replace") as fh:
        first_line = fh.readline().rstrip("\n\r")
    return first_line.split("\t")


# float32 tops out at ~1.18e-38 for normal values; real PCxN p-values/FDRs
# routinely go down to 1e-64 or smaller, so a strict float32 cast collapses
# ~18% of a typical PCxN edge list's p-values to exactly 0.0, destroying both
# the magnitude and the relative ranking of the most-significant edges. To
# avoid that silently, p_value/fdr are (by default) stored as -log10(x)
# float32 - a standard bioinformatics convention that keeps values like 1e-300
# comfortably representable (as 300.0) with full ordering preserved. This is
# a schema-visible choice: the resulting column names (neg_log10_p_value /
# neg_log10_fdr vs. plain p_value / fdr) and `pvalue_mode` are recorded in
# pathway_index.json so downstream consumers know which they got.
NEG_LOG10_FLOOR = 1e-300  # clip before log10 so exact-zero inputs don't yield inf


@dataclass
class MetricFields:
    p_value_col: str
    fdr_col: str
    transform: str  # "neg_log10" | "raw"

    @classmethod
    def for_mode(cls, mode: str) -> "MetricFields":
        if mode == "neg_log10":
            return cls(p_value_col="neg_log10_p_value", fdr_col="neg_log10_fdr", transform="neg_log10")
        if mode == "raw":
            return cls(p_value_col="p_value", fdr_col="fdr", transform="raw")
        raise ValueError(f"Unknown --pvalue-mode value: {mode}")

    @property
    def more_significant_is_larger(self) -> bool:
        return self.transform == "neg_log10"


def build_sort_columns(sort_by: str, metrics: MetricFields) -> Tuple[List[pl.Expr], List[str], List[bool]]:
    """Return (helper exprs, sort keys, descending flags) for a bucket's sort."""
    sig_desc = metrics.more_significant_is_larger  # True => larger value = more significant
    if sort_by == "fdr":
        keys = [metrics.fdr_col, "_abs_correlation"]
        desc = [sig_desc, True]
    elif sort_by == "p_value":
        keys = [metrics.p_value_col, "_abs_correlation"]
        desc = [sig_desc, True]
    elif sort_by == "abs_correlation":
        keys = ["_abs_correlation", metrics.fdr_col]
        desc = [True, sig_desc]
    else:
        raise ValueError(f"Unknown --sort-by value: {sort_by}")
    helper_exprs = [pl.col("correlation").abs().alias("_abs_correlation")]
    return helper_exprs, keys, desc


def _metric_exprs(metrics: MetricFields) -> List[pl.Expr]:
    """Build the p_value/fdr -> final metric column(s) expressions, applied
    right after ingestion while the raw values are still float64."""
    if metrics.transform == "raw":
        return [
            pl.col("_raw_p_value").cast(pl.Float32).alias("p_value"),
            pl.col("_raw_fdr").cast(pl.Float32).alias("fdr"),
        ]
    return [
        (-(pl.col("_raw_p_value").clip(NEG_LOG10_FLOOR, None).log10()))
        .cast(pl.Float32)
        .alias("neg_log10_p_value"),
        (-(pl.col("_raw_fdr").clip(NEG_LOG10_FLOOR, None).log10())).cast(pl.Float32).alias("neg_log10_fdr"),
    ]


def symmetrize_and_cast(batch: pa.RecordBatch, colmap: ColumnMap, metrics: MetricFields) -> pl.DataFrame:
    df = pl.from_arrow(batch)
    rename_map = {
        colmap.pathway_a: "pathway_a",
        colmap.pathway_b: "pathway_b",
        colmap.correlation: "correlation",
        colmap.p_value: "_raw_p_value",
        colmap.fdr: "_raw_fdr",
    }
    df = df.rename(rename_map)
    keep = ["pathway_a", "pathway_b", "correlation", "_raw_p_value", "_raw_fdr"] + colmap.extra
    df = df.select(keep).with_columns(
        [
            pl.col("pathway_a").cast(pl.Utf8),
            pl.col("pathway_b").cast(pl.Utf8),
            pl.col("correlation").cast(pl.Float32),
        ]
        + _metric_exprs(metrics)
    )
    metric_cols = [metrics.p_value_col, metrics.fdr_col]
    df = df.select(["pathway_a", "pathway_b", "correlation"] + metric_cols + colmap.extra)

    swapped = df.select(
        [
            pl.col("pathway_b").alias("pathway_a"),
            pl.col("pathway_a").alias("pathway_b"),
            pl.col("correlation"),
            *[pl.col(c) for c in metric_cols],
            *[pl.col(c) for c in colmap.extra],
        ]
    )
    return pl.concat([df, swapped], how="vertical")


class BucketWriters:
    """Lazily-opened ParquetWriters, one per hash bucket, for the shuffle phase."""

    def __init__(self, tmp_dir: Path, num_buckets: int):
        self.tmp_dir = tmp_dir
        self.num_buckets = num_buckets
        self._writers: Dict[int, pq.ParquetWriter] = {}
        self._schema: Optional[pa.Schema] = None

    def path_for(self, idx: int) -> Path:
        return self.tmp_dir / f"bucket_{idx:05d}.parquet"

    def write(self, idx: int, table: pa.Table) -> None:
        if self._schema is None:
            self._schema = table.schema
        writer = self._writers.get(idx)
        if writer is None:
            writer = pq.ParquetWriter(str(self.path_for(idx)), self._schema, compression="zstd")
            self._writers[idx] = writer
        writer.write_table(table)

    def close(self) -> None:
        for w in self._writers.values():
            w.close()


def ingest_and_shuffle(
    input_path: Path,
    colmap: ColumnMap,
    metrics: MetricFields,
    tmp_dir: Path,
    num_buckets: int,
    batch_size_mb: int,
) -> int:
    """Phase 1: stream the raw TSV(.gz), symmetrize + downcast, hash-partition
    into on-disk buckets by pathway_a. Returns total raw input row count."""
    tmp_dir.mkdir(parents=True, exist_ok=True)
    read_options = pa_csv.ReadOptions(block_size=batch_size_mb * 1024 * 1024)
    parse_options = pa_csv.ParseOptions(delimiter="\t")
    convert_options = pa_csv.ConvertOptions(
        column_types={
            colmap.pathway_a: pa.string(),
            colmap.pathway_b: pa.string(),
            colmap.p_value: pa.float64(),
            colmap.fdr: pa.float64(),
        }
    )

    writers = BucketWriters(tmp_dir, num_buckets)
    total_raw_rows = 0
    batch_num = 0
    t0 = time.time()

    stream = pa.input_stream(str(input_path), compression="detect")
    reader = pa_csv.open_csv(
        stream, read_options=read_options, parse_options=parse_options, convert_options=convert_options
    )
    try:
        for batch in reader:
            batch_num += 1
            total_raw_rows += batch.num_rows
            combined = symmetrize_and_cast(batch, colmap, metrics)
            combined = combined.with_columns((pl.col("pathway_a").hash(seed=42) % num_buckets).alias("_bucket"))
            groups = combined.partition_by(["_bucket"], as_dict=True, maintain_order=False)
            for key, gdf in groups.items():
                idx = int(key[0] if isinstance(key, tuple) else key)
                writers.write(idx, gdf.drop("_bucket").to_arrow())
            if batch_num % 20 == 0:
                logger.info(
                    "  ingest: batch %d, %d raw rows so far (%.1fs elapsed)",
                    batch_num,
                    total_raw_rows,
                    time.time() - t0,
                )
    finally:
        writers.close()

    logger.info(
        "Phase 1 done: %d raw input rows -> %d batches -> %d bucket files in %.1fs",
        total_raw_rows,
        batch_num,
        sum(1 for i in range(num_buckets) if writers.path_for(i).exists()),
        time.time() - t0,
    )
    return total_raw_rows


def write_json_slice(path: Path, gdf: pl.DataFrame, metric_cols: List[str], extra_cols: List[str]) -> None:
    cols = ["pathway_b", "correlation"] + metric_cols + extra_cols
    float_cols = {"correlation"} | set(metric_cols)
    columns = {c: gdf[c].to_list() for c in cols}
    edges = []
    n = gdf.height
    for i in range(n):
        edge = {}
        for c in cols:
            v = columns[c][i]
            edge[c] = _fmt_float(v) if c in float_cols else v
        edges.append(edge)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"edges": edges}, fh, separators=(",", ":"))


def process_buckets(
    tmp_dir: Path,
    num_buckets: int,
    master_writer: pq.ParquetWriter,
    pathways_dir: Path,
    slice_format: str,
    sort_by: str,
    metrics: MetricFields,
    max_edges_per_slice: int,
    extra_cols: List[str],
) -> Tuple[Dict[str, dict], int, int]:
    """Phase 2: for every bucket, dedup direction pairs, sort, append to the
    master writer, and emit one slice file + index entry per pathway."""
    index: Dict[str, dict] = {}
    slug_registry: Dict[str, str] = {}  # slug -> original pathway name
    total_master_rows = 0
    total_dupes_removed = 0
    metric_cols = [metrics.p_value_col, metrics.fdr_col]
    helper_exprs, sort_keys, sort_desc = build_sort_columns(sort_by, metrics)

    for idx in range(num_buckets):
        bpath = tmp_dir / f"bucket_{idx:05d}.parquet"
        if not bpath.exists():
            continue
        table = pq.read_table(bpath)
        df = pl.from_arrow(table)
        before = df.height
        df = df.unique(subset=["pathway_a", "pathway_b"], keep="first")
        total_dupes_removed += before - df.height

        df = df.with_columns(helper_exprs)
        df = df.sort(["pathway_a"] + sort_keys, descending=[False] + sort_desc)
        df = df.drop("_abs_correlation")

        master_writer.write_table(df.to_arrow())
        total_master_rows += df.height

        groups = df.partition_by(["pathway_a"], as_dict=True, maintain_order=True)
        for key, gdf in groups.items():
            pname = key[0] if isinstance(key, tuple) else key
            gdf = gdf.drop("pathway_a")

            full_edge_count = gdf.height
            correlation = gdf["correlation"].to_list()
            min_corr = min(correlation) if correlation else None
            max_corr = max(correlation) if correlation else None
            metric_stats = {}
            for mc in metric_cols:
                vals = gdf[mc].to_list()
                metric_stats[f"min_{mc}"] = _fmt_float(min(vals)) if vals else None
                metric_stats[f"max_{mc}"] = _fmt_float(max(vals)) if vals else None

            sliced = gdf.head(max_edges_per_slice) if max_edges_per_slice > 0 else gdf

            base_slug = slugify(pname)
            slug = base_slug
            if slug in slug_registry and slug_registry[slug] != pname:
                slug = f"{base_slug}~{_short_hash(pname)}"
            slug_registry[slug] = pname

            out_path = pathways_dir / f"{slug}.{ 'json' if slice_format == 'json' else 'parquet' }"
            if slice_format == "json":
                write_json_slice(out_path, sliced, metric_cols, extra_cols)
            else:
                pq.write_table(sliced.to_arrow(), out_path, compression="zstd")

            index[slug] = {
                "pathway_id": slug,
                "pathway_name": pname,
                "edge_count": full_edge_count,
                "edge_count_in_slice": sliced.height,
                "min_correlation": _fmt_float(min_corr),
                "max_correlation": _fmt_float(max_corr),
                **metric_stats,
            }
        logger.info("  bucket %d/%d: %d rows, %d pathways so far", idx + 1, num_buckets, df.height, len(index))

    return index, total_master_rows, total_dupes_removed


def run_validation(
    master_path: Path,
    pathways_dir: Path,
    slice_format: str,
    index: Dict[str, dict],
    input_path: Path,
    total_raw_rows: int,
    total_master_rows: int,
    total_dupes_removed: int,
) -> bool:
    ok = True
    logger.info("=" * 70)
    logger.info("VALIDATION REPORT")
    logger.info("=" * 70)
    logger.info("Raw input rows read:              %d", total_raw_rows)
    logger.info("Master (bidirectional, deduped):  %d", total_master_rows)
    logger.info("Duplicate directed pairs removed: %d", total_dupes_removed)
    logger.info("Unique pathways indexed:          %d", len(index))

    input_size = input_path.stat().st_size
    master_size = master_path.stat().st_size
    slice_files = list(pathways_dir.glob(f"*.{ 'json' if slice_format == 'json' else 'parquet' }"))
    slice_total_size = sum(p.stat().st_size for p in slice_files)
    logger.info("Input file size:    %.2f MB", input_size / 1e6)
    logger.info("Master parquet size: %.2f MB (%.1fx vs raw input)", master_size / 1e6, input_size / max(master_size, 1))
    logger.info("Slice files:        %d files, %.2f MB total, %.2f KB avg", len(slice_files), slice_total_size / 1e6, (slice_total_size / max(len(slice_files), 1)) / 1e3)

    if not index:
        logger.error("No pathways were indexed - validation FAILED")
        return False

    sample_id = max(index, key=lambda k: index[k]["edge_count"])
    entry = index[sample_id]
    logger.info("Sampling pathway %r (name=%r, edge_count=%d) for round-trip check", sample_id, entry["pathway_name"], entry["edge_count"])

    slice_path = pathways_dir / f"{sample_id}.{ 'json' if slice_format == 'json' else 'parquet' }"
    if slice_format == "json":
        with open(slice_path, encoding="utf-8") as fh:
            slice_data = json.load(fh)
        slice_targets = {e["pathway_b"] for e in slice_data["edges"]}
    else:
        slice_targets = set(pl.read_parquet(slice_path)["pathway_b"].to_list())

    master_lf = pl.scan_parquet(master_path).filter(pl.col("pathway_a") == entry["pathway_name"])
    master_df = master_lf.collect()
    master_targets = set(master_df["pathway_b"].to_list())

    expected_in_slice = min(entry["edge_count"], entry["edge_count_in_slice"])
    if master_df.height != entry["edge_count"]:
        logger.error(
            "MISMATCH: master has %d rows for pathway_a=%r but index says edge_count=%d",
            master_df.height, entry["pathway_name"], entry["edge_count"],
        )
        ok = False
    if not slice_targets.issubset(master_targets):
        logger.error("MISMATCH: slice contains pathway_b values not present in master rows for this pathway")
        ok = False
    if len(slice_targets) != expected_in_slice:
        logger.error(
            "MISMATCH: slice file has %d edges, expected %d (edge_count_in_slice)",
            len(slice_targets), expected_in_slice,
        )
        ok = False

    logger.info("Round-trip check: %s", "PASS" if ok else "FAIL")
    logger.info("=" * 70)
    return ok


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="ETL: raw PCxN pairwise TSV(.gz) -> master parquet + per-pathway slices + index.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--input-tsv",
        nargs="+",
        type=Path,
        default=None,
        help="One or more input .tsv/.tsv.gz files, each treated as a separate tissue. "
        "Combine with --input-dir if useful.",
    )
    p.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Directory to scan for tissue TSVs using --glob-pattern (each match is one tissue)",
    )
    p.add_argument(
        "--glob-pattern",
        default="*_pathway_estimates.tsv*",
        help="Glob (relative to --input-dir) used to find tissue files",
    )
    p.add_argument(
        "--tissue-suffix",
        default="_pathway_estimates",
        help="Suffix stripped from each filename stem to derive the tissue id "
        "(e.g. 'combined_pathway_estimates.tsv' -> tissue id 'combined')",
    )
    p.add_argument("--output-dir", required=True, type=Path, help="Root directory to write per-tissue outputs + tissues_index.json into")
    p.add_argument("--slice-format", choices=["json", "parquet"], default="json")

    p.add_argument("--col-pathway-a", default=None, help="Override auto-detected pathway_a column name")
    p.add_argument("--col-pathway-b", default=None, help="Override auto-detected pathway_b column name")
    p.add_argument("--col-correlation", default=None, help="Override auto-detected correlation column name")
    p.add_argument("--col-pvalue", default=None, help="Override auto-detected p_value column name")
    p.add_argument("--col-fdr", default=None, help="Override auto-detected fdr column name")
    p.add_argument("--extra-cols", default="", help="Comma-separated extra columns to carry through untouched (e.g. layer,weight)")

    p.add_argument("--num-buckets", type=int, default=256, help="Hash-partition fan-out for the shuffle phase; raise for larger inputs")
    p.add_argument("--batch-size-mb", type=int, default=64, help="Approx. CSV read block size per streaming batch")
    p.add_argument(
        "--pvalue-mode",
        choices=["neg_log10", "raw"],
        default="neg_log10",
        help=(
            "'neg_log10' stores -log10(p_value)/-log10(fdr) as float32 (columns "
            "neg_log10_p_value/neg_log10_fdr) so very small p-values stay distinguishable "
            "instead of underflowing float32 to 0.0. 'raw' stores p_value/fdr as plain float32."
        ),
    )
    p.add_argument("--sort-by", choices=["fdr", "p_value", "abs_correlation"], default="fdr", help="Primary significance ranking within each pathway slice")
    p.add_argument("--max-edges-per-slice", type=int, default=0, help="Cap edges written per slice file (0 = unlimited; master parquet is never truncated)")

    p.add_argument("--pretty-index", action="store_true", help="Pretty-print pathway_index.json (default: compact)")
    p.add_argument("--keep-temp", action="store_true", help="Keep the intermediate bucket files for debugging")
    p.add_argument("--log-level", default="INFO")
    return p.parse_args(argv)


def resolve_input_paths(args: argparse.Namespace) -> List[Path]:
    """Combine --input-tsv and --input-dir into an ordered, de-duplicated list."""
    paths: List[Path] = []
    if args.input_tsv:
        paths.extend(args.input_tsv)
    if args.input_dir:
        matches = sorted(args.input_dir.glob(args.glob_pattern))
        if not matches:
            raise SystemExit(f"No files under {args.input_dir} matched glob {args.glob_pattern!r}")
        paths.extend(matches)
    if not paths:
        raise SystemExit("No input files given. Pass --input-tsv (one or more) and/or --input-dir.")

    seen = set()
    unique_paths = []
    for p in paths:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            unique_paths.append(p)
    return unique_paths


def build_tissue_inputs(paths: List[Path], suffix: str) -> List[Tuple[str, Path]]:
    """Validate each input exists and derive its tissue id, failing fast on collisions."""
    result: List[Tuple[str, Path]] = []
    seen: Dict[str, Path] = {}
    for p in paths:
        if not p.exists():
            raise SystemExit(f"Input file not found: {p}")
        tissue_id = derive_tissue_id(p, suffix)
        if tissue_id in seen:
            raise SystemExit(
                f"Tissue id collision: both {seen[tissue_id]} and {p} derive tissue_id={tissue_id!r}. "
                "Rename one of the input files or adjust --tissue-suffix."
            )
        seen[tissue_id] = p
        result.append((tissue_id, p))
    return result


def run_pipeline_for_tissue(tissue_id: str, input_path: Path, output_dir: Path, args: argparse.Namespace) -> dict:
    """Run the full 3-phase pipeline for a single tissue's input file, writing
    into its own output_dir. Returns a summary dict for the tissue manifest."""
    pathways_dir = output_dir / "pathways"
    tmp_dir = output_dir / "_tmp_buckets"
    output_dir.mkdir(parents=True, exist_ok=True)
    pathways_dir.mkdir(parents=True, exist_ok=True)
    master_path = output_dir / "pcxn_master.parquet"
    index_path = output_dir / "pathway_index.json"

    header = peek_header(input_path)
    colmap = resolve_columns(header, args)
    metrics = MetricFields.for_mode(args.pvalue_mode)

    t_start = time.time()
    logger.info("Phase 1/3: streaming ingest + symmetrize + hash-partition into %d buckets", args.num_buckets)
    total_raw_rows = ingest_and_shuffle(input_path, colmap, metrics, tmp_dir, args.num_buckets, args.batch_size_mb)

    logger.info("Phase 2/3: per-bucket dedup + sort + master append + slice/index generation")
    bucket_files = [tmp_dir / f"bucket_{i:05d}.parquet" for i in range(args.num_buckets)]
    bucket_files = [p for p in bucket_files if p.exists()]
    if not bucket_files:
        raise RuntimeError(f"No data was read from {input_path} (0 bucket files produced).")
    # Derive the master schema from an actual bucket file rather than hardcoding
    # one, since extra passthrough columns can be of any Arrow type.
    master_schema = pq.read_schema(str(bucket_files[0]))
    master_writer = pq.ParquetWriter(str(master_path), master_schema, compression="zstd")
    try:
        index, total_master_rows, total_dupes_removed = process_buckets(
            tmp_dir,
            args.num_buckets,
            master_writer,
            pathways_dir,
            args.slice_format,
            args.sort_by,
            metrics,
            args.max_edges_per_slice,
            colmap.extra,
        )
    finally:
        master_writer.close()

    if not args.keep_temp:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    logger.info("Phase 3/3: writing pathway_index.json")
    index_payload = {
        "tissue_id": tissue_id,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source_file": str(input_path),
        "slice_format": args.slice_format,
        "sort_by": args.sort_by,
        "pvalue_mode": metrics.transform,
        "columns": {
            "pathway_b": "pathway_b",
            "correlation": "correlation",
            "p_value": metrics.p_value_col,
            "fdr": metrics.fdr_col,
        },
        "total_pathways": len(index),
        "total_edges": total_master_rows,
        "pathways": [index[k] for k in sorted(index.keys())],
    }
    with open(index_path, "w", encoding="utf-8") as fh:
        if args.pretty_index:
            json.dump(index_payload, fh, indent=2)
        else:
            json.dump(index_payload, fh, separators=(",", ":"))

    logger.info("Tissue %r finished in %.1fs", tissue_id, time.time() - t_start)

    ok = run_validation(
        master_path,
        pathways_dir,
        args.slice_format,
        index,
        input_path,
        total_raw_rows,
        total_master_rows,
        total_dupes_removed,
    )

    return {
        "tissue_id": tissue_id,
        "source_file": str(input_path),
        "master_parquet": f"{tissue_id}/pcxn_master.parquet",
        "pathway_index": f"{tissue_id}/pathway_index.json",
        "pathways_dir": f"{tissue_id}/pathways/",
        "slice_format": args.slice_format,
        "total_pathways": len(index),
        "total_edges": total_master_rows,
        "validation_passed": ok,
    }


def write_tissue_manifest(path: Path, entries: List[dict], pretty: bool) -> None:
    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "total_tissues": len(entries),
        "tissues": sorted(entries, key=lambda e: e["tissue_id"]),
    }
    with open(path, "w", encoding="utf-8") as fh:
        if pretty:
            json.dump(payload, fh, indent=2)
        else:
            json.dump(payload, fh, separators=(",", ":"))
    logger.info("Wrote tissue manifest: %s (%d tissues)", path, len(entries))


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(message)s")

    output_root: Path = args.output_dir
    output_root.mkdir(parents=True, exist_ok=True)

    input_paths = resolve_input_paths(args)
    tissue_inputs = build_tissue_inputs(input_paths, args.tissue_suffix)
    logger.info("Discovered %d tissue file(s): %s", len(tissue_inputs), [t for t, _ in tissue_inputs])

    manifest_entries = []
    overall_ok = True
    for tissue_id, input_path in tissue_inputs:
        logger.info("=" * 70)
        logger.info("TISSUE: %s  (source: %s)", tissue_id, input_path)
        logger.info("=" * 70)
        try:
            summary = run_pipeline_for_tissue(tissue_id, input_path, output_root / tissue_id, args)
            manifest_entries.append(summary)
            if not summary["validation_passed"]:
                overall_ok = False
        except Exception:
            logger.exception("Tissue %r failed - skipping", tissue_id)
            overall_ok = False

    write_tissue_manifest(output_root / "tissues_index.json", manifest_entries, args.pretty_index)

    return 0 if overall_ok else 1


if __name__ == "__main__":
    sys.exit(main())
