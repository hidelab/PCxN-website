#!/usr/bin/env python3
"""
Converts a whole `.bin`-format pipeline run into the standard
`<tissue>_pathway_estimates.tsv` files prepare_pcxn_dataset.py already knows
how to ingest via --input-dir.

Background (fully traced against the sibling `pcxn_v2` R pipeline, not
guessed): a run's output directory always has this shape (confirmed against
several real runs, e.g. `pcxn_v2/output/mac_run/<geneset_build>/`):

    <run_dir>/
      combined_estimates/
        pathway_estimates.tsv       <- final combined tissue, already a
                                        complete Pathway.A/B/PathCor/p.value/
                                        p.adj TSV (pcxn_02_combine.R's output)
                                        AND the pairing/row-order key every
                                        .bin file below shares
      mean_pcor2_barcode_tables/
        <tissue>_estimate.bin        <- one triplet of headerless raw
        <tissue>_pvalue.bin             little-endian float64 arrays per
        <tissue>_statistic.bin         tissue (pcxn_01_explevel.R's output)
        <tissue2>_estimate.bin
        ...

The `.bin` files carry no pathway names or header of any kind - just the
raw values, in the same row order as `combined_estimates/pathway_estimates.tsv`.
This script uses that TSV's Pathway.A/Pathway.B columns as the shared
pairing key for every tissue's `.bin` triplet (joined by position), copies
the combined TSV through as-is, and computes each tissue's own BH-FDR from
its own p-values (pcxn_02_combine.R only ever computes p.adj at its later
multi-experiment combine stage, so a single experiment-level `_pvalue.bin`
has no p.adj of its own).

Usage:
    python scripts/convert_bin_to_tsv.py --input-dir /path/to/run_dir --output-dir data
    python scripts/convert_bin_to_tsv.py --input-dir /path/to/run_dir --tissues CervixUteri_CervixEctocervix,Blood_WholeBlood
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import polars as pl

logger = logging.getLogger("convert_bin_to_tsv")

FLOAT64_SIZE = 8
COMBINED_SUBDIR = "combined_estimates"
COMBINED_FILENAME = "pathway_estimates.tsv"
PER_TISSUE_SUBDIR = "mean_pcor2_barcode_tables"

# Same alias style as prepare_pcxn_dataset.py's COLUMN_ALIASES, so a
# combined pathway_estimates.tsv with slightly different header spellings
# still works.
PATHWAY_A_CANDIDATES = ["pathway_a", "pathwaya", "pathway1", "pathwayasource", "sourcepathway", "nodea", "source"]
PATHWAY_B_CANDIDATES = ["pathway_b", "pathwayb", "pathway2", "pathwaybtarget", "targetpathway", "nodeb", "target"]


def _normalize_header(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def resolve_pair_columns(header: List[str]):
    normalized = {_normalize_header(h): h for h in header}

    def find(candidates, label):
        for cand in candidates:
            key = _normalize_header(cand)
            if key in normalized:
                return normalized[key]
        raise SystemExit(f"Could not find a '{label}'-like column in {COMBINED_FILENAME}'s header: {header}")

    return find(PATHWAY_A_CANDIDATES, "pathway_a"), find(PATHWAY_B_CANDIDATES, "pathway_b")


def load_pairs(path: Path) -> pl.DataFrame:
    header = pl.read_csv(path, separator="\t", n_rows=0).columns
    a_col, b_col = resolve_pair_columns(header)
    logger.info("Pairing key %s: using %r as Pathway.A, %r as Pathway.B", path, a_col, b_col)
    lf = pl.scan_csv(path, separator="\t", schema_overrides={a_col: pl.Utf8, b_col: pl.Utf8})
    return lf.select([pl.col(a_col).alias("Pathway.A"), pl.col(b_col).alias("Pathway.B")]).collect()


def read_bin_values(path: Path) -> np.ndarray:
    size = path.stat().st_size
    if size % FLOAT64_SIZE != 0:
        raise ValueError(
            f"{path}: size {size} bytes is not a multiple of {FLOAT64_SIZE} - "
            "not a valid headerless float64 .bin file."
        )
    return np.fromfile(path, dtype="<f8")


def require_matching_length(label: str, values: np.ndarray, n_pairs: int) -> None:
    if len(values) != n_pairs:
        raise ValueError(
            f"Row-count mismatch: the pairing key has {n_pairs} rows, but {label} has "
            f"{len(values)} values. These must come from the exact same run/build - "
            "refusing to silently mismatch every row."
        )


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Standard BH-FDR (same procedure pcxn_02_combine.R's p.adjust(method="fdr")
    uses), computed directly (no scipy/statsmodels dependency): rank ascending,
    p * n / rank, enforce monotonicity via a reverse running minimum, clip to
    [0, 1]."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranks = np.arange(1, n + 1)
    adjusted = pvals[order] * n / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return out


def copy_combined(combined_tsv: Path, output_dir: Path, gzip_output: bool) -> Path:
    suffix = ".tsv.gz" if gzip_output else ".tsv"
    out_path = output_dir / f"combined_pathway_estimates{suffix}"
    if gzip_output:
        with open(combined_tsv, "rb") as src, gzip.open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copyfile(combined_tsv, out_path)
    logger.info("Copied %s -> %s (%.1f MB)", combined_tsv, out_path, out_path.stat().st_size / 1e6)
    return out_path


def discover_tissues(mean_pcor_dir: Path) -> List[str]:
    return sorted(p.name[: -len("_estimate.bin")] for p in mean_pcor_dir.glob("*_estimate.bin"))


def process_tissue(tissue: str, mean_pcor_dir: Path, pairs_df: pl.DataFrame, output_dir: Path, gzip_output: bool) -> dict:
    n_pairs = pairs_df.height
    estimate_path = mean_pcor_dir / f"{tissue}_estimate.bin"
    pvalue_path = mean_pcor_dir / f"{tissue}_pvalue.bin"
    statistic_path = mean_pcor_dir / f"{tissue}_statistic.bin"

    estimate = read_bin_values(estimate_path)
    require_matching_length(str(estimate_path), estimate, n_pairs)

    data = {
        "Pathway.A": pairs_df["Pathway.A"],
        "Pathway.B": pairs_df["Pathway.B"],
        "PathCor": estimate.astype(np.float32),
    }

    has_pvalue = pvalue_path.exists()
    padj = None
    pvalue = None
    if has_pvalue:
        pvalue = read_bin_values(pvalue_path)
        require_matching_length(str(pvalue_path), pvalue, n_pairs)
        padj = benjamini_hochberg(pvalue)
        data["p.value"] = pvalue.astype(np.float32)
        data["p.adj"] = padj.astype(np.float32)
    else:
        logger.warning(
            "%s: no matching %s - output will have Pathway.A/Pathway.B/PathCor only, "
            "no p.value/p.adj. prepare_pcxn_dataset.py requires those columns, so this "
            "tissue can't be ingested by it until a matching _pvalue.bin is supplied.",
            tissue,
            pvalue_path.name,
        )

    if statistic_path.exists():
        statistic = read_bin_values(statistic_path)
        require_matching_length(str(statistic_path), statistic, n_pairs)
        data["statistic"] = statistic.astype(np.float32)

    out_df = pl.DataFrame(data)
    suffix = ".tsv.gz" if gzip_output else ".tsv"
    out_path = output_dir / f"{tissue}_pathway_estimates{suffix}"
    if gzip_output:
        with gzip.open(out_path, "wt") as fh:
            out_df.write_csv(fh, separator="\t")
    else:
        out_df.write_csv(out_path, separator="\t")

    frac_sig = float(np.mean(padj < 0.05)) if padj is not None else None
    logger.info(
        "%s: wrote %s (%d rows, %.1f MB)%s",
        tissue,
        out_path,
        out_df.height,
        out_path.stat().st_size / 1e6,
        f", {frac_sig:.4f} of pairs at FDR<0.05" if frac_sig is not None else "",
    )
    return {
        "tissue": tissue,
        "output": str(out_path),
        "rows": out_df.height,
        "has_pvalue": has_pvalue,
        "frac_significant": frac_sig,
    }


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help=f"Run output directory containing {COMBINED_SUBDIR}/{COMBINED_FILENAME} and {PER_TISSUE_SUBDIR}/*_estimate.bin",
    )
    p.add_argument("--output-dir", type=Path, default=Path("data"), help="Where to write <tissue>_pathway_estimates.tsv files")
    p.add_argument(
        "--tissues",
        type=str,
        default=None,
        help="Comma-separated tissue names to process (default: every tissue found in mean_pcor2_barcode_tables/)",
    )
    p.add_argument("--skip-combined", action="store_true", help="Don't copy the combined_estimates TSV through")
    p.add_argument("--gzip", action="store_true", help="Write .tsv.gz instead of .tsv")
    p.add_argument("--log-level", default="INFO")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(message)s")

    combined_tsv = args.input_dir / COMBINED_SUBDIR / COMBINED_FILENAME
    if not combined_tsv.exists():
        raise SystemExit(f"Expected pairing-key file not found: {combined_tsv}")
    mean_pcor_dir = args.input_dir / PER_TISSUE_SUBDIR
    if not mean_pcor_dir.exists():
        raise SystemExit(f"Expected directory not found: {mean_pcor_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    pairs_df = load_pairs(combined_tsv)
    n_pairs = pairs_df.height
    logger.info("Pairing key has %d rows.", n_pairs)

    if not args.skip_combined:
        copy_combined(combined_tsv, args.output_dir, args.gzip)

    available = discover_tissues(mean_pcor_dir)
    if not available:
        raise SystemExit(f"No *_estimate.bin files found under {mean_pcor_dir}")

    if args.tissues:
        requested = [t.strip() for t in args.tissues.split(",") if t.strip()]
        missing = [t for t in requested if t not in available]
        if missing:
            logger.error("Requested tissue(s) not found in %s: %s", mean_pcor_dir, ", ".join(missing))
        tissues = [t for t in requested if t in available]
    else:
        tissues = available

    logger.info("Processing %d tissue(s): %s", len(tissues), ", ".join(tissues))

    results = []
    failed = []
    for tissue in tissues:
        try:
            results.append(process_tissue(tissue, mean_pcor_dir, pairs_df, args.output_dir, args.gzip))
        except Exception:
            logger.exception("Tissue %r failed - skipping", tissue)
            failed.append(tissue)

    logger.info("=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    logger.info("Succeeded: %d/%d tissues", len(results), len(tissues))
    no_pvalue = [r["tissue"] for r in results if not r["has_pvalue"]]
    if no_pvalue:
        logger.info("Correlation-only (no matching _pvalue.bin): %s", ", ".join(no_pvalue))
    if failed:
        logger.error("Failed: %s", ", ".join(failed))
    logger.info("=" * 70)

    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
