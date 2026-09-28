#!/usr/bin/env python3
"""
Runs the full .bin -> tsv -> parquet/json -> tagged pipeline as one command,
driven by a single YAML config file, instead of invoking
convert_bin_to_tsv.py / prepare_pcxn_dataset.py / tag_pathway_collections.py
by hand in sequence.

Three stages, always in this order (each independently toggleable):
  1. bin_conversion   -> convert_bin_to_tsv.py   (skip if you already have
                         <tissue>_pathway_estimates.tsv files in data/)
  2. website_etl      -> prepare_pcxn_dataset.py (the core ETL; on by default)
  3. tag_collections  -> tag_pathway_collections.py (on by default)

Each stage's script is imported and invoked via its own `main(argv)` - the
exact same code path as running it standalone, just orchestrated from one
config file instead of copy-pasted commands. A stage that returns/raises a
non-zero exit stops the pipeline immediately (later stages depend on earlier
ones' output, so there's no reason to continue).

See pipeline.example.yaml for a documented starting config.

Usage:
    python scripts/run_pipeline.py --config pipeline.yaml
    python scripts/run_pipeline.py --config pipeline.yaml --dry-run
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List

import yaml

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import convert_bin_to_tsv  # noqa: E402
import prepare_pcxn_dataset  # noqa: E402
import tag_pathway_collections  # noqa: E402

logger = logging.getLogger("run_pipeline")


def _flag(args: List[str], name: str, enabled: bool) -> None:
    if enabled:
        args.append(name)


def _as_csv(value: Any) -> str:
    return ",".join(value) if isinstance(value, list) else str(value)


def build_bin_conversion_args(cfg: Dict[str, Any]) -> List[str]:
    args = ["--input-dir", str(cfg["input_dir"])]
    args += ["--output-dir", str(cfg.get("output_dir", "data"))]
    if cfg.get("tissues"):
        args += ["--tissues", _as_csv(cfg["tissues"])]
    _flag(args, "--skip-combined", cfg.get("skip_combined", False))
    _flag(args, "--gzip", cfg.get("gzip", False))
    args += ["--log-level", str(cfg.get("log_level", "INFO"))]
    return args


def build_website_etl_args(cfg: Dict[str, Any]) -> List[str]:
    args: List[str] = []
    if cfg.get("input_tsv"):
        args += ["--input-tsv", *[str(p) for p in cfg["input_tsv"]]]
    if cfg.get("input_dir"):
        args += ["--input-dir", str(cfg["input_dir"])]
    if cfg.get("glob_pattern"):
        args += ["--glob-pattern", str(cfg["glob_pattern"])]
    if cfg.get("tissue_suffix"):
        args += ["--tissue-suffix", str(cfg["tissue_suffix"])]
    args += ["--output-dir", str(cfg["output_dir"])]
    args += ["--slice-format", str(cfg.get("slice_format", "json"))]

    for key, flag in [
        ("col_pathway_a", "--col-pathway-a"),
        ("col_pathway_b", "--col-pathway-b"),
        ("col_correlation", "--col-correlation"),
        ("col_pvalue", "--col-pvalue"),
        ("col_fdr", "--col-fdr"),
    ]:
        if cfg.get(key):
            args += [flag, str(cfg[key])]

    if cfg.get("extra_cols"):
        args += ["--extra-cols", _as_csv(cfg["extra_cols"])]
    if cfg.get("num_buckets"):
        args += ["--num-buckets", str(cfg["num_buckets"])]
    if cfg.get("batch_size_mb"):
        args += ["--batch-size-mb", str(cfg["batch_size_mb"])]

    args += ["--pvalue-mode", str(cfg.get("pvalue_mode", "neg_log10"))]
    args += ["--sort-by", str(cfg.get("sort_by", "fdr"))]
    if cfg.get("max_edges_per_slice") is not None:
        args += ["--max-edges-per-slice", str(cfg["max_edges_per_slice"])]

    _flag(args, "--pretty-index", cfg.get("pretty_index", False))
    _flag(args, "--keep-temp", cfg.get("keep_temp", False))
    args += ["--log-level", str(cfg.get("log_level", "INFO"))]
    return args


def build_tag_collections_args(cfg: Dict[str, Any]) -> List[str]:
    args: List[str] = ["--output-dir", str(cfg["output_dir"])]
    args += ["--log-level", str(cfg.get("log_level", "INFO"))]
    return args


def run_stage(name: str, main_fn, argv: List[str], dry_run: bool) -> None:
    logger.info("=" * 70)
    logger.info("STAGE: %s", name)
    logger.info("=" * 70)
    logger.info("argv: %s", " ".join(argv))
    if dry_run:
        return
    rc = main_fn(argv)
    if rc:
        raise SystemExit(f"Stage '{name}' failed (exit code {rc}) - stopping pipeline.")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True, type=Path, help="Pipeline config YAML file (see pipeline.example.yaml)")
    p.add_argument("--dry-run", action="store_true", help="Print each stage's resolved argv without running anything")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(message)s")

    with open(args.config, encoding="utf-8") as fh:
        config = yaml.safe_load(fh) or {}

    bin_cfg = config.get("bin_conversion") or {}
    etl_cfg = config.get("website_etl") or {}
    tag_cfg = config.get("tag_collections") or {}

    if bin_cfg.get("enabled", False):
        if not bin_cfg.get("input_dir"):
            raise SystemExit(
                "bin_conversion.enabled is true but bin_conversion.input_dir is not set "
                "(must contain combined_estimates/pathway_estimates.tsv and "
                "mean_pcor2_barcode_tables/*.bin - see convert_bin_to_tsv.py's docstring)."
            )
        run_stage("bin_conversion", convert_bin_to_tsv.main, build_bin_conversion_args(bin_cfg), args.dry_run)
    else:
        logger.info("Skipping bin_conversion stage (bin_conversion.enabled is not true).")

    if etl_cfg.get("enabled", True):
        if not etl_cfg.get("output_dir"):
            raise SystemExit("website_etl.output_dir is required.")
        if not etl_cfg.get("input_dir") and not etl_cfg.get("input_tsv"):
            raise SystemExit("website_etl needs input_dir and/or input_tsv.")
        run_stage("website_etl", prepare_pcxn_dataset.main, build_website_etl_args(etl_cfg), args.dry_run)
    else:
        logger.info("Skipping website_etl stage (website_etl.enabled is false).")

    if tag_cfg.get("enabled", True):
        output_dir = tag_cfg.get("output_dir") or etl_cfg.get("output_dir")
        if not output_dir:
            raise SystemExit(
                "tag_collections needs an output_dir - set tag_collections.output_dir, "
                "or rely on website_etl.output_dir."
            )
        run_stage(
            "tag_collections",
            tag_pathway_collections.main,
            build_tag_collections_args({**tag_cfg, "output_dir": output_dir}),
            args.dry_run,
        )
    else:
        logger.info("Skipping tag_collections stage (tag_collections.enabled is false).")

    logger.info("Pipeline finished successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
