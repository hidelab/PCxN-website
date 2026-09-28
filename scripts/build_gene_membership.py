#!/usr/bin/env python3
"""
Builds gene_membership.json for one geneset build, directly from the
sibling pcxn_v2 R pipeline's own pathway_table.csv - no MSigDB GMT files
needed (the previous version of this script depended on GMT files that
only ever existed transiently on a different machine and never persisted
here).

Gene membership is a property of a geneset build, not of any one tissue:
every tissue built from the same pathway_table.csv shares the identical
pathway universe (confirmed empirically - see docs/gene_membership_plan.md),
so one gene_membership.json per build is shared across every tissue that
uses it. tissues_index.json should point each tissue at the right build's
file.

Inputs:
  --pathway-table   A build's pathway_table.csv (set_name,genes,node_type;
                     genes = Entrez IDs), e.g.
                     pcxn_v2/input/gene_sets/<build>/pathway_table.csv
  --entrez-symbol-db  A local Entrez ID -> symbol/description JSON snapshot
                     (see snapshot_entrez_symbols.py), not a SQL database.
  --build-name      Label for this build (default: pathway-table's parent
                     directory name).
  --output          Where to write gene_membership.json.

Usage:
    python scripts/build_gene_membership.py \\
      --pathway-table pcxn_v2/input/gene_sets/MSigDB2024__gtextoil_gfilter/pathway_table.csv \\
      --entrez-symbol-db scripts/resources/entrez_gene_symbols.json \\
      --output output/gene_membership_MSigDB2024__gtextoil_gfilter.json
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import polars as pl

SCRIPTS_DIR = Path(__file__).resolve().parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from tag_pathway_collections import classify  # noqa: E402

logger = logging.getLogger("build_gene_membership")


def load_pathway_table(path: Path) -> pl.DataFrame:
    df = pl.read_csv(path, schema_overrides={"genes": pl.Utf8})
    return df.group_by("set_name").agg(pl.col("genes").alias("entrez_ids"))


def load_symbol_snapshot(path: Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def build_membership(pathway_df: pl.DataFrame, symbol_map: dict) -> dict:
    pathways: dict = {}
    all_entrez_ids: set = set()
    unmatched_entrez_ids: set = set()

    for set_name, entrez_ids in pathway_df.iter_rows():
        genes = []
        for entrez_id in entrez_ids:
            all_entrez_ids.add(entrez_id)
            info = symbol_map.get(entrez_id)
            if info is None:
                unmatched_entrez_ids.add(entrez_id)
                continue
            genes.append({"entrez": entrez_id, "symbol": info["symbol"], "description": info["description"]})
        collection, source_db = classify(set_name)
        pathways[set_name] = {
            "size": len(genes),
            "genes": genes,
            "collection": collection,
            "source_db": source_db,
        }

    return {
        "universe_size": len(all_entrez_ids),
        "unique_entrez_ids_matched": len(all_entrez_ids) - len(unmatched_entrez_ids),
        "unique_entrez_ids_unmatched": len(unmatched_entrez_ids),
        "total_pathways": len(pathways),
        "pathways": pathways,
    }, sorted(unmatched_entrez_ids)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pathway-table", required=True, type=Path, help="A build's pathway_table.csv")
    p.add_argument("--entrez-symbol-db", required=True, type=Path, help="Entrez ID -> symbol/description JSON snapshot")
    p.add_argument("--build-name", default=None, help="Label for this build (default: pathway-table's parent directory name)")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--log-level", default="INFO")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(message)s")

    build_name = args.build_name or args.pathway_table.resolve().parent.name

    logger.info("Loading pathway table: %s", args.pathway_table)
    pathway_df = load_pathway_table(args.pathway_table)
    logger.info("Loading Entrez->symbol snapshot: %s", args.entrez_symbol_db)
    symbol_map = load_symbol_snapshot(args.entrez_symbol_db)

    result, unmatched_entrez_ids = build_membership(pathway_df, symbol_map)

    membership = {
        "build": build_name,
        "pathway_table_source": str(args.pathway_table),
        "entrez_symbol_source": str(args.entrez_symbol_db),
        **result,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(membership, separators=(",", ":")), encoding="utf-8")

    logger.info("=" * 70)
    logger.info("COVERAGE REPORT: %s", build_name)
    logger.info("=" * 70)
    logger.info("Total pathways:                %d", result["total_pathways"])
    logger.info("Unique Entrez IDs in universe:  %d", result["universe_size"])
    logger.info("Matched to a symbol:            %d", result["unique_entrez_ids_matched"])
    logger.info("Unmatched (no symbol found):    %d", result["unique_entrez_ids_unmatched"])
    if unmatched_entrez_ids:
        logger.info("Sample unmatched Entrez IDs: %s", unmatched_entrez_ids[:20])
    logger.info("Wrote %s (%.1f MB)", args.output, args.output.stat().st_size / 1e6)

    return 0


if __name__ == "__main__":
    sys.exit(main())
