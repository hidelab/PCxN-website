#!/usr/bin/env python3
"""
One-off extraction: dump org.Hs.eg.db's Entrez ID -> gene symbol/description
mapping into a plain local JSON file, so downstream scripts (e.g.
build_gene_membership.py) never need a SQL engine or a dependency on
another user's shared R library path at runtime.

This is NOT part of the regular pipeline - it only needs to be re-run if a
fresher annotation snapshot is ever wanted later.

Source: org.Hs.eg.db (Bioconductor) 3.18.0, whose bundled SQLite database
has a `genes` table (_id PK, gene_id = Entrez ID as text) joined to
`gene_info` (_id FK, symbol, gene_name). Entrez ID -> symbol is 1:1 by
construction (both _id columns are unique-constrained).

Usage:
    python scripts/snapshot_entrez_symbols.py \\
      --org-hs-eg-db /data/resources/tools/R/4.3.2/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite \\
      --output scripts/resources/entrez_gene_symbols.json
"""

from __future__ import annotations

import argparse
import json
import logging
import sqlite3
import sys
from pathlib import Path

logger = logging.getLogger("snapshot_entrez_symbols")

DEFAULT_SOURCE = "/data/resources/tools/R/4.3.2/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite"


def extract(db_path: Path) -> dict:
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        rows = con.execute(
            "SELECT g.gene_id, gi.symbol, gi.gene_name FROM genes g JOIN gene_info gi ON g._id = gi._id"
        ).fetchall()
    finally:
        con.close()
    return {gene_id: {"symbol": symbol, "description": gene_name} for gene_id, symbol, gene_name in rows}


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--org-hs-eg-db", type=Path, default=Path(DEFAULT_SOURCE), help="Path to org.Hs.eg.db's bundled SQLite file")
    p.add_argument("--output", type=Path, default=Path("scripts/resources/entrez_gene_symbols.json"))
    p.add_argument("--log-level", default="INFO")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(message)s")

    if not args.org_hs_eg_db.exists():
        raise SystemExit(f"Not found: {args.org_hs_eg_db}")

    mapping = extract(args.org_hs_eg_db)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(mapping, separators=(",", ":")), encoding="utf-8")

    logger.info("Extracted %d Entrez ID -> symbol/description entries from %s", len(mapping), args.org_hs_eg_db)
    logger.info("Wrote %s (%.1f MB)", args.output, args.output.stat().st_size / 1e6)
    return 0


if __name__ == "__main__":
    sys.exit(main())
