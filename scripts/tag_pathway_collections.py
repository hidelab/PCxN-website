#!/usr/bin/env python3
"""
Enrichment pass over ETL output: tags each pathway in a tissue's
pathway_index.json with a "collection" (broad bucket, matching the old
PCxN app's 4-way MSigDB/Pathprint taxonomy) and "source_db" (finer facet -
which specific database within that bucket), derived from the pathway name
itself via regex.

This is a separate, optional pass — not part of prepare_pcxn_dataset.py -
because the classification is a naming-convention heuristic, not something
knowable from the raw correlation TSVs, and its rules may need tuning per
project independently of the core ETL.

Rules (verified against data/combined_pathway_estimates.tsv - see
scripts/README.md and the design doc for the real prefix counts that
grounded these). The "Pathway." prefix is optional because different
geneset builds use different naming conventions for the same underlying
databases - confirmed by inspecting real data: the older blood/brain build
uses "Pathway.BIOCARTA_..." while the newer MSigDB2024__gtextoil_gfilter
build (combined, CervixUteri_CervixEctocervix) drops that prefix entirely
("BIOCARTA_41BB_PATHWAY"), so requiring it caused every pathway in that
build to silently fall through to Pathprint (2730/2730 on a first run of
this script - not a real "unclassified" count, a naming-convention miss):
  ^(?:Pathway\\.)?(REACTOME|BIOCARTA|PID|KEGG|SIG|SA|NABA|WNT)_   -> MSigDB_C2_CP
  ^HALLMARK_                                                -> MSigDB_H_Hallmark
  ^GO_                                                       -> MSigDB_C5_GO_BP
  anything else                                              -> Pathprint

A tissue's actual mix of these varies - a zero count for a bucket on a given
tissue is expected, not a bug (e.g. this sample data has zero HALLMARK_/GO_
names). The per-tissue coverage summary this script prints is how you'd spot
an actual naming-convention mismatch (a large, unexpected "unclassified"
count) rather than a merely-absent category.

Usage:
    python scripts/tag_pathway_collections.py --output-dir output
    python scripts/tag_pathway_collections.py --pathway-index output/combined/pathway_index.json
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from collections import Counter
from pathlib import Path
from typing import List, Tuple

logger = logging.getLogger("tag_pathway_collections")

# (regex, collection, source_db-or-None). None means "use the regex's own
# capture group as source_db"; order matters, first match wins.
MSIGDB_C2_CP_DBS = ("REACTOME", "BIOCARTA", "PID", "KEGG", "SIG", "SA", "NABA", "WNT", "WP")
RULES: List[Tuple[re.Pattern, str]] = [
    (re.compile(r"^(?:Pathway\.)?(" + "|".join(MSIGDB_C2_CP_DBS) + r")_"), "MSigDB_C2_CP"),
    (re.compile(r"^HALLMARK_"), "MSigDB_H_Hallmark"),
    (re.compile(r"^GO_"), "MSigDB_C5_GO_BP"),
]
PATHPRINT_PREFIX = re.compile(r"^([A-Za-z0-9]+)_")


def classify(pathway_name: str) -> Tuple[str, str]:
    """Returns (collection, source_db) for one pathway name."""
    for pattern, collection in RULES:
        m = pattern.match(pathway_name)
        if m:
            source_db = m.group(1) if m.groups() else collection
            return collection, source_db
    m = PATHPRINT_PREFIX.match(pathway_name)
    source_db = m.group(1) if m else "Pathprint"
    return "Pathprint", source_db


def tag_one_index(index_path: Path) -> None:
    with open(index_path, encoding="utf-8") as fh:
        payload = json.load(fh)

    counts = Counter()
    for entry in payload["pathways"]:
        collection, source_db = classify(entry["pathway_name"])
        entry["collection"] = collection
        entry["source_db"] = source_db
        counts[collection] += 1

    payload["collection_summary"] = dict(sorted(counts.items(), key=lambda kv: -kv[1]))

    with open(index_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, separators=(",", ":"))

    total = sum(counts.values())
    logger.info(
        "%s: %d pathways tagged - %s",
        index_path,
        total,
        ", ".join(f"{k}: {v}" for k, v in payload["collection_summary"].items()),
    )


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output-dir", type=Path, default=None, help="ETL output root - tags every tissue listed in its tissues_index.json")
    p.add_argument("--pathway-index", nargs="+", type=Path, default=None, help="One or more specific pathway_index.json files to tag directly")
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()), format="%(asctime)s %(levelname)s %(message)s")

    targets: List[Path] = []
    if args.pathway_index:
        targets.extend(args.pathway_index)
    if args.output_dir:
        manifest_path = args.output_dir / "tissues_index.json"
        if not manifest_path.exists():
            raise SystemExit(f"No tissues_index.json found under {args.output_dir}")
        with open(manifest_path, encoding="utf-8") as fh:
            manifest = json.load(fh)
        for tissue in manifest["tissues"]:
            targets.append(args.output_dir / tissue["pathway_index"])
    if not targets:
        raise SystemExit("Pass --output-dir and/or --pathway-index.")

    for path in targets:
        if not path.exists():
            logger.error("Not found, skipping: %s", path)
            continue
        tag_one_index(path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
