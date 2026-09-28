# prepare_pcxn_dataset.py

ETL pipeline that turns one or more raw PCxN pairwise pathway co-expression
TSV files - one per tissue - into the artifacts needed for the web deployment
architecture (Cloudflare R2 + GitHub Pages + Workers): a canonical
`pcxn_master.parquet` archive, per-pathway slice files for fast R2/Worker
reads, and `pathway_index.json` / `tissues_index.json` directories for
frontend autocomplete and routing.

See the module docstring at the top of [prepare_pcxn_dataset.py](prepare_pcxn_dataset.py)
for the full output layout and design rationale.

## Setup

This dev host (CentOS 7, glibc 2.17) can't install current polars/pyarrow
wheels with the system Python - pip falls back to a source build that fails
(numpy needs GCC >= 9.3; this host has 4.8.5/8.2.0). Use a conda env with a
newer Python and these pinned, older package versions instead:

```bash
conda create -n pcxn-etl python=3.11 -y
conda activate pcxn-etl
pip install "numpy<2" "pyarrow==14.0.2" "polars==0.20.31"
```

On a machine with a newer glibc (>= 2.28), the unpinned constraints in
[requirements.txt](requirements.txt) (`pip install -r scripts/requirements.txt`)
should work fine and are preferred there.

Next time, just `conda activate pcxn-etl` - no need to recreate the env.

## Running

From the project root, with the env active:

```bash
conda activate pcxn-etl

# Scan a directory for <tissue>_pathway_estimates.tsv[.gz] files:
python scripts/prepare_pcxn_dataset.py \
  --input-dir data \
  --output-dir output \
  --slice-format json

# Or name specific files explicitly (mix of plain/gzipped is fine):
python scripts/prepare_pcxn_dataset.py \
  --input-tsv data/brain_pathway_estimates.tsv data/liver_pathway_estimates.tsv.gz \
  --output-dir output

# Recommended flags when the output feeds explore-site/ (see its README):
# --sort-by abs_correlation matches how Explore ranks neighbors, and slices
# are uncapped by default (no --max-edges-per-slice needed) since data
# hosting there is off GitHub, not size-constrained.
python scripts/prepare_pcxn_dataset.py \
  --input-dir data \
  --output-dir output \
  --slice-format json \
  --sort-by abs_correlation
```

After generating `output/`, tag each pathway with a `collection`/`source_db`
facet (see [tag_pathway_collections.py](tag_pathway_collections.py) for the
classifier rules) - `explore-site/`'s collection filter reads these fields:

```bash
python scripts/tag_pathway_collections.py --output-dir output
```

Each input file is treated as one tissue (id derived from the filename by
stripping `_pathway_estimates` - configurable via `--tissue-suffix`) and
processed into its own subtree:

```
output/
  <tissue_id>/pcxn_master.parquet
  <tissue_id>/pathways/<pathway_id>.json   (or .parquet, via --slice-format)
  <tissue_id>/pathway_index.json
  tissues_index.json                        <- manifest of all tissues
```

The script logs progress per phase (ingest -> shuffle, per-bucket
processing, index write) and per tissue, then ends with a validation report
(row counts, compression ratio, slice file stats, and a round-trip check on
a sample pathway). A tissue that fails is logged and skipped rather than
aborting the whole batch; the process exits non-zero if any tissue failed
or failed validation.

Run `python scripts/prepare_pcxn_dataset.py --help` for the full list of
options (column overrides, `--num-buckets` / `--batch-size-mb` for tuning
memory use on very large files, `--pvalue-mode`, `--sort-by`,
`--max-edges-per-slice`, etc).

## Running the whole pipeline in one command

[run_pipeline.py](run_pipeline.py) chains all three scripts above -
`convert_bin_to_tsv.py` (optional, off by default) -> `prepare_pcxn_dataset.py`
-> `tag_pathway_collections.py` - driven by one YAML config instead of
running each by hand:

```bash
cp scripts/pipeline.example.yaml pipeline.yaml   # then edit paths in it
python scripts/run_pipeline.py --config pipeline.yaml

# See the resolved command for each stage without running anything:
python scripts/run_pipeline.py --config pipeline.yaml --dry-run
```

Each stage has its own `enabled` switch in the config - e.g. leave
`bin_conversion.enabled: false` (the default) once `data/` already has the
`<tissue>_pathway_estimates.tsv` files you need, and only `website_etl` +
`tag_collections` will run. A stage that fails stops the pipeline
immediately, since each stage's output feeds the next. See
[pipeline.example.yaml](pipeline.example.yaml) for every available option
and what each one maps to.

## Gene-to-pathway membership (`gene_membership.json`)

[build_gene_membership.py](build_gene_membership.py) builds a
`gene_membership.json` for one geneset **build** (not per tissue - every
tissue built from the same `pathway_table.csv` shares the identical
pathway universe, so one file covers all of them):

```bash
python scripts/build_gene_membership.py \
  --pathway-table /path/to/pcxn_v2/input/gene_sets/<build>/pathway_table.csv \
  --entrez-symbol-db scripts/resources/entrez_gene_symbols.json \
  --output output/gene_membership/<build>.json
```

It reads the sibling `pcxn_v2` R pipeline's own `pathway_table.csv`
directly (no MSigDB GMT files needed - those are what the old version of
this script depended on, and they don't persist on this host), joins each
pathway's Entrez IDs against a local symbol lookup, and tags each pathway
with the same `collection`/`source_db` classifier
[tag_pathway_collections.py](tag_pathway_collections.py) uses.

The symbol lookup (`--entrez-symbol-db`) is a plain JSON file, not a SQL
database - a one-time snapshot of `org.Hs.eg.db`, produced by
[snapshot_entrez_symbols.py](snapshot_entrez_symbols.py):

```bash
python scripts/snapshot_entrez_symbols.py \
  --org-hs-eg-db /data/resources/tools/R/4.3.2/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite \
  --output scripts/resources/entrez_gene_symbols.json
```

This only needs to be re-run if a fresher annotation snapshot is ever
wanted - it's not part of the regular pipeline. See
[../docs/gene_membership_plan.md](../docs/gene_membership_plan.md) for the
full design rationale (why `pathway_table.csv` instead of GMT files, why
Entrez→symbol needs no per-tissue input, and how this compares to the
previous GMT-based approach).
