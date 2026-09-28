# PCxN Explore

A static **Pathway Co-Activity Map**. Given a pathway — or a group of pathways — it shows which other pathways rise and fall with it across public expression experiments, including pairs that share no genes.

The site is plain HTML, CSS, and JavaScript. There is no backend. The UI lives in `explore-site/`. Pathway slices and indexes are served from [Hugging Face](https://huggingface.co/datasets/sgupta15/pcxn-explore-data).

**The question it answers:** which pathways correlate with the one you are studying?

## What you can do

- Pick one or more **tissue-cell types** (54 GTEx tissues plus a `combined` default) and a **collection** (Reactome, WikiPathways, KEGG, BioCarta, PID, and smaller source DBs).
- Add pathways to one or more named **groups**, then **Build network**.
- See the neighborhood three ways: a Cytoscape **network**, an **edge table**, and a **heatmap**.
- Look up pathways **from a gene list** (over-representation with BH-adjusted q-values). Example lists in that dialog are illustrative marker genes, not experimental data.
- Hide or add nodes from the graph table, view **genes** in a pathway, and view **linked pathways** with correlation and p-value.
- Share a result: the URL keeps tissue, groups, pathways, and cutoffs.

Correlation is Pearson-style pathway co-activity (`r`). Red edges are positive, blue negative; color intensity follows `|r|`. Node size is degree in the current graph. Sliders re-filter without another fetch.

A query is capped at **60 slice fetches**: `tissues × (selected pathways + topN) ≤ 60`. Lower topN, fewer tissues, or fewer seeds if you hit the limit. TopN `0` draws only edges among the selected pathways (used by the two-group comparison example).

## Data

| | |
| --- | --- |
| Pathways | 2,730 |
| Tissue-cell types | 54 (plus a `combined` network) |
| Collections | 8, from each pathway's `source_db` |
| Hosted data | [sgupta15/pcxn-explore-data](https://huggingface.co/datasets/sgupta15/pcxn-explore-data) |

Each tissue folder has a `pathway_index.json` and one JSON **slice** per pathway (that pathway's correlations to every other pathway). `tissues_index.json` lists tissues. `gene_membership.json` is tissue-independent and is loaded only when you open From genes or view genes.

`explore-site/config.js` points the site at the Hugging Face dataset's `resolve/main` URL.

## Repository layout

```
explore-site/     # live UI (host this on GitHub Pages)
scripts/          # ETL: TSV → per-pathway JSON slices + indexes
.nojekyll         # allow GitHub Pages to serve this tree as-is
```

## Run the site

Serve `explore-site/` with any static file server. The page fetches data from Hugging Face.

```bash
cd explore-site
python3 -m http.server 8000
```

Open `http://localhost:8000/index.html`.

From the repo root you can instead run `python3 -m http.server 8000` and open `/explore-site/index.html`. Any static file server works (`npx serve`, VS Code Live Server).

## Publish

1. Confirm `explore-site/config.js` `PRODUCTION_DATA_BASE_URL` is the Hugging Face dataset's `resolve/main` URL (no trailing slash).
2. Push this repo and enable **GitHub Pages** on `explore-site/` (or the repo root if you prefer `/explore-site/index.html`).
3. Set `GITHUB_REPO_URL` in `config.js` if you want the header GitHub link. Documentation is always shown and opens the in-site README.

To refresh the hosted data, rebuild with `scripts/` (see [scripts/README.md](scripts/README.md)) and upload the generated tree to the Hugging Face dataset: `tissues_index.json` at the root, one folder per tissue, plus `gene_membership.json`.

## How to cite

A manuscript describing this tool is in preparation. Citation details will be added here when it is available.
