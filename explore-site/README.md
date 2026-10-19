# PCxN Explore (static edition)

A no-backend reimplementation of the old Explore mode (selected pathways →
correlated neighborhood). Plain HTML/CSS/JS, no build step. Data is fetched
at runtime from `config.js`'s `PRODUCTION_DATA_BASE_URL` (the Hugging Face
dataset [sgupta15/pcxn-explore-data](https://huggingface.co/datasets/sgupta15/pcxn-explore-data)).

The question it answers: which pathways correlate with the one you are studying?

See the repo-root [README.md](../README.md) for features, citation, and
publishing.

## Run the site

No backend process is needed — just a static file server:

```bash
cd explore-site
python3 -m http.server 8000
```

Open `http://localhost:8000/index.html`. The page loads pathway indexes and
slices from Hugging Face.

## Publish

1. Confirm `config.js`'s `PRODUCTION_DATA_BASE_URL` is the dataset
   `resolve/main` URL (no trailing slash).
2. Push `explore-site/` to GitHub and enable Pages.
3. Optionally set `GITHUB_REPO_URL` so the header GitHub link appears.
   Documentation is always shown and opens `documentation.html` (the
   project README).
