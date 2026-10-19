// Data is served from this Hugging Face dataset (no trailing slash).
const PRODUCTION_DATA_BASE_URL = "https://huggingface.co/datasets/sgupta15/pcxn-explore-data/resolve/main";

function resolveDataBaseUrl() {
  return PRODUCTION_DATA_BASE_URL.replace(/\/$/, "");
}

const DATA_BASE_URL = resolveDataBaseUrl();

// Optional public repo link shown in the header. Leave empty to hide it.
const GITHUB_REPO_URL = "";

// Header Documentation link. Points at the in-site rendering of README.md.
const DOCUMENTATION_URL = "documentation.html";
