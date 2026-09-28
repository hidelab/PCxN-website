"use strict";
/* PCxN Explore - static edition.
 * Reads config.js's DATA_BASE_URL, then replicates the old app's two-stage
 * "getCorrelation.jsp" query as static-file fetches (see the design doc /
 * scripts/prepare_pcxn_dataset.py for the data shape this consumes).
 *
 * Multi-pathway, multi-group Explore (see the design doc's later addition):
 * you can select pathways into one or more named "groups" (like the legacy
 * app's up/down-regulated seed lists, generalized to N groups). The query
 * algorithm treats the union of every group's pathways as one combined seed
 * set - group identity only ever affects rendering (node color), never the
 * query itself, exactly like the legacy app's up/down split never affected
 * its SQL query, only its frontend coloring.
 */

// ---------- DOM ----------
const el = (id) => document.getElementById(id);
const errorBanner = el("error-banner");
const errorBannerText = el("error-banner-text");
const tissueList = el("tissueList");
const tissueCombobox = el("tissueCombobox");
const tissueSummary = el("tissueSummary");
const tissueMenu = el("tissueMenu");
const tissueFilter = el("tissueFilter");
const tissueTabList = el("tissueTabList");
const collectionSelect = el("collectionSelect");
const groupsContainer = el("groupsContainer");
const addGroupButton = el("addGroupButton");
const fetchCounterEl = el("fetchCounter");
const topNInput = el("topNInput");
const corrCutoffInput = el("corrCutoffInput");
const pCutoffInput = el("pCutoffInput");
const queryForm = el("query-form");
const goButton = el("goButton");
const resetButton = el("resetButton");
const queryStatus = el("query-status");
const indexMeta = el("index-meta");
const emptyState = el("empty-state");
const results = el("results");
const resultMeta = el("result-meta");
const liveFilters = el("live-filters");
const corrSlider = el("corrSlider");
const pSlider = el("pSlider");
const corrSliderValue = el("corrSliderValue");
const pSliderValue = el("pSliderValue");
const networkSection = el("network-section");
const heatmapSection = el("heatmap-section");
const tableSection = el("table-section");
const nodeInfo = el("node-info");
const legendEl = el("legend");
const exportCsvButton = el("exportCsvButton");
const copyPathwaysButton = el("copyPathwaysButton");
const copyGenesButton = el("copyGenesButton");
const genesDialog = el("genesDialog");
const genesDialogTitle = el("genesDialogTitle");
const genesDialogMeta = el("genesDialogMeta");
const genesDialogList = el("genesDialogList");
const copyPathwayGenesButton = el("copyPathwayGenesButton");
const genesDialogClose = el("genesDialogClose");
const linksDialog = el("linksDialog");
const linksDialogTitle = el("linksDialogTitle");
const linksDialogMeta = el("linksDialogMeta");
const linksDialogList = el("linksDialogList");
const copyPathwayLinksButton = el("copyPathwayLinksButton");
const linksDialogClose = el("linksDialogClose");
const layoutSelect = el("layoutSelect");
const fitButton = el("fitButton");
const saveNetworkPngButton = el("saveNetworkPngButton");
const saveHeatmapPngButton = el("saveHeatmapPngButton");
const graphNodeTableBody = el("graphNodeTable").querySelector("tbody");
const graphNodeMeta = el("graphNodeMeta");
const graphShowAll = el("graphShowAll");
const graphAddInput = el("graphAddInput");
const graphAddButton = el("graphAddButton");
const graphAddMenu = el("graphAddMenu");
const copyLinkButton = el("copyLinkButton");
const aboutButton = el("aboutButton");
const aboutDialog = el("aboutDialog");
const bulkLookupDialog = el("bulkLookupDialog");
const bulkLookupTitle = el("bulkLookupTitle");
const bulkLookupTerms = el("bulkLookupTerms");
const bulkLookupCount = el("bulkLookupCount");
const bulkLookupList = el("bulkLookupList");
const bulkLookupSelectAll = el("bulkLookupSelectAll");
const bulkLookupClearVisible = el("bulkLookupClearVisible");
const bulkLookupDone = el("bulkLookupDone");
const fromGenesDialog = el("fromGenesDialog");
const fromGenesTitle = el("fromGenesTitle");
const fromGenesList = el("fromGenesList");
const fromGenesQCutoff = el("fromGenesQCutoff");
const fromGenesStatus = el("fromGenesStatus");
const fromGenesCount = el("fromGenesCount");
const fromGenesTable = el("fromGenesTable");
const fromGenesTableBody = fromGenesTable.querySelector("tbody");
const fromGenesSelectAll = el("fromGenesSelectAll");
const fromGenesClearVisible = el("fromGenesClearVisible");
const fromGenesDone = el("fromGenesDone");
const githubLink = el("githubLink");
const docsLink = el("docsLink");
const citeButton = el("citeButton");
const citeDialog = el("citeDialog");
const homeLink = el("homeLink");
const homeBrand = el("homeBrand");

el("data-base-url-display").textContent = DATA_BASE_URL;
if (typeof GITHUB_REPO_URL === "string" && GITHUB_REPO_URL) {
  githubLink.href = GITHUB_REPO_URL;
  githubLink.hidden = false;
}
if (docsLink) {
  const docsUrl = typeof DOCUMENTATION_URL === "string" && DOCUMENTATION_URL
    ? DOCUMENTATION_URL
    : "documentation.html";
  docsLink.href = docsUrl;
  docsLink.hidden = false;
}

// ---------- App state ----------
let tissuesIndexData = null;
let currentTissueId = null;
let pathwayIndexData = null; // full pathway_index.json for the current tissue
// Edge records inside a slice file identify the OTHER pathway by its raw
// pathway_name (that's all the ETL stores per-edge - see
// scripts/prepare_pcxn_dataset.py's write_json_slice), not by pathway_id.
// Slice filenames, however, are the slugified pathway_id (colons etc.
// replaced - see slugify() in that same script). So fetching a neighbor's
// (or a query member's) own slice requires translating name -> id first;
// this map (keyed by the display-facing pathway_name) does that and doubles
// as the lookup for per-node metadata (collection/source_db) on click.
let pathwayByName = new Map(); // pathway_name -> index entry
let pColumn = null; // e.g. "neg_log10_p_value" or "p_value"
let fdrColumn = null;
let pvalueMode = "raw";

let lastQuery = null; // {groupByName, groups: [{id,name,color,pathways:[name,...]}], corrCutOff, pCutOff}
let cachedEdges = null; // unfiltered induced-subgraph edges for the active tissue
let cachedMaxAbsCorr = 0;
let tissueResults = {}; // tissue_id -> {status, edges, neighborNames, failedQueryNames, failedNames, maxAbsCorr, errorMessage}
let activeTissueId = null;
let selectedTissues = []; // tissue_id[] in tissues_index order
let tissueMenuItems = [];
let tissueMenuIndex = -1;

let cy = null;
let dataTable = null;
let cxChart = null;
let hiddenGraphNames = new Set();
let graphAddItems = [];
let graphAddIndex = -1;

const ROUTING_PARAMS = ["data", "dataBase", "base"];
let activeTab = "network";
// Total fetches per query = (pathways selected across all groups) + topN,
// all parallel. Kept conservative for v1 (see design doc §9 "Validation /
// limits") - both for fetch fan-out and for how dense/unreadable the
// resulting network gets, since stage 2 is deliberately unfiltered by the
// cutoffs (matching the legacy app's spec).
const MAX_TOTAL_FETCHES = 60;

// Node colors for group 1, 2, 3, ... (cycles if you add more groups than
// this has entries). Deliberately avoids the neighbor blue (#2563eb, see
// nodeColor()) so seeds and neighbors never look like the same category.
// Group 1 defaults to the same near-black the app previously used
// unconditionally for "the query pathway", so the common single-group case
// looks exactly as it did before groups existed.
const GROUP_COLORS = ["#111827", "#dc2626", "#059669", "#7c3aed", "#d97706", "#0e7490", "#db2777", "#65a30d"];
const NEIGHBOR_COLOR = "#2563eb";
const ADDED_COLOR = "#5b4a3a";
const SAFE_CSS_COLOR = /^#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})$/;

function cssColor(value, fallback = "#111827") {
  return typeof value === "string" && SAFE_CSS_COLOR.test(value.trim()) ? value.trim() : fallback;
}

// ---------- Groups state ----------
// Each group: {id, name, color, pathways: [pathway_index entry, ...], rowEl,
// combobox: {inputEl, menuEl, containerEl, items, index}}
let groups = [];
let groupIdCounter = 0;
let bulkLookupGroup = null;
let bulkLookupItems = [];
let fromGenesGroup = null;
let fromGenesItems = [];
let fromGenesPrefillKey = null;
let lastFromGenesHits = [];
let fromGenesSortKey = "q";
let fromGenesSortDir = "asc";
let geneMembership = null;
let geneMembershipError = null;

function makeGroup() {
  groupIdCounter += 1;
  const color = GROUP_COLORS[(groupIdCounter - 1) % GROUP_COLORS.length];
  const group = { id: `g${groupIdCounter}`, name: `Group ${groupIdCounter}`, color, pathways: [], rowEl: null, combobox: null };
  groups.push(group);
  return group;
}

function clearAllGroupRows() {
  for (const g of groups) {
    if (g.rowEl) g.rowEl.remove();
  }
  groups = [];
  groupIdCounter = 0;
}

function resetGroups() {
  clearAllGroupRows();
  addGroupRow(makeGroup());
  updateFetchCounter();
}

function loadGroupsFromPayload(payload) {
  clearAllGroupRows();
  for (const gp of payload) {
    const group = makeGroup();
    if (gp && typeof gp.name === "string" && gp.name.trim()) group.name = gp.name.trim();
    if (gp && typeof gp.color === "string" && SAFE_CSS_COLOR.test(gp.color.trim())) group.color = gp.color.trim();
    addGroupRow(group);
    if (gp && Array.isArray(gp.pathways)) {
      for (const name of gp.pathways) {
        const entry = pathwayByName.get(name);
        if (entry && !allChippedNames().has(name)) group.pathways.push(entry);
      }
    }
    renderChips(group);
  }
  if (groups.length === 0) {
    resetGroups();
  } else {
    updateRemoveButtonsVisibility();
    updateFetchCounter();
  }
}

function allChippedNames() {
  const names = new Set();
  for (const g of groups) {
    for (const p of g.pathways) names.add(p.pathway_name);
  }
  return names;
}

function flattenQuerySet() {
  const seen = new Set();
  const result = [];
  for (const g of groups) {
    for (const p of g.pathways) {
      if (!seen.has(p.pathway_name)) {
        seen.add(p.pathway_name);
        result.push({ id: p.pathway_id, name: p.pathway_name });
      }
    }
  }
  return result;
}

function updateRemoveButtonsVisibility() {
  const hide = groups.length <= 1;
  for (const g of groups) {
    const btn = g.rowEl.querySelector(".remove-group-btn");
    btn.hidden = hide;
  }
}

function selectedTissueIds() {
  return selectedTissues.slice();
}

function defaultTissueIds() {
  const ids = (tissuesIndexData && tissuesIndexData.tissues ? tissuesIndexData.tissues : []).map((t) => t.tissue_id);
  if (ids.includes("combined")) return ["combined"];
  return ids.slice(0, 1);
}

function referenceTissueId(ids = selectedTissueIds()) {
  if (ids.includes("combined")) return "combined";
  return ids[0] || null;
}

function orderedTissueIds(ids) {
  const order = (tissuesIndexData && tissuesIndexData.tissues ? tissuesIndexData.tissues : []).map((t) => t.tissue_id);
  const wanted = new Set(ids);
  return order.filter((id) => wanted.has(id));
}

function tissueSummaryText(ids = selectedTissueIds()) {
  if (ids.length === 0) return "Select tissues";
  if (ids.length <= 2) return ids.join(", ");
  return `${ids.length} tissues`;
}

function updateTissueSummary() {
  tissueSummary.textContent = tissueSummaryText();
}

function setTissueChecks(ids) {
  selectedTissues = orderedTissueIds(ids);
  updateTissueSummary();
  syncTissueMenuCheckboxes();
}

function parseUrlTissues() {
  const params = new URLSearchParams(window.location.search);
  const valid = new Set((tissuesIndexData.tissues || []).map((t) => t.tissue_id));
  const fromList = (params.get("tissues") || "")
    .split(",")
    .map((s) => s.trim())
    .filter((id) => valid.has(id));
  if (fromList.length) return fromList;
  const one = params.get("tissue");
  if (one && valid.has(one)) return [one];
  return defaultTissueIds();
}

function plannedFetchCount() {
  const seedCount = flattenQuerySet().length;
  const topN = Math.max(0, parseInt(topNInput.value, 10) || 0);
  const tissueCount = Math.max(1, selectedTissueIds().length);
  return {
    seedCount,
    topN,
    tissueCount,
    perTissue: seedCount + topN,
    total: (seedCount + topN) * tissueCount,
  };
}

function fetchWouldExceedCap() {
  const plan = plannedFetchCount();
  return plan.seedCount > 0 && plan.total > MAX_TOTAL_FETCHES;
}

function updateFetchCounter() {
  const plan = plannedFetchCount();
  const over = fetchWouldExceedCap();
  const maxTissuesOneSeed = Math.max(1, Math.floor(MAX_TOTAL_FETCHES / Math.max(1, 1 + plan.topN)));
  if (plan.seedCount === 0) {
    if (plan.tissueCount > maxTissuesOneSeed) {
      fetchCounterEl.textContent = `At Top N ${plan.topN}, at most ${maxTissuesOneSeed} tissue(s) fit the ${MAX_TOTAL_FETCHES}-file limit with one pathway. ${plan.tissueCount} selected.`;
      fetchCounterEl.classList.add("warn");
    } else {
      fetchCounterEl.textContent = "";
      fetchCounterEl.classList.remove("warn");
    }
  } else if (over) {
    fetchCounterEl.textContent = `${plan.seedCount} pathway(s) × ${plan.tissueCount} tissue(s) would fetch ${plan.total} files, above the ${MAX_TOTAL_FETCHES}-file limit. Reduce tissues, pathways, or Top N.`;
    fetchCounterEl.classList.add("warn");
  } else {
    fetchCounterEl.textContent = `${plan.seedCount} pathway(s) across ${groups.length} group(s) → ${plan.perTissue} file(s) × ${plan.tissueCount} tissue(s) = ${plan.total} (limit ${MAX_TOTAL_FETCHES}).`;
    fetchCounterEl.classList.remove("warn");
  }
  if (!queryBusy) goButton.disabled = over;
}

function removeGroup(id) {
  if (groups.length <= 1) return;
  const idx = groups.findIndex((g) => g.id === id);
  if (idx === -1) return;
  const [group] = groups.splice(idx, 1);
  group.rowEl.remove();
  updateRemoveButtonsVisibility();
  updateFetchCounter();
}

function renderChips(group) {
  const chipList = group.rowEl.querySelector(".chip-list");
  chipList.innerHTML = "";
  for (const p of group.pathways) {
    const chip = document.createElement("span");
    chip.className = "chip";
    const label = document.createElement("span");
    label.className = "chip-label";
    label.textContent = p.pathway_name;
    label.title = p.pathway_name;
    chip.appendChild(label);
    const removeBtn = document.createElement("button");
    removeBtn.type = "button";
    removeBtn.className = "chip-remove";
    removeBtn.setAttribute("aria-label", `Remove ${p.pathway_name}`);
    removeBtn.textContent = "×";
    removeBtn.addEventListener("click", () => {
      group.pathways = group.pathways.filter((x) => x.pathway_name !== p.pathway_name);
      renderChips(group);
      syncMenuCheckbox(group, p.pathway_name, false);
      syncBulkCheckbox(p.pathway_name, false);
      syncFromGenesCheckbox(p.pathway_name, false);
      updateFetchCounter();
    });
    chip.appendChild(removeBtn);
    chipList.appendChild(chip);
  }
  updateFetchCounter();
}

function addGroupRow(group) {
  const row = document.createElement("div");
  row.className = "group-row";
  row.dataset.groupId = group.id;
  row.innerHTML = `
    <div class="group-row-header">
      <span class="group-swatch" style="background:${cssColor(group.color)}"></span>
      <input type="text" class="group-name-input" value="${escapeHtml(group.name)}" aria-label="Group name">
      <button type="button" class="remove-group-btn text-btn" aria-label="Remove group">×</button>
    </div>
    <div class="chip-list"></div>
    <div class="combobox">
      <input type="text" class="group-pathway-input" placeholder="Add a pathway…" autocomplete="off" spellcheck="false" role="combobox" aria-autocomplete="list" aria-controls="menu-${group.id}" aria-expanded="false">
      <ul id="menu-${group.id}" class="combobox-menu" role="listbox" hidden></ul>
    </div>
    <div class="group-lookup">
      <p class="group-lookup-label">Bulk lookup</p>
      <div class="group-lookup-actions">
        <button type="button" class="bulk-lookup-btn text-btn">From pathway</button>
        <button type="button" class="from-genes-btn text-btn">From genes</button>
      </div>
    </div>
  `;
  groupsContainer.appendChild(row);
  group.rowEl = row;

  const nameInput = row.querySelector(".group-name-input");
  nameInput.addEventListener("input", () => {
    group.name = nameInput.value;
    if (bulkLookupDialog.open && bulkLookupGroup === group) {
      bulkLookupTitle.textContent = `From pathway — ${group.name}`;
    }
    if (fromGenesDialog.open && fromGenesGroup === group) {
      fromGenesTitle.textContent = `Pathway lookup from genes — ${group.name}`;
    }
  });

  const removeBtn = row.querySelector(".remove-group-btn");
  removeBtn.addEventListener("click", () => removeGroup(group.id));

  const inputEl = row.querySelector(".group-pathway-input");
  const menuEl = row.querySelector(".combobox-menu");
  const containerEl = row.querySelector(".combobox");
  group.combobox = { inputEl, menuEl, containerEl, items: [], index: -1 };
  attachGroupComboboxEvents(group);

  row.querySelector(".bulk-lookup-btn").addEventListener("click", () => openBulkLookup(group));
  row.querySelector(".from-genes-btn").addEventListener("click", () => openFromGenes(group));

  renderChips(group);
  updateRemoveButtonsVisibility();
}

// ---------- Small helpers ----------
// Additive, not overwrite: a query can fail in more than one independent
// way in the same run (e.g. a neighbor fetch 404 AND a render library
// missing), and both should stay visible rather than the later one
// clobbering the earlier one.
function showError(message) {
  const current = errorBannerText.textContent;
  errorBannerText.textContent = current ? `${current}\n${message}` : message;
  errorBanner.hidden = false;
}
function clearError() {
  errorBanner.hidden = true;
  errorBannerText.textContent = "";
}
function setStatus(message) {
  queryStatus.textContent = message || "";
}
let queryBusy = false;
function setBusy(isBusy) {
  queryBusy = isBusy;
  goButton.disabled = isBusy || fetchWouldExceedCap();
  goButton.textContent = isBusy ? "…" : "Build network";
}
function showResults(visible) {
  emptyState.hidden = visible;
  results.hidden = !visible;
}

async function fetchJSON(path, options) {
  const url = `${DATA_BASE_URL}/${path}`;
  let res;
  try {
    res = await fetch(url, options);
  } catch (networkErr) {
    throw new Error(`Network error fetching ${url}: ${networkErr.message}`);
  }
  if (!res.ok) {
    throw new Error(`Fetch failed (${res.status} ${res.statusText}): ${url}`);
  }
  return res.json();
}

function pathwaySliceSlugs(pathwayId) {
  const slugs = [];
  const add = (s) => {
    if (s && !slugs.includes(s)) slugs.push(s);
  };
  add(pathwayId);
  if (pathwayId && pathwayId.startsWith("Pathway.")) {
    add(pathwayId.slice("Pathway.".length));
  } else if (pathwayId) {
    add(`Pathway.${pathwayId}`);
  }
  return slugs;
}

async function fetchPathwaySlice(tissueId, pathwayId) {
  let lastErr;
  for (const slug of pathwaySliceSlugs(pathwayId)) {
    try {
      return await fetchJSON(`${tissueId}/pathways/${slug}.json`);
    } catch (err) {
      lastErr = err;
    }
  }
  throw lastErr;
}

function passesPCutoff(edge, pCutOff) {
  if (pvalueMode === "neg_log10") {
    return edge[pColumn] > -Math.log10(pCutOff);
  }
  return edge[pColumn] < pCutOff;
}

// Render-friendly metric display: convert neg_log10 storage back to an
// ordinary p-value-looking string, since that's what a human expects to read
// even though the on-disk/filtering representation is log-transformed.
function formatMetric(value) {
  if (value === null || value === undefined) return "";
  if (pvalueMode === "neg_log10") {
    return Math.pow(10, -value).toExponential(2);
  }
  return typeof value === "number" ? value.toExponential(2) : String(value);
}
function formatCorrelation(value) {
  return typeof value === "number" ? value.toFixed(3) : String(value);
}

function routingParams() {
  const current = new URLSearchParams(window.location.search);
  const keep = new URLSearchParams();
  for (const key of ROUTING_PARAMS) {
    if (current.has(key)) keep.set(key, current.get(key));
  }
  return keep;
}

function setHomeLink() {
  const params = routingParams();
  const qs = params.toString();
  const href = qs ? `index.html?${qs}` : "index.html";
  homeLink.href = href;
  if (homeBrand) homeBrand.href = href;
}

// Pathway names in this dataset can contain almost any punctuation (colons,
// commas, parens - see the design doc), so groups/their pathway lists are
// round-tripped as one JSON blob rather than a hand-rolled delimiter scheme.
// URLSearchParams handles percent-encoding on write/decode on read.
function writeUrlState() {
  const params = routingParams();
  const tissues = selectedTissueIds();
  if (tissues.length) params.set("tissues", tissues.join(","));
  if (collectionSelect.value) params.set("collection", collectionSelect.value);
  const groupsPayload = groups
    .filter((g) => g.pathways.length > 0)
    .map((g) => ({ name: g.name, color: g.color, pathways: g.pathways.map((p) => p.pathway_name) }));
  if (groupsPayload.length > 0) {
    params.set("groups", JSON.stringify(groupsPayload));
  }
  params.set("topN", topNInput.value);
  params.set("corr", corrCutoffInput.value);
  params.set("p", pCutoffInput.value);
  const next = `${window.location.pathname}?${params.toString()}`;
  history.replaceState(null, "", next);
}

function applyUrlFields() {
  const params = new URLSearchParams(window.location.search);
  if (params.get("collection") && collectionSelect.querySelector(`option[value="${CSS.escape(params.get("collection"))}"]`)) {
    collectionSelect.value = params.get("collection");
  }
  if (params.get("topN")) topNInput.value = params.get("topN");
  if (params.get("corr")) corrCutoffInput.value = params.get("corr");
  if (params.get("p")) pCutoffInput.value = params.get("p");

  const groupsParam = params.get("groups");
  if (groupsParam) {
    try {
      const parsed = JSON.parse(groupsParam);
      if (Array.isArray(parsed) && parsed.length > 0) {
        loadGroupsFromPayload(parsed);
      }
    } catch (err) {
      console.error("Malformed ?groups= URL parameter, ignoring:", err);
    }
  }
}

// ---------- Bootstrap: load tissues, then pathway index ----------
async function init() {
  try {
    tissuesIndexData = await fetchJSON("tissues_index.json", { cache: "no-store" });
  } catch (err) {
    showError(`Could not load tissues_index.json from ${DATA_BASE_URL}. ${err.message}`);
    return;
  }
  if (!tissuesIndexData.tissues || tissuesIndexData.tissues.length === 0) {
    showError("tissues_index.json has no tissues.");
    return;
  }
  const initialTissues = parseUrlTissues();
  populateTissueList(initialTissues);
  const initialTissue = referenceTissueId(initialTissues);
  await loadTissue(initialTissue);
  applyUrlFields();
  if (flattenQuerySet().length > 0) {
    await onGo();
  }
}

function allTissueEntries() {
  const list = tissuesIndexData && tissuesIndexData.tissues ? tissuesIndexData.tissues.slice() : [];
  list.sort((a, b) => {
    if (a.tissue_id === "combined") return -1;
    if (b.tissue_id === "combined") return 1;
    return a.tissue_id.localeCompare(b.tissue_id);
  });
  return list;
}

function filteredTissues(query) {
  const needle = query.trim().toLowerCase();
  const list = allTissueEntries();
  if (!needle) return list;
  return list.filter((t) => t.tissue_id.toLowerCase().includes(needle));
}

function syncTissueMenuCheckboxes() {
  const selected = new Set(selectedTissues);
  tissueList.querySelectorAll("li[data-tissue]").forEach((li) => {
    const box = li.querySelector("input[type=checkbox]");
    if (box) box.checked = selected.has(li.dataset.tissue);
  });
}

function highlightTissueItem() {
  tissueMenuItems.forEach((item, i) => {
    item.setAttribute("aria-selected", i === tissueMenuIndex ? "true" : "false");
  });
  const active = tissueMenuItems[tissueMenuIndex];
  if (active) {
    active.scrollIntoView({ block: "nearest" });
    tissueFilter.setAttribute("aria-activedescendant", active.id);
  } else {
    tissueFilter.removeAttribute("aria-activedescendant");
  }
}

function renderTissueMenu() {
  const matches = filteredTissues(tissueFilter.value);
  tissueList.innerHTML = "";
  tissueMenuItems = [];
  tissueMenuIndex = -1;
  tissueFilter.removeAttribute("aria-activedescendant");
  if (matches.length === 0) {
    const empty = document.createElement("li");
    empty.className = "combo-empty";
    empty.textContent = "No matching tissues";
    tissueList.appendChild(empty);
    return;
  }
  const selected = new Set(selectedTissues);
  const frag = document.createDocumentFragment();
  for (const t of matches) {
    const li = document.createElement("li");
    li.id = `tissue-option-${t.tissue_id}`;
    li.setAttribute("role", "option");
    li.dataset.tissue = t.tissue_id;
    const checked = selected.has(t.tissue_id);
    li.innerHTML = `<div class="combo-option"><input type="checkbox" tabindex="-1"${
      checked ? " checked" : ""
    }><span class="combo-text"><span class="combo-name">${escapeHtml(t.tissue_id)}</span></span></div>`;
    li.addEventListener("mousedown", (evt) => {
      evt.preventDefault();
      toggleTissue(t.tissue_id);
    });
    frag.appendChild(li);
    tissueMenuItems.push(li);
  }
  tissueList.appendChild(frag);
}

function openTissueMenu() {
  if (!tissuesIndexData) return;
  tissueMenu.hidden = false;
  tissueSummary.setAttribute("aria-expanded", "true");
  tissueFilter.setAttribute("aria-expanded", "true");
  renderTissueMenu();
  tissueFilter.focus();
}

function closeTissueMenu() {
  if (selectedTissues.length === 0) {
    setTissueChecks(defaultTissueIds());
    updateFetchCounter();
    writeUrlState();
    const ref = referenceTissueId();
    if (ref && ref !== currentTissueId) loadTissue(ref);
  }
  tissueMenu.hidden = true;
  tissueSummary.setAttribute("aria-expanded", "false");
  tissueFilter.setAttribute("aria-expanded", "false");
  tissueFilter.removeAttribute("aria-activedescendant");
  tissueFilter.value = "";
  tissueList.innerHTML = "";
  tissueMenuItems = [];
  tissueMenuIndex = -1;
}

function populateTissueList(selectedIds) {
  setTissueChecks(selectedIds);
  updateFetchCounter();
}

async function toggleTissue(tissueId) {
  if (selectedTissues.includes(tissueId)) {
    setTissueChecks(selectedTissues.filter((id) => id !== tissueId));
  } else {
    setTissueChecks([...selectedTissues, tissueId]);
  }
  await onTissueChange();
}

async function onTissueChange() {
  const ids = selectedTissueIds();
  updateFetchCounter();
  writeUrlState();
  const ref = referenceTissueId(ids);
  if (ref && ref !== currentTissueId) {
    await loadTissue(ref);
  }
}

async function loadTissue(tissueId) {
  clearError();
  setStatus("Loading pathway index…");
  const manifestEntry = tissuesIndexData.tissues.find((t) => t.tissue_id === tissueId);
  if (!manifestEntry) {
    showError(`Unknown tissue: ${tissueId}`);
    return;
  }
  try {
    pathwayIndexData = await fetchJSON(manifestEntry.pathway_index, { cache: "no-store" });
  } catch (err) {
    showError(`Could not load pathway_index.json for tissue '${tissueId}'. ${err.message}`);
    return;
  }
  currentTissueId = tissueId;
  pColumn = pathwayIndexData.columns.p_value;
  fdrColumn = pathwayIndexData.columns.fdr;
  pvalueMode = pathwayIndexData.pvalue_mode;

  pathwayByName = new Map();
  for (const p of pathwayIndexData.pathways) {
    pathwayByName.set(p.pathway_name, p);
    if (p.pathway_name.startsWith("Pathway.")) {
      pathwayByName.set(p.pathway_name.slice("Pathway.".length), p);
    } else {
      pathwayByName.set(`Pathway.${p.pathway_name}`, p);
    }
  }

  populateCollectionSelect();
  closeAllGroupMenus();
  indexMeta.textContent = "";
  updateEmptyCatalog();
  setStatus(`${pathwayIndexData.pathways.length} pathways loaded for '${tissueId}'.`);
}

const SOURCE_DB_LABELS = {
  REACTOME: "Reactome",
  WP: "WikiPathways",
  KEGG: "KEGG",
  BIOCARTA: "BioCarta",
  PID: "PID",
  SA: "SA",
  SIG: "Sig",
  WNT: "Wnt",
  NABA: "NABA",
};

function sourceDbLabel(id) {
  return SOURCE_DB_LABELS[id] || id;
}

function pathwayMetaLabel(p) {
  return p && p.source_db ? sourceDbLabel(p.source_db) : "";
}

function sourceDbSummary() {
  const counts = {};
  for (const p of pathwayIndexData.pathways || []) {
    const id = p.source_db;
    if (!id) continue;
    counts[id] = (counts[id] || 0) + 1;
  }
  return Object.entries(counts).sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
}

function pathwayMatchesCollection(p) {
  const selected = collectionSelect.value;
  if (!selected) return true;
  return p.source_db === selected;
}

function tissueCellTypeCount() {
  const list = tissuesIndexData && tissuesIndexData.tissues ? tissuesIndexData.tissues : [];
  return list.filter((t) => t.tissue_id !== "combined").length;
}

function updateEmptyCatalog() {
  const catalog = el("emptyCatalog");
  if (!catalog) return;
  const nPathways = pathwayIndexData && pathwayIndexData.pathways ? pathwayIndexData.pathways.length : 0;
  const nTissues = tissueCellTypeCount();
  const nCollections = sourceDbSummary().length;
  if (!nPathways) {
    catalog.textContent = "";
    return;
  }
  const pathwayBit = `${nPathways.toLocaleString()} pathway${nPathways === 1 ? "" : "s"}`;
  const tissueBit = `${nTissues.toLocaleString()} tissues-celltypes`;
  const collectionBit = `${nCollections.toLocaleString()} collection${nCollections === 1 ? "" : "s"}`;
  catalog.textContent = `${pathwayBit} across ${tissueBit} in ${collectionBit}`;
}

function populateCollectionSelect() {
  const previous = collectionSelect.value;
  collectionSelect.innerHTML = '<option value="">All</option>';
  for (const [id, count] of sourceDbSummary()) {
    const opt = document.createElement("option");
    opt.value = id;
    opt.textContent = `${sourceDbLabel(id)} (${count})`;
    collectionSelect.appendChild(opt);
  }
  if (previous && collectionSelect.querySelector(`option[value="${CSS.escape(previous)}"]`)) {
    collectionSelect.value = previous;
  }
}

function findPathwayByName(name) {
  const needle = name.trim().toLowerCase();
  const aliases = [needle];
  if (needle.startsWith("pathway.")) aliases.push(needle.slice("pathway.".length));
  else aliases.push(`pathway.${needle}`);
  for (const p of pathwayIndexData.pathways) {
    const pname = p.pathway_name.toLowerCase();
    if (aliases.includes(pname)) return p;
  }
  return null;
}

function namesInOtherGroups(group) {
  const names = new Set();
  for (const g of groups) {
    if (g === group) continue;
    for (const p of g.pathways) names.add(p.pathway_name);
  }
  return names;
}

function groupHasPathway(group, name) {
  return group.pathways.some((p) => p.pathway_name === name);
}

function filteredPathways(query, excludeNames) {
  const q = (query || "").trim().toLowerCase();
  const matches = [];
  for (const p of pathwayIndexData.pathways || []) {
    if (!pathwayMatchesCollection(p)) continue;
    if (excludeNames.has(p.pathway_name)) continue;
    if (q && !p.pathway_name.toLowerCase().includes(q)) continue;
    matches.push(p);
  }
  return matches;
}

function syncMenuCheckbox(group, name, checked) {
  const cb = group.combobox;
  if (!cb || cb.menuEl.hidden) return;
  const li = cb.items.find((item) => item.dataset.name === name);
  if (!li) return;
  const box = li.querySelector("input[type=checkbox]");
  if (box) box.checked = checked;
}

// ---------- Per-group combobox (generalizes what used to be one global combobox) ----------
function openGroupMenu(group) {
  if (!pathwayIndexData) return;
  const cb = group.combobox;
  const excluded = namesInOtherGroups(group);
  const selected = new Set(group.pathways.map((p) => p.pathway_name));
  const matches = filteredPathways(cb.inputEl.value, excluded);
  cb.menuEl.innerHTML = "";
  cb.items = [];
  cb.index = -1;
  if (matches.length === 0) {
    const empty = document.createElement("li");
    empty.className = "combo-empty";
    empty.textContent = "No matching pathways";
    cb.menuEl.appendChild(empty);
  } else {
    const frag = document.createDocumentFragment();
    for (const p of matches) {
      const li = document.createElement("li");
      li.setAttribute("role", "option");
      li.dataset.name = p.pathway_name;
      const checked = selected.has(p.pathway_name);
      li.innerHTML = `<div class="combo-option"><input type="checkbox" tabindex="-1"${
        checked ? " checked" : ""
      }><span class="combo-text"><span class="combo-name">${escapeHtml(p.pathway_name)}</span><span class="combo-meta">${escapeHtml(
        pathwayMetaLabel(p)
      )}</span></span></div>`;
      li.addEventListener("mousedown", (evt) => {
        evt.preventDefault();
        togglePathwayInGroup(group, p.pathway_name);
      });
      frag.appendChild(li);
      cb.items.push(li);
    }
    cb.menuEl.appendChild(frag);
  }
  cb.menuEl.hidden = false;
  cb.inputEl.setAttribute("aria-expanded", "true");
}

function closeGroupMenu(group) {
  const cb = group.combobox;
  cb.menuEl.hidden = true;
  cb.menuEl.innerHTML = "";
  cb.items = [];
  cb.index = -1;
  cb.inputEl.setAttribute("aria-expanded", "false");
}

function closeAllGroupMenus() {
  for (const g of groups) closeGroupMenu(g);
}

function highlightGroupItem(group) {
  const cb = group.combobox;
  cb.items.forEach((item, i) => item.setAttribute("aria-selected", i === cb.index ? "true" : "false"));
  const active = cb.items[cb.index];
  if (active) active.scrollIntoView({ block: "nearest" });
}

function togglePathwayInGroup(group, name) {
  const entry = pathwayByName.get(name);
  if (!entry) return;
  if (namesInOtherGroups(group).has(name)) return;
  if (groupHasPathway(group, name)) {
    group.pathways = group.pathways.filter((p) => p.pathway_name !== name);
    renderChips(group);
    syncMenuCheckbox(group, name, false);
    syncBulkCheckbox(name, false);
    syncFromGenesCheckbox(name, false);
  } else {
    group.pathways.push(entry);
    renderChips(group);
    syncMenuCheckbox(group, name, true);
    syncBulkCheckbox(name, true);
    syncFromGenesCheckbox(name, true);
  }
}

function attachGroupComboboxEvents(group) {
  const { inputEl } = group.combobox;
  inputEl.addEventListener("focus", () => openGroupMenu(group));
  inputEl.addEventListener("click", () => openGroupMenu(group));
  inputEl.addEventListener("input", () => openGroupMenu(group));
  inputEl.addEventListener("keydown", (evt) => {
    const cb = group.combobox;
    if (evt.key === "Backspace" && inputEl.value === "" && group.pathways.length > 0) {
      evt.preventDefault();
      const removed = group.pathways[group.pathways.length - 1];
      group.pathways = group.pathways.slice(0, -1);
      renderChips(group);
      if (removed) syncMenuCheckbox(group, removed.pathway_name, false);
      return;
    }
    if (cb.menuEl.hidden && (evt.key === "ArrowDown" || evt.key === "ArrowUp")) {
      openGroupMenu(group);
    }
    if (cb.menuEl.hidden) return;
    if (evt.key === "ArrowDown") {
      evt.preventDefault();
      cb.index = Math.min(cb.items.length - 1, cb.index + 1);
      highlightGroupItem(group);
    } else if (evt.key === "ArrowUp") {
      evt.preventDefault();
      cb.index = Math.max(0, cb.index - 1);
      highlightGroupItem(group);
    } else if (evt.key === "Enter" && cb.index >= 0) {
      evt.preventDefault();
      togglePathwayInGroup(group, cb.items[cb.index].dataset.name);
    } else if (evt.key === "Escape") {
      evt.preventDefault();
      closeGroupMenu(group);
    }
  });
}

// ---------- Bulk pathway lookup (per-group dialog over the same index) ----------
function parseBulkTerms(text) {
  const seen = new Set();
  const terms = [];
  for (const raw of (text || "").split(/[\n,;]+/)) {
    const term = raw.trim().toLowerCase();
    if (!term || seen.has(term)) continue;
    seen.add(term);
    terms.push(term);
  }
  return terms;
}

function bulkLookupMode() {
  const checked = document.querySelector('input[name="bulkLookupMode"]:checked');
  return checked && checked.value === "and" ? "and" : "or";
}

function bulkFilteredPathways(terms, mode, excludeNames) {
  const matches = [];
  for (const p of pathwayIndexData.pathways || []) {
    if (!pathwayMatchesCollection(p)) continue;
    if (excludeNames.has(p.pathway_name)) continue;
    const name = p.pathway_name.toLowerCase();
    const hit = mode === "and" ? terms.every((t) => name.includes(t)) : terms.some((t) => name.includes(t));
    if (hit) matches.push(p);
  }
  return matches;
}

function syncBulkCheckbox(name, checked) {
  if (!bulkLookupDialog.open) return;
  const li = bulkLookupItems.find((item) => item.dataset.name === name);
  if (!li) return;
  const box = li.querySelector("input[type=checkbox]");
  if (box) box.checked = checked;
}

function renderBulkLookupResults() {
  bulkLookupList.innerHTML = "";
  bulkLookupItems = [];
  const terms = parseBulkTerms(bulkLookupTerms.value);
  if (!bulkLookupGroup || !pathwayIndexData || terms.length === 0) {
    bulkLookupCount.textContent = "";
    const empty = document.createElement("li");
    empty.className = "bulk-lookup-empty";
    empty.textContent = "Paste one or more terms to search.";
    bulkLookupList.appendChild(empty);
    return;
  }
  const selected = new Set(bulkLookupGroup.pathways.map((p) => p.pathway_name));
  const matches = bulkFilteredPathways(terms, bulkLookupMode(), namesInOtherGroups(bulkLookupGroup));
  const n = matches.length;
  bulkLookupCount.textContent = n === 1 ? "1 matching pathway" : `${n} matching pathways`;
  if (n === 0) {
    const empty = document.createElement("li");
    empty.className = "bulk-lookup-empty";
    empty.textContent = "No matching pathways";
    bulkLookupList.appendChild(empty);
    return;
  }
  const frag = document.createDocumentFragment();
  for (const p of matches) {
    const li = document.createElement("li");
    li.setAttribute("role", "option");
    li.dataset.name = p.pathway_name;
    const checked = selected.has(p.pathway_name);
    li.innerHTML = `<div class="combo-option"><input type="checkbox" tabindex="-1"${
      checked ? " checked" : ""
    }><span class="combo-text"><span class="combo-name">${escapeHtml(p.pathway_name)}</span><span class="combo-meta">${escapeHtml(
      pathwayMetaLabel(p)
    )}</span></span></div>`;
    li.addEventListener("mousedown", (evt) => {
      evt.preventDefault();
      togglePathwayInGroup(bulkLookupGroup, p.pathway_name);
    });
    frag.appendChild(li);
    bulkLookupItems.push(li);
  }
  bulkLookupList.appendChild(frag);
}

function openBulkLookup(group) {
  if (!pathwayIndexData) return;
  closeAllGroupMenus();
  bulkLookupGroup = group;
  bulkLookupTitle.textContent = `From pathway — ${group.name}`;
  bulkLookupTerms.value = "";
  const orRadio = document.querySelector('input[name="bulkLookupMode"][value="or"]');
  if (orRadio) orRadio.checked = true;
  renderBulkLookupResults();
  bulkLookupDialog.showModal();
  bulkLookupTerms.focus();
}

function setVisibleBulkSelection(checked) {
  if (!bulkLookupGroup) return;
  const other = namesInOtherGroups(bulkLookupGroup);
  for (const li of bulkLookupItems) {
    const name = li.dataset.name;
    const has = groupHasPathway(bulkLookupGroup, name);
    if (checked && !has) {
      const entry = pathwayByName.get(name);
      if (entry && !other.has(name)) bulkLookupGroup.pathways.push(entry);
    } else if (!checked && has) {
      bulkLookupGroup.pathways = bulkLookupGroup.pathways.filter((p) => p.pathway_name !== name);
    }
  }
  renderChips(bulkLookupGroup);
  renderBulkLookupResults();
}

// ---------- From genes: ORA over the current collection, then add as seeds ----------
const logFactCache = [0];
function logFact(n) {
  if (n < 0) return Number.NEGATIVE_INFINITY;
  for (let i = logFactCache.length; i <= n; i++) {
    logFactCache[i] = logFactCache[i - 1] + Math.log(i);
  }
  return logFactCache[n];
}
function logChoose(n, k) {
  if (k < 0 || k > n) return Number.NEGATIVE_INFINITY;
  return logFact(n) - logFact(k) - logFact(n - k);
}
function logAddExp(a, b) {
  if (a === Number.NEGATIVE_INFINITY) return b;
  if (b === Number.NEGATIVE_INFINITY) return a;
  const m = Math.max(a, b);
  return m + Math.log(Math.exp(a - m) + Math.exp(b - m));
}
function hypergeomUpper(k, n, K, N) {
  if (k <= 0) return 1;
  if (n <= 0 || N <= 0 || k > n || k > K) return 0;
  const logDen = logChoose(N, n);
  if (!Number.isFinite(logDen)) return 0;
  let logSum = Number.NEGATIVE_INFINITY;
  const maxI = Math.min(n, K);
  for (let i = k; i <= maxI; i++) {
    const logTerm = logChoose(K, i) + logChoose(N - K, n - i) - logDen;
    if (Number.isFinite(logTerm)) logSum = logAddExp(logSum, logTerm);
  }
  if (logSum === Number.NEGATIVE_INFINITY) return 0;
  const p = Math.exp(logSum);
  return p > 1 ? 1 : p;
}
function bhAdjust(pvals) {
  const m = pvals.length;
  const q = new Array(m);
  if (m === 0) return q;
  const order = pvals.map((_, i) => i).sort((a, b) => pvals[a] - pvals[b]);
  let minQ = 1;
  for (let rank = m; rank >= 1; rank -= 1) {
    const i = order[rank - 1];
    minQ = Math.min(minQ, (pvals[i] * m) / rank);
    q[i] = minQ;
  }
  return q;
}
function parseGeneTokens(text) {
  const tokens = [];
  const seen = new Set();
  for (const raw of (text || "").split(/[\n,;]+/)) {
    const token = raw.trim();
    if (!token) continue;
    const key = token.toUpperCase();
    if (seen.has(key)) continue;
    seen.add(key);
    tokens.push(token);
  }
  return tokens;
}
function formatQValue(value) {
  if (!Number.isFinite(value)) return "";
  if (value === 0) return "0";
  return value.toExponential(2);
}

async function ensureGeneMembership() {
  if (geneMembership) return geneMembership;
  if (geneMembershipError) throw geneMembershipError;
  try {
    geneMembership = await fetchJSON("gene_membership.json");
    return geneMembership;
  } catch (err) {
    geneMembershipError = err;
    throw err;
  }
}

function geneRecordSymbol(g) {
  if (!g) return "";
  if (typeof g === "string") return g;
  return g.symbol || "";
}

function geneRecordEntrez(g) {
  if (!g || typeof g === "string") return "";
  return g.entrez || "";
}

function geneRecordDescription(g) {
  if (!g || typeof g === "string") return "";
  return g.description || "";
}

function membershipRecord(pathwayName) {
  if (!geneMembership) return null;
  const map = geneMembership.pathways || {};
  if (map[pathwayName]) return map[pathwayName];
  if (pathwayName.startsWith("Pathway.")) return map[pathwayName.slice("Pathway.".length)] || null;
  return map[`Pathway.${pathwayName}`] || null;
}

function pathwayGeneRecords(pathwayName) {
  const rec = membershipRecord(pathwayName);
  return rec && Array.isArray(rec.genes) ? rec.genes : [];
}

function pathwayGeneSymbols(pathwayName) {
  const seen = new Set();
  const symbols = [];
  for (const g of pathwayGeneRecords(pathwayName)) {
    const symbol = geneRecordSymbol(g);
    if (!symbol) continue;
    const key = symbol.toUpperCase();
    if (seen.has(key)) continue;
    seen.add(key);
    symbols.push(symbol);
  }
  return symbols;
}

function visiblePathwayNames() {
  const names = new Set();
  for (const e of currentlyVisibleEdges()) {
    names.add(e.pathway_a);
    names.add(e.pathway_b);
  }
  return Array.from(names).sort((a, b) => a.localeCompare(b));
}

async function copyText(text, button, label) {
  try {
    await navigator.clipboard.writeText(text);
    if (button) {
      button.textContent = "Copied";
      setTimeout(() => {
        button.textContent = label;
      }, 1400);
    }
  } catch {
    showError("Could not copy automatically.");
  }
}

let genesDialogPathway = "";

async function openPathwayGenes(pathwayName) {
  genesDialogPathway = pathwayName;
  genesDialogTitle.textContent = pathwayName;
  genesDialogMeta.textContent = "Loading genes…";
  genesDialogList.innerHTML = "";
  copyPathwayGenesButton.hidden = true;
  genesDialog.showModal();
  try {
    await ensureGeneMembership();
  } catch (err) {
    genesDialogMeta.textContent = `Could not load gene membership. ${err.message}`;
    return;
  }
  renderPathwayGenesDialog(pathwayName);
}

function renderPathwayGenesDialog(pathwayName) {
  const records = pathwayGeneRecords(pathwayName);
  if (records.length === 0) {
    genesDialogMeta.textContent = "No gene membership for this pathway.";
    genesDialogList.innerHTML = `<p class="genes-dialog-empty">This pathway is not in the membership table.</p>`;
    copyPathwayGenesButton.hidden = true;
    return;
  }
  genesDialogMeta.textContent = `${records.length} gene${records.length === 1 ? "" : "s"}`;
  copyPathwayGenesButton.hidden = false;
  const rows = records
    .map((g) => {
      const symbol = geneRecordSymbol(g) || geneRecordEntrez(g) || "—";
      const entrez = geneRecordEntrez(g);
      const desc = geneRecordDescription(g);
      return `<tr><td>${escapeHtml(symbol)}</td><td>${escapeHtml(entrez)}</td><td>${escapeHtml(desc)}</td></tr>`;
    })
    .join("");
  genesDialogList.innerHTML = `<table><thead><tr><th>Symbol</th><th>Entrez</th><th>Description</th></tr></thead><tbody>${rows}</tbody></table>`;
}

async function copyVisiblePathways() {
  const names = visiblePathwayNames();
  if (names.length === 0) {
    showError("No pathways in the current edge table.");
    return;
  }
  await copyText(names.join("\n"), copyPathwaysButton, "Copy pathways");
}

async function copyVisiblePathwayGenes() {
  const names = visiblePathwayNames();
  if (names.length === 0) {
    showError("No pathways in the current edge table.");
    return;
  }
  try {
    await ensureGeneMembership();
  } catch (err) {
    showError(`Could not load gene membership. ${err.message}`);
    return;
  }
  const seen = new Set();
  const symbols = [];
  let missing = 0;
  for (const name of names) {
    const genes = pathwayGeneSymbols(name);
    if (genes.length === 0) missing += 1;
    for (const symbol of genes) {
      const key = symbol.toUpperCase();
      if (seen.has(key)) continue;
      seen.add(key);
      symbols.push(symbol);
    }
  }
  symbols.sort((a, b) => a.localeCompare(b));
  if (symbols.length === 0) {
    showError(
      missing === names.length
        ? "None of the pathways in the edge table have gene membership."
        : "No gene symbols to copy."
    );
    return;
  }
  await copyText(symbols.join("\n"), copyGenesButton, "Copy genes");
}

function scoreFromGenes(tokens) {
  const membershipPathways = geneMembership.pathways || {};
  const entrezToSymbol = geneMembership.entrez_to_symbol || {};
  const tested = [];
  const universe = new Set();
  let skippedNoMembership = 0;
  for (const p of pathwayIndexData.pathways || []) {
    if (!pathwayMatchesCollection(p)) continue;
    const rec = membershipPathways[p.pathway_name];
    if (!rec) {
      skippedNoMembership += 1;
      continue;
    }
    const geneSet = new Set(
      (rec.genes || []).map(geneRecordSymbol).filter(Boolean).map((g) => g.toUpperCase())
    );
    geneSet.forEach((g) => universe.add(g));
    tested.push({ p, rec, geneSet });
  }
  const recognized = [];
  const unrecognized = [];
  const seen = new Set();
  for (const token of tokens) {
    let symbol = /^\d+$/.test(token) ? entrezToSymbol[token] : token;
    if (!symbol) {
      unrecognized.push(token);
      continue;
    }
    const key = symbol.toUpperCase();
    if (!universe.has(key)) {
      unrecognized.push(token);
      continue;
    }
    if (seen.has(key)) continue;
    seen.add(key);
    recognized.push(key);
  }
  const N = universe.size;
  const n = recognized.length;
  const querySet = new Set(recognized);
  const rows = [];
  for (const t of tested) {
    let overlap = 0;
    for (const g of querySet) {
      if (t.geneSet.has(g)) overlap += 1;
    }
    const size = t.rec.size || t.geneSet.size;
    const pval = hypergeomUpper(overlap, n, size, N);
    const fold = n && N && size ? overlap / n / (size / N) : 0;
    rows.push({
      name: t.p.pathway_name,
      size,
      overlap,
      fold,
      pval,
    });
  }
  const qvals = bhAdjust(rows.map((r) => r.pval));
  rows.forEach((r, i) => {
    r.q = qvals[i];
  });
  rows.sort((a, b) => a.q - b.q || a.pval - b.pval);
  return { rows, recognized, unrecognized, universeSize: N, skippedNoMembership, tested: tested.length };
}

function syncFromGenesCheckbox(name, checked) {
  if (!fromGenesDialog.open) return;
  const row = fromGenesItems.find((item) => item.dataset.name === name);
  if (!row) return;
  const box = row.querySelector("input[type=checkbox]");
  if (box) box.checked = checked;
}

function fromGenesEmptyRow(message) {
  const tr = document.createElement("tr");
  tr.className = "from-genes-empty";
  tr.innerHTML = `<td class="from-genes-empty" colspan="6">${escapeHtml(message)}</td>`;
  fromGenesTableBody.appendChild(tr);
}

function sortFromGenesHits(hits) {
  const mul = fromGenesSortDir === "asc" ? 1 : -1;
  return [...hits].sort((a, b) => {
    if (fromGenesSortKey === "name") {
      return a.name.localeCompare(b.name) * mul;
    }
    const av = Number(a[fromGenesSortKey]);
    const bv = Number(b[fromGenesSortKey]);
    const aMissing = !Number.isFinite(av);
    const bMissing = !Number.isFinite(bv);
    if (aMissing && bMissing) return 0;
    if (aMissing) return 1;
    if (bMissing) return -1;
    if (av === bv) return 0;
    return (av < bv ? -1 : 1) * mul;
  });
}

function updateFromGenesSortHeaders() {
  fromGenesTable.querySelectorAll("th[data-sort]").forEach((th) => {
    if (th.dataset.sort === fromGenesSortKey) {
      th.setAttribute("aria-sort", fromGenesSortDir === "asc" ? "ascending" : "descending");
    } else {
      th.setAttribute("aria-sort", "none");
    }
  });
}

function onFromGenesSortClick(evt) {
  const th = evt.target.closest("th[data-sort]");
  if (!th) return;
  const key = th.dataset.sort;
  if (fromGenesSortKey === key) {
    fromGenesSortDir = fromGenesSortDir === "asc" ? "desc" : "asc";
  } else {
    fromGenesSortKey = key;
    fromGenesSortDir = "asc";
  }
  updateFromGenesSortHeaders();
  if (lastFromGenesHits.length) renderFromGenesRows(lastFromGenesHits);
}

function renderFromGenesRows(hits, fallbackNote) {
  fromGenesTableBody.innerHTML = "";
  fromGenesItems = [];
  updateFromGenesSortHeaders();
  if (hits.length === 0) {
    fromGenesEmptyRow(fallbackNote || "No pathways to show.");
    return;
  }
  const ordered = sortFromGenesHits(hits);
  const selected = new Set(fromGenesGroup.pathways.map((p) => p.pathway_name));
  const frag = document.createDocumentFragment();
  for (const hit of ordered) {
    const tr = document.createElement("tr");
    tr.dataset.name = hit.name;
    const checked = selected.has(hit.name);
    tr.innerHTML = `<td><input type="checkbox" tabindex="-1"${checked ? " checked" : ""}></td>
      <td>${escapeHtml(hit.name)}</td>
      <td>${hit.size}</td>
      <td>${hit.overlap}</td>
      <td>${Number.isFinite(hit.fold) ? hit.fold.toFixed(1) : ""}</td>
      <td>${formatQValue(hit.q)}</td>`;
    tr.addEventListener("mousedown", (evt) => {
      evt.preventDefault();
      togglePathwayInGroup(fromGenesGroup, hit.name);
    });
    frag.appendChild(tr);
    fromGenesItems.push(tr);
  }
  fromGenesTableBody.appendChild(frag);
}

async function renderFromGenesResults() {
  fromGenesTableBody.innerHTML = "";
  fromGenesItems = [];
  fromGenesCount.textContent = "";
  fromGenesStatus.textContent = "";
  const tokens = parseGeneTokens(fromGenesList.value);
  if (!fromGenesGroup) return;
  if (tokens.length === 0) {
    fromGenesEmptyRow("Paste gene symbols or Entrez IDs to search.");
    return;
  }
  try {
    await ensureGeneMembership();
  } catch (err) {
    fromGenesStatus.textContent = `Could not load gene membership. ${err.message}`;
    fromGenesEmptyRow("Gene membership file is missing.");
    return;
  }
  const qCut = Number(fromGenesQCutoff.value);
  const cutoff = Number.isFinite(qCut) ? qCut : 0.05;
  const result = scoreFromGenes(tokens);
  const collectionLabel = collectionSelect.value || "all collections";
  const skipBit = result.skippedNoMembership
    ? ` ${result.skippedNoMembership} pathway(s) skipped (no gene membership).`
    : "";
  fromGenesStatus.textContent = `${result.recognized.length} recognized gene(s), ${result.unrecognized.length} unrecognized. Tested ${result.tested} pathway(s) in ${collectionLabel} (universe ${result.universeSize}).${skipBit}`;
  if (result.unrecognized.length && result.unrecognized.length <= 12) {
    fromGenesStatus.textContent += ` Unrecognized: ${result.unrecognized.join(", ")}.`;
  } else if (result.unrecognized.length > 12) {
    fromGenesStatus.textContent += ` Unrecognized include: ${result.unrecognized.slice(0, 8).join(", ")}…`;
  }
  if (result.tested === 0) {
    fromGenesEmptyRow("No pathways in this collection have gene membership.");
    return;
  }
  if (result.recognized.length === 0) {
    fromGenesEmptyRow("None of those genes are in the current collection universe.");
    return;
  }
  const other = namesInOtherGroups(fromGenesGroup);
  const ranked = result.rows.filter((r) => !other.has(r.name));
  const significant = ranked.filter((r) => r.q < cutoff);
  let hits = significant;
  let usedFallback = false;
  if (hits.length === 0) {
    hits = ranked.slice(0, 10);
    usedFallback = true;
    fromGenesCount.textContent = `No pathways with q < ${cutoff}. Showing top 10.`;
  } else {
    fromGenesCount.textContent =
      hits.length === 1 ? "1 pathway with q below cutoff" : `${hits.length} pathways with q below cutoff`;
  }
  const prefillKey = `${tokens.join("\n")}||${collectionSelect.value}`;
  if (fromGenesPrefillKey !== prefillKey && !usedFallback) {
    fromGenesPrefillKey = prefillKey;
    for (const hit of hits.slice(0, 3)) {
      if (!groupHasPathway(fromGenesGroup, hit.name)) {
        const entry = pathwayByName.get(hit.name);
        if (entry) fromGenesGroup.pathways.push(entry);
      }
    }
    renderChips(fromGenesGroup);
  }
  lastFromGenesHits = hits;
  renderFromGenesRows(hits, usedFallback ? "No pathways to show." : "No pathways with q below cutoff.");
}

async function openFromGenes(group) {
  if (!pathwayIndexData) return;
  closeAllGroupMenus();
  if (bulkLookupDialog.open) bulkLookupDialog.close();
  fromGenesGroup = group;
  fromGenesPrefillKey = null;
  fromGenesTitle.textContent = `Pathway lookup from genes — ${group.name}`;
  fromGenesList.value = "";
  fromGenesQCutoff.value = "0.05";
  fromGenesStatus.textContent = "";
  fromGenesCount.textContent = "";
  await renderFromGenesResults();
  fromGenesDialog.showModal();
  fromGenesList.focus();
}

function setVisibleFromGenesSelection(checked) {
  if (!fromGenesGroup) return;
  const other = namesInOtherGroups(fromGenesGroup);
  for (const row of fromGenesItems) {
    const name = row.dataset.name;
    if (!name) continue;
    const has = groupHasPathway(fromGenesGroup, name);
    if (checked && !has) {
      const entry = pathwayByName.get(name);
      if (entry && !other.has(name)) fromGenesGroup.pathways.push(entry);
    } else if (!checked && has) {
      fromGenesGroup.pathways = fromGenesGroup.pathways.filter((p) => p.pathway_name !== name);
    }
  }
  renderChips(fromGenesGroup);
  renderFromGenesRows(lastFromGenesHits);
}

// ---------- Query-replication algorithm (see design doc §9's generalization) ----------
// querySet: deduped array of {id, name} - the UNION of every group's members,
// flattened. Group identity plays no role here, only in rendering (below) -
// same as the legacy app's up/down split never affected its SQL query, only
// its frontend coloring.
async function runExploreQuery(querySet, topN, corrCutOff, pCutOff, tissueId) {
  const queryNames = new Set(querySet.map((m) => m.name));

  // Stage 1: fetch every query member's own slice, then rank candidates
  // OUTSIDE the query set by max(abs(correlation)) across whichever members
  // they're connected to. A single-branch scan over each member's own
  // (already bidirectional) slice correctly reproduces the legacy SQL's
  // two-branch UNION...MAX(abs(correlation)) query - see
  // scripts/prepare_pcxn_dataset.py's docstring on why symmetrization makes
  // that equivalence hold.
  const failedQueryNames = [];
  const memberEdgesByName = new Map(); // name -> edges[] (successful fetches only)
  const memberFetches = querySet.map((m) =>
    fetchPathwaySlice(tissueId, m.id)
      .then((d) => ({ member: m, edges: d.edges }))
      .catch((err) => ({ member: m, error: err }))
  );
  const memberResults = await Promise.all(memberFetches);
  for (const r of memberResults) {
    if (r.error) {
      failedQueryNames.push(r.member.name);
      continue;
    }
    memberEdgesByName.set(r.member.name, r.edges);
  }
  if (failedQueryNames.length === querySet.length) {
    throw new Error(`All ${querySet.length} query pathway fetch(es) failed: ${failedQueryNames.join(", ")}`);
  }

  const candidateScore = new Map(); // pathway_b name -> max abs(correlation)
  for (const edges of memberEdgesByName.values()) {
    for (const e of edges) {
      if (queryNames.has(e.pathway_b)) continue; // exclude other seeds from candidacy
      if (Math.abs(e.correlation) > corrCutOff && passesPCutoff(e, pCutOff)) {
        const score = Math.abs(e.correlation);
        const prev = candidateScore.get(e.pathway_b);
        if (prev === undefined || score > prev) candidateScore.set(e.pathway_b, score);
      }
    }
  }
  const neighborNames = Array.from(candidateScore.entries())
    .sort((a, b) => b[1] - a[1])
    .slice(0, topN)
    .map(([name]) => name);

  // Stage 2: induced subgraph among querySet ∪ neighborNames, unfiltered by
  // the cutoffs (matching the legacy app's spec). Query-query edges among
  // the seeds themselves appear automatically since both endpoints are
  // already in `allowedNames` - no special-casing needed.
  const allowedNames = new Set([...queryNames, ...neighborNames]);
  const edgeMap = new Map();
  function addEdge(aName, bName, rec) {
    const key = [aName, bName].sort().join("|");
    if (!edgeMap.has(key)) {
      edgeMap.set(key, { pathway_a: aName, pathway_b: bName, correlation: rec.correlation, [pColumn]: rec[pColumn], [fdrColumn]: rec[fdrColumn] });
    }
  }
  for (const [memberName, edges] of memberEdgesByName) {
    for (const e of edges) {
      if (allowedNames.has(e.pathway_b)) addEdge(memberName, e.pathway_b, e);
    }
  }

  const failedNames = [];
  const neighborFetches = neighborNames.map((name) => {
    const entry = pathwayByName.get(name);
    if (!entry) return Promise.resolve({ name, error: new Error("not found in pathway_index") });
    return fetchPathwaySlice(tissueId, entry.pathway_id)
      .then((d) => ({ name, edges: d.edges }))
      .catch((err) => ({ name, error: err }));
  });
  const neighborResults = await Promise.all(neighborFetches);
  for (const r of neighborResults) {
    if (r.error) {
      failedNames.push(r.name);
      continue;
    }
    for (const e of r.edges) {
      if (allowedNames.has(e.pathway_b)) addEdge(r.name, e.pathway_b, e);
    }
  }

  return { edges: Array.from(edgeMap.values()), neighborNames, failedQueryNames, failedNames };
}

function buildGroupByNameMap() {
  const map = new Map();
  groups.forEach((g, idx) => {
    const name = g.name && g.name.trim() ? g.name.trim() : `Group ${idx + 1}`;
    for (const p of g.pathways) {
      map.set(p.pathway_name, { id: g.id, name, color: g.color });
    }
  });
  return map;
}

function renderResultMeta(selectedCount, groupCount, neighborCount, edgeCount) {
  resultMeta.hidden = false;
  copyLinkButton.hidden = false;
  resultMeta.textContent = `${selectedCount} pathway(s) across ${groupCount} group(s) + ${neighborCount} neighbor(s), ${edgeCount} edges`;
}

function renderLegend() {
  if (!lastQuery) {
    legendEl.innerHTML = "";
    return;
  }
  const items = [];
  for (const g of lastQuery.groups) {
    if (g.pathways.length === 0) continue;
    items.push(
      `<li><span class="legend-swatch" style="background:${cssColor(g.color)}"></span>${escapeHtml(g.name)}</li>`
    );
  }
  items.push(`<li><span class="legend-swatch" style="background:${NEIGHBOR_COLOR}"></span>Neighbor</li>`);
  const extras = activeTissueCache() && activeTissueCache().extraNames ? activeTissueCache().extraNames : [];
  if (extras.length) {
    items.push(`<li><span class="legend-swatch" style="background:${ADDED_COLOR}"></span>Added</li>`);
  }
  legendEl.innerHTML = `<ul class="legend-list">${items.join("")}</ul>`;
}

function renderTissueTabs(tissueIds) {
  tissueTabList.innerHTML = "";
  tissueTabList.hidden = tissueIds.length === 0;
  for (const id of tissueIds) {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "tab tissue-tab";
    btn.setAttribute("role", "tab");
    btn.dataset.tissue = id;
    btn.textContent = id;
    const cache = tissueResults[id];
    if (cache && cache.status === "error") {
      btn.classList.add("is-error");
      btn.title = cache.errorMessage || "Failed to load this tissue";
    } else if (cache && (cache.failedQueryNames.length || cache.failedNames.length)) {
      btn.title = "Some files failed to load for this tissue";
    }
    btn.addEventListener("click", () => setActiveTissue(id));
    tissueTabList.appendChild(btn);
  }
}

function summarizeTissueWarnings(tissueId, result) {
  const bits = [];
  if (result.failedQueryNames.length) {
    bits.push(
      `${result.failedQueryNames.length} query pathway fetch(es) failed and were excluded: ${result.failedQueryNames.join(", ")}`
    );
  }
  if (result.failedNames.length) {
    bits.push(
      `${result.failedNames.length} neighbor fetch(es) failed and were omitted: ${result.failedNames.join(", ")}`
    );
  }
  if (!bits.length) return "";
  return tissueId ? `${tissueId}: ${bits.join(" ")}` : bits.join(" ");
}

function applyTissueCache(tissueId) {
  const cache = tissueResults[tissueId];
  if (!cache || cache.status === "error") {
    cachedEdges = [];
    cachedMaxAbsCorr = 0;
    resultMeta.hidden = false;
    copyLinkButton.hidden = false;
    resultMeta.textContent = cache && cache.errorMessage ? cache.errorMessage : `No results for ${tissueId}`;
    setStatus(cache && cache.errorMessage ? cache.errorMessage : "");
    renderLegend();
    renderAll();
    return;
  }
  cachedEdges = cache.edges;
  cachedMaxAbsCorr = cache.maxAbsCorr;
  if (lastQuery) {
    corrSlider.max = Math.max(1, cachedMaxAbsCorr, lastQuery.corrCutOff);
    pSlider.max = Math.max(0.05, lastQuery.pCutOff);
  }
  renderResultMeta(
    lastQuery ? lastQuery.selectedCount : flattenQuerySet().length,
    lastQuery ? lastQuery.groupCount : groups.length,
    cache.neighborNames.length,
    cache.edges.length
  );
  const nSelected = lastQuery ? lastQuery.selectedCount : flattenQuerySet().length;
  const nGroups = lastQuery ? lastQuery.groupCount : groups.length;
  setStatus(
    `${cache.edges.length} edges among ${nSelected + cache.neighborNames.length} pathways (${nSelected} selected across ${nGroups} group(s) + ${cache.neighborNames.length} neighbors) in ${tissueId}.`
  );
  renderLegend();
  renderAll();
}

function setActiveTissue(tissueId, rerender = true) {
  activeTissueId = tissueId;
  tissueTabList.querySelectorAll(".tissue-tab").forEach((btn) => {
    btn.setAttribute("aria-selected", btn.dataset.tissue === tissueId ? "true" : "false");
  });
  if (rerender) applyTissueCache(tissueId);
}

// ---------- GO / Reset handlers ----------
async function onGo() {
  clearError();
  closeAllGroupMenus();
  closeTissueMenu();
  if (!pathwayIndexData) {
    showError("Pathway index is still loading.");
    return;
  }
  const querySet = flattenQuerySet();
  if (querySet.length === 0) {
    showError("Add at least one pathway to a group before running Explore.");
    return;
  }
  const tissueIds = selectedTissueIds();
  if (tissueIds.length === 0) {
    showError("Select at least one tissue before running Explore.");
    return;
  }
  const topN = Math.max(0, parseInt(topNInput.value, 10) || 0);
  const corrCutOff = parseFloat(corrCutoffInput.value) || 0;
  const pCutOff = parseFloat(pCutoffInput.value) || 1;

  const perTissue = querySet.length + topN;
  const totalFetches = perTissue * tissueIds.length;
  if (totalFetches > MAX_TOTAL_FETCHES) {
    showError(
      `This query would fetch ${totalFetches} files (${tissueIds.length} tissue(s) × (${querySet.length} selected pathway(s) + topN ${topN})), above the ${MAX_TOTAL_FETCHES}-file limit. Reduce the number of tissues, selected pathways, or topN.`
    );
    return;
  }

  setStatus(tissueIds.length === 1 ? "Querying…" : `Querying ${tissueIds.join(", ")}…`);
  setBusy(true);
  try {
    const settled = await Promise.all(
      tissueIds.map(async (tissueId) => {
        try {
          const result = await runExploreQuery(querySet, topN, corrCutOff, pCutOff, tissueId);
          return {
            tissueId,
            status: "ok",
            edges: result.edges,
            neighborNames: result.neighborNames,
            extraNames: [],
            failedQueryNames: result.failedQueryNames,
            failedNames: result.failedNames,
            maxAbsCorr: result.edges.reduce((m, e) => Math.max(m, Math.abs(e.correlation)), 0),
            errorMessage: "",
          };
        } catch (err) {
          return {
            tissueId,
            status: "error",
            edges: [],
            neighborNames: [],
            extraNames: [],
            failedQueryNames: [],
            failedNames: [],
            maxAbsCorr: 0,
            errorMessage: err.message,
          };
        }
      })
    );

    tissueResults = {};
    const warnings = [];
    for (const item of settled) {
      tissueResults[item.tissueId] = item;
      if (item.status === "error") {
        warnings.push(`${item.tissueId}: ${item.errorMessage}`);
      } else {
        const warn = summarizeTissueWarnings(item.tissueId, item);
        if (warn) warnings.push(warn);
      }
    }
    if (warnings.length) showError(warnings.join(" "));

    const firstOk = tissueIds.find((id) => tissueResults[id].status === "ok");
    if (!firstOk) {
      showError(`Query failed for every selected tissue. ${warnings.join(" ")}`);
      setStatus("");
      return;
    }

    hiddenGraphNames = new Set();
    lastQuery = {
      groupByName: buildGroupByNameMap(),
      groups: groups.map((g) => ({
        id: g.id,
        name: g.name && g.name.trim() ? g.name.trim() : g.name,
        color: g.color,
        pathways: g.pathways.map((p) => p.pathway_name),
      })),
      corrCutOff,
      pCutOff,
      selectedCount: querySet.length,
      groupCount: groups.length,
      tissueIds,
    };

    // Sliders default to the old app's ranges (corr 0-1, p 0.001-0.05); widen
    // them if this query's own cutoffs fall outside that, so .value isn't
    // silently clamped to a stricter cutoff than what was actually queried.
    const firstCache = tissueResults[firstOk];
    corrSlider.max = Math.max(1, firstCache.maxAbsCorr, corrCutOff);
    pSlider.max = Math.max(0.05, pCutOff);
    corrSlider.value = corrCutOff;
    pSlider.value = pCutOff;
    updateSliderLabels();
    liveFilters.hidden = false;
    showResults(true);
    renderTissueTabs(tissueIds);
    setActiveTab("network", false);
    setActiveTissue(firstOk, false);
    writeUrlState();
    applyTissueCache(firstOk);
  } catch (err) {
    showError(`Query failed: ${err.message}`);
    setStatus("");
  } finally {
    setBusy(false);
  }
}

async function onReset() {
  clearError();
  closeAllGroupMenus();
  closeTissueMenu();
  resetGroups();
  collectionSelect.value = "";
  topNInput.value = 10;
  corrCutoffInput.value = 0.5;
  pCutoffInput.value = 0.05;
  lastQuery = null;
  hiddenGraphNames = new Set();
  closeGraphAddMenu();
  cachedEdges = null;
  cachedMaxAbsCorr = 0;
  tissueResults = {};
  activeTissueId = null;
  tissueTabList.innerHTML = "";
  tissueTabList.hidden = true;
  const defaults = defaultTissueIds();
  setTissueChecks(defaults);
  const ref = referenceTissueId(defaults);
  if (ref && ref !== currentTissueId) {
    await loadTissue(ref);
  }
  liveFilters.hidden = true;
  resultMeta.hidden = true;
  copyLinkButton.hidden = true;
  resultMeta.textContent = "";
  legendEl.innerHTML = "";
  setStatus("");
  nodeInfo.textContent = "";
  showResults(false);
  networkSection.hidden = true;
  heatmapSection.hidden = true;
  tableSection.hidden = true;
  writeUrlState();
  if (cy) {
    cy.destroy();
    cy = null;
  }
  if (dataTable) {
    dataTable.destroy();
    el("edgeTable").querySelector("tbody").innerHTML = "";
    dataTable = null;
  }
}

function resolveExamplePathway(pathwayName) {
  return (
    findPathwayByName(pathwayName) ||
    (pathwayName.startsWith("Pathway.") ? null : findPathwayByName(`Pathway.${pathwayName}`))
  );
}

function fillExampleGroup(group, pathwayNames, groupName) {
  if (groupName) {
    group.name = groupName;
    const nameInput = group.rowEl && group.rowEl.querySelector(".group-name-input");
    if (nameInput) nameInput.value = groupName;
  }
  group.pathways = [];
  for (const pathwayName of pathwayNames) {
    const entry = resolveExamplePathway(pathwayName);
    if (entry && !allChippedNames().has(entry.pathway_name)) {
      group.pathways.push(entry);
    }
  }
  renderChips(group);
}

async function runExample(spec, collection, groupName) {
  if (Array.isArray(spec) || typeof spec === "string") {
    spec = {
      groups: [{ name: groupName || "", pathways: Array.isArray(spec) ? spec : [spec] }],
      collection,
    };
  }
  if (spec.collection) {
    if (collectionSelect.querySelector(`option[value="${CSS.escape(spec.collection)}"]`)) {
      collectionSelect.value = spec.collection;
    }
  } else {
    collectionSelect.value = "";
  }
  topNInput.value = spec.topN != null && Number.isFinite(Number(spec.topN)) ? String(spec.topN) : "10";
  corrCutoffInput.value =
    spec.corr != null && Number.isFinite(Number(spec.corr)) ? String(spec.corr) : "0.5";
  const payloads = (spec.groups || []).filter((g) => g && (g.pathways || []).length);
  if (payloads.length === 0) return;
  clearAllGroupRows();
  for (const gp of payloads) {
    const group = makeGroup();
    addGroupRow(group);
    fillExampleGroup(group, gp.pathways, gp.name);
  }
  updateRemoveButtonsVisibility();
  updateFetchCounter();
  await onGo();
}

async function copyShareLink() {
  writeUrlState();
  const url = window.location.href;
  try {
    await navigator.clipboard.writeText(url);
    copyLinkButton.textContent = "Copied";
    setTimeout(() => {
      copyLinkButton.textContent = "Copy link";
    }, 1400);
  } catch {
    showError(`Could not copy automatically. URL: ${url}`);
  }
}

// ---------- Live re-filter sliders (no new fetch, matches old app) ----------
function updateSliderLabels() {
  corrSliderValue.textContent = Number(corrSlider.value).toFixed(2);
  pSliderValue.textContent = Number(pSlider.value).toFixed(3);
}
function onSliderChange() {
  updateSliderLabels();
  if (cachedEdges) renderAll();
}

// ---------- Rendering ----------
function sliderPassingEdges() {
  if (!cachedEdges) return [];
  const corrCutOff = parseFloat(corrSlider.value);
  const pCutOff = parseFloat(pSlider.value);
  return cachedEdges.filter((e) => Math.abs(e.correlation) > corrCutOff && passesPCutoff(e, pCutOff));
}

function currentlyVisibleEdges() {
  return sliderPassingEdges().filter(
    (e) => !hiddenGraphNames.has(e.pathway_a) && !hiddenGraphNames.has(e.pathway_b)
  );
}

function activeTissueCache() {
  return activeTissueId ? tissueResults[activeTissueId] : null;
}

function graphRosterEntries() {
  const seen = new Set();
  const rows = [];
  const cache = activeTissueCache();
  const groupByName = lastQuery ? lastQuery.groupByName : new Map();

  function add(name, fallbackRole) {
    if (!name || seen.has(name)) return;
    seen.add(name);
    const entry = pathwayByName.get(name);
    const groupInfo = groupByName.get(name);
    const extra = cache && cache.extraNames && cache.extraNames.includes(name);
    rows.push({
      name,
      role: groupInfo ? groupInfo.name : extra ? "Added" : fallbackRole,
      roleKind: groupInfo ? "selected" : extra ? "added" : "neighbor",
      color: groupInfo ? cssColor(groupInfo.color) : extra ? ADDED_COLOR : NEIGHBOR_COLOR,
      collection: entry && entry.source_db ? sourceDbLabel(entry.source_db) : "",
      hidden: hiddenGraphNames.has(name),
    });
  }

  if (lastQuery) {
    for (const g of lastQuery.groups) {
      for (const name of g.pathways) add(name, g.name);
    }
  }
  if (cache) {
    for (const name of cache.neighborNames || []) add(name, "Neighbor");
    for (const name of cache.extraNames || []) add(name, "Added");
    for (const e of cache.edges || []) {
      add(e.pathway_a, "Neighbor");
      add(e.pathway_b, "Neighbor");
    }
  }
  const order = { selected: 0, added: 1, neighbor: 2 };
  rows.sort((a, b) => order[a.roleKind] - order[b.roleKind] || a.name.localeCompare(b.name));
  return rows;
}

function degreeIfShown(name) {
  return connectionsFor(name).length;
}

function pathwayRoleInfo(name) {
  const groupByName = lastQuery ? lastQuery.groupByName : new Map();
  const groupInfo = groupByName.get(name);
  if (groupInfo) return { role: groupInfo.name, color: cssColor(groupInfo.color) };
  const cache = activeTissueCache();
  if (cache && cache.extraNames && cache.extraNames.includes(name)) {
    return { role: "Added", color: ADDED_COLOR };
  }
  return { role: "Neighbor", color: NEIGHBOR_COLOR };
}

function connectionsFor(name) {
  const rows = [];
  for (const e of sliderPassingEdges()) {
    const other = e.pathway_a === name ? e.pathway_b : e.pathway_b === name ? e.pathway_a : null;
    if (!other || hiddenGraphNames.has(other)) continue;
    const info = pathwayRoleInfo(other);
    rows.push({
      name: other,
      role: info.role,
      color: info.color,
      correlation: e.correlation,
      p: e[pColumn],
      fdr: e[fdrColumn],
    });
  }
  rows.sort((a, b) => Math.abs(b.correlation) - Math.abs(a.correlation));
  return rows;
}

let linksDialogPathway = "";

function openPathwayLinks(pathwayName) {
  linksDialogPathway = pathwayName;
  const rows = connectionsFor(pathwayName);
  linksDialogTitle.textContent = pathwayName;
  if (rows.length === 0) {
    linksDialogMeta.textContent = "No connections in the current view.";
    linksDialogList.innerHTML = `<p class="genes-dialog-empty">This pathway has no edges after the current cutoffs and hidden nodes.</p>`;
    copyPathwayLinksButton.hidden = true;
  } else {
    linksDialogMeta.textContent = `${rows.length} connection${rows.length === 1 ? "" : "s"}`;
    copyPathwayLinksButton.hidden = false;
    const body = rows
      .map(
        (link) => `<tr>
          <td>${escapeHtml(link.name)}</td>
          <td><span class="graph-swatch" style="background:${cssColor(link.color, NEIGHBOR_COLOR)}"></span>${escapeHtml(link.role)}</td>
          <td class="num">${formatCorrelation(link.correlation)}</td>
          <td class="num">${formatMetric(link.p)}</td>
          <td class="num">${formatMetric(link.fdr)}</td>
        </tr>`
      )
      .join("");
    linksDialogList.innerHTML = `<table><thead><tr><th>Pathway</th><th>Role</th><th>Correlation</th><th>p-value</th><th>FDR</th></tr></thead><tbody>${body}</tbody></table>`;
  }
  linksDialog.showModal();
}

function copyConnections(name, button) {
  const rows = connectionsFor(name);
  if (rows.length === 0) {
    showError("No connections to copy for this pathway.");
    return;
  }
  const header = ["pathway", "role", "correlation", "p_value", "fdr"];
  const lines = [
    header.join("\t"),
    ...rows.map((r) =>
      [r.name, r.role, formatCorrelation(r.correlation), formatMetric(r.p), formatMetric(r.fdr)].join("\t")
    ),
  ];
  copyText(lines.join("\n"), button, "Copy");
}

function setGraphNodeHidden(name, hidden) {
  if (hidden) hiddenGraphNames.add(name);
  else hiddenGraphNames.delete(name);
  renderAll();
}

function edgePairKey(a, b) {
  return a < b ? `${a}|${b}` : `${b}|${a}`;
}

function closeGraphAddMenu() {
  graphAddMenu.hidden = true;
  graphAddMenu.innerHTML = "";
  graphAddItems = [];
  graphAddIndex = -1;
  graphAddInput.setAttribute("aria-expanded", "false");
}

function openGraphAddMenu() {
  if (!pathwayIndexData || !lastQuery) return;
  const visible = new Set(graphRosterEntries().filter((row) => !row.hidden).map((row) => row.name));
  const matches = filteredPathways(graphAddInput.value, visible).slice(0, 40);
  graphAddMenu.innerHTML = "";
  graphAddItems = [];
  graphAddIndex = -1;
  if (matches.length === 0) {
    const empty = document.createElement("li");
    empty.className = "combo-empty";
    empty.textContent = "No matching pathways";
    graphAddMenu.appendChild(empty);
  } else {
    const frag = document.createDocumentFragment();
    for (const p of matches) {
      const li = document.createElement("li");
      li.setAttribute("role", "option");
      li.dataset.name = p.pathway_name;
      li.innerHTML = `<div class="combo-option"><span class="combo-text"><span class="combo-name">${escapeHtml(
        p.pathway_name
      )}</span><span class="combo-meta">${escapeHtml(pathwayMetaLabel(p))}</span></span></div>`;
      li.addEventListener("mousedown", (evt) => {
        evt.preventDefault();
        addPathwayToGraph(p);
      });
      frag.appendChild(li);
      graphAddItems.push(li);
    }
    graphAddMenu.appendChild(frag);
  }
  graphAddMenu.hidden = false;
  graphAddInput.setAttribute("aria-expanded", "true");
}

function highlightGraphAddItem() {
  graphAddItems.forEach((item, i) => {
    item.setAttribute("aria-selected", i === graphAddIndex ? "true" : "false");
  });
  const current = graphAddItems[graphAddIndex];
  if (current) current.scrollIntoView({ block: "nearest" });
}

async function addPathwayToGraph(entry) {
  if (!entry || !activeTissueId || !lastQuery) return;
  const name = entry.pathway_name;
  hiddenGraphNames.delete(name);
  const cache = activeTissueCache();
  const roster = new Set(graphRosterEntries().map((row) => row.name));
  if (roster.has(name)) {
    graphAddInput.value = "";
    closeGraphAddMenu();
    renderAll();
    return;
  }
  graphAddButton.disabled = true;
  try {
    const data = await fetchPathwaySlice(activeTissueId, entry.pathway_id);
    const existing = new Set(cache.edges.map((e) => edgePairKey(e.pathway_a, e.pathway_b)));
    for (const rec of data.edges || []) {
      if (!roster.has(rec.pathway_b)) continue;
      const key = edgePairKey(name, rec.pathway_b);
      if (existing.has(key)) continue;
      existing.add(key);
      cache.edges.push({
        pathway_a: name,
        pathway_b: rec.pathway_b,
        correlation: rec.correlation,
        [pColumn]: rec[pColumn],
        [fdrColumn]: rec[fdrColumn],
      });
    }
    cache.extraNames = cache.extraNames || [];
    if (!cache.extraNames.includes(name)) cache.extraNames.push(name);
    cache.maxAbsCorr = cache.edges.reduce((m, e) => Math.max(m, Math.abs(e.correlation)), 0);
    cachedEdges = cache.edges;
    cachedMaxAbsCorr = cache.maxAbsCorr;
    graphAddInput.value = "";
    closeGraphAddMenu();
    renderAll();
  } catch (err) {
    showError(`Could not add ${name} to the graph. ${err.message}`);
  } finally {
    graphAddButton.disabled = false;
  }
}

function submitGraphAddFromInput() {
  if (!pathwayIndexData) return;
  const visible = new Set(graphRosterEntries().filter((row) => !row.hidden).map((row) => row.name));
  const matches = filteredPathways(graphAddInput.value, visible);
  const typed = graphAddInput.value.trim().toLowerCase();
  const exact = matches.find((p) => p.pathway_name.toLowerCase() === typed);
  const pick = exact || (graphAddItems[graphAddIndex] && resolveExamplePathway(graphAddItems[graphAddIndex].dataset.name)) || matches[0];
  if (pick) addPathwayToGraph(pick);
}

function renderGraphNodeTable() {
  const rows = graphRosterEntries();
  const hiddenCount = rows.filter((row) => row.hidden).length;
  graphNodeMeta.textContent = rows.length
    ? `${rows.length - hiddenCount} of ${rows.length} shown`
    : "";
  graphShowAll.hidden = hiddenCount === 0;
  graphNodeTableBody.innerHTML = "";
  rows.forEach((row, i) => {
    const tr = document.createElement("tr");
    if (row.hidden) tr.classList.add("is-hidden");
    const id = `graph-node-${i}`;
    tr.innerHTML = `<td><input id="${id}" type="checkbox"${row.hidden ? "" : " checked"}></td>
      <td class="graph-name"><label for="${id}">${escapeHtml(row.name)}</label></td>
      <td class="graph-role"><span class="graph-swatch" style="background:${cssColor(row.color, NEIGHBOR_COLOR)}"></span>${escapeHtml(row.role)}${
        row.collection ? ` · ${escapeHtml(row.collection)}` : ""
      }</td>
      <td class="graph-edges">${connectionsFor(row.name).length}</td>
      <td><button type="button" class="genes-view-btn">View</button></td>
      <td><button type="button" class="links-view-btn">View</button></td>`;
    const box = tr.querySelector("input");
    box.addEventListener("change", () => setGraphNodeHidden(row.name, !box.checked));
    tr.querySelector(".genes-view-btn").addEventListener("click", () => openPathwayGenes(row.name));
    tr.querySelector(".links-view-btn").addEventListener("click", () => openPathwayLinks(row.name));
    graphNodeTableBody.appendChild(tr);
  });
}

// Each panel is rendered independently: one panel failing (e.g. a CDN
// library that didn't load) must not prevent the others from rendering, and
// must not look like the query itself failed (see onGo, which no longer
// wraps this call in its own try/catch for exactly that reason).
function setActiveTab(tab, rerender = true) {
  activeTab = tab;
  document.querySelectorAll(".view-tab").forEach((btn) => {
    btn.setAttribute("aria-selected", btn.dataset.tab === tab ? "true" : "false");
  });
  if (rerender && cachedEdges) renderAll();
}

function renderAll() {
  const visible = currentlyVisibleEdges();
  networkSection.hidden = activeTab !== "network";
  tableSection.hidden = activeTab !== "table";
  heatmapSection.hidden = activeTab !== "heatmap";

  const renderErrors = [];
  function tryRender(label, fn) {
    try {
      fn(visible);
    } catch (err) {
      console.error(`${label} failed to render:`, err);
      renderErrors.push(`${label} failed to render (${err.message}).`);
    }
  }
  if (activeTab === "network") {
    tryRender("Network", renderNetwork);
    tryRender("Pathway table", renderGraphNodeTable);
    if (isCyLive()) {
      requestAnimationFrame(() => {
        if (!isCyLive()) return;
        cy.resize();
        fitNetwork();
      });
    }
  }
  if (activeTab === "table") tryRender("Edge table", renderTable);
  if (activeTab === "heatmap") tryRender("Heatmap", renderHeatmap);
  tryRender("Legend", renderLegend);

  if (renderErrors.length > 0) {
    showError(
      renderErrors.join(" ") +
        " If this mentions a library (Cytoscape/DataTable/CanvasXpress) being undefined, its CDN script " +
        "didn't load - check your browser's console/network tab for a blocked or failed request."
    );
  }
}

function edgeColor(correlation) {
  const ratio = cachedMaxAbsCorr > 0 ? Math.min(1, Math.abs(correlation) / cachedMaxAbsCorr) : 0;
  const intensity = Math.round(80 + ratio * 175);
  return correlation >= 0 ? `rgb(${intensity},40,40)` : `rgb(40,40,${intensity})`;
}

// Query-set members use their group color; everything else is a neighbor.
function nodeColor(node) {
  const groupColor = node.data("groupColor");
  if (groupColor) return cssColor(groupColor, NEIGHBOR_COLOR);
  if (node.data("role") === "added") return ADDED_COLOR;
  return NEIGHBOR_COLOR;
}

const NODE_LABEL_MAX = 15;
const NODE_LABEL_PREFIXES = [
  [/^Pathway\./, ""],
  [/^REACTOME_/, "R_"],
  [/^BIOCARTA_/, "B_"],
  [/^WIKIPATHWAYS_/, "WP_"],
  [/^KEGG_/, "K_"],
  [/^PID_/, "P_"],
];

function formatNodeLabel(name) {
  let short = String(name || "");
  for (const [pattern, repl] of NODE_LABEL_PREFIXES) {
    if (pattern.test(short)) short = short.replace(pattern, repl);
  }
  if (short.length <= NODE_LABEL_MAX) return short;
  return `${short.slice(0, NODE_LABEL_MAX - 1)}…`;
}

function networkLayoutOptions() {
  const name = layoutSelect.value;
  const shared = { name, fit: false, animate: false, padding: 40, stop: fitNetwork };
  if (name === "concentric") {
    return {
      ...shared,
      minNodeSpacing: 80,
      equidistant: true,
      spacingFactor: 1.75,
      concentric: (n) => (n.data("groupColor") ? 2 : 1),
      levelWidth: () => 1,
    };
  }
  if (name === "circle" || name === "grid") {
    return { ...shared, spacingFactor: 1.35, avoidOverlap: true };
  }
  return {
    ...shared,
    idealEdgeLength: 240,
    nodeOverlap: 60,
    nodeRepulsion: () => 20000,
    componentSpacing: 160,
    nestingFactor: 1.2,
    gravity: 0.05,
  };
}

function nodeBodySize(node) {
  return 6 + Math.min(8, node.data("degree") || 1);
}

function zoomNetwork(factor) {
  if (!isCyLive()) return;
  cy.zoom({
    level: cy.zoom() * factor,
    renderedPosition: { x: cy.width() / 2, y: cy.height() / 2 },
  });
}

function isCyLive() {
  return Boolean(cy && !cy.destroyed());
}

function expandCompactPositions() {
  if (!isCyLive()) return;
  let minX = Infinity;
  let maxX = -Infinity;
  let minY = Infinity;
  let maxY = -Infinity;
  cy.nodes().forEach((n) => {
    const p = n.position();
    minX = Math.min(minX, p.x);
    maxX = Math.max(maxX, p.x);
    minY = Math.min(minY, p.y);
    maxY = Math.max(maxY, p.y);
  });
  const span = Math.max(maxX - minX, maxY - minY, 1);
  const minSpan = Math.min(cy.width(), cy.height()) * 0.82;
  if (span >= minSpan) return;
  const scale = minSpan / span;
  const cx = (minX + maxX) / 2;
  const cy0 = (minY + maxY) / 2;
  cy.nodes().positions((n) => {
    const p = n.position();
    return { x: cx + (p.x - cx) * scale, y: cy0 + (p.y - cy0) * scale };
  });
}

function fitNetwork() {
  if (!isCyLive() || cy.nodes().empty()) return;
  if (cy.width() < 40 || cy.height() < 40) {
    requestAnimationFrame(() => {
      cy.resize();
      fitNetwork();
    });
    return;
  }
  if (layoutSelect.value === "cose") expandCompactPositions();
  cy.fit(undefined, 40);
  let maxRendered = 0;
  cy.nodes().forEach((n) => {
    maxRendered = Math.max(maxRendered, n.renderedWidth());
  });
  const cap = 16;
  if (maxRendered > cap) {
    cy.zoom({
      level: cy.zoom() * (cap / maxRendered),
      renderedPosition: { x: cy.width() / 2, y: cy.height() / 2 },
    });
  }
}

function hideNodeTooltip() {
  const tip = el("cyTooltip");
  if (tip) tip.hidden = true;
}

function showNodeTooltip(node) {
  const tip = el("cyTooltip");
  if (!tip) return;
  const fullName = node.data("fullName") || node.id();
  const groupName = node.data("groupName") || "Neighbor";
  tip.innerHTML = `<strong>${escapeHtml(fullName)}</strong>${escapeHtml(groupName)}`;
  const pos = node.renderedPosition();
  tip.style.left = `${Math.round(pos.x + 14)}px`;
  tip.style.top = `${Math.round(pos.y + 14)}px`;
  tip.hidden = false;
}

function applyLayout() {
  if (!isCyLive()) return;
  cy.layout(networkLayoutOptions()).run();
}

function renderNetwork(edges) {
  const nodeIds = new Set();
  for (const e of edges) {
    nodeIds.add(e.pathway_a);
    nodeIds.add(e.pathway_b);
  }
  for (const row of graphRosterEntries()) {
    if (!row.hidden) nodeIds.add(row.name);
  }
  const degree = new Map();
  for (const e of edges) {
    degree.set(e.pathway_a, (degree.get(e.pathway_a) || 0) + 1);
    degree.set(e.pathway_b, (degree.get(e.pathway_b) || 0) + 1);
  }
  const groupByName = lastQuery.groupByName;
  const extraNames = new Set((activeTissueCache() && activeTissueCache().extraNames) || []);
  const nodes = Array.from(nodeIds).map((name) => {
    const entry = pathwayByName.get(name);
    const groupInfo = groupByName.get(name);
    const added = !groupInfo && extraNames.has(name);
    return {
      data: {
        id: name,
        fullName: name,
        label: formatNodeLabel(name),
        degree: degree.get(name) || 1,
        role: groupInfo ? groupInfo.id : added ? "added" : "neighbor",
        groupName: groupInfo ? groupInfo.name : added ? "Added" : "Neighbor",
        groupColor: groupInfo ? cssColor(groupInfo.color) : added ? ADDED_COLOR : null,
        collection: entry ? pathwayMetaLabel(entry) : "",
      },
    };
  });
  const edgeEls = edges.map((e, i) => ({
    data: {
      id: `e${i}`,
      source: e.pathway_a,
      target: e.pathway_b,
      correlation: e.correlation,
      color: edgeColor(e.correlation),
    },
  }));

  hideNodeTooltip();
  if (cy) {
    try {
      cy.stop();
      cy.destroy();
    } catch {
      // A layout callback can fire after destroy; ignore the stale instance.
    }
    cy = null;
  }
  cy = cytoscape({
    container: el("cy"),
    elements: { nodes, edges: edgeEls },
    style: [
      {
        selector: "node",
        style: {
          label: "data(label)",
          "font-size": "11px",
          width: nodeBodySize,
          height: nodeBodySize,
          "background-color": nodeColor,
          "text-wrap": "ellipsis",
          "text-max-width": "80px",
          "text-valign": "bottom",
          "text-margin-y": 6,
          color: "#1a1f2b",
        },
      },
      {
        selector: "edge",
        style: {
          width: 1.6,
          "line-color": "data(color)",
          "curve-style": "haystack",
          "haystack-radius": 0,
        },
      },
      { selector: "node.faded", style: { opacity: 0.15 } },
      { selector: "edge.faded", style: { opacity: 0.1 } },
    ],
    layout: networkLayoutOptions(),
  });

  cy.on("tap", "node", (evt) => {
    const node = evt.target;
    const neighborhood = node.closedNeighborhood();
    cy.elements().addClass("faded");
    neighborhood.removeClass("faded");
    const name = node.data("fullName") || node.id();
    const deg = node.data("degree");
    const entry = pathwayByName.get(name);
    const collectionInfo = entry && entry.source_db ? ` · ${sourceDbLabel(entry.source_db)}` : "";
    nodeInfo.innerHTML = `${escapeHtml(name)} — ${deg} connection(s) in this view${escapeHtml(
      collectionInfo
    )} <button type="button" class="text-btn" data-action="genes">View genes</button>
    <button type="button" class="text-btn" data-action="links">View connections</button>`;
    const genesBtn = nodeInfo.querySelector('[data-action="genes"]');
    const linksBtn = nodeInfo.querySelector('[data-action="links"]');
    if (genesBtn) genesBtn.addEventListener("click", () => openPathwayGenes(name));
    if (linksBtn) linksBtn.addEventListener("click", () => openPathwayLinks(name));
    showNodeTooltip(node);
  });
  cy.on("tap", (evt) => {
    if (evt.target === cy) {
      cy.elements().removeClass("faded");
      nodeInfo.textContent = "";
      hideNodeTooltip();
    }
  });
  cy.on("mouseover", "node", (evt) => {
    showNodeTooltip(evt.target);
  });
  cy.on("mouseout", "node", hideNodeTooltip);
  cy.on("pan zoom", hideNodeTooltip);
}

function renderTable(edges) {
  if (dataTable) {
    dataTable.destroy();
    dataTable = null;
  }
  const tbody = el("edgeTable").querySelector("tbody");
  tbody.innerHTML = "";
  for (const e of edges) {
    const tr = document.createElement("tr");
    tr.innerHTML = `<td><button type="button" class="pathway-link" data-pathway="${escapeHtml(
      e.pathway_a
    )}">${escapeHtml(e.pathway_a)}</button></td><td><button type="button" class="pathway-link" data-pathway="${escapeHtml(
      e.pathway_b
    )}">${escapeHtml(e.pathway_b)}</button></td><td>${formatCorrelation(
      e.correlation
    )}</td><td>${formatMetric(e[pColumn])}</td><td>${formatMetric(e[fdrColumn])}</td>`;
    tbody.appendChild(tr);
  }
  dataTable = new DataTable("#edgeTable", {
    order: [[2, "desc"]],
  });
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

// Hand-rolled CSV export (no library dependency): builds a CSV blob from the
// currently-visible edges and triggers a download via a throwaway <a>.
function csvField(value) {
  const s = String(value);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}
function downloadHref(href, filename) {
  const a = document.createElement("a");
  a.href = href;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
}

function figureFilename(kind) {
  const tissue = activeTissueId || currentTissueId || "tissue";
  return `pcxn_${kind}_${tissue}_${Date.now()}.png`;
}

function saveNetworkPng() {
  if (!isCyLive() || cy.elements().empty()) {
    showError("No network to save.");
    return;
  }
  const png = cy.png({
    output: "blob-promise",
    bg: "#fffdf8",
    full: true,
    scale: 2,
    maxWidth: 4096,
    maxHeight: 4096,
  });
  Promise.resolve(png)
    .then((blob) => {
      if (!blob) {
        showError("Could not save the network image.");
        return;
      }
      if (typeof blob === "string") {
        downloadHref(blob, figureFilename("network"));
        return;
      }
      const url = URL.createObjectURL(blob);
      downloadHref(url, figureFilename("network"));
      setTimeout(() => URL.revokeObjectURL(url), 1500);
    })
    .catch((err) => {
      showError(`Could not save the network image. ${err && err.message ? err.message : err}`);
    });
}

function saveHeatmapPng() {
  const wrap = el("heatmap");
  const canvases = [];
  if (wrap) {
    if (wrap.tagName === "CANVAS") canvases.push(wrap);
    canvases.push(...wrap.querySelectorAll("canvas"));
    if (wrap.parentElement) {
      for (const c of wrap.parentElement.querySelectorAll("canvas")) {
        if (!canvases.includes(c)) canvases.push(c);
      }
    }
  }
  const usable = canvases.filter((c) => c.width && c.height);
  if (!cxChart || usable.length === 0) {
    showError("No heatmap to save.");
    return;
  }
  let canvas = usable[0];
  if (usable.length > 1) {
    const w = Math.max(...usable.map((c) => c.width));
    const h = Math.max(...usable.map((c) => c.height));
    canvas = document.createElement("canvas");
    canvas.width = w;
    canvas.height = h;
    const ctx = canvas.getContext("2d");
    ctx.fillStyle = "#fffdf8";
    ctx.fillRect(0, 0, w, h);
    for (const c of usable) ctx.drawImage(c, 0, 0);
  }
  if (typeof canvas.toBlob !== "function") {
    downloadHref(canvas.toDataURL("image/png"), figureFilename("heatmap"));
    return;
  }
  canvas.toBlob((blob) => {
    if (!blob) {
      showError("Could not save the heatmap image.");
      return;
    }
    const url = URL.createObjectURL(blob);
    downloadHref(url, figureFilename("heatmap"));
    setTimeout(() => URL.revokeObjectURL(url), 1500);
  }, "image/png");
}

function exportCsv() {
  const edges = currentlyVisibleEdges();
  if (edges.length === 0) {
    showError("No edges to export.");
    return;
  }
  const header = ["pathway_a", "pathway_b", "correlation", pColumn, fdrColumn];
  const rows = edges.map((e) => [e.pathway_a, e.pathway_b, e.correlation, e[pColumn], e[fdrColumn]]);
  const csv = [header, ...rows].map((row) => row.map(csvField).join(",")).join("\r\n");
  const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  const groupCount = lastQuery ? lastQuery.groups.length : groups.length;
  const pathwayCount = lastQuery
    ? lastQuery.groups.reduce((n, g) => n + g.pathways.length, 0)
    : flattenQuerySet().length;
  a.download = `pcxn_${activeTissueId || currentTissueId}_${groupCount}groups_${pathwayCount}pathways_${Date.now()}.csv`;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

function renderHeatmap(edges) {
  const labels = Array.from(new Set(edges.flatMap((e) => [e.pathway_a, e.pathway_b])));
  const corrLookup = new Map();
  for (const e of edges) {
    corrLookup.set(`${e.pathway_a}|${e.pathway_b}`, e.correlation);
    corrLookup.set(`${e.pathway_b}|${e.pathway_a}`, e.correlation);
  }
  const matrix = labels.map((a) =>
    labels.map((b) => (a === b ? 1 : corrLookup.get(`${a}|${b}`) || 0))
  );
  const cxData = {
    y: {
      vars: labels,
      smps: labels,
      data: matrix,
    },
  };
  if (cxChart) {
    el("heatmap").innerHTML = "";
  }
  cxChart = new CanvasXpress("heatmap", cxData, {
    graphType: "Heatmap",
    heatmapIndicatorHistogram: false,
    showTransition: false,
    showToolbar: false,
  });
}

// ---------- Wire up events ----------
collectionSelect.addEventListener("change", () => {
  closeAllGroupMenus();
  if (bulkLookupDialog.open) renderBulkLookupResults();
  if (fromGenesDialog.open) {
    fromGenesPrefillKey = null;
    renderFromGenesResults();
  }
});
addGroupButton.addEventListener("click", () => {
  addGroupRow(makeGroup());
});
topNInput.addEventListener("input", updateFetchCounter);
queryForm.addEventListener("submit", (e) => {
  e.preventDefault();
  onGo();
});
resetButton.addEventListener("click", onReset);
saveNetworkPngButton.addEventListener("click", saveNetworkPng);
saveHeatmapPngButton.addEventListener("click", saveHeatmapPng);
exportCsvButton.addEventListener("click", exportCsv);
copyPathwaysButton.addEventListener("click", copyVisiblePathways);
copyGenesButton.addEventListener("click", copyVisiblePathwayGenes);
copyPathwayGenesButton.addEventListener("click", () => {
  const symbols = pathwayGeneSymbols(genesDialogPathway);
  if (symbols.length === 0) {
    showError("No genes to copy for this pathway.");
    return;
  }
  copyText(symbols.join("\n"), copyPathwayGenesButton, "Copy genes");
});
genesDialogClose.addEventListener("click", () => genesDialog.close());
copyPathwayLinksButton.addEventListener("click", () => {
  copyConnections(linksDialogPathway, copyPathwayLinksButton);
});
linksDialogClose.addEventListener("click", () => linksDialog.close());
el("edgeTable").addEventListener("click", (evt) => {
  const btn = evt.target.closest("[data-pathway]");
  if (btn) openPathwayGenes(btn.dataset.pathway);
});
copyLinkButton.addEventListener("click", copyShareLink);
corrSlider.addEventListener("input", onSliderChange);
pSlider.addEventListener("input", onSliderChange);
layoutSelect.addEventListener("change", applyLayout);
fitButton.addEventListener("click", fitNetwork);
graphShowAll.addEventListener("click", () => {
  hiddenGraphNames = new Set();
  renderAll();
});
graphAddInput.addEventListener("focus", openGraphAddMenu);
graphAddInput.addEventListener("input", openGraphAddMenu);
graphAddInput.addEventListener("keydown", (evt) => {
  if (evt.key === "ArrowDown") {
    evt.preventDefault();
    if (graphAddMenu.hidden) openGraphAddMenu();
    if (graphAddItems.length) {
      graphAddIndex = (graphAddIndex + 1) % graphAddItems.length;
      highlightGraphAddItem();
    }
  } else if (evt.key === "ArrowUp") {
    evt.preventDefault();
    if (graphAddItems.length) {
      graphAddIndex = (graphAddIndex - 1 + graphAddItems.length) % graphAddItems.length;
      highlightGraphAddItem();
    }
  } else if (evt.key === "Enter") {
    evt.preventDefault();
    submitGraphAddFromInput();
  } else if (evt.key === "Escape") {
    evt.preventDefault();
    closeGraphAddMenu();
  }
});
graphAddButton.addEventListener("click", submitGraphAddFromInput);
el("cyZoomIn").addEventListener("click", () => zoomNetwork(1.25));
el("cyZoomOut").addEventListener("click", () => zoomNetwork(0.8));
el("cyFit").addEventListener("click", fitNetwork);
if (citeButton && citeDialog) {
  citeButton.addEventListener("click", () => citeDialog.showModal());
}
aboutButton.addEventListener("click", () => aboutDialog.showModal());
bulkLookupTerms.addEventListener("input", renderBulkLookupResults);
document.querySelectorAll('input[name="bulkLookupMode"]').forEach((radio) => {
  radio.addEventListener("change", renderBulkLookupResults);
});
bulkLookupSelectAll.addEventListener("click", () => setVisibleBulkSelection(true));
bulkLookupClearVisible.addEventListener("click", () => setVisibleBulkSelection(false));
bulkLookupDone.addEventListener("click", () => bulkLookupDialog.close());
bulkLookupDialog.addEventListener("close", () => {
  bulkLookupGroup = null;
  bulkLookupItems = [];
  bulkLookupList.innerHTML = "";
});
fromGenesList.addEventListener("input", () => {
  renderFromGenesResults();
});
document.querySelectorAll(".from-genes-example").forEach((btn) => {
  btn.addEventListener("click", () => {
    const genes = (btn.dataset.genes || "").trim().split(/\s+/).filter(Boolean);
    fromGenesList.value = genes.join("\n");
    renderFromGenesResults();
    fromGenesList.focus();
  });
});
fromGenesQCutoff.addEventListener("input", () => {
  renderFromGenesResults();
});
fromGenesSelectAll.addEventListener("click", () => setVisibleFromGenesSelection(true));
fromGenesClearVisible.addEventListener("click", () => setVisibleFromGenesSelection(false));
fromGenesDone.addEventListener("click", () => fromGenesDialog.close());
fromGenesDialog.addEventListener("close", () => {
  fromGenesGroup = null;
  fromGenesItems = [];
  lastFromGenesHits = [];
  fromGenesPrefillKey = null;
  fromGenesSortKey = "q";
  fromGenesSortDir = "asc";
  fromGenesTableBody.innerHTML = "";
});
fromGenesTable.querySelector("thead").addEventListener("click", onFromGenesSortClick);
el("errorDismiss").addEventListener("click", clearError);
document.querySelectorAll(".view-tab").forEach((btn) => {
  btn.addEventListener("click", () => setActiveTab(btn.dataset.tab));
});

document.addEventListener("click", (evt) => {
  for (const g of groups) {
    if (g.combobox && !g.combobox.containerEl.contains(evt.target)) closeGroupMenu(g);
  }
  if (tissueCombobox && !tissueCombobox.contains(evt.target)) closeTissueMenu();
  if (el("graphAddInput") && !el("graphAddInput").closest(".graph-add").contains(evt.target)) {
    closeGraphAddMenu();
  }
});

tissueSummary.addEventListener("click", (evt) => {
  evt.preventDefault();
  if (tissueMenu.hidden) openTissueMenu();
  else closeTissueMenu();
});
tissueSummary.addEventListener("keydown", (evt) => {
  if (evt.key === "ArrowDown" || evt.key === "ArrowUp") {
    evt.preventDefault();
    if (tissueMenu.hidden) openTissueMenu();
  }
});
tissueFilter.addEventListener("input", renderTissueMenu);
tissueFilter.addEventListener("keydown", (evt) => {
  if (tissueMenu.hidden && (evt.key === "ArrowDown" || evt.key === "ArrowUp")) {
    openTissueMenu();
  }
  if (evt.key === "Escape") {
    evt.preventDefault();
    closeTissueMenu();
    tissueSummary.focus();
    return;
  }
  if (tissueMenu.hidden) return;
  if (evt.key === "ArrowDown") {
    evt.preventDefault();
    tissueMenuIndex = Math.min(tissueMenuItems.length - 1, tissueMenuIndex + 1);
    highlightTissueItem();
  } else if (evt.key === "ArrowUp") {
    evt.preventDefault();
    tissueMenuIndex = Math.max(0, tissueMenuIndex - 1);
    highlightTissueItem();
  } else if (evt.key === "Enter") {
    evt.preventDefault();
    if (tissueMenuItems.length === 0) return;
    const idx = tissueMenuIndex >= 0 ? tissueMenuIndex : 0;
    const item = tissueMenuItems[idx];
    if (item) toggleTissue(item.dataset.tissue);
  } else if (evt.key === " ") {
    if (tissueMenuIndex >= 0 && tissueMenuItems[tissueMenuIndex]) {
      evt.preventDefault();
      toggleTissue(tissueMenuItems[tissueMenuIndex].dataset.tissue);
    }
  }
});

function splitExamplePathways(value) {
  return (value || "")
    .split(",")
    .map((s) => s.trim())
    .filter(Boolean);
}

document.querySelectorAll(".example-chip[data-pathway], .example-chip[data-pathways]").forEach((btn) => {
  btn.addEventListener("click", () => {
    const groups = [];
    const first = btn.dataset.pathways
      ? splitExamplePathways(btn.dataset.pathways)
      : btn.dataset.pathway
        ? [btn.dataset.pathway]
        : [];
    if (first.length) groups.push({ name: btn.dataset.groupName || "", pathways: first });
    const second = splitExamplePathways(btn.dataset.pathways2 || btn.dataset["pathways-2"]);
    if (second.length) {
      groups.push({
        name: btn.dataset.groupName2 || btn.dataset["groupName-2"] || "",
        pathways: second,
      });
    }
    const topN = btn.dataset.topN === undefined || btn.dataset.topN === "" ? null : Number(btn.dataset.topN);
    const corr = btn.dataset.corr === undefined || btn.dataset.corr === "" ? null : Number(btn.dataset.corr);
    runExample({ groups, collection: btn.dataset.collection, topN, corr });
  });
});
const emptyFromGenes = el("emptyFromGenes");
if (emptyFromGenes) {
  emptyFromGenes.addEventListener("click", () => {
    if (groups[0]) openFromGenes(groups[0]);
  });
}

setHomeLink();
resetGroups();
init();
