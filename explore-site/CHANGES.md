# Changes

## 2026-09-22 — Multi-pathway, multi-group Explore

Adds the ability to select more than one pathway, organized into one or more
named/colored **groups**, and compare them in a single combined network
(design doc: `/home/sumeet/.claude/plans/i-have-a-local-serialized-whisper.md`
§9). Previously Explore only supported a single seed pathway.

**How it behaves**: the default state is exactly one group ("Group 1") with
zero pathways — functionally identical to the old single-pathway UI once you
add one pathway to it. Click "+ Add group" to create additional named,
colored groups; add pathways to whichever group's input you're typing in.
Hitting GO runs one combined topN query over the union of every group's
pathways, then colors each node by which group it came from (or blue/orange
for a newly-discovered neighbor, unchanged from before).

### Why the query algorithm didn't need to become more complex

In the legacy PCxN app, "up-regulated" and "down-regulated" gene sets were
**never queried separately** — they were concatenated into one combined seed
list before being sent to the same two-stage SQL query
(`pcxn-working/pcxn_reference_spec.md` §3). The up/down split only ever
affected *frontend coloring*, never the query itself. Generalizing "2 fixed
groups" to "N user-named groups" follows that same precedent exactly: the
query algorithm treats all groups' pathways as one flat set; group identity
is purely a rendering concern (node color + a legend).

### `explore-site/app.js`

- **`runExploreQuery(queryPathwayId, queryPathwayName, ...)` → `runExploreQuery(querySet, topN, corrCutOff, pCutOff)`.**
  `querySet` is a deduplicated array of `{id, name}` — the flattened union of
  every group's pathways. A query set of size 1 degenerates exactly to the
  old single-pathway behavior, so there's only one code path, not two.
  - **Stage 1** (neighbor ranking) now fetches *every* query member's own
    slice in parallel (previously just one), and scores each candidate
    outside the query set by `max(abs(correlation))` across whichever
    members it's connected to — reproducing the legacy SQL's
    `UNION ... GROUP BY ... MAX(abs(correlation))` query with a single-branch
    scan, which works because the ETL's slices are already fully
    bidirectional (see the code comment in `runExploreQuery`).
  - **Stage 2** (induced subgraph) now includes edges from *every* query
    member's slice, not just one. Query-query edges among the seeds
    themselves appear automatically (both endpoints are already in the
    allowed set) — no special-casing needed.
  - **New failure category**: `failedQueryNames`, distinct from the existing
    `failedNames` (failed neighbor fetches). A failed *query*-member fetch is
    worse — it silently removes that member's votes from every candidate's
    ranking score, potentially hiding a real neighbor entirely, not just
    dropping one already-known node. If *every* query member's fetch fails,
    the function throws (there's nothing left to query) rather than
    returning a hollow, misleadingly-empty result.
- **New groups state and UI** (`groups` array, `makeGroup`/`addGroupRow`/
  `removeGroup`/`renderChips`/`resetGroups`/`loadGroupsFromPayload`,
  `GROUP_COLORS` palette): each group is `{id, name, color, pathways, rowEl,
  combobox}`. The existing single-select combobox
  (`openPathwayMenu`/`closePathwayMenu`/`selectPathway`/`highlightPathwayItem`
  from the prior revision) is generalized into
  `openGroupMenu`/`closeGroupMenu`/`selectPathwayForGroup`/
  `highlightGroupItem`/`attachGroupComboboxEvents`, parametrized per group
  instead of operating on one global `pathwayInput`/`pathwayMenu`. Adding a
  pathway to a group's chip list only touches that group's own `.chip-list`
  DOM (never the whole groups container), so selecting a pathway keeps the
  input focused and ready to type the next name immediately.
- **Cross-group exclusion**: `filteredPathways()` now also excludes any name
  already chipped *in any group*, not just the current one — a pathway can
  only belong to one group at a time, since it can only render in one color.
- **`nodeColor()`** kept its existing precedence chain and only had its first
  branch generalized: previously `role === "query" → black`; now
  `groupByName.has(name) → that group's color`. The existing
  Pathprint-orange-accent-for-neighbors and default-blue-for-neighbors
  branches are untouched — they apply to a disjoint node category (discovered
  neighbors, never query-set members), so there's no conflict between group
  coloring and the collection accent.
- **Dynamic legend** (`renderLegend()`): the old static legend sentence
  ("Black = query pathway; blue = neighbor; orange = Pathprint...") is
  replaced by a `<ul>` generated after each query — one swatch+label per
  active group (from the snapshot captured in `lastQuery.groups` at
  query-time, not live state), plus fixed "Neighbor" / "Neighbor (Pathprint)"
  entries. A static sentence can't describe an arbitrary number of
  user-named groups.
- **URL state / shareable links**: `?pathway=<one name>` is replaced by
  `?groups=<JSON>`, where the JSON is
  `[{name, color, pathways: [name, ...]}, ...]`. Pathway names in this
  dataset can contain colons, commas, and parentheses, so groups are
  round-tripped as one JSON blob (percent-encoded automatically by
  `URLSearchParams`) rather than a hand-rolled delimiter scheme.
  `loadGroupsFromPayload()` reconstructs the groups UI from a parsed payload,
  falling back to a single empty default group on missing/malformed input
  (wrapped in `try/catch`, matching the existing defensive style).
  `tissue`/`collection`/`topN`/`corr`/`p` and the existing
  `data`/`dataBase`/`base` routing passthrough are unaffected.
- **Fetch-count guardrail**: total fetches per query = (pathways selected
  across all groups) + `topN`, all parallel. `MAX_TOTAL_FETCHES = 60` is
  enforced as a blocking validation message in `onGo()` (not silent
  truncation), and a live counter (`updateFetchCounter()`, wired to chip
  add/remove and the `topN` input) shows the cost before you click GO. This
  bounds both fetch fan-out and how dense/unreadable the resulting network
  can get, since stage 2 is deliberately unfiltered by the cutoffs (matching
  the legacy app's spec).
- `onGo`/`onReset`/`runExample`/`renderResultMeta`/`exportCsv` generalized
  from "one pathway" to "N pathways across M groups" (status text, reset
  clears to a single default group, example chips add into the first group,
  result summary and CSV filename reflect group/pathway counts instead of one
  name).

### `explore-site/index.html`

- Replaced the single `.pathway-field` combobox with a `groups-field`
  section: a "+ Add group" button, an (initially empty, JS-populated)
  `#groupsContainer`, and a `#fetchCounter` hint line.
- Replaced the static legend `<p>` with `<div id="legend">`, populated by
  `renderLegend()`; added a separate static `.legend-note` for the
  correlation-sign edge-color explanation (unchanged, not group-specific).
- Updated the empty-state helper text and the About dialog's description to
  describe groups instead of a single seed pathway.

### `explore-site/style.css`

- Added `.groups-header`, `.group-row`, `.group-row-header`, `.group-swatch`,
  `.group-name-input`, `.remove-group-btn`, `.chip-list`, `.chip`,
  `.chip-remove`, `.hint`/`.hint.warn`, `.legend-list`, `.legend-swatch`,
  `.legend-note`. No existing rules were changed (`.combobox`/
  `.combobox-menu`/`.combo-name`/`.combo-meta` are class-based already, so
  they apply to every group's combobox instance without modification).

### Verification performed

No live browser was available in this sandbox (same limitation as earlier
work on this project), so verification focused on what could actually be
executed:

1. **Syntax check** — `node --check app.js` (Node 20) after every edit.
2. **Algorithm cross-validation against real data** — copied the exact
   `runExploreQuery` logic into a standalone Node script, ran it against the
   real local data server for a 3-pathway query set spanning two collections
   (`Pathway.KEGG_CYSTEINE_AND_METHIONINE_METABOLISM`,
   `Pathway.BIOCARTA_RHO_PATHWAY`, `AICHI001_OCILY10_24H:J07_placenta_down`),
   and independently cross-checked both stage 1 (neighbor ranking) and stage 2
   (induced subgraph) against direct `polars` queries of
   `output/combined/pcxn_master.parquet`. Result: **exact match** on the
   neighbor list and the full edge set (153/153 edges identical), including
   all 3 query-query edges appearing automatically with no special-casing.
3. **`topN=0` multi-seed edge case** — confirmed it correctly degenerates to
   the pure query×query induced subgraph (exactly `C(3,2) = 3` edges, zero
   discovered neighbors), matching the legacy spec's documented behavior.
4. **Not verified**: actual visual rendering (group swatch colors, chip UI
   interactions, legend layout, DataTables/CanvasXpress/Cytoscape rendering)
   — please exercise this in a real browser and report anything that looks
   wrong.
