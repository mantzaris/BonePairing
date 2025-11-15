/* eslint-disable no-undef */
'use strict';

/**
 * Bones Association Explorer (two-pane)
 * - Left: expert/consensus network (as before)
 * - Right: overlay chemistry edges when an expert rating for a pair does not exist
 *   Data sources:
 *     - JSON:  data/Mathew-data-new.json
 *     - CSV:   data/blindtestaveragedspectra.csv  (falls back to data/blinetestaveragedspectra.csv)
 *
 * Chemistry policy (right pane):
 *   - Always connect isolates (degree 0 in the left/base view) to their nearest chemical neighbor.
 *   - Optionally, connect each non-isolate to its single best chemical neighbor if:
 *       (a) that pair has no expert rating in the current view (regardless of threshold),
 *       (b) cosine similarity >= CHEM_MIN_COS, and
 *       (c) that chemistry edge is not already present.
 *   - Chemistry score = 2 * cosineSimilarity => in [-2, 2] to match expert scale.
 */

// ------------------------------ UI elements ------------------------------

const expertSelect = document.getElementById('expertSelect');
const aggSelect = document.getElementById('aggSelect');
const thresholdSlider = document.getElementById('thresholdSlider');
const thresholdValue = document.getElementById('thresholdValue');
const showNegativeCheckbox = document.getElementById('showNegative');
const mustLinkSlider = document.getElementById('mustLink');
const mustLinkValue = document.getElementById('mustLinkValue');

const pairsTableBody = document.querySelector('#pairsTable tbody');
const clustersDiv = document.getElementById('clusters');

// ------------------------------ Layout constants ------------------------------

const FORCE_BASE_DISTANCE = 80;
const FORCE_ALPHA = 0.6;
const CHARGE_STRENGTH = -120;

// Chemistry overlay knobs (sensible defaults for 100 bones)
const CHEM_MIN_COS = 0.88;               // require fairly strong cosine similarity
const CHEM_ADD_FOR_NON_ISOLATES = true;  // also add one best chem neighbor for non-isolates
const CHEM_K_FOR_NON_ISOLATES = 1;       // at most one chem edge per non-isolate

// ------------------------------ Utility functions ------------------------------

function pairKey(a, b) {
  return a < b ? a + '|' + b : b + '|' + a;
}

function groupBy(array, keyFn) {
  const out = new Map();
  for (const item of array) {
    const k = keyFn(item);
    if (!out.has(k)) out.set(k, []);
    out.get(k).push(item);
  }
  return out;
}

function mean(values) {
  if (values.length === 0) return 0;
  return values.reduce((a, b) => a + b, 0) / values.length;
}

function median(values) {
  if (values.length === 0) return 0;
  const sorted = values.slice().sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[mid] : 0.5 * (sorted[mid - 1] + sorted[mid]);
}

function trimmedMean(values, frac = 0.1) {
  if (values.length === 0) return 0;
  const sorted = values.slice().sort((a, b) => a - b);
  const cut = Math.floor(frac * sorted.length);
  const trimmed = sorted.slice(cut, sorted.length - cut);
  return mean(trimmed.length ? trimmed : sorted);
}

function sign(x) {
  return x > 0 ? 1 : x < 0 ? -1 : 0;
}

function weightedConsensus(edgeMapsByExpert, pairKeys) {
  const expertIds = Array.from(edgeMapsByExpert.keys());
  const perExpertScores = new Map();
  const weights = new Map();

  for (const expertId of expertIds) {
    const m = edgeMapsByExpert.get(expertId);
    const view = new Map();
    for (const k of pairKeys) {
      if (m.has(k)) view.set(k, m.get(k));
    }
    perExpertScores.set(expertId, view);
  }

  for (const e of expertIds) {
    let agree = 0;
    let total = 0;
    for (const k of pairKeys) {
      let sum = 0;
      let count = 0;
      for (const f of expertIds) {
        if (f === e) continue;
        const mv = perExpertScores.get(f).get(k);
        if (mv !== undefined) {
          sum += mv;
          count += 1;
        }
      }
      if (count === 0) continue;
      const looSign = sign(sum / count);
      const ev = perExpertScores.get(e).get(k);
      if (ev === undefined) continue;
      total += 1;
      if (sign(ev) === looSign) agree += 1;
    }
    weights.set(e, total > 0 ? agree / total : 1);
  }

  const wsum = Array.from(weights.values()).reduce((a, b) => a + b, 0) || 1;
  for (const e of expertIds) {
    weights.set(e, weights.get(e) / wsum);
  }

  const scores = new Map();
  for (const k of pairKeys) {
    let s = 0;
    for (const e of expertIds) {
      const m = edgeMapsByExpert.get(e);
      if (m.has(k)) s += weights.get(e) * m.get(k);
    }
    scores.set(k, s);
  }
  return { scores, weights };
}

function distanceFromScore(score, base = FORCE_BASE_DISTANCE, alpha = FORCE_ALPHA) {
  return base * Math.exp(-alpha * score);
}

// ------------------------------ Two D3 scenes with zoom ------------------------------

function makeScene(svgSelector) {
  const svg = d3.select(svgSelector);

  const zoomPane = svg.append('rect')
    .attr('class', 'zoom-pane')
    .attr('fill', 'transparent')
    .attr('pointer-events', 'all');

  const group = svg.append('g');
  const linkLayer = group.append('g');
  const labelLayer = group.append('g');
  const nodeLayer = group.append('g');

  let currentTransform = d3.zoomIdentity;
  const zoomBehavior = d3.zoom()
    .scaleExtent([0.25, 4])
    .on('zoom', (event) => {
      currentTransform = event.transform;
      group.attr('transform', currentTransform);
    });

  return {
    svg, zoomPane, group, linkLayer, labelLayer, nodeLayer,
    zoomBehavior, get currentTransform() { return currentTransform; },
    set currentTransform(t) { currentTransform = t; },
    simulation: null
  };
}

const leftScene = makeScene('#graphLeft');
const rightScene = makeScene('#graphRight');

function sizeSceneToContainer(scene, elementId) {
  const el = document.getElementById(elementId);
  const rect = el.getBoundingClientRect();

  scene.svg.attr('width', rect.width).attr('height', rect.height);

  scene.zoomPane
    .attr('x', 0)
    .attr('y', 0)
    .attr('width', rect.width)
    .attr('height', rect.height);

  scene.zoomBehavior.extent([[0, 0], [rect.width, rect.height]]);
  scene.svg.call(scene.zoomBehavior).call(scene.zoomBehavior.transform, scene.currentTransform);
}

function sizeAllScenes() {
  sizeSceneToContainer(leftScene, 'graphLeft');
  sizeSceneToContainer(rightScene, 'graphRight');
}

// ------------------------------ Data loading and indices ------------------------------

/**
 * @typedef {Object} Dataset
 * @property {string} collectionId
 * @property {{id: string, type?: string}[]} bones
 * @property {{id: string, name: string}[]} experts
 * @property {{expertId: string, boneA: string, boneB: string, score: number}[]} ratings
 * @property {{vectors?: Record<string, number[]>, wavelengths?: number[]}} [chemistry]
 */

let dataset = null;
let bones = [];     // array of bone ids as strings
let expertsById = new Map(); // expertId -> expert object
let expertEdgeMaps = new Map(); // expertId -> Map(pairKey -> score)

// Parse JSON to indices
function buildIndices(raw) {
  const bonesSet = new Set(raw.bones.map(b => b.id));
  for (const r of raw.ratings) {
    bonesSet.add(r.boneA);
    bonesSet.add(r.boneB);
  }
  bones = Array.from(bonesSet).sort();

  expertsById = new Map(raw.experts.map(e => [e.id, e]));

  // Canonicalize ratings
  const canonicalRatings = raw.ratings.map(r => {
    const a = r.boneA;
    const b = r.boneB;
    const [x, y] = a < b ? [a, b] : [b, a];
    return { expertId: r.expertId, boneA: x, boneB: y, score: r.score, key: pairKey(x, y) };
  });

  const byExpert = groupBy(canonicalRatings, row => row.expertId);
  expertEdgeMaps = new Map();
  for (const [expertId, rows] of byExpert.entries()) {
    const byKey = groupBy(rows, row => row.key);
    const edgeMap = new Map();
    for (const [k, arr] of byKey.entries()) {
      edgeMap.set(k, mean(arr.map(a => a.score)));
    }
    expertEdgeMaps.set(expertId, edgeMap);
  }

  populateExpertSelector(raw.experts);
}

function populateExpertSelector(experts) {
  for (const e of experts) {
    const opt = document.createElement('option');
    opt.value = e.id;
    opt.textContent = e.name;
    expertSelect.appendChild(opt);
  }
}

// ------------------------------ Chemistry CSV loading ------------------------------

/**
 * Try multiple candidate paths for the spectra CSV, returning parsed rows or null.
 */
async function loadFirstAvailableCsv(candidates) {
  for (const path of candidates) {
    try {
      const res = await fetch(path, { cache: 'no-store' });
      if (res && res.ok) {
        const text = await res.text();
        const rows = d3.csvParse(text, d3.autoType);
        console.info(`[chem] loaded spectra: ${rows.length} rows from ${path}`);
        return { rows, path };
      }
    } catch (err) {
      // try next candidate
    }
  }
  return { rows: null, path: null };
}

/**
 * Build dataset.chemistry.vectors (L2-normalized Float64Array) and wavelengths[].
 * Accepts the array of row objects returned by d3.csvParse(autoType).
 */
function buildChemistryFromCsvRows(rows) {
  if (!rows || rows.length === 0) return;

  // Detect the id column name robustly: "Bone #", "Bone#", "Bone", "bone", "id"
  const sample = rows[0] || {};
  const keys = Object.keys(sample);
  const idKey = keys.find(k => /bone\s*#?$/i.test(k)) || keys.find(k => /^id$/i.test(k)) || keys[0];

  // Spectral columns = everything except the id column; sort numerically by wavelength
  const spectralKeys = keys.filter(k => k !== idKey)
    .filter(k => Number.isFinite(+k)) // keep numeric wavelength columns only
    .sort((a, b) => (+a) - (+b));

  const vectors = {};
  const wavelengths = spectralKeys.map(k => +k);

  function normalizeL2(arr) {
    let sumsq = 0;
    for (let i = 0; i < arr.length; i++) sumsq += arr[i] * arr[i];
    const norm = Math.sqrt(sumsq);
    if (!isFinite(norm) || norm === 0) return null;
    const out = new Float64Array(arr.length);
    const inv = 1 / norm;
    for (let i = 0; i < arr.length; i++) out[i] = arr[i] * inv;
    return out;
  }

  for (const row of rows) {
    const id = String(row[idKey]);
    // Some CSVs might have blank lines; guard against undefined ids
    if (!id || id === 'undefined' || id === 'null') continue;
    const raw = new Float64Array(spectralKeys.length);
    for (let i = 0; i < spectralKeys.length; i++) {
      const v = row[spectralKeys[i]];
      raw[i] = Number.isFinite(v) ? +v : 0;
    }
    const vec = normalizeL2(raw);
    if (vec) vectors[id] = vec;
  }

  if (!dataset.chemistry) dataset.chemistry = {};
  dataset.chemistry.vectors = vectors;
  dataset.chemistry.wavelengths = wavelengths;

  console.info(`[chem] vectors built for ${Object.keys(vectors).length} bones; dims=${wavelengths.length}`);
}

function getChemVector(boneId) {
  if (dataset && dataset.chemistry && dataset.chemistry.vectors) {
    const v = dataset.chemistry.vectors[boneId];
    return v || null;
  }
  return null;
}

function cosineSimilarity(vecA, vecB) {
  // Vectors are L2-normalized; cosine similarity is their dot product
  let s = 0;
  const n = vecA.length;
  for (let i = 0; i < n; i++) s += vecA[i] * vecB[i];
  return s; // in [-1, 1]
}

// ------------------------------ Build edges for current view ------------------------------

/**
 * Returns:
 *  - edges: thresholded edges for the left/base view
 *  - fromLabel: label for the table
 *  - ratedPairKeys: Set of pair keys that DO have an expert rating in the current view
 */
function buildEdgesForCurrentView() {
  const view = expertSelect.value;
  const agg = aggSelect.value;
  const threshold = parseFloat(thresholdSlider.value);
  const showNegative = showNegativeCheckbox.checked;

  thresholdValue.textContent = threshold.toFixed(1);
  aggSelect.disabled = view !== 'consensus';

  let edgeMap = new Map();
  let fromLabel = 'Consensus';

  if (view === 'consensus') {
    // union of all pairs across experts
    const allKeys = new Set();
    for (const m of expertEdgeMaps.values()) for (const k of m.keys()) allKeys.add(k);
    const allPairKeys = Array.from(allKeys);
    const { scores } = weightedConsensus(expertEdgeMaps, allPairKeys);
    edgeMap = scores;
  } else {
    edgeMap = expertEdgeMaps.get(view) || new Map();
    const e = expertsById.get(view);
    fromLabel = e ? e.name : view;
  }

  // ratedPairKeys = all keys present in the chosen view (regardless of threshold)
  const ratedPairKeys = new Set(edgeMap.keys());

  const edges = [];
  for (const [k, s] of edgeMap.entries()) {
    if (Math.abs(s) < threshold) continue;
    if (!showNegative && s < 0) continue;
    const [a, b] = k.split('|');
    edges.push({ source: a, target: b, score: s });
  }
  return { edges, fromLabel, ratedPairKeys };
}

// ------------------------------ Chemistry overlay edges ------------------------------

/**
 * Build chemistry edges for: isolates (always) and optionally for non-isolates.
 * Only add when that pair has NO expert rating (use ratedPairKeys).
 * Uses cosine similarity on normalized spectra; chemistry score = 2 * cosine.
 */
function buildChemAugmentEdges(baseEdges, ratedPairKeys) {
  // degree from baseEdges (strings before simulation mutates)
  const degree = new Map(bones.map(id => [id, 0]));
  const drawnPairKeys = new Set();
  for (const e of baseEdges) {
    const a = typeof e.source === 'string' ? e.source : e.source.id;
    const b = typeof e.target === 'string' ? e.target : e.target.id;
    degree.set(a, (degree.get(a) || 0) + 1);
    degree.set(b, (degree.get(b) || 0) + 1);
    drawnPairKeys.add(pairKey(a, b));
  }

  // Precompute best chemical neighbors for all bones
  const chemNeighbors = new Map(); // id -> sorted array [{id, cos}, ...] best first
  const idsWithVectors = bones.filter(id => !!getChemVector(id));
  if (idsWithVectors.length < 2) {
    return { edges: [], status: 'No or insufficient chemistry vectors' };
  }

  for (const id of idsWithVectors) {
    const vA = getChemVector(id);
    let best = null;
    let bestCos = -Infinity;
    let best2 = null;
    let best2Cos = -Infinity;

    for (const other of idsWithVectors) {
      if (other === id) continue;
      const vB = getChemVector(other);
      const cos = cosineSimilarity(vA, vB);
      if (cos > bestCos) {
        best2 = best; best2Cos = bestCos;
        best = { id: other, cos };
        bestCos = cos;
      } else if (cos > best2Cos) {
        best2 = { id: other, cos };
        best2Cos = cos;
      }
    }
    const arr = [];
    if (best) arr.push(best);
    if (best2) arr.push(best2);
    chemNeighbors.set(id, arr);
  }

  const chemEdges = [];
  const appendedPairKeys = new Set();

  function maybeAddChemEdge(a, b, cos) {
    const k = pairKey(a, b);
    if (ratedPairKeys.has(k)) return;        // expert already rated this pair (even if filtered out)
    if (drawnPairKeys.has(k)) return;        // already drawn by base edges
    if (appendedPairKeys.has(k)) return;     // already added by chemistry
    if (!(cos >= CHEM_MIN_COS)) return;      // not strong enough

    const scoreChem = 2 * cos; // map cosine [-1,1] -> [-2,2]
    chemEdges.push({
      source: a,
      target: b,
      score: scoreChem,
      isChem: true,
      chemCos: cos
    });
    appendedPairKeys.add(k);
  }

  // 1) Always connect isolates to their nearest chemical neighbor (if eligible)
  const isolates = bones.filter(id => (degree.get(id) || 0) === 0);
  for (const id of isolates) {
    const neigh = chemNeighbors.get(id);
    if (!neigh || neigh.length === 0) continue;
    // try best, then second-best, to avoid conflicts with rated pairs
    for (const cand of neigh) {
      maybeAddChemEdge(id, cand.id, cand.cos);
      if (appendedPairKeys.has(pairKey(id, cand.id))) break;
    }
  }

  // 2) Optionally add one chem neighbor for non-isolates (to surface strong non-rated affinities)
  if (CHEM_ADD_FOR_NON_ISOLATES) {
    const nonIsolates = bones.filter(id => (degree.get(id) || 0) > 0);
    for (const id of nonIsolates) {
      const neigh = chemNeighbors.get(id);
      if (!neigh || neigh.length === 0) continue;
      let added = 0;
      for (const cand of neigh) {
        maybeAddChemEdge(id, cand.id, cand.cos);
        if (appendedPairKeys.has(pairKey(id, cand.id))) {
          added += 1;
          if (added >= CHEM_K_FOR_NON_ISOLATES) break;
        }
      }
    }
  }

  const status = chemEdges.length
    ? `Chem overlay: ${chemEdges.length} edges (min cos ${CHEM_MIN_COS.toFixed(2)})`
    : 'No chemistry edges qualified (threshold too high or all pairs rated)';
  return { edges: chemEdges, status };
}

// ------------------------------ Rendering ------------------------------

function drawScene(scene, nodes, edgeList) {
  if (scene.simulation) scene.simulation.stop();

  const widthScale = d3.scaleLinear().domain([0, 2]).range([1, 4]);

  const linkKey = (d) => {
    const a = d.source.id || d.source;
    const b = d.target.id || d.target;
    return (a < b ? a + '|' + b : b + '|' + a) + (d.isChem ? '|chem' : '|human');
  };

  const linksSel = scene.linkLayer.selectAll('line').data(edgeList, linkKey);
  linksSel.exit().remove();
  const linksEnter = linksSel.enter().append('line');
  const links = linksEnter.merge(linksSel)
    .attr('stroke', d => d.isChem ? '#8e24aa' : (d.score >= 0 ? '#1e9e45' : '#c62828'))
    .attr('stroke-width', d => d.isChem ? 2 : widthScale(Math.min(2, Math.abs(d.score))))
    .attr('stroke-opacity', d => d.isChem ? 0.8 : 0.9)
    .attr('stroke-dasharray', d => d.isChem ? '2,2' : (d.score < 0 ? '4,3' : '0'));

  const labelSel = scene.labelLayer.selectAll('text').data(edgeList, linkKey);
  labelSel.exit().remove();
  const labelEnter = labelSel.enter().append('text').attr('class', 'edge-label');
  const labels = labelEnter.merge(labelSel)
    .text(d => d.isChem ? ('chem s=' + (d.chemCos !== undefined ? d.chemCos.toFixed(2) : '')) : d.score.toFixed(2))
    .attr('font-size', 10)
    .attr('fill', d => d.isChem ? '#6a1b9a' : '#111827')
    .attr('opacity', 0.9);

  const nodeSel = scene.nodeLayer.selectAll('g.node').data(nodes, d => d.id);
  nodeSel.exit().remove();
  const nodeEnter = nodeSel.enter().append('g').attr('class', 'node');
  nodeEnter.append('circle')
    .attr('r', 10)
    .attr('fill', '#3b82f6')
    .attr('stroke', '#1e40af')
    .attr('stroke-width', 1.2);
  nodeEnter.append('text')
    .attr('y', -14)
    .attr('text-anchor', 'middle')
    .attr('font-size', 11)
    .text(d => d.id);

  const nodesMerged = nodeEnter.merge(nodeSel);

  nodesMerged.on('mouseover', function (_evt, d) {
    const id = d.id;
    links.attr('stroke-opacity', l => {
      const a = l.source.id || l.source; const b = l.target.id || l.target;
      return (a === id || b === id) ? 1.0 : 0.15;
    });
    nodesMerged.selectAll('circle').attr('opacity', nd => nd.id === id ? 1.0 : 0.6);
  }).on('mouseout', function () {
    links.attr('stroke-opacity', l => l.isChem ? 0.8 : 0.9);
    nodesMerged.selectAll('circle').attr('opacity', 1.0);
  });

  const svgW = +scene.svg.attr('width');
  const svgH = +scene.svg.attr('height');
  scene.simulation = d3.forceSimulation(nodes)
    .force('charge', d3.forceManyBody().strength(CHARGE_STRENGTH))
    .force('center', d3.forceCenter(svgW / 2, svgH / 2))
    .force('link', d3.forceLink(edgeList)
      .id(d => d.id)
      .distance(d => {
        if (d.isChem) return distanceFromScore(d.score, FORCE_BASE_DISTANCE * 1.1, FORCE_ALPHA * 0.7);
        return distanceFromScore(d.score);
      })
      .strength(d => d.isChem ? 0.35 : 0.7))
    .on('tick', () => {
      links
        .attr('x1', d => d.source.x)
        .attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x)
        .attr('y2', d => d.target.y);

      labels
        .attr('x', d => (d.source.x + d.target.x) / 2)
        .attr('y', d => (d.source.y + d.target.y) / 2);

      nodesMerged.attr('transform', d => 'translate(' + d.x + ',' + d.y + ')');
    });
}

// Show a small status message inside the right pane
function updateChemStatus(scene, text) {
  const sel = scene.group.selectAll('text.chem-status').data(text ? [text] : []);
  sel.exit().remove();
  sel.enter().append('text')
    .attr('class', 'chem-status')
    .attr('x', 12)
    .attr('y', 20)
    .attr('font-size', 12)
    .attr('fill', '#6b7280')
    .attr('opacity', 0.9)
    .text(d => d)
    .merge(sel)
    .text(d => d);
}

// ------------------------------ Panels and clustering (left view) ------------------------------

function updatePairsTable(edges, fromLabel) {
  pairsTableBody.innerHTML = '';
  const sorted = edges.slice().sort((a, b) => Math.abs(b.score) - Math.abs(a.score)).slice(0, 30);
  for (const e of sorted) {
    const tr = document.createElement('tr');
    tr.innerHTML =
      '<td>' + (e.source.id || e.source) + '</td>' +
      '<td>' + (e.target.id || e.target) + '</td>' +
      '<td>' + e.score.toFixed(2) + '</td>' +
      '<td>' + fromLabel + '</td>';
    pairsTableBody.appendChild(tr);
  }
}

function buildConsensus(method) {
  // This is only used by updateClusters which needs scores >= threshold
  const keys = new Set();
  for (const m of expertEdgeMaps.values()) for (const k of m.keys()) keys.add(k);
  const allPairKeys = Array.from(keys);

  if (method === 'weighted') {
    return weightedConsensus(expertEdgeMaps, allPairKeys);
  }

  const scores = new Map();
  for (const k of allPairKeys) {
    const vals = [];
    for (const m of expertEdgeMaps.values()) {
      if (m.has(k)) vals.push(m.get(k));
    }
    let v = 0;
    if (method === 'mean') v = mean(vals);
    else if (method === 'median') v = median(vals);
    else if (method === 'trimmed') v = trimmedMean(vals, 0.1);
    scores.set(k, v);
  }
  return { scores, weights: new Map() };
}

function updateClusters(view, agg) {
  const tPos = parseFloat(mustLinkSlider.value);
  mustLinkValue.textContent = tPos.toFixed(1);

  let edgeMap = new Map();
  if (view === 'consensus') {
    const { scores } = buildConsensus(agg);
    edgeMap = scores;
  } else {
    edgeMap = expertEdgeMaps.get(view) || new Map();
  }

  const adj = new Map(bones.map(b => [b, new Set()]));
  for (const [k, s] of edgeMap.entries()) {
    if (s >= tPos) {
      const [a, b] = k.split('|');
      adj.get(a).add(b);
      adj.get(b).add(a);
    }
  }

  const seen = new Set();
  const clusters = [];
  for (const b of bones) {
    if (seen.has(b)) continue;
    const comp = [];
    const stack = [b];
    seen.add(b);
    while (stack.length) {
      const u = stack.pop();
      comp.push(u);
      for (const v of adj.get(u)) {
        if (!seen.has(v)) {
          seen.add(v);
          stack.push(v);
        }
      }
    }
    clusters.push(comp.sort());
  }
  clusters.sort((a, b) => b.length - a.length);

  clustersDiv.innerHTML = '';
  for (const comp of clusters) {
    const line = document.createElement('div');
    for (const id of comp) {
      const span = document.createElement('span');
      span.className = 'cluster-badge';
      span.textContent = id;
      line.appendChild(span);
    }
    clustersDiv.appendChild(line);
  }
}

// ------------------------------ Top-level redraw ------------------------------

function redrawAll() {
  sizeAllScenes();

  const { edges: baseEdges, fromLabel, ratedPairKeys } = buildEdgesForCurrentView();
  const nodes = bones.map(id => ({ id }));

  // Left: base network
  drawScene(leftScene, nodes.map(n => ({ ...n })), baseEdges.map(e => ({ ...e })));

  // Right: base + chemistry overlay
  const chemResult = buildChemAugmentEdges(baseEdges, ratedPairKeys);
  const chemEdges = chemResult.edges || [];
  const statusText = chemResult.status || null;

  const rightEdges = baseEdges.concat(chemEdges);
  drawScene(rightScene, nodes.map(n => ({ ...n })), rightEdges);
  updateChemStatus(rightScene, statusText);

  // Side panels reflect the left/base view
  updatePairsTable(baseEdges, fromLabel);
  updateClusters(expertSelect.value, aggSelect.value);

  // Console diagnostics
  const uniqueBaseVertices = new Set(baseEdges.flatMap(e => [e.source, e.target]).map(x => (typeof x === 'string' ? x : x.id)));
  console.info(`[chem] candidate bones with vectors: ${Object.keys((dataset.chemistry && dataset.chemistry.vectors) || {}).length}`);
  console.info(`[chem] isolates in base view: ${bones.length - uniqueBaseVertices.size}`);
  console.info(`[chem] chem edges drawn: ${chemEdges.length} (min cos ${CHEM_MIN_COS})`);
}

// ------------------------------ Wiring and init ------------------------------

async function init() {
  try {
    // Load JSON first
    const jsonRes = await fetch('data/Mathew-data-new.json', { cache: 'no-store' });
    const json = await jsonRes.json();
    dataset = /** @type {Dataset} */ (json);
    buildIndices(dataset);

    // Load spectra CSV (try two spellings)
    const csv = await loadFirstAvailableCsv([
      'data/blindtestaveragedspectra.csv',
      'data/blinetestaveragedspectra.csv'
    ]);
    if (csv.rows && csv.rows.length) {
      buildChemistryFromCsvRows(csv.rows);
    } else {
      console.warn('[chem] spectra CSV not found; chemistry overlay will be empty.');
      dataset.chemistry = dataset.chemistry || {};
      dataset.chemistry.vectors = dataset.chemistry.vectors || {};
    }

    sizeAllScenes();
    redrawAll();
  } catch (err) {
    console.error('Initialization failed:', err);
  }
}

window.addEventListener('resize', () => {
  sizeAllScenes();
  for (const scene of [leftScene, rightScene]) {
    if (scene.simulation) {
      const w = +scene.svg.attr('width'), h = +scene.svg.attr('height');
      scene.simulation.force('center', d3.forceCenter(w / 2, h / 2));
      scene.simulation.alpha(0.3).restart();
    }
  }
});

expertSelect.addEventListener('change', redrawAll);
aggSelect.addEventListener('change', redrawAll);
thresholdSlider.addEventListener('input', redrawAll);
showNegativeCheckbox.addEventListener('change', redrawAll);
mustLinkSlider.addEventListener('input', () => updateClusters(expertSelect.value, aggSelect.value));

init();
