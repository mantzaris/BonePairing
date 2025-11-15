/* eslint-disable no-undef */
'use strict';

/**
 * Bones Association Explorer MVP (two-pane)
 * - Left: original expert/consensus network
 * - Right: augmented network where isolated bones get a chemistry nearest-neighbor edge
 */

// ------------------------------ DOM Elements ------------------------------

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

// ------------------------------ Utility functions ------------------------------

function countChemVectors() {
  let c = 0;
  for (const id of bones) if (getChemVector(id)) c++;
  return c;
}

/** Add: small in-graph status message for the right pane */
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

// ------------------------------ Two D3 scene helpers ------------------------------

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

// ------------------------------ Data loading and indexing ------------------------------

/**
 * @typedef {Object} Dataset
 * @property {string} collectionId
 * @property {{id: string, type?: string, chemistry?: {vector?: number[]}}[]} bones
 * @property {{id: string, name: string}[]} experts
 * @property {{expertId: string, boneA: string, boneB: string, score: number}[]} ratings
 * @property {{vectors?: Record<string, number[]>}} chemistry
 */

let dataset = null;
let bones = [];     // string[]
let expertsById = new Map(); // id -> expert
let expertEdgeMaps = new Map(); // expertId -> Map(pairKey -> score)

/** Canonicalize rating entries and build indices. */
function buildIndices(raw) {
  const bonesSet = new Set(raw.bones.map(b => b.id));
  for (const r of raw.ratings) {
    bonesSet.add(r.boneA);
    bonesSet.add(r.boneB);
  }
  bones = Array.from(bonesSet).sort();

  expertsById = new Map(raw.experts.map(e => [e.id, e]));

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

function allPairKeys() {
  const keys = new Set();
  for (const m of expertEdgeMaps.values()) {
    for (const k of m.keys()) keys.add(k);
  }
  return Array.from(keys);
}

function buildConsensus(method) {
  const keys = allPairKeys();
  if (method === 'weighted') {
    return weightedConsensus(expertEdgeMaps, keys);
  }
  const scores = new Map();
  for (const k of keys) {
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

// ------------------------------ Chemistry helpers ------------------------------

function getChemVector(boneId) {
  if (dataset && dataset.chemistry && dataset.chemistry.vectors && dataset.chemistry.vectors[boneId]) {
    return dataset.chemistry.vectors[boneId];
  }
  const bone = dataset.bones.find(b => b.id === boneId);
  if (bone && bone.chemistry && Array.isArray(bone.chemistry.vector)) {
    return bone.chemistry.vector;
  }
  return null;
}

function euclidean(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i++) {
    const d = a[i] - b[i];
    s += d * d;
  }
  return Math.sqrt(s);
}


function computeMaxChemDistance() {
  const ids = bones.filter(id => getChemVector(id));
  if (ids.length < 2) return 0; // not enough data for a max distance

  let maxd = 0;
  for (let i = 0; i < ids.length; i++) {
    const ai = getChemVector(ids[i]);
    if (!ai) continue;
    for (let j = i + 1; j < ids.length; j++) {
      const bj = getChemVector(ids[j]);
      if (!bj) continue;
      const d = euclidean(ai, bj);
      if (d > maxd) maxd = d;
    }
  }
  return maxd; // may be 0 if all vectors are identical
}




function nearestByChem(boneId) {
  const a = getChemVector(boneId);
  if (!a) return null;
  let bestId = null;
  let bestDist = Infinity;
  for (const id of bones) {
    if (id === boneId) continue;
    const b = getChemVector(id);
    if (!b) continue;
    const d = euclidean(a, b);
    if (d < bestDist) {
      bestDist = d;
      bestId = id;
    }
  }
  return bestId ? { id: bestId, dist: bestDist } : null;
}


// ---------- Chemistry coverage helpers (fallback + diagnostics) ----------
function countChemCoverage() {
  let withVec = 0;
  for (const id of bones) if (getChemVector(id)) withVec += 1;
  return { withVec, total: bones.length };
}

/**
 * If the JSON lacks chemistry vectors, create deterministic placeholders
 * so the augmented graph still works. No-op when vectors already exist.
 */
function ensureChemVectorsIfMissing() {
  const cov = countChemCoverage();
  if (cov.withVec > 0) {
    console.info(`[chem] vectors present: ${cov.withVec}/${cov.total}`);
    return;
  }
  console.warn('[chem] no vectors found; synthesizing placeholders');

  // Simple seeded LCG so vectors are reproducible per bone id
  function lcg(seed) {
    let s = seed >>> 0;
    return () => {
      s = (1664525 * s + 1013904223) >>> 0;
      return s / 4294967296;
    };
  }

  for (const b of dataset.bones) {
    const idNum = parseInt(b.id, 10) || 0;
    const rnd = lcg(12345 + idNum * 2654435761);
    const raw = [rnd(), rnd(), rnd(), rnd(), rnd()];
    const sum = raw.reduce((a, c) => a + c, 0) || 1;
    const vec = raw.map(v => v / sum); // unit-sum 5D signature
    if (!b.chemistry) b.chemistry = {};
    b.chemistry.vector = vec;
  }
  const cov2 = countChemCoverage();
  console.info(`[chem] synthesized vectors: ${cov2.withVec}/${cov2.total}`);
}



// ------------------------------ Edge construction ------------------------------

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
    const { scores } = buildConsensus(agg);
    edgeMap = scores;
  } else {
    edgeMap = expertEdgeMaps.get(view) || new Map();
    const e = expertsById.get(view);
    fromLabel = e ? e.name : view;
  }

  const edges = [];
  for (const [k, s] of edgeMap.entries()) {
    if (Math.abs(s) < threshold) continue;
    if (!showNegative && s < 0) continue;
    const [a, b] = k.split('|');
    edges.push({ source: a, target: b, score: s });
  }
  return { edges, fromLabel };
}

function buildChemAugmentEdges(edgesBase) {
  // Build degree from base edges (strings before simulation mutates)
  const deg = new Map(bones.map(id => [id, 0]));
  const haveEdgeKey = new Set();
  for (const e of edgesBase) {
    const a = typeof e.source === 'string' ? e.source : e.source.id;
    const b = typeof e.target === 'string' ? e.target : e.target.id;
    deg.set(a, (deg.get(a) || 0) + 1);
    deg.set(b, (deg.get(b) || 0) + 1);
    haveEdgeKey.add(pairKey(a, b));
  }

  const isolates = bones.filter(id => (deg.get(id) || 0) === 0);

  // If no chemistry vectors exist at all, or we cannot compute a usable scale, bail early.
  const chemCount = countChemVectors();
  if (chemCount === 0) {
    console.warn('[chem] No chemistry vectors found in dataset.');
    return { edges: [], status: 'No chemistry vectors in data' };
  }

  if (isolates.length === 0) {
    // Not an error: you currently have no isolates under these filters.
    return { edges: [], status: 'No isolates under current filters' };
  }

  const dmax = computeMaxChemDistance();
  if (dmax <= 0) {
    console.warn('[chem] Chemistry vectors have zero spread; cannot form similarity.');
    return { edges: [], status: 'Chemistry vectors have zero spread' };
  }

  const chemEdges = [];
  const chemKeys = new Set();

  for (const id of isolates) {
    const nn = nearestByChem(id);
    if (!nn) continue; // this isolate lacks a vector or neighbors lack vectors
    const k = pairKey(id, nn.id);
    if (haveEdgeKey.has(k) || chemKeys.has(k)) continue;

    const sim = Math.max(0, Math.min(1, 1 - nn.dist / dmax)); // normalize to [0,1]
    const scoreChem = 0.8 + 1.2 * sim; // gentle pull in [0.8, 2.0]

    chemEdges.push({
      source: id,
      target: nn.id,
      score: scoreChem,
      isChem: true,
      chemSim: sim,
      chemDist: nn.dist
    });
    chemKeys.add(k);
  }
  return { edges: chemEdges, status: chemEdges.length ? null : 'Isolates lacked usable chemistry' };
}


// ------------------------------ Rendering ------------------------------

function drawScene(scene, nodes, edgeList) {
  // stop previous sim
  if (scene.simulation) scene.simulation.stop();

  // Scales
  const widthScale = d3.scaleLinear().domain([0, 2]).range([1, 4]);

  // Links
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
    .attr('stroke-opacity', d => d.isChem ? 0.75 : 0.9)
    .attr('stroke-dasharray', d => d.isChem ? '2,2' : (d.score < 0 ? '4,3' : '0'));

  // Edge labels
  const labelSel = scene.labelLayer.selectAll('text').data(edgeList, linkKey);
  labelSel.exit().remove();
  const labelEnter = labelSel.enter().append('text').attr('class', 'edge-label');
  const labels = labelEnter.merge(labelSel)
    .text(d => d.isChem ? ('chem ' + (d.chemSim !== undefined ? d.chemSim.toFixed(2) : '')) : d.score.toFixed(2))
    .attr('font-size', 10)
    .attr('fill', d => d.isChem ? '#6a1b9a' : '#111827')
    .attr('opacity', 0.85);

  // Nodes
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
    links.attr('stroke-opacity', l => l.isChem ? 0.75 : 0.9);
    nodesMerged.selectAll('circle').attr('opacity', 1.0);
  });

  // Force simulation
  const svgW = +scene.svg.attr('width');
  const svgH = +scene.svg.attr('height');
  scene.simulation = d3.forceSimulation(nodes)
    .force('charge', d3.forceManyBody().strength(CHARGE_STRENGTH))
    .force('center', d3.forceCenter(svgW / 2, svgH / 2))
    .force('link', d3.forceLink(edgeList)
      .id(d => d.id)
      .distance(d => {
        if (d.isChem) {
          // Slightly longer baseline and gentler decay for chemistry edges
          return distanceFromScore(d.score, FORCE_BASE_DISTANCE * 1.1, FORCE_ALPHA * 0.7);
        }
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

// ------------------------------ Top-level redraw ------------------------------



function redrawAll() {
  sizeAllScenes();

  const { edges: baseEdges, fromLabel } = buildEdgesForCurrentView();
  const nodes = bones.map(id => ({ id }));

  // Left: base network
  drawScene(leftScene, nodes.map(n => ({ ...n })), baseEdges.map(e => ({ ...e })));

  // Right: base + chemistry augmentation for isolates
  const chemResult = buildChemAugmentEdges(baseEdges);
  const chemEdges = Array.isArray(chemResult) ? chemResult : (chemResult.edges || []);
  const statusText = Array.isArray(chemResult) ? null : (chemResult.status || null);

  const rightEdges = baseEdges.concat(chemEdges);
  drawScene(rightScene, nodes.map(n => ({ ...n })), rightEdges);

  // Show why chemistry edges may not be visible
  updateChemStatus(rightScene, statusText);

  // Side panels reflect the left/base view
  updatePairsTable(baseEdges, fromLabel);
  updateClusters(expertSelect.value, aggSelect.value);

  // Helpful logs for quick inspection
  console.info(`[chem] vectors: ${countChemVectors()}, isolates: ${bones.length - new Set(baseEdges.flatMap(e => [e.source, e.target]).map(x => (typeof x === 'string' ? x : x.id))).size}, chemEdges drawn: ${chemEdges.length}`);
}





/** Update the top pairs table (left view). */
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

/** Very simple must-link clustering by connected components (left view). */
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

// ------------------------------ Wiring ------------------------------

async function init() {
  try {
    const res = await fetch('data/Mathew-data.json');
    const json = await res.json();
    dataset = /** @type {Dataset} */ (json);

    buildIndices(dataset);
    ensureChemVectorsIfMissing();
    sizeAllScenes();
    redrawAll();

  } catch (err) {
    console.error('Failed to load data:', err);
  }
}

window.addEventListener('resize', () => {
  sizeAllScenes();
  // recenter both simulations
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
