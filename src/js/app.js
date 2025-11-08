/* eslint-disable no-undef */
'use strict';

/**
 * Bones Association Explorer MVP
 * - Data model: bones, experts, ratings (Likert in [-2, 2])
 * - Views: per-expert and consensus
 * - Consensus: mean, median, trimmed mean, weighted by inter-rater agreement
 * - UI: threshold, negative edge toggle, must-link clustering
 *
 * Notes for contributors:
 * - Keep variable names descriptive. No slang or cryptic abbreviations.
 * - Treat scores as symmetric and store canonical pairs (a|b with a < b).
 */

// ------------------------------ DOM Elements ------------------------------

const svg = d3.select('#graph');

// new: transparent full-size rect to capture zoom/pan gestures
const zoomPane = svg.append('rect')
  .attr('class', 'zoom-pane')
  .attr('fill', 'transparent')
  .attr('pointer-events', 'all');

const graphGroup = svg.append('g');
const linkLayer = graphGroup.append('g');
const labelLayer = graphGroup.append('g');
const nodeLayer = graphGroup.append('g');

// new: shared zoom state and behavior
let currentTransform = d3.zoomIdentity;
const zoomBehavior = d3.zoom()
  .scaleExtent([0.25, 4]) // adjust min/max zoom
  .on('zoom', (event) => {
    currentTransform = event.transform;
    graphGroup.attr('transform', currentTransform);
  });


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

/** Canonical pair key "a|b" with a < b */
function pairKey(a, b) {
  return a < b ? a + '|' + b : b + '|' + a;
}

/** Shallow groupBy into Map(key -> array) */
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

/**
 * Weighted consensus by simple leave-one-out sign agreement.
 * Returns { scores: Map(pairKey -> score), weights: Map(expertId -> weight) }.
 */
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

  // Compute agreement ratio per expert
  for (const e of expertIds) {
    let agree = 0;
    let total = 0;
    for (const k of pairKeys) {
      // leave-one-out mean sign
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

  // Normalize to sum to 1, avoid division by zero
  const wsum = Array.from(weights.values()).reduce((a, b) => a + b, 0) || 1;
  for (const e of expertIds) {
    weights.set(e, weights.get(e) / wsum);
  }

  // Weighted mean per key
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

/** Distance mapping: positive scores pull together, negative push apart. */
function distanceFromScore(score, base = FORCE_BASE_DISTANCE, alpha = FORCE_ALPHA) {
  return base * Math.exp(-alpha * score);
}

/** Resize the SVG to the container size. */
function sizeSvgToContainer() {
  const el = document.getElementById('graph');
  const rect = el.getBoundingClientRect();

  svg.attr('width', rect.width).attr('height', rect.height);

  // new: keep the zoom-capture rect in sync with the SVG size
  zoomPane
    .attr('x', 0)
    .attr('y', 0)
    .attr('width', rect.width)
    .attr('height', rect.height);

  // new: update zoom extents for the new size and reapply the current transform
  zoomBehavior.extent([[0, 0], [rect.width, rect.height]]);
  svg.call(zoomBehavior).call(zoomBehavior.transform, currentTransform);
}


// ------------------------------ Data loading and indexing ------------------------------

/**
 * @typedef {Object} Dataset
 * @property {string} collectionId
 * @property {{id: string, type?: string}[]} bones
 * @property {{id: string, name: string}[]} experts
 * @property {{expertId: string, boneA: string, boneB: string, score: number}[]} ratings
 */

let dataset = null; // will hold the loaded JSON
let bones = [];     // string[]
let expertsById = new Map(); // id -> expert
let expertEdgeMaps = new Map(); // expertId -> Map(pairKey -> score)
let simulation = null;

/** Canonicalize rating entries and build indices. */
function buildIndices(raw) {
  // Copy bones and experts maps
  const bonesSet = new Set(raw.bones.map(b => b.id));
  for (const r of raw.ratings) {
    bonesSet.add(r.boneA);
    bonesSet.add(r.boneB);
  }
  bones = Array.from(bonesSet).sort();

  expertsById = new Map(raw.experts.map(e => [e.id, e]));

  // Canonicalize ratings (ensure boneA < boneB and add .key)
  const canonicalRatings = raw.ratings.map(r => {
    const a = r.boneA;
    const b = r.boneB;
    const [x, y] = a < b ? [a, b] : [b, a];
    return { expertId: r.expertId, boneA: x, boneB: y, score: r.score, key: pairKey(x, y) };
  });

  // Build expert -> Map(pairKey -> mean score if duplicates exist)
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

  // Populate expert selector
  populateExpertSelector(raw.experts);
}

/** Add experts to the dropdown after the default consensus. */
function populateExpertSelector(experts) {
  for (const e of experts) {
    const opt = document.createElement('option');
    opt.value = e.id;
    opt.textContent = e.name;
    expertSelect.appendChild(opt);
  }
}

/** All pair keys present in any expert map. */
function allPairKeys() {
  const keys = new Set();
  for (const m of expertEdgeMaps.values()) {
    for (const k of m.keys()) keys.add(k);
  }
  return Array.from(keys);
}

/** Build consensus map with the selected method. */
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

// ------------------------------ Rendering ------------------------------

/** Redraw the whole view based on UI state. */
function redraw() {
  sizeSvgToContainer();

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

  const nodes = bones.map(id => ({ id }));

  // Stop any previous simulation
  if (simulation) simulation.stop();

  // Scales
  const widthScale = d3.scaleLinear().domain([0, 2]).range([1, 4]);

  // Links
  const linksSel = linkLayer.selectAll('line').data(edges, d => d.source + '|' + d.target);
  linksSel.exit().remove();
  const linksEnter = linksSel.enter().append('line');
  const links = linksEnter.merge(linksSel)
    .attr('stroke', d => d.score >= 0 ? '#1e9e45' : '#c62828')
    .attr('stroke-width', d => widthScale(Math.min(2, Math.abs(d.score))))
    .attr('stroke-opacity', 0.9)
    .attr('stroke-dasharray', d => d.score < 0 ? '4,3' : '0');

  // Edge labels
  const labelSel = labelLayer.selectAll('text').data(edges, d => d.source + '|' + d.target);
  labelSel.exit().remove();
  const labelEnter = labelSel.enter().append('text').attr('class', 'edge-label');
  const labels = labelEnter.merge(labelSel).text(d => d.score.toFixed(2));

  // Nodes
  const nodeSel = nodeLayer.selectAll('g.node').data(nodes, d => d.id);
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
    links.attr('stroke-opacity', l => (l.source.id === id || l.target.id === id) ? 1.0 : 0.15);
    nodesMerged.selectAll('circle').attr('opacity', nd => nd.id === id ? 1.0 : 0.6);
  }).on('mouseout', function () {
    links.attr('stroke-opacity', 0.9);
    nodesMerged.selectAll('circle').attr('opacity', 1.0);
  });

  // Force simulation
  simulation = d3.forceSimulation(nodes)
    .force('charge', d3.forceManyBody().strength(CHARGE_STRENGTH))
    .force('center', d3.forceCenter(svg.attr('width') / 2, svg.attr('height') / 2))
    .force('link', d3.forceLink(edges)
      .id(d => d.id)
      .distance(d => distanceFromScore(d.score))
      .strength(0.7))
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

  // Update table and clusters
  updatePairsTable(edges, fromLabel);
  updateClusters(view, agg);
}

/** Update the top pairs table. */
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

/** Very simple must-link clustering by connected components. */
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

  // Build adjacency for edges with score >= must-link threshold
  const adj = new Map(bones.map(b => [b, new Set()]));
  for (const [k, s] of edgeMap.entries()) {
    if (s >= tPos) {
      const [a, b] = k.split('|');
      adj.get(a).add(b);
      adj.get(b).add(a);
    }
  }

  // Connected components DFS
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

  // Render clusters
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
    sizeSvgToContainer();
    redraw();
  } catch (err) {
    console.error('Failed to load data:', err);
  }
}

window.addEventListener('resize', () => {
  sizeSvgToContainer();
  // When resizing, recenter the force layout target
  if (simulation) {
    simulation.force('center', d3.forceCenter(svg.attr('width') / 2, svg.attr('height') / 2));
    simulation.alpha(0.3).restart();
  }
});

expertSelect.addEventListener('change', redraw);
aggSelect.addEventListener('change', redraw);
thresholdSlider.addEventListener('input', redraw);
showNegativeCheckbox.addEventListener('change', redraw);
mustLinkSlider.addEventListener('input', () => updateClusters(expertSelect.value, aggSelect.value));

// Kick off
init();
