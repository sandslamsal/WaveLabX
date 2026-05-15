/* =========================================================================
 * Wave Reflection Analysis — three-probe method
 *   Goda & Suzuki (1976); Kobayashi, Cox & Wurjanto (1990)
 * Pure client-side. No data leaves the browser.
 * ========================================================================= */

"use strict";

const G = 9.81;            // gravitational acceleration [m/s^2]
const COND_LIMIT = 30.0;   // 3-probe -> 2-probe fallback threshold
const DL_MIN = 0.05;       // valid spacing range  0.05 <= dl/L <= 0.45
const DL_MAX = 0.45;

/* ---------------------------------------------------------------------------
 * Complex arithmetic — complex numbers are {re, im}
 * ------------------------------------------------------------------------- */
const cAdd = (a, b) => ({ re: a.re + b.re, im: a.im + b.im });
const cSub = (a, b) => ({ re: a.re - b.re, im: a.im - b.im });
const cMul = (a, b) => ({
  re: a.re * b.re - a.im * b.im,
  im: a.re * b.im + a.im * b.re,
});
const cConj = (a) => ({ re: a.re, im: -a.im });
const cAbs = (a) => Math.hypot(a.re, a.im);
const cExp = (theta) => ({ re: Math.cos(theta), im: Math.sin(theta) }); // e^{i*theta}

/* Solve a 2x2 complex system  [[a,b],[c,d]] x = [p,q]. */
function solve2x2(a, b, c, d, p, q) {
  const det = cSub(cMul(a, d), cMul(b, c));
  const dd = det.re * det.re + det.im * det.im;
  if (dd < 1e-30) return [{ re: 0, im: 0 }, { re: 0, im: 0 }];
  const invDet = { re: det.re / dd, im: -det.im / dd };
  const x0 = cMul(cSub(cMul(d, p), cMul(b, q)), invDet);
  const x1 = cMul(cSub(cMul(a, q), cMul(c, p)), invDet);
  return [x0, x1];
}

/* Condition number of a 2x2 Hermitian PD matrix [[a,b],[b*,d]], a,d real. */
function condHermitian2(a, b, d) {
  const tr = a + d;
  const disc = Math.sqrt(((a - d) / 2) ** 2 + (b.re ** 2 + b.im ** 2));
  const lMax = tr / 2 + disc;
  const lMin = tr / 2 - disc;
  if (lMin <= 1e-30) return Infinity;
  return Math.sqrt(lMax / lMin);
}

/* ---------------------------------------------------------------------------
 * Linear dispersion relation  omega^2 = g k tanh(k d)  -> wave number k
 * ------------------------------------------------------------------------- */
function dispersion(freq, depth) {
  const omega = 2 * Math.PI * freq;
  const target = omega * omega;
  let k = target / G; // deep-water guess
  if (!isFinite(k) || k <= 0) return target / G;
  for (let it = 0; it < 100; it++) {
    const th = Math.tanh(k * depth);
    const F = G * k * th - target;
    const dF = G * th + G * k * depth * (1 - th * th);
    if (Math.abs(dF) < 1e-30) break;
    const kNew = k - F / dF;
    if (kNew <= 0) { k = k / 2; continue; }
    if (Math.abs(kNew - k) < 1e-12) { k = kNew; break; }
    k = kNew;
  }
  return k;
}
const wavelength = (freq, depth) => (2 * Math.PI) / dispersion(freq, depth);

/* ---------------------------------------------------------------------------
 * Hann window of length N (cached).
 * ------------------------------------------------------------------------- */
const _hannCache = {};
function hann(N) {
  if (_hannCache[N]) return _hannCache[N];
  const w = new Float64Array(N);
  let s = 0;
  for (let n = 0; n < N; n++) {
    w[n] = 0.5 * (1 - Math.cos((2 * Math.PI * n) / (N - 1)));
    s += w[n];
  }
  w.sum = s;
  _hannCache[N] = w;
  return w;
}

/* Per-column mean (DC offset). */
function colMean(col) {
  let m = 0;
  for (let n = 0; n < col.length; n++) m += col[n];
  return m / col.length;
}

/* ---------------------------------------------------------------------------
 * Detect the dominant wave frequency from a signal, within [fmin,fmax].
 * Uses the gauge with the largest variance; band-limited DFT magnitude scan
 * with a twiddle-factor recurrence for speed.
 * ------------------------------------------------------------------------- */
function detectFrequency(columns, fs, fmin, fmax) {
  // pick the strongest gauge
  let sig = columns[0], bestVar = -1;
  for (const c of columns) {
    const m = colMean(c);
    let v = 0;
    for (let n = 0; n < c.length; n++) v += (c[n] - m) * (c[n] - m);
    if (v > bestVar) { bestVar = v; sig = c; }
  }
  const N = sig.length;
  const mean = colMean(sig);
  const w = hann(N);
  const df = fs / N;
  const iLo = Math.max(1, Math.floor(fmin / df));
  const iHi = Math.min(Math.floor(N / 2), Math.ceil(fmax / df));
  if (iHi <= iLo) return null;

  let bestMag = -1, bestIdx = iLo;
  for (let idx = iLo; idx <= iHi; idx++) {
    const ang = (-2 * Math.PI * idx) / N;
    const cs = Math.cos(ang), sn = Math.sin(ang);
    let tr = 1, ti = 0;        // twiddle^n
    let re = 0, im = 0;
    for (let n = 0; n < N; n++) {
      const v = (sig[n] - mean) * w[n];
      re += v * tr;
      im += v * ti;
      const ntr = tr * cs - ti * sn;
      ti = tr * sn + ti * cs;
      tr = ntr;
    }
    const mag = re * re + im * im;
    if (mag > bestMag) { bestMag = mag; bestIdx = idx; }
  }
  return bestIdx * df;
}

/* ---------------------------------------------------------------------------
 * Single-bin DFT amplitudes of three Hann-windowed gauge records.
 * ------------------------------------------------------------------------- */
function gaugeAmplitudes(columns, fs, fTarget) {
  const N = columns[0].length;
  const w = hann(N);
  const df = fs / N;
  let idx = Math.round(fTarget / df);
  idx = Math.max(1, Math.min(idx, Math.floor(N / 2)));
  const binFreq = idx * df;

  const A = [];
  for (let j = 0; j < 3; j++) {
    const col = columns[j];
    const mean = colMean(col);
    let re = 0, im = 0;
    for (let n = 0; n < N; n++) {
      const ang = (-2 * Math.PI * idx * n) / N;
      const v = (col[n] - mean) * w[n];
      re += v * Math.cos(ang);
      im += v * Math.sin(ang);
    }
    A.push({ re: (2 * re) / w.sum, im: (2 * im) / w.sum });
  }
  return { A, binFreq };
}

/* ---------------------------------------------------------------------------
 * Goda-Suzuki three-probe separation for one gauge array.
 * ------------------------------------------------------------------------- */
function threeProbe(columns, fs, fTarget, depth, pos) {
  const { A, binFreq } = gaugeAmplitudes(columns, fs, fTarget);
  const k = dispersion(binFreq, depth);

  const inc = pos.map((x) => cExp(-k * x));   // incident   e^{-i k x}
  const ref = pos.map((x) => cExp(k * x));    // reflected  e^{+i k x}

  // (M^H M)[0][1] = sum_j e^{+2 i k x_j};  diagonal entries = 3
  let off = { re: 0, im: 0 };
  for (let j = 0; j < 3; j++) off = cAdd(off, cExp(2 * k * pos[j]));
  const cond3 = condHermitian2(3, off, 3);

  let Hi, Hr, fallback = false;

  if (cond3 < COND_LIMIT) {
    let p = { re: 0, im: 0 }, q = { re: 0, im: 0 };
    for (let j = 0; j < 3; j++) {
      p = cAdd(p, cMul(cConj(inc[j]), A[j]));
      q = cAdd(q, cMul(cConj(ref[j]), A[j]));
    }
    const sol = solve2x2(
      { re: 3, im: 0 }, off,
      cConj(off), { re: 3, im: 0 },
      p, q
    );
    Hi = 2 * cAbs(sol[0]);
    Hr = 2 * cAbs(sol[1]);
  } else {
    fallback = true;
    const pairs = [[0, 1], [1, 2], [0, 2]];
    let best = pairs[0], bestCond = Infinity;
    for (const [i, j] of pairs) {
      const o = cAdd(cExp(2 * k * pos[i]), cExp(2 * k * pos[j]));
      const c = condHermitian2(2, o, 2);
      if (c < bestCond) { bestCond = c; best = [i, j]; }
    }
    const [i, j] = best;
    const sol = solve2x2(inc[i], ref[i], inc[j], ref[j], A[i], A[j]);
    Hi = 2 * cAbs(sol[0]);
    Hr = 2 * cAbs(sol[1]);
  }
  return { Hi, Hr, fallback, k, binFreq };
}

/* ---------------------------------------------------------------------------
 * CSV parsing — 6 numeric columns, optional header row.
 * ------------------------------------------------------------------------- */
function parseCSV(text) {
  const lines = text.split(/\r?\n/);
  const cols = [[], [], [], [], [], []];
  let started = false;
  for (let li = 0; li < lines.length; li++) {
    const line = lines[li].trim();
    if (!line) continue;
    const parts = line.split(/[,;\t]/);
    if (parts.length < 6) continue;
    const vals = parts.slice(0, 6).map((p) => parseFloat(p));
    if (vals.some((v) => Number.isNaN(v))) { continue; } // header / junk
    started = true;
    for (let c = 0; c < 6; c++) cols[c].push(vals[c]);
  }
  return started ? cols : null;
}

/* Water depth from file name (frequency is detected from the signal). */
function parseDepth(name) {
  const d = name.match(/Depth\s*=\s*([0-9]*\.?[0-9]+)/i);
  return d ? parseFloat(d[1]) : null;
}

/* =========================================================================
 * APPLICATION
 * ========================================================================= */
const state = { records: [] };
let nextId = 1;
const $ = (id) => document.getElementById(id);

function getSettings() {
  const fs = parseFloat($("fs").value) || 100;
  const fmin = parseFloat($("fmin").value) || 0.1;
  const fmax = parseFloat($("fmax").value) || 2.0;
  const depth = parseFloat($("depth").value) || 0.25;
  const pos1 = [...document.querySelectorAll(".pos1")]
    .sort((a, b) => a.dataset.i - b.dataset.i)
    .map((el) => parseFloat(el.value) || 0);
  const pos2 = [...document.querySelectorAll(".pos2")]
    .sort((a, b) => a.dataset.i - b.dataset.i)
    .map((el) => parseFloat(el.value) || 0);
  return { fs, fmin, fmax, depth, pos1, pos2 };
}

/* Analyse one record. redetect = re-run frequency detection. */
function analyzeRecord(rec, redetect) {
  rec.error = null;
  const s = getSettings();

  if (!rec.cols || rec.cols[0].length < 16) {
    rec.error = "Not a valid 6-column time series";
    rec.result = null;
    return;
  }
  if (redetect && !rec.freqManual) {
    const f = detectFrequency(rec.cols.slice(0, 3), s.fs, s.fmin, s.fmax);
    if (f) rec.freq = f;
  }
  if (!(rec.depth > 0) || !(rec.freq > 0)) {
    rec.error = "Set depth / frequency";
    rec.result = null;
    return;
  }
  try {
    const a1 = threeProbe(rec.cols.slice(0, 3), s.fs, rec.freq, rec.depth, s.pos1);
    const a2 = threeProbe(rec.cols.slice(3, 6), s.fs, rec.freq, rec.depth, s.pos2);
    const L = wavelength(rec.freq, rec.depth);
    const dl1 = (Math.max(...s.pos1) - Math.min(...s.pos1)) / L;
    const dl2 = (Math.max(...s.pos2) - Math.min(...s.pos2)) / L;
    rec.result = {
      Hi1: a1.Hi, Hr1: a1.Hr, Kr1: a1.Hi > 0 ? a1.Hr / a1.Hi : NaN,
      Hi2: a2.Hi, Hr2: a2.Hr, Kr2: a2.Hi > 0 ? a2.Hr / a2.Hi : NaN,
      Kt: a1.Hi > 0 ? a2.Hi / a1.Hi : NaN,
      period: 1 / rec.freq,
      fallback: a1.fallback || a2.fallback,
      L, dl1, dl2,
      ratioWarn: dl1 < DL_MIN || dl1 > DL_MAX || dl2 < DL_MIN || dl2 > DL_MAX,
    };
  } catch (e) {
    rec.error = "Computation failed: " + e.message;
    rec.result = null;
  }
}

function analyzeAll(redetect) {
  state.records.forEach((r) => analyzeRecord(r, redetect));
  renderTable();
}

/* ---------------------------------------------------------------------------
 * Rendering
 * ------------------------------------------------------------------------- */
const fmt = (x, p = 4) =>
  x == null || Number.isNaN(x) ? "&mdash;" : Number(x).toFixed(p);

function renderTable() {
  const body = $("resultsBody");
  body.innerHTML = "";
  let anyFallback = false, anyRatio = false;

  state.records.forEach((rec) => {
    const tr = document.createElement("tr");
    const r = rec.result;
    if (r && r.fallback) anyFallback = true;
    if (r && r.ratioWarn) anyRatio = true;

    const depthCell = `<td class="editable">
        <input type="number" step="0.01" min="0" value="${rec.depth ?? ""}"
               data-id="${rec.id}" data-field="depth" /></td>`;
    const freqCell = `<td class="editable">
        <input type="number" step="0.001" min="0"
               value="${rec.freq != null ? rec.freq.toFixed(3) : ""}"
               data-id="${rec.id}" data-field="freq" /></td>`;

    if (rec.error) {
      tr.innerHTML = `
        <td title="${rec.name}">${rec.name}</td>
        ${depthCell}${freqCell}
        <td colspan="8" class="badge-err">${rec.error}</td>
        <td><button class="row-del" data-del="${rec.id}">&times;</button></td>`;
    } else {
      const cls = r.fallback || r.ratioWarn ? ' class="fallback"' : "";
      const flag = r.fallback ? " &#9888;" : r.ratioWarn ? " &#9888;" : "";
      tr.innerHTML = `
        <td title="${rec.name}">${rec.name}${flag}</td>
        ${depthCell}${freqCell}
        <td>${fmt(r.period, 3)}</td>
        <td${cls}>${fmt(r.Hi1)}</td>
        <td${cls}>${fmt(r.Hr1)}</td>
        <td${cls}>${fmt(r.Kr1, 3)}</td>
        <td${cls}>${fmt(r.Hi2)}</td>
        <td${cls}>${fmt(r.Hr2)}</td>
        <td${cls}>${fmt(r.Kr2, 3)}</td>
        <td>${fmt(r.Kt, 3)}</td>
        <td><button class="row-del" data-del="${rec.id}">&times;</button></td>`;
    }
    body.appendChild(tr);
  });

  $("rowCount").innerHTML =
    state.records.length ? `&middot; ${state.records.length} file(s)` : "";
  $("resultsCard").hidden = state.records.length === 0;

  const warn = $("warnNote");
  const msgs = [];
  if (anyFallback)
    msgs.push("&#9888; A two-gauge fallback was used where the three-gauge " +
      "system was ill-conditioned (spacing near a multiple of L/2).");
  if (anyRatio)
    msgs.push("&#9888; Some rows have gauge spacing outside the valid range " +
      "0.05 &le; &Delta;l/L &le; 0.45 — interpret those with caution.");
  if (msgs.length) { warn.hidden = false; warn.innerHTML = msgs.join("<br>"); }
  else warn.hidden = true;

  body.querySelectorAll("input[data-field]").forEach((inp) => {
    inp.addEventListener("change", () => {
      const rec = state.records.find((x) => x.id === +inp.dataset.id);
      if (!rec) return;
      const val = parseFloat(inp.value);
      rec[inp.dataset.field] = Number.isNaN(val) ? null : val;
      if (inp.dataset.field === "freq") rec.freqManual = true;
      analyzeRecord(rec, false);
      renderTable();
    });
  });
  body.querySelectorAll("[data-del]").forEach((btn) => {
    btn.addEventListener("click", () => {
      state.records = state.records.filter((x) => x.id !== +btn.dataset.del);
      renderTable();
    });
  });
}

function updateSpacingReadout() {
  const s = getSettings();
  const sp = (p) =>
    `Spacing: 1&ndash;2 = ${(p[1] - p[0]).toFixed(2)} m, ` +
    `2&ndash;3 = ${(p[2] - p[1]).toFixed(2)} m, ` +
    `1&ndash;3 = ${(p[2] - p[0]).toFixed(2)} m`;
  $("spacing1").innerHTML = sp(s.pos1);
  $("spacing2").innerHTML = sp(s.pos2);
}

/* ---------------------------------------------------------------------------
 * File ingestion
 * ------------------------------------------------------------------------- */
function readFiles(fileList) {
  const files = [...fileList].filter((f) => /\.csv$/i.test(f.name));
  if (!files.length) return;
  const s = getSettings();
  let pending = files.length;

  files.forEach((file) => {
    const reader = new FileReader();
    reader.onload = (e) => {
      const cols = parseCSV(e.target.result);
      const rec = {
        id: nextId++,
        name: file.name,
        depth: parseDepth(file.name) ?? s.depth,
        freq: null,
        freqManual: false,
        cols,
        result: null,
        error: null,
      };
      analyzeRecord(rec, true); // detect frequency on load
      state.records.push(rec);
      if (--pending === 0) {
        state.records.sort((a, b) => a.name.localeCompare(b.name));
        renderTable();
      }
    };
    reader.onerror = () => {
      state.records.push({
        id: nextId++, name: file.name, depth: null, freq: null,
        freqManual: false, cols: null, result: null, error: "Could not read file",
      });
      if (--pending === 0) renderTable();
    };
    reader.readAsText(file);
  });
}

/* ---------------------------------------------------------------------------
 * Export
 * ------------------------------------------------------------------------- */
function exportCSV() {
  const head = [
    "File", "Depth_d_m", "Freq_f_Hz", "Period_T_s",
    "Hi1_m", "Hr1_m", "Kr1", "Hi2_m", "Hr2_m", "Kr2",
    "Kt_Hi2overHi1", "Wavelength_L_m", "Fallback2probe", "SpacingRatioWarn",
  ];
  const rows = [head.join(",")];
  state.records.forEach((rec) => {
    const r = rec.result;
    if (!r) {
      rows.push([`"${rec.name}"`, rec.depth ?? "", rec.freq ?? "",
        "", "", "", "", "", "", "", "", "", "", rec.error || "error"].join(","));
      return;
    }
    const g = (x, p) => (x == null || Number.isNaN(x) ? "" : x.toFixed(p));
    rows.push([
      `"${rec.name}"`, g(rec.depth, 4), g(rec.freq, 5), g(r.period, 4),
      g(r.Hi1, 6), g(r.Hr1, 6), g(r.Kr1, 4),
      g(r.Hi2, 6), g(r.Hr2, 6), g(r.Kr2, 4),
      g(r.Kt, 4), g(r.L, 4),
      r.fallback ? "yes" : "no", r.ratioWarn ? "yes" : "no",
    ].join(","));
  });
  const blob = new Blob([rows.join("\n")], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "reflection_analysis_results.csv";
  a.click();
  URL.revokeObjectURL(url);
}

/* ---------------------------------------------------------------------------
 * Wiring
 * ------------------------------------------------------------------------- */
function init() {
  const dz = $("dropzone");
  const fileInput = $("fileInput");

  dz.addEventListener("click", () => fileInput.click());
  fileInput.addEventListener("change", (e) => {
    readFiles(e.target.files);
    fileInput.value = "";
  });
  ["dragenter", "dragover"].forEach((ev) =>
    dz.addEventListener(ev, (e) => { e.preventDefault(); dz.classList.add("drag"); })
  );
  ["dragleave", "drop"].forEach((ev) =>
    dz.addEventListener(ev, (e) => { e.preventDefault(); dz.classList.remove("drag"); })
  );
  dz.addEventListener("drop", (e) => {
    if (e.dataTransfer?.files) readFiles(e.dataTransfer.files);
  });

  $("recomputeBtn").addEventListener("click", () => analyzeAll(true));
  $("exportBtn").addEventListener("click", exportCSV);
  $("clearBtn").addEventListener("click", () => {
    state.records = [];
    renderTable();
  });
  $("applyDepth").addEventListener("click", () => {
    const d = parseFloat($("depth").value);
    if (!(d > 0)) return;
    state.records.forEach((r) => { r.depth = d; });
    analyzeAll(false);
  });

  // gauge layout / sampling changes -> recompute (no re-detect needed)
  document.querySelectorAll("#fs, .pos1, .pos2").forEach((el) =>
    el.addEventListener("change", () => {
      updateSpacingReadout();
      if (state.records.length) analyzeAll(false);
    })
  );
  // detection-band changes -> recompute WITH re-detection
  document.querySelectorAll("#fmin, #fmax").forEach((el) =>
    el.addEventListener("change", () => {
      if (state.records.length) analyzeAll(true);
    })
  );

  updateSpacingReadout();
}

if (typeof document !== "undefined") {
  document.addEventListener("DOMContentLoaded", init);
}

/* expose core functions for Node-based testing */
if (typeof module !== "undefined" && module.exports) {
  module.exports = { threeProbe, dispersion, wavelength, detectFrequency, parseCSV, parseDepth };
}
