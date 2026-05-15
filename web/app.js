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
  // Each array is defined by two spacings (X12, X23); gauge positions
  // are 0, X12, X12+X23.
  const spacings = (cls) => {
    const s = [...document.querySelectorAll(cls)]
      .sort((a, b) => a.dataset.i - b.dataset.i)
      .map((el) => parseFloat(el.value) || 0);
    return [0, s[0], s[0] + s[1]];
  };
  const pos1 = spacings(".sp1");
  const pos2 = spacings(".sp2");
  const skipWaves = Math.max(0, parseInt($("skipWaves").value, 10) || 0);
  const numWaves = Math.max(0, parseInt($("numWaves").value, 10) || 0);
  return { fs, fmin, fmax, depth, skipWaves, numWaves, pos1, pos2 };
}

/* Current wave-type mode: "regular" or "irregular". */
function getMode() {
  const r = document.querySelector('input[name="wavemode"]:checked');
  return r ? r.value : "regular";
}

/* Analyse one record. redetect = re-run frequency detection (regular mode). */
function analyzeRecord(rec, redetect) {
  rec.error = null;
  const s = getSettings();
  const mode = getMode();

  if (!rec.cols || rec.cols[0].length < 16) {
    rec.error = "Not a valid 6-column time series";
    rec.result = null;
    return;
  }

  // ===== IRREGULAR (spectral) mode ====================================
  if (mode === "irregular") {
    if (!(rec.depth > 0)) {
      rec.error = "Set water depth";
      rec.result = null;
      return;
    }
    const SP = typeof WaveLabXSpectral !== "undefined" ? WaveLabXSpectral
      : typeof window !== "undefined" ? window.WaveLabXSpectral : null;
    if (!SP) {
      rec.error = "Spectral module not loaded";
      rec.result = null;
      return;
    }
    try {
      const a1 = SP.threeProbeArray(rec.cols.slice(0, 3), s.fs, rec.depth, s.pos1);
      const a2 = SP.threeProbeArray(rec.cols.slice(3, 6), s.fs, rec.depth, s.pos2);
      rec.freq = a1.fp;                       // spectral-peak frequency
      const ret = Math.min(a1.retained, a2.retained);
      rec.result = {
        Hi1: a1.Hi, Hr1: a1.Hr, Kr1: a1.Kr,
        Hi2: a2.Hi, Hr2: a2.Hr, Kr2: a2.Kr,
        Kt: a1.Hi > 0 ? a2.Hi / a1.Hi : NaN,
        period: a1.Tp,
        fallback: false, ratioWarn: false,
        retained: ret, retainedWarn: !(ret >= 0.8),
        spectral: true,
      };
    } catch (e) {
      rec.error = "Computation failed: " + e.message;
      rec.result = null;
    }
    return;
  }

  // ===== REGULAR (single-frequency) mode ==============================
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
    // Optional analysis window: skip the first M waves, then use N waves.
    // Frequency is always detected on the full record; the window only
    // restricts the three-probe computation.
    let cols = rec.cols;
    rec.windowInfo = null;
    if ((s.skipWaves > 0 || s.numWaves > 0) && rec.freq > 0) {
      const N = rec.cols[0].length;
      const spw = s.fs / rec.freq;                 // samples per wave
      let start = Math.round(s.skipWaves * spw);
      let len = s.numWaves > 0 ? Math.round(s.numWaves * spw) : N - start;
      start = Math.min(Math.max(start, 0), N);
      len = Math.min(Math.max(len, 0), N - start);
      if (len < 16) {
        rec.error = "Analysis window too short — reduce skip or increase waves";
        rec.result = null;
        return;
      }
      cols = rec.cols.map((c) => c.slice(start, start + len));
      rec.windowInfo = { start, len, waves: len / spw };
    }
    const a1 = threeProbe(cols.slice(0, 3), s.fs, rec.freq, rec.depth, s.pos1);
    const a2 = threeProbe(cols.slice(3, 6), s.fs, rec.freq, rec.depth, s.pos2);
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
      spectral: false,
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
  const irregular = getMode() === "irregular";
  let anyFallback = false, anyRatio = false, anyRetained = false;

  state.records.forEach((rec) => {
    const tr = document.createElement("tr");
    const r = rec.result;
    if (r && r.fallback) anyFallback = true;
    if (r && r.ratioWarn) anyRatio = true;
    if (r && r.retainedWarn) anyRetained = true;

    const depthCell = `<td class="editable">
        <input type="number" step="0.01" min="0" value="${rec.depth ?? ""}"
               data-id="${rec.id}" data-field="depth" /></td>`;
    const freqVal = rec.freq != null ? rec.freq.toFixed(3) : "";
    // frequency is user-editable in regular mode; derived (peak) in irregular
    const freqCell = irregular
      ? `<td>${freqVal || "&mdash;"}</td>`
      : `<td class="editable">
        <input type="number" step="0.001" min="0" value="${freqVal}"
               data-id="${rec.id}" data-field="freq" /></td>`;

    if (rec.error) {
      tr.innerHTML = `
        <td title="${rec.name}">${rec.name}</td>
        ${depthCell}${freqCell}
        <td colspan="8" class="badge-err">${rec.error}</td>
        <td><button class="row-del" data-del="${rec.id}">&times;</button></td>`;
    } else {
      const warnRow = r.fallback || r.ratioWarn || r.retainedWarn;
      const cls = warnRow ? ' class="fallback"' : "";
      const flag = warnRow ? " &#9888;" : "";
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
  if (anyRetained)
    msgs.push("&#9888; Some rows retained less than 80% of spectral energy " +
      "within the valid frequency band — interpret those with caution.");
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

  refreshVizFiles();
}

/* ===========================================================================
 * VISUALIZATION — time-series and energy-spectrum plots (canvas)
 * ========================================================================= */
const VIZ_COLORS = ["#1f5fa6", "#c0392b", "#1f7a4d", "#b9591a", "#6a4ca8", "#0e8a8a"];
const VIZ_SPEC = { inc: "#1f5fa6", ref: "#c0392b", tra: "#1f7a4d" };

let vizView = null;   // visible x-window {x0,x1}; null = full range
let vizYView = null;  // visible y-window {y0,y1}; null = autoscale
let vizMode = "none"; // active drag tool: "none" | "box" | "pan"
let vizDrag = null;   // in-progress drag state
let vizGeom = null;   // geometry of the last draw, for pixel<->data inversion
let vizSpecCache = null; // { key, a1, a2 } — cached spectral analysis

const vizType = () => ($("vizType") ? $("vizType").value : "series");

/* Repopulate the file selector from the current records (selection kept). */
function refreshVizFiles() {
  const sel = $("vizFile");
  if (!sel) return;
  const prev = sel.value;
  sel.innerHTML = "";
  state.records.forEach((rec) => {
    if (!rec.cols) return;
    const o = document.createElement("option");
    o.value = String(rec.id);
    o.textContent = rec.name;
    sel.appendChild(o);
  });
  if ([...sel.options].some((o) => o.value === prev)) sel.value = prev;
  if ($("vizBox") && $("vizBox").open) drawViz();
}

/* Adaptive numeric format for an axis, given the range it spans. */
function fmtAxis(range) {
  const a = Math.abs(range);
  if (a >= 200) return (v) => v.toFixed(0);
  if (a >= 20) return (v) => v.toFixed(1);
  if (a >= 2) return (v) => v.toFixed(2);
  if (a >= 0.2) return (v) => v.toFixed(3);
  if (a >= 0.002) return (v) => v.toFixed(5);
  return (v) => v.toExponential(1);
}

/* Build a time-series plot object for the time-series mode. */
function buildSeriesPlot(rec) {
  const probes = [...document.querySelectorAll(".viz-probe")]
    .filter((c) => c.checked)
    .map((c) => +c.value);
  if (!probes.length) return { empty: "Select at least one probe to plot." };

  const fs = parseFloat($("fs").value) || 100;
  const N = rec.cols[0].length;
  const dt = 1 / fs;
  const x = new Float64Array(N);
  for (let i = 0; i < N; i++) x[i] = i * dt;

  const series = probes.map((p) => ({
    label: "Probe " + (p + 1),
    color: VIZ_COLORS[p % VIZ_COLORS.length],
    x,
    y: rec.cols[p],
  }));
  return {
    xLabel: "Time (s)",
    yLabel: "Surface elevation",
    yUnit: "m",
    xMin: 0,
    xMax: (N - 1) * dt,
    series,
    info: `${rec.name} — ${N} samples, ${((N - 1) * dt).toFixed(1)} s at ${fs} Hz`,
  };
}

/* Build an energy-spectrum plot object: incident/reflected from probes
 * 1-3, transmitted from probes 4-6, via the WaveLabX spectral method. */
function buildSpectrumPlot(rec) {
  const SP = typeof WaveLabXSpectral !== "undefined" ? WaveLabXSpectral
    : typeof window !== "undefined" ? window.WaveLabXSpectral : null;
  if (!SP) return { empty: "Spectral module not loaded." };
  if (!(rec.depth > 0)) return { empty: "Set a water depth for this file first." };

  const fs = parseFloat($("fs").value) || 100;
  const s = getSettings();
  const key = [rec.id, fs, rec.depth, s.pos1.join(","), s.pos2.join(",")].join("|");
  if (!vizSpecCache || vizSpecCache.key !== key) {
    vizSpecCache = {
      key,
      a1: SP.threeProbeArray(rec.cols.slice(0, 3), fs, rec.depth, s.pos1),
      a2: SP.threeProbeArray(rec.cols.slice(3, 6), fs, rec.depth, s.pos2),
    };
  }
  const { a1, a2 } = vizSpecCache;

  const want = [...document.querySelectorAll(".viz-curve")]
    .filter((c) => c.checked)
    .map((c) => c.value);
  if (!want.length) return { empty: "Select at least one spectrum to plot." };

  const series = [];
  if (want.includes("inc"))
    series.push({ label: "Incident", color: VIZ_SPEC.inc, x: a1.spectra.f, y: a1.spectra.Si });
  if (want.includes("ref"))
    series.push({ label: "Reflected", color: VIZ_SPEC.ref, x: a1.spectra.f, y: a1.spectra.Sr });
  if (want.includes("tra"))
    series.push({ label: "Transmitted", color: VIZ_SPEC.tra, x: a2.spectra.f, y: a2.spectra.Si });

  let xMin = Infinity, xMax = -Infinity;
  for (const ser of series)
    for (const xv of ser.x) { if (xv < xMin) xMin = xv; if (xv > xMax) xMax = xv; }
  if (!isFinite(xMin)) return { empty: "No spectral data in the valid frequency band." };

  return {
    xLabel: "Frequency (Hz)",
    yLabel: "Spectral density S(f)",
    yUnit: "m²·s",
    xMin, xMax, series,
    info: `${rec.name} — incident/reflected from probes 1–3, ` +
      `transmitted from probes 4–6`,
  };
}

/* Draw the active plot (time series or energy spectrum). */
function drawViz() {
  const canvas = $("vizCanvas");
  if (!canvas) return;
  const hint = $("vizHint");
  const rec = state.records.find((r) => String(r.id) === $("vizFile").value);

  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || 1000;
  const cssH = 380;
  canvas.width = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, cssH);

  if (!rec || !rec.cols) {
    vizGeom = null;
    hint.textContent = "No file selected — load CSV files first.";
    return;
  }
  let plot;
  try {
    plot = vizType() === "spectrum" ? buildSpectrumPlot(rec) : buildSeriesPlot(rec);
  } catch (e) {
    vizGeom = null;
    hint.textContent = "Plot failed: " + e.message;
    return;
  }
  if (!plot || plot.empty || !plot.series || !plot.series.length) {
    vizGeom = null;
    hint.textContent = plot && plot.empty ? plot.empty : "Nothing to plot.";
    return;
  }
  renderPlot(ctx, cssW, cssH, plot);
}

/* Shared renderer: axes, gridlines, line+dot series, legend, zoom window. */
function renderPlot(ctx, cssW, cssH, plot) {
  const mL = 62, mR = 14, mT = 26, mB = 34;
  const pW = cssW - mL - mR, pH = cssH - mT - mB;

  // visible x-window
  let x0 = plot.xMin, x1 = plot.xMax;
  if (x1 <= x0) x1 = x0 + 1;
  const xFull = x1 - x0;
  if (vizView) {
    x0 = Math.max(plot.xMin, Math.min(vizView.x0, plot.xMax));
    x1 = Math.min(plot.xMax, Math.max(vizView.x1, x0 + xFull * 1e-4));
  }

  // visible y-window: explicit if box-zoomed, else autoscale within x-window
  let yMin, yMax;
  if (vizYView) {
    yMin = vizYView.y0; yMax = vizYView.y1;
  } else {
    yMin = Infinity; yMax = -Infinity;
    for (const ser of plot.series) {
      const n = ser.y.length;
      for (let i = 0; i < n; i++) {
        const xv = ser.x[i];
        if (xv < x0 || xv > x1) continue;
        const yv = ser.y[i];
        if (Number.isNaN(yv)) continue;
        if (yv < yMin) yMin = yv;
        if (yv > yMax) yMax = yv;
      }
    }
    if (!isFinite(yMin)) { yMin = 0; yMax = 1; }
    if (yMin === yMax) { yMin -= 1; yMax += 1; }
    const pad = (yMax - yMin) * 0.06;
    yMin -= pad; yMax += pad;
  }

  const xOf = (t) => mL + ((t - x0) / (x1 - x0)) * pW;
  const yOf = (v) => mT + pH - ((v - yMin) / (yMax - yMin)) * pH;
  const xFmt = fmtAxis(x1 - x0);
  // factor the y values into a "×10^k" label multiplier so the tick
  // numbers stay in a readable mantissa range (applies to both plots)
  let yScaleExp = 0;
  const yMaxAbs = Math.max(Math.abs(yMin), Math.abs(yMax));
  if (yMaxAbs > 0) {
    yScaleExp = -Math.floor(Math.log10(yMaxAbs));
  }
  const yScale = Math.pow(10, yScaleExp);
  const yFmt = fmtAxis((yMax - yMin) * yScale);

  // grid + ticks
  ctx.strokeStyle = "#e2e7ee";
  ctx.lineWidth = 1;
  ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
  ctx.fillStyle = "#6b7787";
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (let i = 0; i <= 5; i++) {
    const v = yMin + (i / 5) * (yMax - yMin);
    const y = yOf(v);
    ctx.beginPath(); ctx.moveTo(mL, y); ctx.lineTo(mL + pW, y); ctx.stroke();
    ctx.fillText(yFmt(v * yScale), mL - 6, y);
  }
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  for (let i = 0; i <= 6; i++) {
    const t = x0 + (i / 6) * (x1 - x0);
    const x = xOf(t);
    ctx.beginPath(); ctx.moveTo(x, mT); ctx.lineTo(x, mT + pH); ctx.stroke();
    ctx.fillText(xFmt(t), x, mT + pH + 6);
  }

  // axis titles
  ctx.fillStyle = "#1a2433";
  ctx.fillText(plot.xLabel, mL + pW / 2, cssH - 13);

  // rotated y-axis title, with a "×10^k" factor (k as a superscript)
  const baseFont = "11px -apple-system, BlinkMacSystemFont, sans-serif";
  const supFont = "8px -apple-system, BlinkMacSystemFont, sans-serif";
  const segs = [{ t: plot.yLabel, sup: false }];
  if (yScaleExp !== 0) {
    // axis reads value × 10^labelExp, so labelExp = -yScaleExp
    const labelExp = -yScaleExp;
    segs.push({ t: "  ×10", sup: false });
    segs.push({ t: (labelExp < 0 ? "−" : "") + Math.abs(labelExp), sup: true });
  }
  segs.push({ t: " (" + (plot.yUnit || "") + ")", sup: false });
  ctx.save();
  ctx.translate(13, mT + pH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "left";
  ctx.textBaseline = "middle";
  let segTotal = 0;
  for (const sg of segs) {
    ctx.font = sg.sup ? supFont : baseFont;
    segTotal += ctx.measureText(sg.t).width;
  }
  let segX = -segTotal / 2;
  for (const sg of segs) {
    ctx.font = sg.sup ? supFont : baseFont;
    ctx.fillText(sg.t, segX, sg.sup ? -4 : 0);
    segX += ctx.measureText(sg.t).width;
  }
  ctx.font = baseFont;
  ctx.restore();

  // series — line + dot at each plotted point, clipped to the plot area
  ctx.save();
  ctx.beginPath();
  ctx.rect(mL, mT, pW, pH);
  ctx.clip();
  for (const ser of plot.series) {
    const n = ser.y.length;
    let lo = 0, hi = n - 1;            // visible index range (x ascending)
    while (lo < n && ser.x[lo] < x0) lo++;
    while (hi >= 0 && ser.x[hi] > x1) hi--;
    lo = Math.max(0, lo - 1);
    hi = Math.min(n - 1, hi + 1);
    if (hi < lo) continue;
    const nVis = hi - lo + 1;
    const stride = Math.max(1, Math.floor(nVis / Math.max(800, pW)));
    const dotR = nVis / stride < pW / 6 ? 2.6 : 1.7;

    ctx.strokeStyle = ser.color;
    ctx.lineWidth = 0.9;
    ctx.beginPath();
    let first = true;
    for (let i = lo; i <= hi; i += stride) {
      const yv = ser.y[i];
      if (Number.isNaN(yv)) { first = true; continue; }
      const x = xOf(ser.x[i]), y = yOf(yv);
      if (first) { ctx.moveTo(x, y); first = false; }
      else ctx.lineTo(x, y);
    }
    ctx.stroke();

    ctx.fillStyle = ser.color;
    for (let i = lo; i <= hi; i += stride) {
      const yv = ser.y[i];
      if (Number.isNaN(yv)) continue;
      ctx.beginPath();
      ctx.arc(xOf(ser.x[i]), yOf(yv), dotR, 0, 2 * Math.PI);
      ctx.fill();
    }
  }
  ctx.restore();

  // legend along the top
  ctx.textAlign = "left";
  ctx.textBaseline = "middle";
  let lx = mL;
  for (const ser of plot.series) {
    ctx.fillStyle = ser.color;
    ctx.fillRect(lx, mT - 15, 14, 3);
    ctx.fillStyle = "#1a2433";
    ctx.fillText(ser.label, lx + 18, mT - 14);
    lx += ctx.measureText(ser.label).width + 32;
  }

  // geometry for the drag tools
  vizGeom = { mL, mT, pW, pH, vx0: x0, vx1: x1, vy0: yMin, vy1: yMax,
              xMin: plot.xMin, xMax: plot.xMax };

  // rubber-band rectangle while box-zooming
  if (vizDrag && vizDrag.mode === "box") {
    const rx = Math.min(vizDrag.x0, vizDrag.x1);
    const ry = Math.min(vizDrag.y0, vizDrag.y1);
    const rw = Math.abs(vizDrag.x1 - vizDrag.x0);
    const rh = Math.abs(vizDrag.y1 - vizDrag.y0);
    ctx.fillStyle = "rgba(31,95,166,0.12)";
    ctx.fillRect(rx, ry, rw, rh);
    ctx.strokeStyle = "rgba(31,95,166,0.85)";
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 3]);
    ctx.strokeRect(rx, ry, rw, rh);
    ctx.setLineDash([]);
  }

  $("vizHint").textContent =
    plot.info + (vizView || vizYView ? "  ·  zoomed (Reset to clear)" : "");
}

/* ---- drag tools: box-zoom and pan -------------------------------------- */

/* Toggle a drag tool; clicking the active tool turns it off. */
function vizSetMode(mode) {
  vizMode = vizMode === mode ? "none" : mode;
  $("vizBoxZoom").classList.toggle("active", vizMode === "box");
  $("vizPan").classList.toggle("active", vizMode === "pan");
  const cv = $("vizCanvas");
  cv.classList.toggle("mode-box", vizMode === "box");
  cv.classList.toggle("mode-pan", vizMode === "pan");
}

function vizPointer(ev) {
  const r = $("vizCanvas").getBoundingClientRect();
  return { x: ev.clientX - r.left, y: ev.clientY - r.top };
}

function vizOnDown(ev) {
  if (vizMode === "none" || !vizGeom) return;
  const p = vizPointer(ev);
  vizDrag = {
    mode: vizMode,
    x0: p.x, y0: p.y, x1: p.x, y1: p.y,
    g: { ...vizGeom },
    yStart: vizYView ? { ...vizYView } : null,
  };
  ev.preventDefault();
  window.addEventListener("mousemove", vizOnMove);
  window.addEventListener("mouseup", vizOnUp);
}

function vizOnMove(ev) {
  if (!vizDrag) return;
  const p = vizPointer(ev);
  vizDrag.x1 = p.x;
  vizDrag.y1 = p.y;
  if (vizDrag.mode === "pan") {
    const g = vizDrag.g;
    const span = g.vx1 - g.vx0;
    const full = g.xMax - g.xMin;
    let a = g.vx0 - ((p.x - vizDrag.x0) / g.pW) * span;
    let b = a + span;
    if (a < g.xMin) { b += g.xMin - a; a = g.xMin; }
    if (b > g.xMax) { a -= b - g.xMax; b = g.xMax; }
    a = Math.max(g.xMin, a);
    vizView = b - a >= full - 1e-9 ? null : { x0: a, x1: b };
    if (vizDrag.yStart) {
      const shift = ((p.y - vizDrag.y0) / g.pH) * (g.vy1 - g.vy0);
      vizYView = { y0: vizDrag.yStart.y0 + shift, y1: vizDrag.yStart.y1 + shift };
    }
  }
  drawViz();
}

function vizOnUp() {
  window.removeEventListener("mousemove", vizOnMove);
  window.removeEventListener("mouseup", vizOnUp);
  if (!vizDrag) return;
  if (vizDrag.mode === "box") {
    const g = vizDrag.g;
    const xa = Math.max(g.mL, Math.min(vizDrag.x0, vizDrag.x1));
    const xb = Math.min(g.mL + g.pW, Math.max(vizDrag.x0, vizDrag.x1));
    const ya = Math.max(g.mT, Math.min(vizDrag.y0, vizDrag.y1));
    const yb = Math.min(g.mT + g.pH, Math.max(vizDrag.y0, vizDrag.y1));
    if (xb - xa > 6 && yb - ya > 6) {       // ignore tiny accidental drags
      vizView = {
        x0: g.vx0 + ((xa - g.mL) / g.pW) * (g.vx1 - g.vx0),
        x1: g.vx0 + ((xb - g.mL) / g.pW) * (g.vx1 - g.vx0),
      };
      // screen y is inverted: the top edge (ya) maps to the larger value
      vizYView = {
        y0: g.vy0 + ((g.mT + g.pH - yb) / g.pH) * (g.vy1 - g.vy0),
        y1: g.vy0 + ((g.mT + g.pH - ya) / g.pH) * (g.vy1 - g.vy0),
      };
    }
  }
  vizDrag = null;
  drawViz();
}

/* Zoom the visualization x-axis by a factor about the window centre. */
function vizZoom(factor) {
  if (!vizGeom) return;
  const g = vizGeom;
  const full = g.xMax - g.xMin;
  let a = vizView ? vizView.x0 : g.xMin;
  let b = vizView ? vizView.x1 : g.xMax;
  const c = (a + b) / 2;
  let span = (b - a) * factor;
  span = Math.min(span, full);
  span = Math.max(span, full * 0.002);
  a = c - span / 2;
  b = c + span / 2;
  if (a < g.xMin) { b += g.xMin - a; a = g.xMin; }
  if (b > g.xMax) { a -= b - g.xMax; b = g.xMax; }
  a = Math.max(g.xMin, a);
  vizView = b - a >= full - 1e-9 ? null : { x0: a, x1: b };
  drawViz();
}

function updateSpacingReadout() {
  const s = getSettings();
  const sp = (p) =>
    `Gauge positions from gauge&nbsp;1: ` +
    `${p[0].toFixed(2)}, ${p[1].toFixed(2)}, ${p[2].toFixed(2)} m`;
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
  document.querySelectorAll("#fs, #skipWaves, #numWaves, .sp1, .sp2").forEach((el) =>
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

  // wave-type toggle -> adapt the UI and recompute everything
  document.querySelectorAll('input[name="wavemode"]').forEach((el) =>
    el.addEventListener("change", () => {
      applyMode();
      if (state.records.length) analyzeAll(true);
    })
  );

  // visualization controls
  $("vizFile").addEventListener("change", () => {
    vizView = null; vizYView = null; drawViz();
  });
  $("vizType").addEventListener("change", () => {
    vizView = null; vizYView = null;
    applyVizType();
    drawViz();
  });
  $("vizZoomIn").addEventListener("click", () => vizZoom(0.6));
  $("vizZoomOut").addEventListener("click", () => vizZoom(1 / 0.6));
  $("vizBoxZoom").addEventListener("click", () => vizSetMode("box"));
  $("vizPan").addEventListener("click", () => vizSetMode("pan"));
  $("vizCanvas").addEventListener("mousedown", vizOnDown);
  $("vizReset").addEventListener("click", () => {
    vizView = null; vizYView = null; drawViz();
  });
  document.querySelectorAll(".viz-probe").forEach((cb) =>
    cb.addEventListener("change", () => {
      $("vizAll").checked =
        [...document.querySelectorAll(".viz-probe")].every((c) => c.checked);
      drawViz();
    })
  );
  $("vizAll").addEventListener("change", () => {
    const on = $("vizAll").checked;
    document.querySelectorAll(".viz-probe").forEach((c) => { c.checked = on; });
    drawViz();
  });
  document.querySelectorAll(".viz-curve").forEach((cb) =>
    cb.addEventListener("change", () => {
      $("vizCurveAll").checked =
        [...document.querySelectorAll(".viz-curve")].every((c) => c.checked);
      drawViz();
    })
  );
  $("vizCurveAll").addEventListener("change", () => {
    const on = $("vizCurveAll").checked;
    document.querySelectorAll(".viz-curve").forEach((c) => { c.checked = on; });
    drawViz();
  });
  $("vizBox").addEventListener("toggle", () => {
    if ($("vizBox").open) drawViz();
  });
  window.addEventListener("resize", () => {
    if ($("vizBox").open) drawViz();
  });

  applyMode();
  applyVizType();
  updateSpacingReadout();
}

/* Show the probe checkboxes for the time-series plot, or the
 * spectrum-curve checkboxes for the energy-spectrum plot — never both.
 * (display is set inline because the .viz-probes class would otherwise
 * override the [hidden] attribute.) */
function applyVizType() {
  const spec = vizType() === "spectrum";
  if ($("vizProbes")) $("vizProbes").style.display = spec ? "none" : "flex";
  if ($("vizCurves")) $("vizCurves").style.display = spec ? "flex" : "none";
}

/* Show/hide regular-only settings and update the mode note + table headers. */
function applyMode() {
  const irregular = getMode() === "irregular";
  ["fld-fmin", "fld-fmax", "fld-skip", "fld-num"].forEach((id) => {
    const el = $(id);
    if (el) el.style.display = irregular ? "none" : "";
  });
  const note = $("modeNote");
  if (note) {
    note.innerHTML = irregular
      ? "Spectral method — incident/reflected separated at every frequency bin; "
        + "H<sub>i</sub>, H<sub>r</sub> are H<sub>m0</sub> significant wave heights."
      : "Single-frequency method — separation at the detected dominant wave frequency.";
  }
  // relabel the f / T column headers for the active mode
  const thF = $("th-f"), thT = $("th-T");
  if (thF) thF.innerHTML = irregular ? "<i>f</i><sub>p</sub> (Hz)" : "<i>f</i> (Hz)";
  if (thT) thT.innerHTML = irregular ? "<i>T</i><sub>p</sub> (s)" : "<i>T</i> (s)";
}

if (typeof document !== "undefined") {
  document.addEventListener("DOMContentLoaded", init);
}

/* expose core functions for Node-based testing */
if (typeof module !== "undefined" && module.exports) {
  module.exports = { threeProbe, dispersion, wavelength, detectFrequency, parseCSV, parseDepth };
}
