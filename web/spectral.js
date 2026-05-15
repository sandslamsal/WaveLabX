/* =========================================================================
 * spectral.js — irregular-wave (spectral) three-probe reflection analysis.
 *
 * Faithful JavaScript port of WaveLabX's wavelabx/three_probe_array
 * (Goda-style frequency-domain redundant-array method for random waves).
 *
 * The incident/reflected separation is performed at every frequency bin;
 * incident and reflected spectra are integrated to spectral wave heights
 * Hm0 = 4*sqrt(m0).
 * ========================================================================= */

"use strict";

const SP_G = 9.81;
const SP_COND_MAX = 1e6;   // mask a pair's frequency if cond exceeds this

/* ---------------------------------------------------------------------------
 * Radix-2 iterative FFT (length must be a power of two), in place.
 * ------------------------------------------------------------------------- */
function fftPow2(re, im, inverse) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      let t = re[i]; re[i] = re[j]; re[j] = t;
      t = im[i]; im[i] = im[j]; im[j] = t;
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = ((inverse ? 2 : -2) * Math.PI) / len;
    const wlr = Math.cos(ang), wli = Math.sin(ang);
    const half = len >> 1;
    for (let i = 0; i < n; i += len) {
      let wr = 1, wi = 0;
      for (let k = 0; k < half; k++) {
        const ur = re[i + k], ui = im[i + k];
        const xr = re[i + k + half], xi = im[i + k + half];
        const vr = xr * wr - xi * wi;
        const vi = xr * wi + xi * wr;
        re[i + k] = ur + vr; im[i + k] = ui + vi;
        re[i + k + half] = ur - vr; im[i + k + half] = ui - vi;
        const nwr = wr * wlr - wi * wli;
        wi = wr * wli + wi * wlr;
        wr = nwr;
      }
    }
  }
  if (inverse) {
    for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  }
}

/* ---------------------------------------------------------------------------
 * Bluestein DFT of a real signal of arbitrary length N.
 * Returns the full complex spectrum {re, im} of length N — identical to
 * numpy.fft.fft for a real input.
 * ------------------------------------------------------------------------- */
function dftReal(x) {
  const N = x.length;
  const cos = new Float64Array(N), sin = new Float64Array(N);
  for (let n = 0; n < N; n++) {
    const m = (n * n) % (2 * N);          // keep the angle argument small
    const ang = (Math.PI * m) / N;
    cos[n] = Math.cos(ang);
    sin[n] = Math.sin(ang);
  }
  let M = 1;
  while (M < 2 * N - 1) M <<= 1;

  // a[n] = x[n] * exp(-i*pi*n^2/N)
  const are = new Float64Array(M), aim = new Float64Array(M);
  for (let n = 0; n < N; n++) {
    are[n] = x[n] * cos[n];
    aim[n] = -x[n] * sin[n];
  }
  // b[n] = exp(+i*pi*n^2/N), symmetric: b[M-n] = b[n]
  const bre = new Float64Array(M), bim = new Float64Array(M);
  bre[0] = cos[0]; bim[0] = sin[0];
  for (let n = 1; n < N; n++) {
    bre[n] = cos[n]; bim[n] = sin[n];
    bre[M - n] = cos[n]; bim[M - n] = sin[n];
  }
  fftPow2(are, aim, false);
  fftPow2(bre, bim, false);
  for (let i = 0; i < M; i++) {
    const r = are[i] * bre[i] - aim[i] * bim[i];
    const ii = are[i] * bim[i] + aim[i] * bre[i];
    are[i] = r; aim[i] = ii;
  }
  fftPow2(are, aim, true);

  // X[k] = c[k] * exp(-i*pi*k^2/N)
  const Xre = new Float64Array(N), Xim = new Float64Array(N);
  for (let k = 0; k < N; k++) {
    Xre[k] = are[k] * cos[k] + aim[k] * sin[k];
    Xim[k] = -are[k] * sin[k] + aim[k] * cos[k];
  }
  return { re: Xre, im: Xim };
}

/* ---------------------------------------------------------------------------
 * Linear wavelength — faithful port of wavelabx/core.compute_wavelength
 * (fixed-point iteration L_{n+1} = (g T^2 / 2pi) tanh(k h)).
 * ------------------------------------------------------------------------- */
function spComputeWavelength(h, T) {
  const L0 = (SP_G * T * T) / (2 * Math.PI);
  let L = L0;
  for (let it = 0; it < 50; it++) {
    const k = (2 * Math.PI) / L;
    const Lnext = L0 * Math.tanh(k * h);
    if (Math.abs(Lnext - L) < 1e-5) { L = Lnext; break; }
    L = Lnext;
  }
  return L;
}

/* mean of up to three values, ignoring NaN (numpy nanmean over a pair axis) */
function nanmean3(a, b, c) {
  let s = 0, n = 0;
  if (!Number.isNaN(a)) { s += a; n++; }
  if (!Number.isNaN(b)) { s += b; n++; }
  if (!Number.isNaN(c)) { s += c; n++; }
  return n ? s / n : NaN;
}

/* ---------------------------------------------------------------------------
 * three_probe_array — spectral incident/reflected separation for one
 * three-gauge array.
 *
 *   cols : array of three number[]  (gauge time series, equal length)
 *   fs   : sampling frequency [Hz]
 *   h    : water depth [m]
 *   pos  : [x1,x2,x3] gauge positions [m]
 *
 * Returns { Hi, Hr, Kr, Htot, retained, Tp, fp, validRange, badCond }.
 * Hi, Hr are spectral Hm0 wave heights.
 * ------------------------------------------------------------------------- */
function threeProbeArray(cols, fs, h, pos) {
  const N = cols[0].length;
  const dt = 1.0 / fs;
  const df = fs / N;
  const half = N >> 1;
  const nf = half - 1;                 // frequency bins 1 .. half-1
  if (nf < 2) throw new Error("record too short for spectral analysis");

  // Detrend (remove mean) and FFT each gauge
  const An = [new Float64Array(nf), new Float64Array(nf), new Float64Array(nf)];
  const Bn = [new Float64Array(nf), new Float64Array(nf), new Float64Array(nf)];
  const Sn = [new Float64Array(nf), new Float64Array(nf), new Float64Array(nf)];
  for (let j = 0; j < 3; j++) {
    const col = cols[j];
    let mean = 0;
    for (let n = 0; n < N; n++) mean += col[n];
    mean /= N;
    const z = new Float64Array(N);
    for (let n = 0; n < N; n++) z[n] = col[n] - mean;
    const F = dftReal(z);
    for (let i = 0; i < nf; i++) {
      const re = F.re[i + 1], im = F.im[i + 1];
      An[j][i] = (2 * re) / N;
      Bn[j][i] = (-2 * im) / N;
      Sn[j][i] = (dt * 2 * (re * re + im * im)) / N;
    }
  }

  // Frequency vector and wavenumber per bin
  const f = new Float64Array(nf);
  const k = new Float64Array(nf);
  for (let i = 0; i < nf; i++) {
    f[i] = df * (i + 1);
    const L = spComputeWavelength(h, 1.0 / f[i]);
    k[i] = L > 0 ? (2 * Math.PI) / L : 0;
  }

  // Per-pair incident/reflected Fourier coefficients
  const g1 = [0, 0, 1], g2 = [1, 2, 2];
  const Ainc = [], Binc = [], Aref = [], Bref = [];
  for (let j = 0; j < 3; j++) {
    Ainc.push(new Float64Array(nf)); Binc.push(new Float64Array(nf));
    Aref.push(new Float64Array(nf)); Bref.push(new Float64Array(nf));
  }
  const nmin = [0, 0, 0], nmax = [nf - 1, nf - 1, nf - 1];

  for (let j = 0; j < 3; j++) {
    const A1 = An[g1[j]], A2 = An[g2[j]];
    const B1 = Bn[g1[j]], B2 = Bn[g2[j]];
    const pos1 = pos[g1[j]], pos2 = pos[g2[j]];
    const dx = Math.abs(pos2 - pos1);

    for (let i = 0; i < nf; i++) {
      const ki = k[i];
      const s1 = Math.sin(ki * pos1), c1 = Math.cos(ki * pos1);
      const s2 = Math.sin(ki * pos2), c2 = Math.cos(ki * pos2);
      const t1 = -A2[i] * s1 + A1[i] * s2 + B2[i] * c1 - B1[i] * c2;
      const t2 =  A2[i] * c1 - A1[i] * c2 + B2[i] * s1 - B1[i] * s2;
      const t3 = -A2[i] * s1 + A1[i] * s2 - B2[i] * c1 + B1[i] * c2;
      const t4 =  A2[i] * c1 - A1[i] * c2 - B2[i] * s1 + B1[i] * s2;
      const denom = 2.0 * Math.sin(ki * dx + 1e-16);
      Ainc[j][i] = t1 / denom;
      Binc[j][i] = t2 / denom;
      Aref[j][i] = t3 / denom;
      Bref[j][i] = t4 / denom;

      // conditioning of the 2x2 pair matrix; mask if extreme
      const sr = Math.cos(2 * ki * pos1) + Math.cos(2 * ki * pos2);
      const si = Math.sin(2 * ki * pos1) + Math.sin(2 * ki * pos2);
      const sAbs = Math.hypot(sr, si);
      const cond = sAbs >= 2 ? Infinity
        : Math.sqrt((2 + sAbs) / (2 - sAbs));
      if (cond > SP_COND_MAX) {
        Ainc[j][i] = Binc[j][i] = Aref[j][i] = Bref[j][i] = NaN;
      }
    }

    // valid wavelength band (Goda guideline)  0.05 <= dx/L <= 0.45
    const kmax = (2 * Math.PI) / (dx / 0.45);
    const kmin = (2 * Math.PI) / (dx / 0.05);
    let lo = -1, hi = -1;
    for (let i = 0; i < nf; i++) {
      if (k[i] >= kmin && k[i] <= kmax) { if (lo < 0) lo = i; hi = i; }
    }
    if (lo < 0) { nmin[j] = 0; nmax[j] = nf - 1; }
    else { nmin[j] = lo; nmax[j] = hi; }

    // mask bins outside this pair's valid range
    for (let i = 0; i < nf; i++) {
      if (i < nmin[j] || i > nmax[j]) {
        Ainc[j][i] = Binc[j][i] = Aref[j][i] = Bref[j][i] = NaN;
      }
    }
  }

  const gMin = Math.min(nmin[0], nmin[1], nmin[2]);
  const gMax = Math.max(nmax[0], nmax[1], nmax[2]);

  // retained-energy diagnostic (gauge-1 spectrum; first bin zeroed)
  let totalE = 0, retainedE = 0;
  for (let i = 1; i < nf; i++) {
    totalE += Sn[0][i];
    if (i >= gMin && i <= gMax) retainedE += Sn[0][i];
  }
  const retained = totalE > 0 ? retainedE / totalE : NaN;

  // pair-averaged spectra, integrated energy
  let Ei = 0, Er = 0, mo = 0;
  let peakS = -1, peakIdx = gMin;
  for (let i = gMin; i <= gMax; i++) {
    const ai = nanmean3(Ainc[0][i], Ainc[1][i], Ainc[2][i]);
    const bi = nanmean3(Binc[0][i], Binc[1][i], Binc[2][i]);
    const ar = nanmean3(Aref[0][i], Aref[1][i], Aref[2][i]);
    const br = nanmean3(Bref[0][i], Bref[1][i], Bref[2][i]);
    const Si = (ai * ai + bi * bi) / (2 * df);
    const Sr = (ar * ar + br * br) / (2 * df);
    if (!Number.isNaN(Si)) { Ei += Si; if (Si > peakS) { peakS = Si; peakIdx = i; } }
    if (!Number.isNaN(Sr)) Er += Sr;
    mo += Sn[0][i];
  }
  Ei *= df; Er *= df; mo *= df;

  const Hi = 4 * Math.sqrt(Math.max(Ei, 0));
  const Hr = 4 * Math.sqrt(Math.max(Er, 0));
  const Htot = 4 * Math.sqrt(Math.max(mo, 0));
  const Kr = Hi > 0 ? Hr / (Hi + 1e-16) : NaN;
  const fp = f[peakIdx];

  return {
    Hi, Hr, Kr, Htot, retained,
    fp, Tp: fp > 0 ? 1 / fp : NaN,
    validRange: [gMin, gMax],
  };
}

/* exports for browser + Node */
if (typeof window !== "undefined") {
  window.WaveLabXSpectral = { threeProbeArray, spComputeWavelength, dftReal };
}
if (typeof module !== "undefined" && module.exports) {
  module.exports = { threeProbeArray, spComputeWavelength, dftReal, fftPow2 };
}
