// Deformation experiment replayer for darboux-lock-simulator.
// Loads run artifacts produced by `darboux-igph` (manifest.json +
// per-run web_bundle.json) and plays them back on a 2D canvas.

const ARTIFACTS_ROOT = "./artifacts_out";
const MANIFEST_URL = `${ARTIFACTS_ROOT}/manifest.json`;
const CURVE_SAMPLES = 256;
const KNOT_STACK_THRESHOLD = -3.0; // count gamma_i < this as "stacked"

const COLORS = {
  background: "#07111d",
  spline: "#66aefc",
  control: "#ff7272",
  target: "rgba(255, 156, 67, 0.55)",
  targetCore: "#ff9c43",
  grid: "rgba(193, 216, 248, 0.08)",
  axis: "rgba(193, 216, 248, 0.18)",
  // Continuity coloring: green C² → amber C¹ → red C⁰. Score s ∈ [0, 1] is
  // produced by the Python continuity_score metric (1.0 ≈ smooth).
  contHigh: [93, 214, 167],   // #5DD6A7
  contMid:  [245, 192, 97],   // #F5C061
  contLow:  [255, 106, 106],  // #FF6A6A
  // Knot heatmap palette: bright cyan = full interval, deep purple = collapsed.
  knotFull: [180, 222, 255],
  knotCollapsed: [124, 58, 237],
};

const els = {
  canvas: document.getElementById("deformation-canvas"),
  renderState: document.getElementById("deformation-state"),
  badgeState: document.getElementById("def-badge-state"),
  badgeStep: document.getElementById("def-badge-step"),
  badgeKnots: document.getElementById("def-badge-knots"),
  sinkhornValue: document.getElementById("def-sinkhorn-value"),
  kineticValue: document.getElementById("def-kinetic-value"),
  hamiltonianValue: document.getElementById("def-hamiltonian-value"),
  gradValue: document.getElementById("def-grad-value"),
  progressBar: document.getElementById("def-progress-bar"),
  progressPercent: document.getElementById("def-progress-percent"),
  runSelect: document.getElementById("def-run-select"),
  methodLabel: document.getElementById("def-method-label"),
  slider: document.getElementById("def-progress-slider"),
  playBtn: document.getElementById("def-play-btn"),
  resetBtn: document.getElementById("def-reset-btn"),
  runInfo: document.getElementById("def-run-info"),
  flowEq: document.getElementById("def-flow-equation"),
  curveEq: document.getElementById("def-curve-equation"),
  knotStrip: document.getElementById("def-knot-strip"),
};

const state = {
  manifest: null,
  bundle: null,
  frameIndex: 0,
  playing: false,
  lastStamp: 0,
  framesPerSecond: 30,
  // accumulator across rAF ticks so we can advance one frame at a time
  frameAccumulator: 0,
};

const ctx = els.canvas ? els.canvas.getContext("2d") : null;
const knotCtx = els.knotStrip ? els.knotStrip.getContext("2d") : null;

// --- Closed cubic uniform B-spline (port of evaluate_curve from
// darboux_igph/splines/bspline.py). controlPts is an array of [x, y].
// Returns an array of [x, y] samples evenly spaced in parameter [0, n).
function evaluateClosedCubicBSpline(controlPts, nSamples) {
  const n = controlPts.length;
  if (n < 4) {
    throw new Error(`closed cubic B-spline needs >= 4 control points; got ${n}`);
  }
  const out = new Array(nSamples);
  for (let s = 0; s < nSamples; s += 1) {
    const u = (s / nSamples) * n;
    const iFloor = ((Math.floor(u) % n) + n) % n;
    const t = u - Math.floor(u);
    const oneMinusT = 1.0 - t;
    const t2 = t * t;
    const t3 = t2 * t;
    const b0 = (oneMinusT * oneMinusT * oneMinusT) / 6.0;
    const b1 = (3.0 * t3 - 6.0 * t2 + 4.0) / 6.0;
    const b2 = (-3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0) / 6.0;
    const b3 = t3 / 6.0;
    const i0 = ((iFloor - 1) % n + n) % n;
    const i1 = iFloor;
    const i2 = (iFloor + 1) % n;
    const i3 = (iFloor + 2) % n;
    const p0 = controlPts[i0];
    const p1 = controlPts[i1];
    const p2 = controlPts[i2];
    const p3 = controlPts[i3];
    out[s] = [
      b0 * p0[0] + b1 * p1[0] + b2 * p2[0] + b3 * p3[0],
      b0 * p0[1] + b1 * p1[1] + b2 * p2[1] + b3 * p3[1],
    ];
  }
  return out;
}

// --- Resizing the canvas backing store to its CSS size with devicePixelRatio.
function resizeCanvasToDisplaySize() {
  if (!els.canvas) return;
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const width = els.canvas.clientWidth;
  const height = els.canvas.clientHeight;
  if (!width || !height) return;
  const desiredW = Math.round(width * dpr);
  const desiredH = Math.round(height * dpr);
  if (els.canvas.width !== desiredW || els.canvas.height !== desiredH) {
    els.canvas.width = desiredW;
    els.canvas.height = desiredH;
  }
}

function computeBoundsFromBundle(bundle) {
  let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
  const target = bundle.target;
  for (let i = 0; i < target.length; i += 1) {
    const [x, y] = target[i];
    if (x < minX) minX = x;
    if (y < minY) minY = y;
    if (x > maxX) maxX = x;
    if (y > maxY) maxY = y;
  }
  const q = bundle.trajectory.q;
  for (let f = 0; f < q.length; f += 1) {
    const frame = q[f];
    for (let i = 0; i < frame.length; i += 1) {
      const [x, y] = frame[i];
      if (x < minX) minX = x;
      if (y < minY) minY = y;
      if (x > maxX) maxX = x;
      if (y > maxY) maxY = y;
    }
  }
  // pad a bit
  const padX = (maxX - minX) * 0.12 + 0.05;
  const padY = (maxY - minY) * 0.12 + 0.05;
  return {
    minX: minX - padX,
    maxX: maxX + padX,
    minY: minY - padY,
    maxY: maxY + padY,
  };
}

function makeProjector(bounds, width, height) {
  const dataW = bounds.maxX - bounds.minX;
  const dataH = bounds.maxY - bounds.minY;
  if (dataW <= 0 || dataH <= 0) {
    return (x, y) => [width / 2, height / 2];
  }
  // Preserve aspect ratio.
  const scale = Math.min(width / dataW, height / dataH);
  const offsetX = (width - dataW * scale) / 2;
  const offsetY = (height - dataH * scale) / 2;
  return (x, y) => {
    const px = offsetX + (x - bounds.minX) * scale;
    // Flip y because canvas y grows downward; data y grows upward.
    const py = height - (offsetY + (y - bounds.minY) * scale);
    return [px, py];
  };
}

function drawBackground(width, height, project, bounds) {
  ctx.fillStyle = COLORS.background;
  ctx.fillRect(0, 0, width, height);

  // subtle grid lines aligned to integer data coords if range is reasonable
  ctx.strokeStyle = COLORS.grid;
  ctx.lineWidth = 1;
  ctx.beginPath();
  const stepX = niceStep(bounds.maxX - bounds.minX);
  const stepY = niceStep(bounds.maxY - bounds.minY);
  for (let x = Math.ceil(bounds.minX / stepX) * stepX; x < bounds.maxX; x += stepX) {
    const [px0] = project(x, bounds.minY);
    ctx.moveTo(px0, 0);
    ctx.lineTo(px0, height);
  }
  for (let y = Math.ceil(bounds.minY / stepY) * stepY; y < bounds.maxY; y += stepY) {
    const [, py0] = project(bounds.minX, y);
    ctx.moveTo(0, py0);
    ctx.lineTo(width, py0);
  }
  ctx.stroke();

  // axes through origin if visible
  ctx.strokeStyle = COLORS.axis;
  ctx.beginPath();
  if (bounds.minX < 0 && bounds.maxX > 0) {
    const [px0] = project(0, bounds.minY);
    ctx.moveTo(px0, 0);
    ctx.lineTo(px0, height);
  }
  if (bounds.minY < 0 && bounds.maxY > 0) {
    const [, py0] = project(bounds.minX, 0);
    ctx.moveTo(0, py0);
    ctx.lineTo(width, py0);
  }
  ctx.stroke();
}

function niceStep(span) {
  if (span <= 0 || !Number.isFinite(span)) return 1;
  const rough = span / 6;
  const pow10 = Math.pow(10, Math.floor(Math.log10(rough)));
  const norm = rough / pow10;
  let nice;
  if (norm < 1.5) nice = 1;
  else if (norm < 3) nice = 2;
  else if (norm < 7) nice = 5;
  else nice = 10;
  return nice * pow10;
}

// Linear blend between two [r,g,b] colors at t ∈ [0, 1].
function lerpRGB(a, b, t) {
  const u = Math.max(0, Math.min(1, t));
  return [
    Math.round(a[0] + (b[0] - a[0]) * u),
    Math.round(a[1] + (b[1] - a[1]) * u),
    Math.round(a[2] + (b[2] - a[2]) * u),
  ];
}

// Score s in [0, 1]: 1.0 ≈ C²-clean, 0.0 ≈ C⁰ break. Map green → amber → red.
function continuityColor(score) {
  const s = Math.max(0, Math.min(1, score));
  let rgb;
  if (s >= 0.5) {
    // High to mid: green → amber
    rgb = lerpRGB(COLORS.contMid, COLORS.contHigh, (s - 0.5) * 2);
  } else {
    // Mid to low: amber → red
    rgb = lerpRGB(COLORS.contLow, COLORS.contMid, s * 2);
  }
  return `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
}

function drawFrame(bundle, frameIndex) {
  if (!ctx) return;
  resizeCanvasToDisplaySize();
  const width = els.canvas.width;
  const height = els.canvas.height;

  const bounds = bundle._bounds;
  const project = makeProjector(bounds, width, height);

  drawBackground(width, height, project, bounds);

  // Target cloud
  ctx.fillStyle = COLORS.target;
  for (let i = 0; i < bundle.target.length; i += 1) {
    const [x, y] = bundle.target[i];
    const [px, py] = project(x, y);
    ctx.beginPath();
    ctx.arc(px, py, 1.6, 0, Math.PI * 2);
    ctx.fill();
  }

  // Spline curve — prefer Python-rendered samples from frames.curve (uses
  // the non-uniform evaluator when knot_adaptive is on). Fallback to the JS
  // uniform formula for older bundles without frames.
  const framesCurve = bundle.frames?.curve;
  const samples = framesCurve
    ? framesCurve[frameIndex]
    : evaluateClosedCubicBSpline(bundle.trajectory.q[frameIndex], CURVE_SAMPLES);

  const contScores = bundle.frames?.continuity_score?.[frameIndex];

  ctx.lineWidth = 2.6;
  ctx.lineJoin = "round";
  ctx.lineCap = "round";

  if (contScores && contScores.length === samples.length) {
    // Per-segment coloring by local continuity score.
    for (let i = 0; i < samples.length; i += 1) {
      const a = samples[i];
      const b = samples[(i + 1) % samples.length];
      const score = 0.5 * (contScores[i] + contScores[(i + 1) % samples.length]);
      const [pax, pay] = project(a[0], a[1]);
      const [pbx, pby] = project(b[0], b[1]);
      ctx.strokeStyle = continuityColor(score);
      ctx.beginPath();
      ctx.moveTo(pax, pay);
      ctx.lineTo(pbx, pby);
      ctx.stroke();
    }
  } else {
    // Single-color fallback for older bundles.
    ctx.strokeStyle = COLORS.spline;
    ctx.beginPath();
    for (let i = 0; i <= samples.length; i += 1) {
      const [x, y] = samples[i % samples.length];
      const [px, py] = project(x, y);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    }
    ctx.stroke();
  }

  // Control polygon (faint connecting line + red dots)
  const controlPts = bundle.trajectory.q[frameIndex];
  ctx.strokeStyle = "rgba(255, 114, 114, 0.32)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  for (let i = 0; i <= controlPts.length; i += 1) {
    const [x, y] = controlPts[i % controlPts.length];
    const [px, py] = project(x, y);
    if (i === 0) ctx.moveTo(px, py);
    else ctx.lineTo(px, py);
  }
  ctx.stroke();
  ctx.fillStyle = COLORS.control;
  for (let i = 0; i < controlPts.length; i += 1) {
    const [x, y] = controlPts[i];
    const [px, py] = project(x, y);
    ctx.beginPath();
    ctx.arc(px, py, 3.2, 0, Math.PI * 2);
    ctx.fill();
  }
}

// --- Knot strip heatmap ---------------------------------------------------
// X axis = parameter index 0..n_ctrl-1; Y axis = time (top = step 0, bottom =
// step T). Color brightness = Δ_i value (bright = full interval, dark =
// collapsed). A horizontal pin marks the current frame.
function resizeKnotStripToDisplaySize() {
  if (!els.knotStrip) return;
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = els.knotStrip.clientWidth;
  const h = els.knotStrip.clientHeight;
  if (!w || !h) return;
  const desiredW = Math.round(w * dpr);
  const desiredH = Math.round(h * dpr);
  if (els.knotStrip.width !== desiredW || els.knotStrip.height !== desiredH) {
    els.knotStrip.width = desiredW;
    els.knotStrip.height = desiredH;
  }
}

function drawKnotStrip(bundle, frameIndex) {
  if (!knotCtx || !els.knotStrip) return;
  resizeKnotStripToDisplaySize();
  const width = els.knotStrip.width;
  const height = els.knotStrip.height;

  knotCtx.fillStyle = COLORS.background;
  knotCtx.fillRect(0, 0, width, height);

  const intervals = bundle.frames?.knot_intervals;
  const gammaFrames = bundle.trajectory.gamma;
  if (!intervals && !gammaFrames) return;

  // Source heatmap data: prefer pre-computed intervals; fallback to gamma.
  const data = intervals
    || gammaFrames.map((g) => g.map((v) => Math.log1p(Math.exp(-Math.abs(v))) + Math.max(v, 0)));

  const nFrames = data.length;
  const nCtrl = data[0].length;

  // Normalize Δ values per-bundle so the heatmap has good contrast: map
  // [Δ_min, Δ_max] → [0, 1] for color interpolation.
  let dMin = Infinity, dMax = -Infinity;
  for (let f = 0; f < nFrames; f += 1) {
    for (let i = 0; i < nCtrl; i += 1) {
      const v = data[f][i];
      if (v < dMin) dMin = v;
      if (v > dMax) dMax = v;
    }
  }
  const dRange = Math.max(dMax - dMin, 1e-6);

  const cellW = width / nCtrl;
  const cellH = height / nFrames;

  for (let f = 0; f < nFrames; f += 1) {
    for (let i = 0; i < nCtrl; i += 1) {
      const v = data[f][i];
      const t = (v - dMin) / dRange; // 0 = collapsed, 1 = full
      const rgb = lerpRGB(COLORS.knotCollapsed, COLORS.knotFull, t);
      knotCtx.fillStyle = `rgb(${rgb[0]}, ${rgb[1]}, ${rgb[2]})`;
      knotCtx.fillRect(i * cellW, f * cellH, Math.ceil(cellW) + 1, Math.ceil(cellH) + 1);
    }
  }

  // Stack markers: highlight cells where γ_i crosses the stack threshold.
  if (gammaFrames) {
    for (let f = 0; f < nFrames; f += 1) {
      for (let i = 0; i < nCtrl; i += 1) {
        if (gammaFrames[f][i] < KNOT_STACK_THRESHOLD) {
          knotCtx.strokeStyle = "rgba(255, 255, 255, 0.85)";
          knotCtx.lineWidth = 1;
          knotCtx.strokeRect(i * cellW + 0.5, f * cellH + 0.5, cellW - 1, cellH - 1);
        }
      }
    }
  }

  // Current-frame indicator: bright horizontal line.
  const y = (frameIndex + 0.5) * cellH;
  knotCtx.strokeStyle = "rgba(255, 200, 80, 0.95)";
  knotCtx.lineWidth = 2;
  knotCtx.beginPath();
  knotCtx.moveTo(0, y);
  knotCtx.lineTo(width, y);
  knotCtx.stroke();

  // Subtle axis ticks: vertical grid lines every 4 control points.
  knotCtx.strokeStyle = "rgba(193, 216, 248, 0.10)";
  knotCtx.lineWidth = 1;
  for (let i = 0; i < nCtrl; i += 4) {
    knotCtx.beginPath();
    knotCtx.moveTo(i * cellW + 0.5, 0);
    knotCtx.lineTo(i * cellW + 0.5, height);
    knotCtx.stroke();
  }
}

function formatFloat(v) {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  const abs = Math.abs(v);
  if (abs === 0) return "0.000";
  if (abs >= 1000 || abs < 1e-3) return v.toExponential(2);
  return v.toFixed(4);
}

function countKnotStacks(gammaFrame) {
  let n = 0;
  for (let i = 0; i < gammaFrame.length; i += 1) {
    if (gammaFrame[i] < KNOT_STACK_THRESHOLD) n += 1;
  }
  return n;
}

function updateMetrics(bundle, frameIndex) {
  const scalars = bundle.scalars;
  const total = bundle.trajectory.q.length;

  els.sinkhornValue.textContent = formatFloat(scalars.sinkhorn_divergence?.[frameIndex]);
  els.kineticValue.textContent = formatFloat(scalars.T_kinetic?.[frameIndex]);
  els.hamiltonianValue.textContent = formatFloat(scalars.H_d?.[frameIndex]);
  els.gradValue.textContent = formatFloat(scalars.grad_q_norm?.[frameIndex]);

  // Prefer pre-computed stack counts from the bundle; fall back to live count.
  const knots = bundle.frames?.knot_stack_count
    ? bundle.frames.knot_stack_count[frameIndex]
    : countKnotStacks(bundle.trajectory.gamma[frameIndex]);
  els.badgeKnots.textContent = String(knots);
  els.badgeStep.textContent = `${frameIndex + 1} / ${total}`;
  els.badgeState.textContent = state.playing
    ? "animating"
    : frameIndex === 0
      ? "idle"
      : frameIndex >= total - 1
        ? "settled"
        : "paused";

  const pct = total > 1 ? frameIndex / (total - 1) : 0;
  els.progressBar.style.transform = `scaleX(${pct})`;
  els.progressPercent.textContent = `${Math.round(100 * pct)}%`;
}

function setFrame(index) {
  if (!state.bundle) return;
  const total = state.bundle.trajectory.q.length;
  const clamped = Math.max(0, Math.min(total - 1, index | 0));
  state.frameIndex = clamped;
  els.slider.value = String(clamped);
  drawFrame(state.bundle, clamped);
  drawKnotStrip(state.bundle, clamped);
  updateMetrics(state.bundle, clamped);
}

function animationLoop(stamp) {
  if (!state.bundle) return;
  if (!state.playing) {
    state.lastStamp = 0;
    return;
  }
  if (!state.lastStamp) state.lastStamp = stamp;
  const delta = (stamp - state.lastStamp) / 1000;
  state.lastStamp = stamp;
  state.frameAccumulator += delta * state.framesPerSecond;
  let advanced = false;
  while (state.frameAccumulator >= 1) {
    state.frameAccumulator -= 1;
    let next = state.frameIndex + 1;
    if (next >= state.bundle.trajectory.q.length) {
      state.playing = false;
      els.playBtn.textContent = "Play";
      next = state.bundle.trajectory.q.length - 1;
      setFrame(next);
      return;
    }
    state.frameIndex = next;
    advanced = true;
  }
  if (advanced) setFrame(state.frameIndex);
  // Always schedule the next tick while playing; otherwise the loop dies
  // on the first frame (delta = 0, no advance, no reschedule).
  requestAnimationFrame(animationLoop);
}

function play() {
  if (!state.bundle) return;
  if (state.frameIndex >= state.bundle.trajectory.q.length - 1) {
    state.frameIndex = 0;
  }
  state.playing = true;
  state.lastStamp = 0;
  state.frameAccumulator = 0;
  els.playBtn.textContent = "Pause";
  requestAnimationFrame(animationLoop);
}

function pause() {
  state.playing = false;
  els.playBtn.textContent = "Play";
}

function togglePlay() {
  if (state.playing) pause();
  else play();
}

function resetPlayback() {
  pause();
  setFrame(0);
}

async function loadManifest() {
  const resp = await fetch(MANIFEST_URL, { cache: "no-store" });
  if (!resp.ok) throw new Error(`manifest fetch failed: ${resp.status}`);
  state.manifest = await resp.json();
  els.runSelect.innerHTML = "";
  const runs = (state.manifest.runs || []).filter((r) => r.web_bundle_ref);
  for (const run of runs) {
    const option = document.createElement("option");
    option.value = run.href;
    option.textContent = `${run.example}  ·  ${run.method}  ·  seed ${run.seed}`;
    option.dataset.bundleRef = run.web_bundle_ref;
    option.dataset.method = run.method;
    option.dataset.runId = run.run_id;
    els.runSelect.appendChild(option);
  }
  if (els.runSelect.options.length === 0) {
    setRenderState(
      "No runs in manifest.json have a web_bundle_ref. Re-run the experiment with --emit-web-bundle.",
      "error",
    );
  }
}

async function loadRun(href, bundleRef) {
  setRenderState(`Loading ${href}…`);
  const url = `${ARTIFACTS_ROOT}/${href}/${bundleRef}`;
  const resp = await fetch(url, { cache: "no-store" });
  if (!resp.ok) throw new Error(`bundle fetch failed: ${resp.status}`);
  const bundle = await resp.json();
  bundle._bounds = computeBoundsFromBundle(bundle);
  state.bundle = bundle;
  const total = bundle.trajectory.q.length;
  els.slider.min = "0";
  els.slider.max = String(total - 1);
  els.slider.value = "0";
  els.methodLabel.textContent = bundle.config?.method ?? "—";
  els.runInfo.innerHTML = formatRunInfo(bundle);
  setRenderState("", "ready");
  // Allow pre-positioning the timeline from a URL parameter (useful for
  // headless screenshots and deep-links).
  const params = new URLSearchParams(window.location.search);
  const requestedFrame = Number.parseInt(params.get("frame") ?? "0", 10);
  const startFrame = Number.isFinite(requestedFrame)
    ? Math.max(0, Math.min(total - 1, requestedFrame))
    : 0;
  setFrame(startFrame);
}

function formatRunInfo(bundle) {
  const cfg = bundle.config || {};
  const lines = [
    `run_id: ${bundle.run_id}`,
    `example: ${cfg.example}`,
    `method:  ${cfg.method}`,
    `seed:    ${cfg.seed}`,
    `n_ctrl:  ${cfg.spline?.n_ctrl}`,
    `n_steps: ${cfg.n_steps}`,
    `dt:      ${cfg.dt}`,
    `frames:  ${bundle.trajectory.q.length}`,
  ];
  return lines.map((l) => `<div>${l}</div>`).join("");
}

function setRenderState(message, mode) {
  if (!els.renderState) return;
  if (mode === "ready") {
    els.renderState.classList.add("ready");
    els.renderState.classList.remove("error");
    return;
  }
  els.renderState.classList.remove("ready");
  els.renderState.classList.toggle("error", mode === "error");
  els.renderState.textContent = message;
}

function bindUI() {
  els.runSelect.addEventListener("change", async () => {
    const option = els.runSelect.selectedOptions[0];
    if (!option) return;
    pause();
    try {
      await loadRun(option.value, option.dataset.bundleRef);
    } catch (err) {
      console.error(err);
      setRenderState(`Failed to load bundle: ${err.message}`, "error");
    }
  });

  els.slider.addEventListener("input", () => {
    pause();
    setFrame(Number.parseInt(els.slider.value, 10));
  });

  els.playBtn.addEventListener("click", togglePlay);
  els.resetBtn.addEventListener("click", resetPlayback);

  window.addEventListener("resize", () => {
    if (state.bundle) drawFrame(state.bundle, state.frameIndex);
  });

  // Redraw when this tab becomes active (canvas may have had zero size while hidden).
  window.addEventListener("darboux:tab-changed", (event) => {
    if (event.detail?.tab !== "deformation") return;
    // Wait a frame so layout settles before measuring canvas size.
    requestAnimationFrame(() => {
      if (state.bundle) drawFrame(state.bundle, state.frameIndex);
    });
  });
}

function renderStaticMath() {
  const katexApi = window.katex;
  const inlineNodes = document.querySelectorAll(
    '.deformation-layout [data-tex]',
  );
  const tryRender = () => {
    if (!window.katex) return false;
    inlineNodes.forEach((node) => {
      window.katex.render(node.dataset.tex, node, {
        throwOnError: false,
        displayMode: false,
      });
    });
    if (els.flowEq) {
      window.katex.render(
        String.raw`\dot q = M_q^{-1} p,\quad \dot p = -\nabla_q V - D_q\, M_q^{-1} p`,
        els.flowEq,
        { throwOnError: false, displayMode: true },
      );
    }
    if (els.curveEq) {
      window.katex.render(
        String.raw`c(u) = \sum_{i=0}^{n-1} B_{i,3}(u)\, q_i,\quad u \in [0, n)`,
        els.curveEq,
        { throwOnError: false, displayMode: true },
      );
    }
    return true;
  };
  if (!tryRender()) {
    // KaTeX script is deferred; retry shortly.
    const interval = setInterval(() => {
      if (tryRender()) clearInterval(interval);
    }, 200);
  }
}

async function init() {
  if (!ctx) return;
  bindUI();
  renderStaticMath();
  try {
    await loadManifest();
    const first = els.runSelect.options[0];
    if (first) {
      await loadRun(first.value, first.dataset.bundleRef);
    }
  } catch (err) {
    console.error(err);
    setRenderState(`Failed to load manifest: ${err.message}`, "error");
  }
}

init();
