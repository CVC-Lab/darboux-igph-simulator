# darboux-igph-simulator

Browser-based deformation explorer for the
[darboux-igph](https://github.com/CVC-Lab/darboux-igph) Python package.
Companion to *An Information-Geometric Port-Hamiltonian Framework for
Operator-Constrained Inverse Problems* (Gokhale & Bajaj, SIAM MDS26).

Live: **https://cvc-lab.github.io/darboux-igph-simulator/**

## What it shows

Loads versioned run artifacts produced by `darboux-igph` and replays them
on a 2D canvas. Per frame you get:

- **Spline + target canvas** with per-segment **continuity coloring**
  (green C² → amber C¹ → red C⁰). The curve is the Python-evaluated
  rendering (uniform or non-uniform Cox-de Boor as appropriate), not a
  client-side approximation.
- **Latent knot heatmap** below: Δ_i = softplus(γ_i) per parameter index
  (x-axis) over time (y-axis). White outlines mark stacked cells; an
  amber line tracks the current scrub position.
- **Live metric chips**: S_ε (Sinkhorn divergence), T(p) (kinetic energy),
  H_d (damped Hamiltonian), ‖∇_q V‖ (control-point gradient norm).
- **Run picker** populated from `manifest.json`, **timeline slider**,
  Play / Reset.
- **Top-right badges**: STATE / STEP / KNOT STACKS.
- **Side rail**: machine-readable run header + KaTeX-rendered flow law.

This repo intentionally does NOT include the SO(3) Casimir-leaf demo —
that's a separate concern at
[shubham0704/darboux-lock-simulator](https://github.com/shubham0704/darboux-lock-simulator).

## Run locally

```bash
git clone git@github.com:CVC-Lab/darboux-igph-simulator.git
cd darboux-igph-simulator
python3 -m http.server 8000
```

Open <http://localhost:8000/>.

The committed `artifacts_out/` directory contains the lightweight JSON
bundles needed to drive the page. To regenerate them from a sibling
`darboux-igph` checkout (after running new experiments):

```bash
./scripts/sync_artifacts.sh           # default source: ../darboux-igph
./scripts/sync_artifacts.sh /elsewhere/darboux-igph
```

The script pulls only `manifest.json` + per-run `web_bundle.json`. Heavy
files (`trajectory.npz`, `target.npz`, `preview.gif`, `run.json`) stay in
`darboux-igph` and are gitignored here.

## Repo layout

```
index.html        Static markup; loads deformation.js + KaTeX
styles.css        Dark theme, IBM Plex typography, full layout
deformation.js    Vanilla ES module — fetch / render / scrub / animate
artifacts_out/    Committed JSON bundles per run (+ manifest.json)
scripts/
  sync_artifacts.sh   Pull fresh bundles from a sibling darboux-igph checkout
```

No build step. No npm. Pure static site.

## Adding a new run

In a `darboux-igph` venv:

```bash
darboux-run darboux_igph/experiments/configs/<config>.yaml --animate
darboux-export-web artifacts_out
```

Then back here:

```bash
./scripts/sync_artifacts.sh
git add artifacts_out
git commit -m "Add <new-run> bundle"
git push
```

GitHub Pages auto-rebuilds within a minute.
