# Plan: Model Averaging via Stacking Weights

## Context

The project fits 12 competing Bayesian ODE models for sea urchin foraging preference and currently uses pseudo-BMA weights in `compare_models.R` **only for ranking** — no model averaging is performed on posteriors or predictions. Each model's ODE simulations and plots are produced independently. The goals are:

1. Replace pseudo-BMA with **stacking weights** (`method = "stacking"` in `loo::loo_model_weights()`), which are more robust when multiple similar models compete.
2. **Model-average the preference function** by pooling `(w, q)` draws proportionally to stacking weights, respecting each draw's model family (Logistic vs. vanLeeuwen).
3. **Model-average the ODE simulation plots** — produce the same 9-panel `ODE_simulation_*.pdf` format that exists for each specific model, but for the model-averaged posterior.
4. Add a **model-averaged row** to `Summary_preference.tex` (switch point, baseline preference, log-odds).
5. Add a new **preference function curve plot** (`preference_model_average.pdf`).

---

## Key Technical Insight: Mixing at the ODE Output Level

Different models use different ODE functions and parameter sets — they cannot share draws. However, each model's `outs_parms[[parm_idx]][[init_idx]]` in the cached `.RDA` files already contains **integrated ODE trajectories** — opaque outputs comparable across models. Model averaging works by sampling from these proportionally to stacking weights and pooling, then passing the result to the existing `process_model_sim()` → `plot_model_sim()` pipeline with `model_name = "model_average"`.

For **preference averaging**: pool `(w, q)` draws from each model proportionally, applying `preference(x, w, q, model_name_k)` with the correct family per draw.

Only models in `compare_names` = `setdiff(model_names, exclude_from_comparison)` are used — the four `vanLeeuwen_q` variants are excluded throughout.

---

## Files to Modify

### 1. `.gitignore` *(already done)*
Replace three specific entries with one glob:
```
results/**/ODE_kelp_*.RDA
```
This unblocks `results/tmp/` from being entirely ignored (so posteriors, hashes, processed summaries are trackable) while still ignoring the large raw ODE outputs.

### 2. `code/RunMe.R`
- **After the path variables block** (line 18): add `dir.create(tmp, recursive = TRUE, showWarnings = FALSE)`
- **After the four `run_stage()` calls** (line 179): add:
  ```r
  source('model_average.R')
  if (exists("stacking_wts") && length(compare_names) >= 2) {
    compute_model_average(stacking_wts, compare_names, fits)
  }
  ```

### 3. `code/compare_models.R`
Stacking swap only — no model-averaging code goes here:
- Line 156: `method = "pseudobma"` → `method = "stacking"`
- Rename `pbma_wts` → `stacking_wts` (lines 156, 163, 173, 174, 175, 184)
- Line 204 caption: `"pseudo-BMA method"` → `"stacking method"`
- Line 347 caption: `"pseudo-BMA weight"` → `"stacking weight"`

---

## Files to Create

### 4. `code/model_average.R`
Single function `compute_model_average(stacking_wts, compare_names, fits, n_avg = n_sim_draws)` with four sections:

**A. Pool ODE outputs** — for each kelp level, load each model's `outs_parms`, sample `round(w_k × n_avg)` indices (with rounding correction on the largest-weight model), concatenate, save as `ODE_kelp_<level>_model_average.RDA`, then call `process_model_sim("model_average")` and `plot_model_sim("model_average")`.

**B. Pool `(w, q)` draws** — same proportional sampling (`n_avg_pref = 10000` draws total), deriving `w = q - 4` for models without a sampled `w`, tagged with `model_name` for family dispatch.

**C. Model-averaged preference statistics** — compute `preference(0, w, q, model_name)`, `log_switch_point(0.5, ...)`, and log-odds across pooled draws; format with `fmt_med_ci()`; append "Model average" row to `Summary_preference.tex`.

**D. Preference curve plot** — evaluate `preference(x_grid, w, q, model_name)` for each pooled draw over `x_grid = seq(-5, 5, length.out = 200)`; summarize as median + 95% CI; plot with `ggplot2`; save as `figs/preference_model_average.pdf`.

**Scope notes:** `fmt_med_ci` (local to `compare_models.R`) and `my.theme` (local to `plot_model_sim()`) must each be either redefined inline in `model_average.R` or moved to a shared location. Simplest: redefine both directly in the new script.

---

## New Outputs

| Output | Description |
|---|---|
| `figs/ODE_simulation_model_average.pdf` | 9-panel dynamics plot, same format as per-model plots |
| `figs/preference_model_average.pdf` | Preference curve (median + 95% CI) vs. log(drift/kelp) |
| `tables/Summary_model_comparison.tex` | Updated captions + stacking weights |
| `tables/Summary_preference.tex` | Updated caption + "Model average" row |
| `results/tmp/ODE_kelp_{low,high}_model_average.RDA` | Gitignored (raw pooled ODE outputs) |
| `results/tmp/ODE_toPlot_kelp_{low,high}_model_average.RDA` | Tracked (processed summaries) |

---

## Verification

1. Run `RunMe.R` with `reuse_existing_fits = TRUE`, `reuse_existing_sims = TRUE`.
2. `results/tmp/` created if absent; `ODE_kelp_*.RDA` files gitignored.
3. `Summary_model_comparison.tex` caption says "stacking"; weights differ from pseudo-BMA.
4. `figs/ODE_simulation_model_average.pdf` — 9-panel layout with credible bands.
5. `figs/preference_model_average.pdf` — smooth curve with CI ribbon.
6. `Summary_preference.tex` — "Model average" row present.
