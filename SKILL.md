---
name: recoverlite
description: >
  Run a pre-data recovery test on a planned study design with the
  recoverlite R package. Use when a user asks whether a planned study can
  detect or recover its effect, requests a power analysis, sample-size
  check, or design feasibility review for a two-arm or cluster-randomized
  trial, or wants a PASS/RISK/FAIL design verdict with simulation
  diagnosands (power, bias, coverage, Type S/M, estimand drift, model
  failure) before data collection.
---

# Running recovery tests with recoverlite

recoverlite simulates a *declared* design-analysis pair over a crossed
scenario grid (null and target effects x declared and pessimistic
nuisance assumptions) and applies a fixed, versioned PASS/RISK/FAIL rule.
It answers "can this planned study recover its target estimand?", which
is strictly more than "does it have 80% power?".

## Setup

```r
# R >= 4.1; imports DeclareDesign/fabricatr/randomizr
# remotes::install_github("heidihelena/recoverlite")
library(recoverlite)
```

Mixed-model answer strategies additionally need `lme4` + `lmerTest`
(`pbkrtest` for Kenward-Roger inference).

## Core workflow (always this shape)

```r
design <- declare_recovery(
  target = target_estimand(
    estimand = "<one sentence: quantity, population, scale>",
    scale    = "latent-outcome standardized mean difference",
    sesoi    = 0.40                       # smallest effect size of interest
  ),
  data_strategy   = two_arm_trial(n_per_arm = 115),   # RECRUITED per arm
  measurement     = measured_outcome(reliability = 0.70),
  missingness     = attrition_model(rate = 0.15, mechanism = "differential"),
  answer_strategy = planned_analysis(
    estimator = "linear_model",
    formula   = y_observed ~ treatment
  )
)

result <- recovery_test(design, sims = 2000,
                        scenarios = "confirmatory_grid", seed = <int>)
verdict(result)    # PASS / RISK / FAIL + strict/lenient recomputation
report(result)     # standalone report — ALWAYS show/attach it with the verdict
```

Cluster designs: `cluster_trial(n_clusters, n_per_cluster, icc,
icc_pessimistic)` with `planned_analysis("lmm_random_intercept",
y_observed ~ treatment + (1 | cluster), inference = "satterthwaite" |
"kenward_roger" | "wald_z")` or `planned_analysis("cluster_mean_ttest",
y_observed ~ treatment)`. Simulated data columns available to formulas:
`y_observed`, `treatment`, `baseline`, `cluster`.

## Rules an agent must not break

1. **Fix everything before simulating.** Thresholds, scenario values, and
   the SESOI are chosen *before* looking at results. Never adjust
   `recovery_thresholds()` after seeing a verdict to make it pass; every
   deviation from a shipped profile is echoed in the report. If the user
   wants different thresholds, set them first and rerun.
2. **Never shrink the effect to build the pessimistic scenario.** The
   pessimistic rows perturb nuisance assumptions only (attrition x1.5,
   reliability -0.10, ICC to its declared upper bound, noncompliance
   +50%). Effect-size sensitivity is `effect_fragility()`, reported
   separately, outside the verdict. Target rows always run at the SESOI.
3. **The verdict never travels alone.** Always give the user
   `report(result)` (or its key numbers) with the PASS/RISK/FAIL label.
   The verdict is a decision convention, not a validity classification.
4. **A PASS is evidence about the instrument, not about the world.** Do
   not tell a user a PASS means the study will find an effect or that the
   assumptions are true. A FAIL means the design cannot recover the
   target *even under the user's own declared assumptions* — that is the
   strong claim.
5. **Ask for real values; label defaults.** Prefer the user's empirical
   evidence (prior attrition rates, published ICC ranges, validated
   reliabilities) via the `evidence` arguments; package defaults are the
   lowest evidence tier and the report says so. Declare
   `icc_pessimistic` explicitly rather than accepting the icc + 0.10
   fallback.
6. **Pre-specify degenerate-fit handling for mixed models.**
   `planned_analysis(..., degenerate_counts = FALSE)` (default) records
   singular/boundary fits and their coverage effect without counting
   them as failures; set `TRUE` only as a deliberate pre-specification.
   Fatal errors and nonconvergence always count.
7. **Set a seed** (any integer) or the run is not reproducible, and the
   report will say so.

## Reading the output

- `rejection_rate` is power on target rows; on null rows it is *test
  size* if the null world leaves the analyzable contrast at zero, and
  *target-null rejection rate* (false claims about the target) when
  selection displaces it — the report labels this automatically.
- `target_bias = estimator_bias + estimand_drift`, all in SESOI units.
  Drift is what the data strategy does to the question (a design
  problem, repaired only by changing design/analysis); estimator bias is
  what the estimator does to the answer. Resources (bigger n, better
  measures) repair precision, **not** drift.
- Type S / Type M are conditional on significance; below 200
  contributing simulations they are flagged UNSTABLE and cap the verdict
  at RISK.
- RISK also means: passes declared rows but fails pessimistic rows, or
  any margin within 2 Monte Carlo SEs of a threshold — the fix for the
  latter is more `sims`, not a different threshold.
- Verdicts are recomputed under shipped strict and lenient profiles; a
  flip across profiles is a finding to surface, not noise to hide.

## Typical agent tasks

- **Feasibility review**: declare the user's design faithfully (do not
  idealize: include reliability, attrition, clustering — silence is
  recorded as an omission in the report), run the confirmatory grid,
  present verdict + report + binding failure mode.
- **Design repair**: change ONE thing per candidate (reliability, n,
  clusters, estimator such as `"mi_baseline_adjusted"` for MAR-given-
  baseline dropout), rerun the full test per candidate, compare which
  criterion each repair fixes. Match the repair to the failure: MI-type
  estimators fix drift-induced bias; n fixes power; nothing fixes both
  automatically.
- **Fragility mapping**: `nuisance_fragility(design, "attrition_rate" |
  "reliability" | "icc", values = ...)` for the binding parameter;
  `effect_fragility(design)` for the effect curve.
- **Estimation-focused designs**: `recovery_test(...,
  scenarios = "target_grid", thresholds =
  recovery_thresholds("estimation"))`, optionally with
  `target_estimand(..., max_width = ...)`.

## Cost expectations

`lm`-based designs: seconds for 2000 sims x 4 rows. `lmer`-based:
minutes (Kenward-Roger slowest). 2000 sims/row is the working minimum;
raise it when any margin is within 2 MCSE. Runs are parallelizable at
the whole-`recovery_test()` level (one seed per call).

## Scope (do not overclaim)

Prototype: two-arm parallel trials (baseline, additive classical
measurement error, MCAR / baseline-dependent MAR attrition, one-sided
noncompliance) and cluster-randomized parallel trials (random-intercept
LMM with z/Satterthwaite/KR inference, cluster-mean t-test). Not
implemented: longitudinal/crossed random effects, Bayesian answer
strategies, prediction models, latent-variable measurement models. The
package does not judge whether the declared model is plausible — that
remains the researcher's job, made auditable.
