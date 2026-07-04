#!/usr/bin/env Rscript
# Reproduces the worked examples of the methods paper (sections 3.1 and
# 3.2) on the crossed confirmatory scenario grid and writes every
# numerical field of the manuscript's results tables, with Monte Carlo
# standard errors, to a markdown file.
#
# Every job is run under the algorithmic doubling stopping rule
# (recovery_test_stable, protocol section 2.4): S starts at `start_sims`
# and doubles up to `max_sims` while any required margin is within 2 MCSE
# of its threshold; a stable declared failure locks FAIL immediately.
# Both worked examples are then re-run across several seeds and a verdict-
# robustness table records whether any verdict flips (section 2.2).
#
# Usage:
#   Rscript run-examples.R [output.md] [start_sims] [cores] [max_sims]
# Defaults: example-results.md in the working directory, start 2000 sims,
# 6 cores, max 16000 sims. Every job carries its own seed, so results are
# reproducible regardless of scheduling.

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[1] else "example-results.md"
SIMS <- if (length(args) >= 2) as.integer(args[2]) else 2000L
CORES <- if (length(args) >= 3) as.integer(args[3]) else 6L
MAX_SIMS <- if (length(args) >= 4) as.integer(args[4]) else 16000L
BASE_SEED <- 20260703L
# Base seed plus additional seeds for the multi-seed robustness sweep.
SEEDS <- c(BASE_SEED, 20260704L, 4711L, 12345L)

suppressPackageStartupMessages({
  library(recoverlite)
  library(parallel)
})

# Run one job under the stopping rule; returns the result, its verdict,
# and the stopping record. Each job's own seed makes it order-independent.
run_stable <- function(design, seed) {
  r <- recovery_test_stable(design, start_sims = SIMS, max_sims = MAX_SIMS,
                            scenarios = "confirmatory_grid", seed = seed)
  list(res = r, v = verdict(r), stop = attr(r, "stopping"))
}

t_start <- Sys.time()
L <- character(0)
add <- function(...) L <<- c(L, paste0(...))

g <- function(run, name, col = "value") {
  run$diagnosands[run$diagnosands$diagnosand == name, col]
}
vm <- function(run, name, digits = 3) {
  v <- g(run, name); m <- g(run, name, "mcse")
  if (is.na(v)) return("n/a")
  sprintf("%.*f (%.4f)", digits, v, m)
}
vm_n <- function(run, name) {
  v <- g(run, name)
  if (is.na(v)) return("n/a")
  sprintf("%.3f (%.4f, n = %d)", v, g(run, name, "mcse"),
          g(run, name, "n_contributing"))
}
# Type S formatter: for a zero sign-flip count show the one-sided 95%
# Wilson upper bound (which the verdict checks against the threshold).
vm_ts <- function(run) {
  v <- g(run, "type_s")
  if (is.na(v)) return("n/a")
  n <- g(run, "type_s", "n_contributing")
  up <- g(run, "type_s", "upper")
  if (!is.null(up) && length(up) && is.finite(up) && v == 0) {
    sprintf("0.000 (n = %d; 95%% Wilson upper %.4f)", n, up)
  } else {
    sprintf("%.3f (%.4f, n = %d)", v, g(run, "type_s", "mcse"), n)
  }
}
# One-line description of where the stopping rule landed.
stop_line <- function(res) {
  st <- attr(res, "stopping")
  if (is.null(st)) return("")
  if (st$final_sims == st$start_sims && st$resolved) {
    sprintf("Stopping rule: resolved at S = %d.", st$final_sims)
  } else if (st$hit_ceiling) {
    sprintf("Stopping rule: doubled %d -> %d (ceiling); %d required margin(s) still within 2 MCSE, so RISK by the rule.",
            st$start_sims, st$final_sims, nrow(st$unresolved))
  } else {
    sprintf("Stopping rule: doubled %d -> %d (verdict determined).",
            st$start_sims, st$final_sims)
  }
}
fndw <- function(run) {
  fc <- attr(run$diagnosands, "failure_classes")
  sprintf("%d/%d/%d/%d", fc["fatal"], fc["nonconverged"], fc["degenerate"],
          fc["warned"])
}
verdict_line <- function(v) {
  agree <- length(unique(c(v$verdict_strict, v$verdict, v$verdict_lenient))) == 1
  sprintf("**%s** (strict %s / default %s / lenient %s — %s). Smallest signed margin: %s.",
          v$verdict, v$verdict_strict, v$verdict, v$verdict_lenient,
          if (agree) "profiles agree" else "profiles DISAGREE (a finding in itself)",
          v$smallest_margin)
}

## ------------------------------------------------------------------
## Job definitions (run in parallel; each job carries its own seed)
## ------------------------------------------------------------------
target_31 <- target_estimand(
  estimand = paste("ITT mean difference in the latent outcome at",
                   "post-test among all randomized participants"),
  scale = "latent-outcome standardized mean difference",
  sesoi = 0.40
)
strategy_31 <- two_arm_trial(n_per_arm = 115, allocation = 0.5,
                             baseline_outcome_cor = 0.5)
attrition_31 <- attrition_model(
  rate = 0.15, mechanism = "differential",
  baseline_slope_treated = -0.5, baseline_slope_control = 0,
  max_rate = 0.6,
  evidence = "illustrative declaration for the worked example")
cc_analysis <- planned_analysis("linear_model", y_observed ~ treatment)

design_31 <- declare_recovery(
  target = target_31, data_strategy = strategy_31,
  measurement = measured_outcome(0.70), missingness = attrition_31,
  answer_strategy = cc_analysis)

design_31_mcar <- declare_recovery(
  target = target_31, data_strategy = strategy_31,
  measurement = measured_outcome(0.70),
  missingness = attrition_model(rate = 0.15, mechanism = "mcar"),
  answer_strategy = cc_analysis)

rep_rel <- declare_recovery(
  target = target_31, data_strategy = strategy_31,
  measurement = measured_outcome(0.90), missingness = attrition_31,
  answer_strategy = cc_analysis)

rep_n <- declare_recovery(
  target = target_31,
  data_strategy = two_arm_trial(n_per_arm = 230, allocation = 0.5,
                                baseline_outcome_cor = 0.5),
  measurement = measured_outcome(0.70), missingness = attrition_31,
  answer_strategy = cc_analysis)

rep_mi <- declare_recovery(
  target = target_31, data_strategy = strategy_31,
  measurement = measured_outcome(0.70), missingness = attrition_31,
  answer_strategy = planned_analysis(
    "mi_baseline_adjusted", y_observed ~ treatment + baseline,
    m_imputations = 20))

# Complete-case ANCOVA repair (reviewer maj 4): under MAR given the
# observed baseline B, a baseline-adjusted complete-case regression
# lm(y_observed ~ treatment + B) targets the ITT treatment coefficient
# without multiple imputation, provided there is no treatment-by-baseline
# interaction. Compare against the MI + ANCOVA arm (rep_mi).
rep_ancova <- declare_recovery(
  target = target_31, data_strategy = strategy_31,
  measurement = measured_outcome(0.70), missingness = attrition_31,
  answer_strategy = planned_analysis(
    "linear_model", y_observed ~ treatment + baseline))

target_32 <- target_estimand(
  estimand = paste("ITT mean difference in pupil outcome between",
                   "intervention and control schools"),
  scale = "student-level standardized mean difference",
  sesoi = 0.40
)
schools_16 <- cluster_trial(n_clusters = 16, n_per_cluster = 30,
                            icc = 0.05, icc_pessimistic = 0.15,
                            evidence = "plausible high-ICC value for school outcomes (illustrative)")
lmm_formula <- y_observed ~ treatment + (1 | cluster)
cluster_design <- function(strategy, analysis) {
  suppressWarnings(declare_recovery(
    target = target_32, data_strategy = strategy,
    answer_strategy = analysis))
}
analyses_32 <- list(
  `Naive pupil-level OLS (negative control)` =
    planned_analysis("linear_model", y_observed ~ treatment),
  `LMM Wald z` =
    planned_analysis("lmm_random_intercept", lmm_formula,
                     inference = "wald_z"),
  `LMM Satterthwaite` =
    planned_analysis("lmm_random_intercept", lmm_formula,
                     inference = "satterthwaite"),
  `LMM Kenward-Roger` =
    planned_analysis("lmm_random_intercept", lmm_formula,
                     inference = "kenward_roger"),
  `Cluster-level t-test (14 df)` =
    planned_analysis("cluster_mean_ttest", y_observed ~ treatment)
)

# Canonical job order. rep_ancova is APPENDED last so that every other
# job keeps the per-job seed offset it had in the v0.1.0 archive (job
# seed = seed_base + offset), which keeps the headline numbers comparable
# across releases. Job offset = position in this list.
jobs <- c(
  list(main_31 = design_31, mcar_31 = design_31_mcar, rep_rel = rep_rel,
       rep_n = rep_n, rep_mi = rep_mi),
  setNames(lapply(analyses_32, function(a) cluster_design(schools_16, a)),
           paste0("cl_", seq_along(analyses_32))),
  list(rep_pupils = cluster_design(
         cluster_trial(16, 60, icc = 0.05, icc_pessimistic = 0.15),
         analyses_32$`LMM Satterthwaite`),
       rep_schools = cluster_design(
         cluster_trial(32, 30, icc = 0.05, icc_pessimistic = 0.15),
         analyses_32$`LMM Satterthwaite`),
       rep_ancova = rep_ancova)
)

# ------- Primary run (base seed) under the stopping rule ---------------
message("Running ", length(jobs), " recovery tests under the doubling ",
        "stopping rule (start ", SIMS, ", max ", MAX_SIMS,
        " sims x 4 rows) on ", CORES, " cores ...")
results <- mclapply(seq_along(jobs), function(i) {
  run_stable(jobs[[i]], BASE_SEED + i)
}, mc.cores = CORES)
names(results) <- names(jobs)
stopifnot(!any(vapply(results, inherits, TRUE, "try-error")))

# ------- Multi-seed robustness sweep (all jobs, all seeds) -------------
# Each (job, seed-family) task carries seed = seed_base + job offset and
# runs the full stopping rule; we record the verdict per seed so that any
# flip across seeds is visible (section 2.2). Run as one flat mclapply so
# the core pool stays saturated.
alt_seeds <- SEEDS[-1]
sweep_tasks <- expand.grid(job = seq_along(jobs), s = seq_along(alt_seeds),
                           KEEP.OUT.ATTRS = FALSE)
message("Multi-seed robustness sweep: ", nrow(sweep_tasks),
        " additional (job x seed) runs ...")
sweep <- mclapply(seq_len(nrow(sweep_tasks)), function(k) {
  i <- sweep_tasks$job[k]; j <- sweep_tasks$s[k]
  out <- run_stable(jobs[[i]], alt_seeds[j] + i)
  list(job = names(jobs)[i], seed = alt_seeds[j],
       verdict = out$v$verdict, verdict_strict = out$v$verdict_strict,
       verdict_lenient = out$v$verdict_lenient,
       final_sims = out$stop$final_sims)
}, mc.cores = CORES)
stopifnot(!any(vapply(sweep, inherits, TRUE, "try-error")))

## ------------------------------------------------------------------
## Markdown
## ------------------------------------------------------------------
add("# Worked-example results for the manuscript, section 3")
add("")
add("> Generated by `recoverlite/inst/paper/run-examples.R` — do not edit the")
add("> numbers by hand. ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
    "; base seed ", BASE_SEED, ".")
add("> Every job runs under the algorithmic doubling stopping rule ",
    "(`recovery_test_stable`, section 2.4): S starts at ", SIMS,
    " and doubles up to ", MAX_SIMS, " while any required margin is within")
add("> 2 MCSE of its threshold; a stable declared failure locks FAIL. The ",
    "final S reached is reported per job. Confirmatory grid rows: ",
    "Null-declared, Null-pessimistic, Target-declared, Target-pessimistic.")
add("> Type S zero-count rows carry a one-sided 95% Wilson upper bound, ",
    "which the verdict checks against the threshold (technical comment C).")
add("> Multi-seed robustness sweep over seeds {",
    paste(SEEDS, collapse = ", "), "} (section 2.2; table at the end).")
add("> ", R.version.string, "; recoverlite ",
    as.character(packageVersion("recoverlite")),
    "; threshold profiles ", recovery_thresholds()$version,
    " (default profile, shipped values; strict/lenient recomputed).")
add("")

scenario_labels <- c(null_declared = "Null-declared",
                     null_pessimistic = "Null-pessimistic",
                     target_declared = "Target-declared",
                     target_pessimistic = "Target-pessimistic")

## ----- Example 3.1 --------------------------------------------------
r31 <- results$main_31$res; v31 <- results$main_31$v
nd <- r31$runs$null_declared

add("## Example 3.1 — powered trial with attrition and measurement error")
add("")
add("Design as declared in the manuscript: two-arm RCT, 230 recruited")
add("(115 per arm), SESOI d = 0.40 on the latent scale, baseline-outcome")
add("correlation rho = 0.5, additive classical measurement error with")
add("reliability 0.70 (Var(e) = 1/0.70 - 1; raw contrast unbiased in")
add("expectation, residual variance inflated), attrition")
add("logit Pr(dropout | Z, B) = alpha_Z + gamma_Z B with gamma_0 = 0,")
add("gamma_1 = -0.5 and both marginal arm rates calibrated to 0.15")
add("(overall 0.15). Planned analysis: complete-case")
add("`lm(y_observed ~ treatment)`, Wald 95% interval, alpha = .05.")
add("")

a0 <- recoverlite:::calibrate_dropout_intercept(0.15, 0)
a1 <- recoverlite:::calibrate_dropout_intercept(0.15, -0.5)
add("Parameter-table fill-ins (declared column): alpha_0 = ",
    sprintf("%.4f", a0), ", alpha_1 = ", sprintf("%.4f", a1),
    ", gamma_0 = 0, gamma_1 = -0.5; marginal attrition control 0.15,")
add("intervention 0.15, overall 0.15. Pessimistic column: rates x1.5 ->")
add("0.225/0.225/0.225 (alpha_Z recalibrated: alpha_0 = ",
    sprintf("%.4f", recoverlite:::calibrate_dropout_intercept(0.225, 0)),
    ", alpha_1 = ",
    sprintf("%.4f", recoverlite:::calibrate_dropout_intercept(0.225, -0.5)),
    "), slopes unchanged, reliability 0.70 -> 0.60.")
add("")
drift_nd <- g(nd, "estimand_drift") * 0.40
add("Null-world fill-in: E[theta_obs] under Null-declared = **",
    sprintf("%.4f", drift_nd), "** latent-SD units (drift ",
    vm(nd, "estimand_drift", 4), " in Delta units) — positive, so null-row")
add("rejections are false claims about the target, partly induced by")
add("selection, and are reported as the target-null rejection rate.")
mcar_nd <- results$mcar_31$res$runs$null_declared
add("Pure test size under the MCAR reference null (footnote †): ",
    vm(mcar_nd, "rejection_rate", 4), ".")
add("")
add("### Recovery results (fills the section 3.1 table)")
add("")
add("| Scenario | Rejection rate / Power | Target bias | Estimator bias | Drift | Coverage (theta) | Analyzable coverage | Type S | Type M (n) | Precision | Failures F/N/D/W |")
add("|---|---|---|---|---|---|---|---|---|---|---|")
for (nm in names(scenario_labels)) {
  run <- r31$runs[[nm]]
  add(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
              scenario_labels[nm], vm(run, "rejection_rate"),
              vm(run, "target_bias"), vm(run, "estimator_bias"),
              vm(run, "estimand_drift"), vm(run, "coverage"),
              vm(run, "analyzable_coverage"), vm_ts(run),
              vm_n(run, "type_m"), vm(run, "precision", 3), fndw(run)))
}
add("")
add("Mean analyzed n per simulation (declared rows): ",
    sprintf("%.1f", attr(r31$runs$target_declared$diagnosands,
                         "mean_n_analyzed")),
    " of 230 recruited; pessimistic rows: ",
    sprintf("%.1f", attr(r31$runs$target_pessimistic$diagnosands,
                         "mean_n_analyzed")), ".")
add("")
add("### Verdict")
add("")
add(verdict_line(v31))
add(stop_line(r31))
if (!is.na(v31$binding)) add("", "Binding failure mode: ", v31$binding)
add("")
add("### Design repairs (each simulated as a full declaration)")
add("")
repair_31 <- function(key, label) {
  r <- results[[key]]$res; v <- results[[key]]$v
  td <- r$runs$target_declared; tp <- r$runs$target_pessimistic
  ndr <- r$runs$null_declared
  add("**", label, "** — Target-declared: power ", vm(td, "rejection_rate"),
      ", target bias ", vm(td, "target_bias"),
      ", drift ", vm(td, "estimand_drift"),
      ", coverage ", vm(td, "coverage"),
      ", Type M ", vm_n(td, "type_m"),
      ", mean analyzed n ", sprintf("%.0f", attr(td$diagnosands, "mean_n_analyzed")),
      "; Null-declared rejection ", vm(ndr, "rejection_rate"),
      "; Target-pessimistic: power ", vm(tp, "rejection_rate"),
      ", target bias ", vm(tp, "target_bias"),
      ", drift ", vm(tp, "estimand_drift"),
      ", coverage ", vm(tp, "coverage"),
      ". Verdict: ", verdict_line(v), " ", stop_line(r))
  add("")
}
repair_31("rep_ancova",
          "Complete-case ANCOVA (lm(y_observed ~ treatment + baseline))")
repair_31("rep_mi",
          "MI + ANCOVA analysis model (m = 20, y ~ treatment + baseline, rho = 0.5)")
repair_31("rep_rel", "Reliability 0.70 -> 0.90")
repair_31("rep_n", "Recruitment 230 -> 460 (230 per arm)")
add("Note for the repair paragraph: the two estimand-matched repairs ")
add("(complete-case ANCOVA and MI + ANCOVA) are evaluated through bias ")
add("against theta, not drift: under dropout that is MAR given the ")
add("observed baseline B, a baseline-adjusted treatment coefficient ")
add("targets the full-data ITT estimand (no treatment-by-baseline ")
add("interaction is declared). Compare the ANCOVA and MI drift/bias/power ")
add("columns to see whether ANCOVA removes the estimand drift as cheaply ")
add("as MI. The reliability and larger-recruitment repairs keep the ")
add("complete-case `lm(y ~ treatment)` estimator, so their drift is ")
add("essentially unchanged — a larger or cleaner complete-case study ")
add("estimates the displaced contrast more precisely, not the ITT target.")
add("")

## ----- Example 3.2 --------------------------------------------------
add("## Example 3.2 — the cluster trap: few clusters and fragile inference")
add("")
add("Design: 16 schools (8/8) x 30 pupils (sizes fixed), d = 0.40")
add("student-SD units (SESOI), ICC 0.05 declared / 0.15 pessimistic. No")
add("attrition or exclusions: theta_obs = theta in every scenario, drift")
add("identically zero, and null-row rejection rates are pure test size")
add("(drift columns omitted, per the manuscript). Convergence handling")
add("pre-specified: fatal errors and nonconvergence count as model")
add("failures; singular fits are recorded as degenerate, retained, and")
add("their marginal effect on coverage is reported in the failure-class")
add("notes. Five candidate answer strategies, each fully specified.")
add("")
for (i in seq_along(analyses_32)) {
  key <- paste0("cl_", i)
  r <- results[[key]]$res; v <- results[[key]]$v
  add("### Analysis ", i, ": ", names(analyses_32)[i])
  add("")
  add("| Scenario | Rejection rate / Power | Coverage | Target bias | Type S | Type M (n) | Precision | Failures F/N/D/W |")
  add("|---|---|---|---|---|---|---|---|")
  for (nm in names(scenario_labels)) {
    run <- r$runs[[nm]]
    add(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
                scenario_labels[nm], vm(run, "rejection_rate"),
                vm(run, "coverage"), vm(run, "target_bias"),
                vm_ts(run), vm_n(run, "type_m"),
                vm(run, "precision", 3), fndw(run)))
  }
  td <- r$runs$target_declared
  mf_note <- g(td, "model_failure", "note")
  if (grepl("coverage among degenerate", mf_note)) {
    add("")
    add("Degenerate-fit detail (Target-declared): ",
        sub(".*(coverage among degenerate[^)]*).*", "\\1", mf_note), ".")
  }
  add("")
  add("Verdict: ", verdict_line(v), " ", stop_line(r))
  if (!is.na(v$binding)) add("", "Binding failure mode: ", v$binding)
  add("")
}
add("### Repairs (pupils vs schools, LMM Satterthwaite)")
add("")
for (key in c("rep_pupils", "rep_schools")) {
  r <- results[[key]]$res; v <- results[[key]]$v
  label <- if (key == "rep_pupils")
    "Pupils per school 30 -> 60 (16 schools)" else
    "Schools 16 -> 32 (30 pupils each)"
  td <- r$runs$target_declared; tp <- r$runs$target_pessimistic
  ndr <- r$runs$null_declared
  add("**", label, "** — Null-declared test size ",
      vm(ndr, "rejection_rate"),
      "; Target-declared: power ", vm(td, "rejection_rate"),
      ", coverage ", vm(td, "coverage"),
      "; Target-pessimistic: power ", vm(tp, "rejection_rate"),
      ", coverage ", vm(tp, "coverage"),
      ". Verdict: ", verdict_line(v), " ", stop_line(r))
  add("")
}

## ----- Multi-seed robustness (section 2.2) --------------------------
add("## Multi-seed robustness")
add("")
add("Every job re-run under the stopping rule at each seed; the cell is ",
    "the default-profile verdict (final S). A job whose verdict is not ",
    "constant across seeds is **seed-sensitive**, which is itself ",
    "reportable (section 2.2).")
add("")
# Primary (base-seed) verdicts from `results`, plus the sweep.
primary <- lapply(names(jobs), function(nm) {
  list(job = nm, seed = BASE_SEED, verdict = results[[nm]]$v$verdict,
       final_sims = attr(results[[nm]]$res, "stopping")$final_sims)
})
all_runs <- c(primary, sweep)
cell <- function(nm, sd) {
  hit <- Filter(function(x) x$job == nm && x$seed == sd, all_runs)[[1]]
  sprintf("%s (%d)", hit$verdict, hit$final_sims)
}
job_label <- c(main_31 = "3.1 main (CC lm)",
               rep_ancova = "3.1 CC ANCOVA", rep_mi = "3.1 MI+ANCOVA",
               rep_rel = "3.1 reliability 0.90", rep_n = "3.1 n=460",
               cl_1 = "3.2 naive OLS", cl_2 = "3.2 LMM Wald z",
               cl_3 = "3.2 LMM Satterthwaite", cl_4 = "3.2 LMM Kenward-Roger",
               cl_5 = "3.2 cluster-t", rep_pupils = "3.2 pupils 30->60",
               rep_schools = "3.2 schools 16->32", mcar_31 = "3.1 MCAR ref")
add(paste0("| Job | ", paste(sprintf("seed %d", SEEDS), collapse = " | "),
           " | Stable? |"))
add(paste0("|---|", paste(rep("---", length(SEEDS) + 1), collapse = "|"), "|"))
flips <- character(0)
for (nm in names(jobs)) {
  cells <- vapply(SEEDS, function(sd) cell(nm, sd), character(1))
  verds <- vapply(SEEDS, function(sd) {
    hit <- Filter(function(x) x$job == nm && x$seed == sd, all_runs)[[1]]
    hit$verdict
  }, character(1))
  stable <- length(unique(verds)) == 1
  if (!stable) flips <- c(flips, nm)
  add(sprintf("| %s | %s | %s |", job_label[nm] %||% nm,
              paste(cells, collapse = " | "),
              if (stable) "yes" else "**NO (flips)**"))
}
add("")
if (length(flips)) {
  add("Seed-sensitive verdicts: **",
      paste(vapply(flips, function(nm) job_label[nm] %||% nm, character(1)),
            collapse = ", "),
      "**. These are reportable per section 2.2: the design sits close ",
      "enough to a bright line that Monte Carlo seed choice moves the ",
      "verdict, so the honest report is the flip itself, not either side.")
} else {
  add("No verdict is seed-sensitive: every job returns the same ",
      "default-profile verdict across all seeds.")
}
add("")

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
final_sims_vec <- vapply(names(jobs),
  function(nm) attr(results[[nm]]$res, "stopping")$final_sims, integer(1))
add("---")
add("")
add(sprintf(paste0("Total computation: %.1f minutes wall clock on %d cores ",
                   "(%d jobs x %d seeds under the doubling stopping rule, ",
                   "start %d / max %d sims x 4 scenario rows; primary-run ",
                   "final S ranged %d-%d)."),
            elapsed, CORES, length(jobs), length(SEEDS), SIMS, MAX_SIMS,
            min(final_sims_vec), max(final_sims_vec)))
add("")
add("Full standalone recovery reports: `report(recovery_test_stable(<design>,",
    " seed = <job seed>))`; job seed = seed_base + job offset (offset = ",
    "position in the script's `jobs` list), seed_base in {",
    paste(SEEDS, collapse = ", "), "}.")

writeLines(L, out_path)
message("Wrote ", out_path, " (", sprintf("%.1f", elapsed), " min)")
