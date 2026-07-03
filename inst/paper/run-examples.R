#!/usr/bin/env Rscript
# Reproduces the worked examples of the methods paper (sections 3.1 and
# 3.2) on the crossed confirmatory scenario grid and writes every
# numerical field of the manuscript's results tables, with Monte Carlo
# standard errors, to a markdown file.
#
# Usage:
#   Rscript run-examples.R [output.md] [sims] [cores]
# Defaults: example-results.md in the working directory, 2000 sims,
# 6 cores. Every recovery_test() call carries its own seed, so results
# are reproducible regardless of scheduling.

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[1] else "example-results.md"
SIMS <- if (length(args) >= 2) as.integer(args[2]) else 2000L
CORES <- if (length(args) >= 3) as.integer(args[3]) else 6L
BASE_SEED <- 20260703L

suppressPackageStartupMessages({
  library(recoverlite)
  library(parallel)
})

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
         analyses_32$`LMM Satterthwaite`))
)

message("Running ", length(jobs), " recovery tests (", SIMS,
        " sims x 4 scenario rows each) on ", CORES, " cores ...")
results <- mclapply(seq_along(jobs), function(i) {
  r <- recovery_test(jobs[[i]], sims = SIMS,
                     scenarios = "confirmatory_grid",
                     seed = BASE_SEED + i)
  list(res = r, v = verdict(r))
}, mc.cores = CORES)
names(results) <- names(jobs)
stopifnot(!any(vapply(results, inherits, TRUE, "try-error")))

## ------------------------------------------------------------------
## Markdown
## ------------------------------------------------------------------
add("# Worked-example results for the manuscript, section 3")
add("")
add("> Generated by `recoverlite/inst/paper/run-examples.R` — do not edit the")
add("> numbers by hand. ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"),
    "; base seed ", BASE_SEED, "; ", SIMS,
    " simulations per scenario row (confirmatory grid: Null-declared,")
add("> Null-pessimistic, Target-declared, Target-pessimistic).")
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
              vm(run, "analyzable_coverage"), vm_n(run, "type_s"),
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
      ", Type M ", vm_n(td, "type_m"),
      ", mean analyzed n ", sprintf("%.0f", attr(td$diagnosands, "mean_n_analyzed")),
      "; Null-declared rejection ", vm(ndr, "rejection_rate"),
      "; Target-pessimistic: power ", vm(tp, "rejection_rate"),
      ", target bias ", vm(tp, "target_bias"),
      ". Verdict: ", verdict_line(v))
  add("")
}
repair_31("rep_rel", "Reliability 0.70 -> 0.90")
repair_31("rep_n", "Recruitment 230 -> 460 (230 per arm)")
repair_31("rep_mi",
          "MI baseline-adjusted estimator (m = 20, y ~ treatment + baseline, rho = 0.5)")
add("Note for the repair paragraph: the MI repair is evaluated through")
add("bias against theta, not drift (its target is the full-data ITT")
add("estimand; identifiable because dropout is MAR given the observed")
add("baseline). The larger-recruitment repair shows drift essentially")
add("unchanged — a larger complete-case study estimates the displaced")
add("contrast more precisely.")
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
                vm_n(run, "type_s"), vm_n(run, "type_m"),
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
  add("Verdict: ", verdict_line(v))
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
      ". Verdict: ", verdict_line(v))
  add("")
}

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
add("---")
add("")
add(sprintf(paste0("Total computation: %.1f minutes wall clock on %d cores ",
                   "(%d recovery tests x 4 scenario rows x %d simulations)."),
            elapsed, CORES, length(jobs), SIMS))
add("")
add("Full standalone recovery reports: `report(recovery_test(<design>,",
    " sims = ", SIMS, ", seed = <job seed>))`; job seeds are ",
    BASE_SEED, " + job index in the order defined in this script.")

writeLines(L, out_path)
message("Wrote ", out_path, " (", sprintf("%.1f", elapsed), " min)")
