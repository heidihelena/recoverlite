#' Declare the target estimand
#'
#' Step 1 of the recovery-test protocol: state, in one sentence, the
#' quantity the study is designed to estimate, the scale on which it is
#' expressed, and the smallest effect size of interest (SESOI). The SESOI
#' should be justified on substantive grounds, not reverse-engineered from
#' the available sample size.
#'
#' @param estimand Character. One-sentence statement of the target
#'   estimand (quantity, population, scale).
#' @param scale Character. The scale on which the estimand is expressed,
#'   e.g. `"latent-outcome standardized mean difference"`.
#' @param sesoi Numeric. The smallest effect size of interest, on the
#'   estimand scale. Must be strictly positive. Target scenario rows are
#'   simulated at the SESOI.
#' @param bias_scale_unit Numeric or `NULL`. The declared substantive unit
#'   Delta used to scale target bias, estimator bias, and estimand drift:
#'   each is (expectation - target) / `bias_scale_unit`. Defaults to the
#'   SESOI.
#' @param max_width Numeric or `NULL`. Maximum acceptable confidence
#'   interval width, if one is declared. Used as the precision criterion
#'   by the estimation threshold profile; reported (not thresholded)
#'   otherwise.
#'
#' @return An object of class `recovery_target`.
#' @seealso [declare_recovery()]
#' @export
#' @examples
#' target_estimand(
#'   estimand = "ITT mean difference at 12 weeks",
#'   scale = "latent-outcome standardized mean difference",
#'   sesoi = 0.40
#' )
target_estimand <- function(estimand, scale, sesoi, bias_scale_unit = NULL,
                            max_width = NULL) {
  stopifnot(
    "`estimand` must be a single character string" =
      is.character(estimand) && length(estimand) == 1L,
    "`scale` must be a single character string" =
      is.character(scale) && length(scale) == 1L,
    "`sesoi` must be a single positive number" =
      is.numeric(sesoi) && length(sesoi) == 1L && is.finite(sesoi) && sesoi > 0
  )
  bias_scale_unit <- bias_scale_unit %||% sesoi
  stopifnot(
    "`bias_scale_unit` must be a single positive number" =
      is.numeric(bias_scale_unit) && length(bias_scale_unit) == 1L &&
        bias_scale_unit > 0
  )
  if (!is.null(max_width)) {
    stopifnot("`max_width` must be a single positive number" =
                is.numeric(max_width) && length(max_width) == 1L &&
                  max_width > 0)
  }
  structure(
    list(estimand = estimand, scale = scale, sesoi = sesoi,
         bias_scale_unit = bias_scale_unit, max_width = max_width),
    class = "recovery_target"
  )
}

#' Declare a two-arm randomized trial data strategy
#'
#' Individual random assignment to two arms, with an observed baseline
#' measure. The data-generating model follows the manuscript (section
#' 3.1): for participant i with assignment Z,
#' \deqn{B_i ~ N(0,1)}
#' \deqn{Y*_i = tau Z_i + rho B_i + sqrt(1 - rho^2) eps_i}
#' so the true (latent) outcome has unit SD given Z and the declared
#' effect is expressed in latent standard-deviation units. The baseline is
#' observed at randomization for every participant and is available to the
#' answer strategy as the column `baseline`.
#'
#' `n_per_arm` is the number of participants *recruited* per arm (under
#' equal allocation); attrition, if declared via [attrition_model()], is
#' applied after recruitment.
#'
#' @param n_per_arm Integer. Participants recruited per arm.
#' @param allocation Numeric in (0, 1). Share assigned to treatment.
#' @param baseline_outcome_cor Numeric in \[0, 1). Correlation rho between
#'   the observed baseline and the true outcome (default 0.5).
#' @param noncompliance Numeric in \[0, 1). Share of treated participants
#'   who do not receive treatment (one-sided never-takers). The target
#'   estimand remains the declared ITT effect.
#'
#' @return An object of class `recovery_data_strategy`.
#' @export
two_arm_trial <- function(n_per_arm, allocation = 0.5,
                          baseline_outcome_cor = 0.5, noncompliance = 0) {
  stopifnot(
    "`n_per_arm` must be a single positive whole number" =
      is.numeric(n_per_arm) && length(n_per_arm) == 1L && n_per_arm >= 2 &&
        n_per_arm == round(n_per_arm),
    "`allocation` must be strictly between 0 and 1" =
      is.numeric(allocation) && length(allocation) == 1L &&
        allocation > 0 && allocation < 1,
    "`baseline_outcome_cor` must be in [0, 1)" =
      is.numeric(baseline_outcome_cor) && length(baseline_outcome_cor) == 1L &&
        baseline_outcome_cor >= 0 && baseline_outcome_cor < 1,
    "`noncompliance` must be in [0, 1)" =
      is.numeric(noncompliance) && length(noncompliance) == 1L &&
        noncompliance >= 0 && noncompliance < 1
  )
  structure(
    list(type = "two_arm_trial", n_per_arm = as.integer(n_per_arm),
         allocation = allocation,
         baseline_outcome_cor = baseline_outcome_cor,
         noncompliance = noncompliance),
    class = c("recovery_two_arm", "recovery_data_strategy")
  )
}

#' Declare a cluster-randomized trial data strategy
#'
#' Cluster-level random assignment with individual outcomes and cluster
#' sizes fixed by design. For pupil i in cluster j,
#' \deqn{Y_ij = tau Z_j + u_j + e_ij}
#' with total individual-level variance standardized to 1 and
#' ICC = var(u) / (var(u) + var(e)). The declared effect is expressed in
#' individual-level (total) standard-deviation units. The simulated data
#' contain a `cluster` id column for use in the analysis formula, e.g.
#' `y_observed ~ treatment + (1 | cluster)`.
#'
#' @param n_clusters Integer. Total number of clusters randomized.
#' @param n_per_cluster Integer. Individuals measured per cluster.
#' @param icc Numeric in \[0, 1). Declared intraclass correlation.
#' @param icc_pessimistic Numeric or `NULL`. The upper end of the
#'   empirically documented ICC range for the outcome domain, used in the
#'   pessimistic scenario rows. If `NULL`, the pessimistic rows default to
#'   `icc + 0.10` and the report records the value as a package default
#'   (the lowest tier of the evidence hierarchy). Declare it explicitly.
#' @param allocation Numeric in (0, 1). Share of clusters assigned to
#'   treatment.
#' @param evidence Character or `NULL`. Evidence tier / source for the
#'   declared and pessimistic ICC values (e.g. a published ICC range for
#'   the outcome domain). Echoed in the report.
#'
#' @return An object of class `recovery_data_strategy`.
#' @export
cluster_trial <- function(n_clusters, n_per_cluster, icc,
                          icc_pessimistic = NULL, allocation = 0.5,
                          evidence = NULL) {
  stopifnot(
    "`n_clusters` must be a whole number >= 4" =
      is.numeric(n_clusters) && length(n_clusters) == 1L && n_clusters >= 4 &&
        n_clusters == round(n_clusters),
    "`n_per_cluster` must be a positive whole number" =
      is.numeric(n_per_cluster) && length(n_per_cluster) == 1L &&
        n_per_cluster >= 1 && n_per_cluster == round(n_per_cluster),
    "`icc` must be in [0, 1)" =
      is.numeric(icc) && length(icc) == 1L && icc >= 0 && icc < 1,
    "`allocation` must be strictly between 0 and 1" =
      is.numeric(allocation) && length(allocation) == 1L &&
        allocation > 0 && allocation < 1
  )
  if (!is.null(icc_pessimistic)) {
    stopifnot(
      "`icc_pessimistic` must be in [0, 1) and >= `icc`" =
        is.numeric(icc_pessimistic) && length(icc_pessimistic) == 1L &&
          icc_pessimistic >= icc && icc_pessimistic < 1
    )
  }
  structure(
    list(type = "cluster_trial", n_clusters = as.integer(n_clusters),
         n_per_cluster = as.integer(n_per_cluster), icc = icc,
         icc_pessimistic = icc_pessimistic, allocation = allocation,
         evidence = evidence),
    class = c("recovery_cluster", "recovery_data_strategy")
  )
}

#' Declare the measurement model for the observed outcome
#'
#' Classical additive measurement error, per the manuscript (section 3.1):
#' the observed outcome is `y_observed = y_true + e` with
#' `Var(e) = 1/reliability - 1`, so that reliability =
#' Var(y_true)/Var(y_observed). Under this model the raw treatment
#' contrast is *not attenuated in expectation*: measurement error inflates
#' residual variance (costing precision, power, and — downstream of low
#' power — Type M exaggeration), and is charged to the variance account,
#' not the bias account. (Attenuation by sqrt(reliability) would apply
#' only if the estimator standardized by the observed-outcome SD.)
#'
#' If measurement is omitted from [declare_recovery()], the outcome is
#' treated as perfectly reliable and the report states this explicitly
#' (silence must not imply ideality).
#'
#' @param reliability Numeric in (0, 1\]. Declared reliability (e.g.
#'   test-retest) of the outcome measure.
#' @param evidence Character or `NULL`. Evidence tier / source for the
#'   declared reliability (e.g. a validation study). Echoed in the report.
#'
#' @return An object of class `recovery_measurement`.
#' @export
measured_outcome <- function(reliability, evidence = NULL) {
  stopifnot(
    "`reliability` must be a single number in (0, 1]" =
      is.numeric(reliability) && length(reliability) == 1L &&
        reliability > 0 && reliability <= 1
  )
  structure(list(reliability = reliability, evidence = evidence),
            class = "recovery_measurement")
}

#' Declare the attrition (missing-outcome) model
#'
#' Dropout follows a declared response model rather than an unstated one.
#' Under `mechanism = "differential"`, the probability of dropout is
#' \deqn{logit Pr(dropout | Z, B) = alpha_Z + gamma_Z B}
#' depending on assignment and the *observed baseline* — not on the
#' unobserved post-test outcome. Because dropout depends only on observed
#' quantities (Z, B), missingness is MAR given the baseline; this keeps
#' the drift mechanism transparent and makes baseline-based repairs
#' (e.g. [planned_analysis()] with `estimator = "mi_baseline_adjusted"`)
#' identifiable for the ITT estimand. The arm-specific intercepts alpha_Z
#' are calibrated numerically so that each arm's *marginal* attrition rate
#' equals its declared value: arm imbalance in retention is a declared
#' quantity, not a by-product.
#'
#' With a negative treated-arm slope and higher scores indicating better
#' functioning, intervention participants with poorer baseline prognosis
#' are more likely to drop out. This mechanism persists when the treatment
#' effect is zero, so under null scenarios the analyzable-data contrast is
#' displaced from 0 and rejections are false claims about the target.
#'
#' `mechanism = "mcar"` drops observations completely at random at the
#' declared rate(s).
#'
#' @param rate Numeric in \[0, 1). Anticipated overall marginal attrition
#'   rate; used for both arms unless arm-specific rates are given.
#' @param mechanism Character. `"differential"` (default) or `"mcar"`.
#' @param rate_control,rate_treated Numeric or `NULL`. Arm-specific
#'   marginal attrition rates; default to `rate`.
#' @param baseline_slope_treated Numeric. gamma_1, the log-odds change in
#'   treated-arm dropout per baseline SD (default -0.5: poorer baseline
#'   prognosis, higher dropout). Ignored under `"mcar"`.
#' @param baseline_slope_control Numeric. gamma_0 (default 0). Ignored
#'   under `"mcar"`.
#' @param max_rate Numeric. Field-relevant cap applied when the
#'   pessimistic scenario multiplies the declared rates by 1.5; must be
#'   stated and justified.
#' @param evidence Character or `NULL`. Evidence tier / source for the
#'   declared attrition values (e.g. observed rates in comparable trials).
#'   Echoed in the report.
#'
#' @return An object of class `recovery_attrition`.
#' @export
attrition_model <- function(rate, mechanism = c("differential", "mcar"),
                            rate_control = NULL, rate_treated = NULL,
                            baseline_slope_treated = -0.5,
                            baseline_slope_control = 0,
                            max_rate = 0.6, evidence = NULL) {
  mechanism <- match.arg(mechanism)
  stopifnot(
    "`rate` must be in [0, 1)" =
      is.numeric(rate) && length(rate) == 1L && rate >= 0 && rate < 1,
    "`baseline_slope_treated` must be a single finite number" =
      is.numeric(baseline_slope_treated) &&
        length(baseline_slope_treated) == 1L &&
        is.finite(baseline_slope_treated),
    "`baseline_slope_control` must be a single finite number" =
      is.numeric(baseline_slope_control) &&
        length(baseline_slope_control) == 1L &&
        is.finite(baseline_slope_control),
    "`max_rate` must be in (0, 1)" =
      is.numeric(max_rate) && length(max_rate) == 1L &&
        max_rate > 0 && max_rate < 1
  )
  rate_control <- rate_control %||% rate
  rate_treated <- rate_treated %||% rate
  for (r in c(rate_control, rate_treated)) {
    stopifnot("arm-specific rates must be in [0, 1)" =
                is.numeric(r) && r >= 0 && r < 1)
  }
  if (mechanism == "mcar") {
    baseline_slope_treated <- 0
    baseline_slope_control <- 0
  }
  structure(
    list(rate = rate, mechanism = mechanism,
         rate_control = rate_control, rate_treated = rate_treated,
         baseline_slope_treated = baseline_slope_treated,
         baseline_slope_control = baseline_slope_control,
         max_rate = max_rate, evidence = evidence),
    class = "recovery_attrition"
  )
}

#' Declare the planned analysis (answer strategy)
#'
#' The analysis that will actually be applied to the observed data,
#' stated exactly — including the specific inference method, since
#' finite-sample behavior can differ sharply across methods that share a
#' model formula (manuscript, section 3.2).
#'
#' Available estimators:
#' * `"linear_model"` — [stats::lm()] on the retained sample
#'   (complete-case), Wald t interval. Baseline adjustment is expressed
#'   through the formula, e.g. `y_observed ~ treatment + baseline`.
#' * `"lmm_random_intercept"` — linear mixed model via `lmerTest::lmer()`
#'   on the retained sample; `inference` selects `"satterthwaite"`
#'   (default), `"kenward_roger"` (requires `pbkrtest`), or `"wald_z"`
#'   (large-sample z interval).
#' * `"cluster_mean_ttest"` — two-sample t-test on equal-weighted cluster
#'   means (with cluster sizes fixed by design, equal and pupil-count
#'   weighting coincide), `n_clusters - 2` degrees of freedom.
#' * `"mi_baseline_adjusted"` — multiple imputation of missing outcomes
#'   from a normal model drawing on (treatment, baseline), followed by the
#'   baseline-adjusted linear model on each completed dataset and Rubin's
#'   rules with Barnard-Rubin degrees of freedom. Identifiable for the ITT
#'   estimand when dropout is MAR given the observed baseline; per the
#'   manuscript (section 2.3), evaluated through bias against the target,
#'   not through drift.
#'
#' @param estimator Character; one of the estimators above.
#' @param formula The planned model formula. The simulated data provide
#'   `y_observed`, `treatment`, `baseline`, and (for cluster designs)
#'   `cluster`. The treatment coefficient is the one whose name starts
#'   with `"treatment"`.
#' @param alpha Numeric. Two-sided significance level of the planned
#'   decision rule (default 0.05); also fixes the nominal `1 - alpha`
#'   confidence interval whose coverage is diagnosed.
#' @param inference Character. Inference method for
#'   `"lmm_random_intercept"`; ignored otherwise.
#' @param m_imputations Integer. Number of imputations for
#'   `"mi_baseline_adjusted"` (default 20).
#' @param degenerate_counts Logical. Whether degenerate fits (singular or
#'   boundary variance estimates; failure class (c) of the manuscript's
#'   taxonomy) count against the model-failure threshold. Design-specific
#'   and must be pre-specified; the choice is echoed in the report. Fatal
#'   errors and nonconvergence always count; diagnostic warnings are
#'   reported but never count. Default `FALSE` (degenerate fits are
#'   recorded, retained, and their marginal effect on coverage reported
#'   separately).
#'
#' @return An object of class `recovery_analysis`.
#' @export
planned_analysis <- function(estimator = c("linear_model",
                                           "lmm_random_intercept",
                                           "cluster_mean_ttest",
                                           "mi_baseline_adjusted"),
                             formula, alpha = 0.05,
                             inference = c("satterthwaite", "kenward_roger",
                                           "wald_z"),
                             m_imputations = 20,
                             degenerate_counts = FALSE) {
  estimator <- match.arg(estimator)
  inference <- match.arg(inference)
  stopifnot(
    "`formula` must be a formula" = inherits(formula, "formula"),
    "`alpha` must be in (0, 1)" =
      is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1,
    "`m_imputations` must be a whole number >= 2" =
      is.numeric(m_imputations) && length(m_imputations) == 1L &&
        m_imputations >= 2 && m_imputations == round(m_imputations),
    "`degenerate_counts` must be TRUE or FALSE" =
      is.logical(degenerate_counts) && length(degenerate_counts) == 1L &&
        !is.na(degenerate_counts)
  )
  structure(
    list(estimator = estimator, formula = formula, alpha = alpha,
         inference = inference, m_imputations = as.integer(m_imputations),
         degenerate_counts = degenerate_counts),
    class = "recovery_analysis"
  )
}

#' Declare a recovery design
#'
#' Step 2 of the protocol: assemble the target estimand, the model and
#' data strategy, and the answer strategy into a single executable
#' declaration. Every major design feature is either declared here or
#' recorded as omitted in the report.
#'
#' @param target A [target_estimand()] object.
#' @param data_strategy A [two_arm_trial()] or [cluster_trial()] object.
#' @param measurement A [measured_outcome()] object, or `NULL` to treat
#'   the outcome as perfectly reliable (recorded as an omission).
#' @param missingness An [attrition_model()] object, or `NULL` to model no
#'   attrition (recorded as an omission).
#' @param answer_strategy A [planned_analysis()] object.
#' @param effect Numeric or `NULL`. The declared (expected) effect on the
#'   target scale. Defaults to the SESOI. Target scenario rows are always
#'   simulated at the SESOI; if the declared effect exceeds the SESOI, an
#'   informational expected-effect scenario is added outside the verdict
#'   (a study justified by an optimistic expected effect but unable to
#'   recover the smallest meaningful effect should not receive a clean
#'   PASS). Pessimistic rows never shrink the effect; fragility to the
#'   effect size itself is a separate curve, [effect_fragility()].
#'
#' @return An object of class `recovery_design`.
#' @export
#' @examples
#' design <- declare_recovery(
#'   target = target_estimand(
#'     estimand = "ITT mean difference at 12 weeks",
#'     scale = "latent-outcome standardized mean difference",
#'     sesoi = 0.40
#'   ),
#'   data_strategy = two_arm_trial(n_per_arm = 100),
#'   measurement = measured_outcome(reliability = 0.70),
#'   missingness = attrition_model(rate = 0.15, mechanism = "differential"),
#'   answer_strategy = planned_analysis(
#'     estimator = "linear_model",
#'     formula = y_observed ~ treatment
#'   )
#' )
declare_recovery <- function(target, data_strategy, measurement = NULL,
                             missingness = NULL, answer_strategy,
                             effect = NULL) {
  stopifnot(
    "`target` must be a target_estimand()" =
      inherits(target, "recovery_target"),
    "`data_strategy` must be a two_arm_trial() or cluster_trial()" =
      inherits(data_strategy, "recovery_data_strategy"),
    "`answer_strategy` must be a planned_analysis()" =
      inherits(answer_strategy, "recovery_analysis")
  )
  if (!is.null(measurement)) {
    stopifnot("`measurement` must be a measured_outcome()" =
                inherits(measurement, "recovery_measurement"))
  }
  if (!is.null(missingness)) {
    stopifnot("`missingness` must be an attrition_model()" =
                inherits(missingness, "recovery_attrition"))
  }
  effect <- effect %||% target$sesoi
  stopifnot(
    "`effect` must be a single non-zero finite number" =
      is.numeric(effect) && length(effect) == 1L && is.finite(effect) &&
        effect != 0
  )
  if (inherits(data_strategy, "recovery_cluster") &&
      answer_strategy$estimator %in% c("linear_model",
                                       "mi_baseline_adjusted")) {
    warning("Cluster-randomized data strategy with an ",
            "independent-observations estimator: observations will be ",
            "treated as independent. This is recorded in the report.",
            call. = FALSE)
  }
  if (inherits(data_strategy, "recovery_two_arm") &&
      answer_strategy$estimator %in% c("lmm_random_intercept",
                                       "cluster_mean_ttest")) {
    stop("Cluster-based estimators require a cluster_trial() data strategy.",
         call. = FALSE)
  }

  omissions <- character(0)
  if (is.null(measurement)) {
    omissions <- c(omissions,
      "Measurement reliability was not declared; outcomes are treated as perfectly reliable.")
  }
  if (is.null(missingness) || missingness$rate == 0) {
    omissions <- c(omissions,
      "Attrition was not modeled; all recruited observations are treated as analyzed.")
  }
  if (inherits(data_strategy, "recovery_two_arm")) {
    omissions <- c(omissions,
      "Observations are treated as independent (no clustering declared).")
  }

  structure(
    list(target = target, effect = effect, data_strategy = data_strategy,
         measurement = measurement, missingness = missingness,
         answer_strategy = answer_strategy, omissions = omissions),
    class = "recovery_design"
  )
}
