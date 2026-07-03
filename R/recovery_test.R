#' Run the recovery test
#'
#' Simulates the declared design over the crossed scenario grid
#' (manuscript, section 2.2, Step 3): null and target effects, each under
#' declared and pessimistically perturbed nuisance assumptions, plus an
#' informational expected-effect row when the declared effect exceeds the
#' SESOI. The planned analysis is applied to every simulated dataset and
#' the diagnosands of section 2.3 are computed with Monte Carlo standard
#' errors and explicit inclusion rules.
#'
#' @param design A [declare_recovery()] object.
#' @param sims Integer. Simulations per scenario row. The default 2000 is
#'   an initial working number, not a standard: the relevant stopping rule
#'   is whether the MCSE is small enough to support the verdict, and
#'   verdicts near a threshold may require many more.
#' @param scenarios `"confirmatory_grid"` (default: the four verdict rows
#'   Null-declared, Null-pessimistic, Target-declared,
#'   Target-pessimistic) or `"target_grid"` (the two target rows only,
#'   for estimation-focused designs evaluated with the estimation
#'   profile).
#' @param thresholds A [recovery_thresholds()] object. The profile and any
#'   deviations from it are echoed in the report.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class `recovery_result`. Use [verdict()] and
#'   [report()] on it.
#' @export
recovery_test <- function(design, sims = 2000,
                          scenarios = "confirmatory_grid",
                          thresholds = recovery_thresholds(),
                          seed = NULL) {
  stopifnot(
    "`design` must be a declare_recovery() object" =
      inherits(design, "recovery_design"),
    "`thresholds` must be a recovery_thresholds() object" =
      inherits(thresholds, "recovery_thresholds"),
    "`sims` must be a single whole number >= 2" =
      is.numeric(sims) && length(sims) == 1L && sims >= 2 &&
        sims == round(sims)
  )
  if (is.null(thresholds$max_width) && !is.null(design$target$max_width)) {
    thresholds$max_width <- design$target$max_width
  }
  if (!is.null(seed)) set.seed(seed)
  scs <- build_scenarios(design, scenarios)
  has_attrition <- !is.null(design$missingness) && design$missingness$rate > 0

  t0 <- Sys.time()
  runs <- lapply(scs, function(sc) {
    sim_df <- run_scenario(design, sc$params, sims)
    diag <- compute_diagnosands(
      sim_df,
      theta = sc$params$effect,
      unit = design$target$bias_scale_unit,
      thresholds = thresholds,
      row_type = sc$row_type,
      has_attrition = has_attrition,
      alpha = design$answer_strategy$alpha
    )
    list(scenario = sc, sim_data = sim_df, diagnosands = diag,
         theta = sc$params$effect)
  })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  structure(
    list(design = design, runs = runs, sims = as.integer(sims),
         thresholds = thresholds, seed = seed,
         alpha = design$answer_strategy$alpha,
         scenario_request = scenarios,
         evidence_tiers = attr(scs, "tiers"),
         elapsed_secs = elapsed,
         session = paste(R.version.string, "| recoverlite",
                         as.character(utils::packageVersion("recoverlite")))),
    class = "recovery_result"
  )
}
