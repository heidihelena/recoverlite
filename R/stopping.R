#' Run the recovery test under the algorithmic doubling stopping rule
#'
#' Wraps [recovery_test()] with the pre-declared simulation stopping rule
#' of the protocol (manuscript, section 2.4). The number of simulations
#' `S` starts at `start_sims` and is doubled whenever a required threshold
#' margin is still within `mcse_margin` (default 2) Monte Carlo standard
#' errors of its threshold — i.e. not yet resolvable at the current `S` —
#' up to a pre-declared ceiling `max_sims`. The rule stops as soon as the
#' verdict is *determined*:
#'
#' * a **stable declared failure** (a declared-nuisance margin failing by
#'   more than `mcse_margin` MCSEs) locks the verdict at FAIL — doubling
#'   cannot un-fail it, so the rule stops immediately; or
#' * **all required margins are resolved** (every required-row criterion
#'   is more than `mcse_margin` MCSEs from its threshold and no required
#'   conditional diagnosand is unstable), giving a clean PASS or a
#'   determined RISK.
#'
#' Anything still within `mcse_margin` MCSEs of a threshold at `max_sims`
#' is, by the protocol, reported as RISK: the run cannot resolve whether
#' the design clears the bright line, and that irresolution is itself the
#' finding.
#'
#' Doubling reuses the same `seed`; each `(seed, S)` run is fully
#' reproducible.
#'
#' @inheritParams recovery_test
#' @param start_sims Initial simulations per scenario row (default 2000).
#' @param max_sims Pre-declared maximum simulations per scenario row
#'   (default 16000, i.e. up to three doublings from 2000).
#' @param verbose Logical; if `TRUE`, report each doubling to the console.
#'
#' @return A `recovery_result` (as from [recovery_test()]) carrying an
#'   additional `"stopping"` attribute: a list with `start_sims`,
#'   `final_sims`, `max_sims`, `resolved` (whether the verdict was
#'   determined before the ceiling), and `unresolved` (a data frame of any
#'   required margins still within the MCSE band at the final `S`).
#' @export
recovery_test_stable <- function(design, start_sims = 2000L,
                                 max_sims = 16000L,
                                 scenarios = "confirmatory_grid",
                                 thresholds = recovery_thresholds(),
                                 seed = NULL, verbose = FALSE) {
  stopifnot(
    "`start_sims` must be a single whole number >= 2" =
      is.numeric(start_sims) && length(start_sims) == 1L && start_sims >= 2,
    "`max_sims` must be >= `start_sims`" =
      is.numeric(max_sims) && length(max_sims) == 1L && max_sims >= start_sims
  )
  s <- as.integer(start_sims)
  repeat {
    res <- recovery_test(design, sims = s, scenarios = scenarios,
                         thresholds = thresholds, seed = seed)
    v <- verdict(res)
    unres <- unresolved_required(res, thresholds)
    # FAIL is locked by a stable declared failure; otherwise the verdict
    # is determined once no required margin remains within the MCSE band.
    determined <- v$verdict == "FAIL" || nrow(unres) == 0L
    if (verbose) {
      message(sprintf("S = %d: verdict %s; %d required margin(s) within %g MCSE",
                      s, v$verdict, nrow(unres), thresholds$mcse_margin))
    }
    if (determined || s >= max_sims) break
    s <- min(2L * s, as.integer(max_sims))
  }
  attr(res, "stopping") <- list(
    start_sims = as.integer(start_sims), final_sims = s,
    max_sims = as.integer(max_sims),
    resolved = nrow(unres) == 0L || v$verdict == "FAIL",
    hit_ceiling = s >= max_sims && nrow(unres) > 0L && v$verdict != "FAIL",
    unresolved = unres)
  res
}

# Required-row criteria whose margins are still within the MCSE band
# (either side of the threshold) or whose conditional diagnosand is
# unstable — the margins that doubling `S` can still resolve. Mirrors the
# required-row selection of verdict_under().
unresolved_required <- function(result, thr) {
  estimation <- thr$profile == "estimation"
  runs <- Filter(function(r) r$scenario$counts_for != "informational" &&
                   (!estimation || r$scenario$row_type == "target"),
                 result$runs)
  required_rows <- if (estimation) {
    c("target_declared", "target_pessimistic")
  } else {
    c("null_declared", "null_pessimistic", "target_declared",
      "target_pessimistic")
  }
  out <- list()
  for (nm in intersect(required_rows, names(runs))) {
    ev <- evaluate_criteria(runs[[nm]]$diagnosands,
                            runs[[nm]]$scenario$row_type, thr, result$alpha)
    near <- !is.na(ev$pass) & is.finite(ev$mcse) &
      abs(ev$margin) <= thr$mcse_margin * ev$mcse
    keep <- near | ev$unstable
    if (any(keep)) out[[nm]] <- cbind(scenario = nm, ev[keep, , drop = FALSE])
  }
  res <- do.call(rbind, c(list(NULL), out))
  if (is.null(res)) {
    data.frame(scenario = character(0), criterion = character(0),
               value = numeric(0), mcse = numeric(0),
               requirement = character(0), margin = numeric(0),
               pass = logical(0), stable = logical(0), unstable = logical(0),
               n_contributing = numeric(0), note = character(0),
               upper = numeric(0), stringsAsFactors = FALSE)
  } else {
    rownames(res) <- NULL
    res
  }
}
