#' Apply the PASS/RISK/FAIL verdict rule
#'
#' Applies the fixed verdict rule (manuscript, section 2.2, Step 5)
#' mechanically to a [recovery_test()] result. Once thresholds and
#' scenarios are fixed, the computation is mechanical — the judgment
#' lives in the declaration and the scenario grid, not here.
#'
#' * **PASS** — all required thresholds met under all scenario rows the
#'   selected profile requires, with every margin exceeding
#'   `mcse_margin` (default 2) Monte Carlo standard errors.
#' * **RISK** — thresholds met under the declared-nuisance rows but not
#'   under a pessimistic row; or any margin within `mcse_margin` MCSEs of
#'   a threshold; or any required conditional diagnosand too unstable to
#'   evaluate; or required scenario rows missing from the run.
#' * **FAIL** — one or more required thresholds fail under a
#'   declared-nuisance row — including an inflated target-null rejection
#'   rate under Null-declared. The planned study cannot recover the
#'   target even if the declared assumptions are correct.
#'
#' The verdict is recomputed under the shipped strict and lenient
#' profiles (section 2.5): agreement across profiles carries more
#' decision weight than any one profile, and a verdict that flips is
#' itself a finding — the RISK category exists to hold it. The verdict is
#' a decision convention, not a validity classification, and the full
#' [report()] always travels with it.
#'
#' @param result A `recovery_result` from [recovery_test()].
#'
#' @return An object of class `recovery_verdict`.
#' @export
verdict <- function(result) {
  stopifnot(
    "`result` must be a recovery_result from recovery_test()" =
      inherits(result, "recovery_result")
  )
  thr <- result$thresholds
  main <- verdict_under(result, thr)
  strict <- lenient <- NULL
  if (thr$profile != "estimation") {
    strict <- verdict_under(result, shipped_variant(thr, "strict"))
    lenient <- verdict_under(result, shipped_variant(thr, "lenient"))
  }
  structure(
    list(verdict = main$verdict, binding = main$binding,
         evaluations = main$evaluations, notes = main$notes,
         smallest_margin = main$smallest_margin,
         verdict_strict = strict$verdict %||% NA_character_,
         verdict_lenient = lenient$verdict %||% NA_character_,
         thresholds = thr),
    class = "recovery_verdict"
  )
}

# A shipped profile carrying over the run's MCSE / stability settings.
shipped_variant <- function(thr, profile) {
  recovery_thresholds(profile = profile,
                      overcoverage_flag = thr$overcoverage_flag,
                      mcse_margin = thr$mcse_margin,
                      min_conditional_n = thr$min_conditional_n,
                      max_width = thr$max_width)
}

# Compute the verdict under one threshold profile.
verdict_under <- function(result, thr) {
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
  missing_rows <- setdiff(required_rows, names(runs))

  evaluations <- lapply(runs, function(run) {
    evaluate_criteria(run$diagnosands, run$scenario$row_type, thr,
                      result$alpha)
  })

  # counts_for = NULL selects every evaluated row.
  pick <- function(counts_for, cond) {
    do.call(rbind, c(list(NULL), lapply(names(evaluations), function(nm) {
      if (!is.null(counts_for) &&
          runs[[nm]]$scenario$counts_for != counts_for) return(NULL)
      ev <- evaluations[[nm]]
      ev <- ev[cond(ev), , drop = FALSE]
      if (nrow(ev)) cbind(scenario = nm, ev) else NULL
    })))
  }
  failed_declared <- pick("declared",
                          function(ev) !is.na(ev$pass) & !ev$pass)
  failed_pess <- pick("pessimistic",
                      function(ev) !is.na(ev$pass) & !ev$pass)
  narrow <- pick(NULL, function(ev) !is.na(ev$pass) & ev$pass & !ev$stable)
  unstable <- pick(NULL, function(ev) ev$unstable)

  all_ev <- do.call(rbind, c(list(NULL), lapply(names(evaluations),
    function(nm) cbind(scenario = nm, evaluations[[nm]]))))
  smallest_margin <- if (!is.null(all_ev) && any(is.finite(all_ev$margin))) {
    i <- which.min(all_ev$margin)
    sprintf("%s under '%s' (signed margin %+.4f)", all_ev$criterion[i],
            all_ev$scenario[i], all_ev$margin[i])
  } else NA_character_

  notes <- character(0)
  if (!is.null(failed_declared) && nrow(failed_declared)) {
    v <- "FAIL"
    binding <- sprintf(
      "Under declared-nuisance rows, %s.",
      paste(sprintf("%s = %.3g under '%s' violates %s",
                    failed_declared$criterion, failed_declared$value,
                    failed_declared$scenario, failed_declared$requirement),
            collapse = "; "))
  } else if (!is.null(failed_pess) && nrow(failed_pess)) {
    v <- "RISK"
    binding <- sprintf(
      "Thresholds hold under declared-nuisance rows but fail under pessimistic rows: %s.",
      paste(sprintf("%s = %.3g under '%s' violates %s", failed_pess$criterion,
                    failed_pess$value, failed_pess$scenario,
                    failed_pess$requirement),
            collapse = "; "))
  } else if (length(missing_rows)) {
    v <- "RISK"
    binding <- sprintf(
      "Required scenario rows not evaluated: %s. A PASS requires all rows the '%s' profile requires; rerun with scenarios = \"%s\".",
      paste(missing_rows, collapse = ", "), thr$profile,
      if (estimation) "target_grid" else "confirmatory_grid")
  } else if ((!is.null(narrow) && nrow(narrow)) ||
             (!is.null(unstable) && nrow(unstable))) {
    v <- "RISK"
    parts <- character(0)
    if (!is.null(narrow) && nrow(narrow)) {
      parts <- c(parts, paste(sprintf(
        "%s margin under '%s' is within %g MCSE of its threshold",
        narrow$criterion, narrow$scenario, thr$mcse_margin),
        collapse = "; "))
    }
    if (!is.null(unstable) && nrow(unstable)) {
      parts <- c(parts, paste(sprintf(
        "%s under '%s' is unstable (%d contributing simulations)",
        unstable$criterion, unstable$scenario, unstable$n_contributing),
        collapse = "; "))
    }
    binding <- paste0(
      "All point estimates meet the thresholds, but simulation precision ",
      "is insufficient to confirm a PASS: ", paste(parts, collapse = ". "),
      ". Increase `sims`.")
  } else {
    v <- "PASS"
    binding <- NA_character_
  }
  if (!is.null(unstable) && nrow(unstable) && v == "FAIL") {
    notes <- c(notes, paste(
      "Some conditional diagnosands were unstable and excluded from the",
      "threshold evaluation; see the diagnosand tables."))
  }

  list(verdict = v, binding = binding, evaluations = evaluations,
       notes = notes, smallest_margin = smallest_margin)
}

# One scenario row's diagnosands -> per-criterion evaluation with signed
# margins (positive = passing). Unstable conditional diagnosands get
# pass = NA: they cannot cause a FAIL, but they block a clean PASS.
evaluate_criteria <- function(diag, row_type, thr, alpha) {
  g <- function(name, col = "value") diag[diag$diagnosand == name, col]
  estimation <- thr$profile == "estimation"
  crit <- list()

  if (!estimation) {
    if (row_type == "null") {
      lim <- thr$null_rejection_mult * alpha
      crit <- c(crit, list(list(
        criterion = "target_null_rejection", value = g("rejection_rate"),
        mcse = g("rejection_rate", "mcse"),
        requirement = sprintf("<= %.4g (%.3g x alpha)", lim,
                              thr$null_rejection_mult),
        margin = lim - g("rejection_rate"),
        n_contributing = g("rejection_rate", "n_contributing"),
        unstable = FALSE, note = "")))
    } else {
      crit <- c(crit, list(list(
        criterion = "power", value = g("rejection_rate"),
        mcse = g("rejection_rate", "mcse"),
        requirement = sprintf(">= %.3g", thr$power),
        margin = g("rejection_rate") - thr$power,
        n_contributing = g("rejection_rate", "n_contributing"),
        unstable = FALSE, note = "")))
    }
  }

  crit <- c(crit, list(list(
    criterion = "target_bias", value = g("target_bias"),
    mcse = g("target_bias", "mcse"),
    requirement = sprintf("|value| <= %.3g Delta", thr$target_bias),
    margin = thr$target_bias - abs(g("target_bias")),
    n_contributing = g("target_bias", "n_contributing"),
    unstable = FALSE, note = "")))

  cov <- g("coverage")
  cov_note <- if (!is.na(cov) && cov > thr$overcoverage_flag) {
    sprintf("overcoverage (> %.3g): flagged as inefficiency, evaluated through precision, not failure",
            thr$overcoverage_flag)
  } else ""
  crit <- c(crit, list(list(
    criterion = "coverage", value = cov, mcse = g("coverage", "mcse"),
    requirement = sprintf(">= %.3g (lower bound)", thr$coverage),
    margin = cov - thr$coverage,
    n_contributing = g("coverage", "n_contributing"),
    unstable = FALSE, note = cov_note)))

  if (!estimation && row_type == "target") {
    crit <- c(crit, list(
      list(criterion = "type_s", value = g("type_s"),
           mcse = g("type_s", "mcse"),
           requirement = sprintf("<= %.3g", thr$type_s),
           margin = thr$type_s - g("type_s"),
           n_contributing = g("type_s", "n_contributing"),
           unstable = isTRUE(g("type_s", "unstable")), note = ""),
      list(criterion = "type_m", value = g("type_m"),
           mcse = g("type_m", "mcse"),
           requirement = sprintf("<= %.3g", thr$type_m),
           margin = thr$type_m - g("type_m"),
           n_contributing = g("type_m", "n_contributing"),
           unstable = isTRUE(g("type_m", "unstable")), note = "")))
  }

  if (estimation) {
    crit <- c(crit, list(list(
      criterion = "estimand_drift", value = g("estimand_drift"),
      mcse = g("estimand_drift", "mcse"),
      requirement = sprintf("|value| <= %.3g Delta", thr$drift),
      margin = thr$drift - abs(g("estimand_drift")),
      n_contributing = g("estimand_drift", "n_contributing"),
      unstable = FALSE, note = "")))
    if (!is.null(thr$max_width)) {
      crit <- c(crit, list(list(
        criterion = "precision", value = g("precision"),
        mcse = g("precision", "mcse"),
        requirement = sprintf("mean width <= declared max %.3g",
                              thr$max_width),
        margin = thr$max_width - g("precision"),
        n_contributing = g("precision", "n_contributing"),
        unstable = FALSE, note = "")))
    }
  }

  crit <- c(crit, list(list(
    criterion = "model_failure", value = g("model_failure"),
    mcse = g("model_failure", "mcse"),
    requirement = sprintf("<= %.3g (counted classes)", thr$model_failure),
    margin = thr$model_failure - g("model_failure"),
    n_contributing = g("model_failure", "n_contributing"),
    unstable = FALSE, note = "")))

  out <- do.call(rbind, lapply(crit, function(cr) {
    unstable <- cr$unstable || is.na(cr$value)
    pass <- if (unstable) NA else cr$margin >= 0
    stable <- if (unstable || is.na(cr$mcse)) FALSE else
      cr$margin > thr$mcse_margin * cr$mcse
    data.frame(criterion = cr$criterion, value = cr$value, mcse = cr$mcse,
               requirement = cr$requirement, margin = cr$margin,
               pass = pass, stable = stable, unstable = unstable,
               n_contributing = cr$n_contributing, note = cr$note,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' @export
print.recovery_verdict <- function(x, ...) {
  cat("Recovery-test verdict:", x$verdict,
      sprintf("(profile '%s')\n", x$thresholds$profile))
  if (!is.na(x$verdict_strict)) {
    cat(sprintf("Under shipped profiles: strict %s | default-family %s | lenient %s%s\n",
                x$verdict_strict, x$verdict, x$verdict_lenient,
                if (length(unique(c(x$verdict_strict, x$verdict,
                                    x$verdict_lenient))) > 1)
                  "  <- profile disagreement is itself a finding" else ""))
  }
  cat("Threshold set:", x$thresholds$version,
      if (length(x$thresholds$modified))
        paste0("(deviations: ", paste(x$thresholds$modified, collapse = ", "), ")")
      else "(shipped profile)", "\n")
  if (!is.na(x$binding)) cat("Binding failure mode:", x$binding, "\n")
  if (!is.na(x$smallest_margin)) {
    cat("Smallest signed margin:", x$smallest_margin, "\n")
  }
  for (nm in names(x$evaluations)) {
    cat("\n--", nm, "--\n")
    ev <- x$evaluations[[nm]]
    ev$value <- signif(ev$value, 3)
    ev$mcse <- signif(ev$mcse, 2)
    ev$margin <- signif(ev$margin, 2)
    print(ev[, c("criterion", "value", "mcse", "requirement", "margin",
                 "pass", "stable", "unstable")], row.names = FALSE)
  }
  if (length(x$notes)) cat("\nNotes:", paste(x$notes, collapse = "\n"), "\n")
  cat("\nThe verdict is a decision convention, not a validity",
      "classification; the full report() must travel with it.\n")
  invisible(x)
}
