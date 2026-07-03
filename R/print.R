#' @export
print.recovery_target <- function(x, ...) {
  cat("Target estimand:", x$estimand, "\n")
  cat("  Scale:", x$scale, "\n")
  cat(sprintf("  SESOI: %.3g   Delta (bias/drift unit): %.3g\n",
              x$sesoi, x$bias_scale_unit))
  if (!is.null(x$max_width)) {
    cat(sprintf("  Maximum acceptable CI width: %.3g\n", x$max_width))
  }
  invisible(x)
}

#' @export
print.recovery_design <- function(x, ...) {
  cat("Recovery design declaration\n")
  cat(strrep("-", 27), "\n", sep = "")
  cat("Target:  ", x$target$estimand, "\n")
  cat(sprintf("         SESOI %.3g on scale '%s'; declared expected effect %.3g\n",
              x$target$sesoi, x$target$scale, x$effect))
  for (ln in describe_design(x)) cat(ln, "\n")
  if (length(x$omissions)) {
    cat("Not modeled (silence must not imply ideality):\n")
    for (om in x$omissions) cat("  -", om, "\n")
  }
  cat("Run recovery_test() to simulate the scenario grid; verdict() and",
      "report() to evaluate.\n")
  invisible(x)
}

#' @export
print.recovery_result <- function(x, ...) {
  cat("Recovery-test result:", x$sims, "simulations per scenario row;",
      length(x$runs), "row(s);",
      sprintf("%.1f s elapsed.\n", x$elapsed_secs))
  for (nm in names(x$runs)) {
    run <- x$runs[[nm]]
    cat("\n--", run$scenario$label, "--\n")
    diag <- run$diagnosands
    diag$value <- signif(diag$value, 4)
    diag$mcse <- signif(diag$mcse, 2)
    print(diag[, c("diagnosand", "value", "mcse", "n_contributing",
                   "unstable")], row.names = FALSE)
  }
  cat("\nUse verdict() for the PASS/RISK/FAIL rule and report() for the",
      "standalone report.\n")
  invisible(x)
}
