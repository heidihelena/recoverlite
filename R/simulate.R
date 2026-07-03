# Simulation engine (protocol Step 4).
#
# Data generation is declared with DeclareDesign / fabricatr / randomizr;
# the planned analysis is applied to each simulated dataset by
# fit_analysis(), which records the estimate, interval, decision, and a
# CLASSIFIED failure indicator (fatal / nonconvergence / degenerate /
# warnings — manuscript, section 2.3). Each simulation also records the
# analyzable-data contrast theta_obs_s: the contrast obtained by applying
# the analysis-population rule to that dataset's full data before
# estimation, evaluated on true (pre-measurement-error) outcomes. This is
# the raw material for the estimand-drift diagnosand and the exact
# decomposition target bias = estimator bias + drift.

# Calibrate the dropout intercept alpha so that the marginal rate
# E_B[plogis(alpha + slope * B)] equals `rate` for B ~ N(0, 1).
calibrate_dropout_intercept <- function(rate, slope) {
  # Finite stand-in for "never drops out": plogis(-30) ~ 1e-13, and it
  # keeps the arm-interaction arithmetic finite when the other arm has a
  # nonzero rate (-Inf would propagate NaN through (a1 - a0) * treatment).
  if (rate <= 0) return(-30)
  if (slope == 0) return(stats::qlogis(rate))
  f <- function(a) {
    stats::integrate(function(b) stats::plogis(a + slope * b) * stats::dnorm(b),
                     -Inf, Inf)$value - rate
  }
  stats::uniroot(f, c(-25, 25))$root
}

# Build the DeclareDesign object (model + assignment + reveal +
# measurement + attrition; no estimator) for one scenario row.
build_dd_design <- function(p) {
  if (p$type == "two_arm_trial") {
    n_total <- 2L * p$n_per_arm
    m_treat <- round(n_total * p$allocation)
    eff <- p$effect
    rho <- p$rho
    nc <- p$noncompliance %||% 0

    dd <- DeclareDesign::declare_model(
      N = n_total,
      baseline = rnorm(N),
      eps = rnorm(N),
      complier = rbinom(N, 1, 1 - nc),
      fabricatr::potential_outcomes(
        y_true ~ eff * complier * treatment + rho * baseline +
          sqrt(1 - rho^2) * eps,
        conditions = list(treatment = c(0, 1))
      )
    ) +
      DeclareDesign::declare_assignment(
        treatment = randomizr::complete_ra(N = N, m = m_treat)
      ) +
      DeclareDesign::declare_reveal(y_true, treatment)
  } else {
    n_cl <- p$n_clusters
    n_pp <- p$n_per_cluster
    m_cl <- round(n_cl * p$allocation)
    eff <- p$effect
    icc <- p$icc

    dd <- DeclareDesign::declare_model(
      cluster = fabricatr::add_level(
        N = n_cl,
        u_c = rnorm(N, 0, sqrt(icc))
      ),
      pupil = fabricatr::add_level(
        N = n_pp,
        e = rnorm(N, 0, sqrt(1 - icc)),
        baseline = rnorm(N),
        fabricatr::potential_outcomes(
          y_true ~ eff * treatment + u_c + e,
          conditions = list(treatment = c(0, 1))
        )
      )
    ) +
      DeclareDesign::declare_assignment(
        treatment = randomizr::cluster_ra(clusters = cluster, m = m_cl)
      ) +
      DeclareDesign::declare_reveal(y_true, treatment)
  }

  # Measurement: classical ADDITIVE error, Var(e) = 1/r - 1, so that
  # reliability = Var(y_true) / Var(y_observed). The raw treatment
  # contrast is unbiased in expectation; the error inflates residual
  # variance (a precision/power failure, not a bias failure).
  rel <- p$reliability
  if (!is.null(rel) && rel < 1) {
    err_sd <- sqrt(1 / rel - 1)
    dd <- dd + DeclareDesign::declare_measurement(
      y_observed = y_true + err_sd * rnorm(N)
    )
  } else {
    dd <- dd + DeclareDesign::declare_measurement(y_observed = y_true)
  }

  # Attrition: logit Pr(dropout | Z, B) = alpha_Z + gamma_Z * B, with
  # alpha_Z calibrated to the declared arm-specific marginal rates.
  # retained == TRUE enters the analysis population.
  att <- p$attrition
  if (!is.null(att)) {
    a0 <- calibrate_dropout_intercept(att$rate_control, att$slope_control)
    a1 <- calibrate_dropout_intercept(att$rate_treated, att$slope_treated)
    g0 <- att$slope_control
    g1 <- att$slope_treated
    dd <- dd + DeclareDesign::declare_measurement(
      p_drop = plogis((a0 + (a1 - a0) * treatment) +
                        (g0 + (g1 - g0) * treatment) * baseline),
      retained = runif(N) > p_drop
    )
  } else {
    dd <- dd + DeclareDesign::declare_measurement(retained = rep(TRUE, N))
  }

  dd
}

# Run one scenario row: draw `sims` datasets and apply the planned
# analysis to each.
run_scenario <- function(design, params, sims) {
  dd <- build_dd_design(params)
  analysis <- design$answer_strategy
  theta <- params$effect

  rows <- vector("list", sims)
  for (i in seq_len(sims)) {
    dat <- DeclareDesign::draw_data(dd)
    fit <- fit_analysis(analysis, dat)
    ret <- dat$retained
    z <- dat$treatment
    # Analyzable-data contrast on true (pre-measurement-error) outcomes.
    theta_obs <- mean(dat$y_true[ret & z == 1]) -
      mean(dat$y_true[ret & z == 0])
    counted <- fit$fatal || fit$nonconverged ||
      (fit$degenerate && analysis$degenerate_counts)
    rows[[i]] <- data.frame(
      sim = i,
      estimate = fit$estimate, se = fit$se,
      ci_lo = fit$ci_lo, ci_hi = fit$ci_hi, p = fit$p,
      sig = if (is.na(fit$p)) NA else fit$p < analysis$alpha,
      covered = if (is.na(fit$ci_lo)) NA else
        (fit$ci_lo <= theta && theta <= fit$ci_hi),
      covered_obs = if (is.na(fit$ci_lo)) NA else
        (fit$ci_lo <= theta_obs && theta_obs <= fit$ci_hi),
      fatal = fit$fatal, nonconverged = fit$nonconverged,
      degenerate = fit$degenerate, warned = fit$warned,
      counted_failure = counted,
      theta_obs = theta_obs,
      n_analyzed = sum(ret),
      attrition_realized = 1 - mean(ret)
    )
  }
  out <- do.call(rbind, rows)
  attr(out, "degenerate_counts") <- analysis$degenerate_counts
  out
}

# ---------------------------------------------------------------------
# Answer strategies. Each returns estimate, CI, p, and the four-class
# failure record. Hard errors are recorded, not raised: a fit that cannot
# be produced is itself a diagnosand.
# ---------------------------------------------------------------------
fit_analysis <- function(analysis, dat) {
  out <- list(estimate = NA_real_, se = NA_real_, ci_lo = NA_real_,
              ci_hi = NA_real_, p = NA_real_, fatal = FALSE,
              nonconverged = FALSE, degenerate = FALSE, warned = FALSE)
  adat <- dat[dat$retained, , drop = FALSE]

  res <- switch(analysis$estimator,
    linear_model = fit_lm(analysis, adat),
    lmm_random_intercept = fit_lmm(analysis, adat),
    cluster_mean_ttest = fit_cluster_ttest(analysis, adat),
    mi_baseline_adjusted = fit_mi(analysis, dat)
  )
  utils::modifyList(out, res)
}

extract_treatment_row <- function(co) {
  row <- grep("^treatment", rownames(co))[1]
  if (is.na(row)) NULL else row
}

fit_lm <- function(analysis, adat) {
  alpha <- analysis$alpha
  fit <- tryCatch(stats::lm(analysis$formula, data = adat),
                  error = function(e) e)
  if (inherits(fit, "error")) return(list(fatal = TRUE))
  co <- stats::coef(summary(fit))
  row <- extract_treatment_row(co)
  if (is.null(row)) return(list(fatal = TRUE))
  crit <- stats::qt(1 - alpha / 2, fit$df.residual)
  est <- co[row, 1]; se <- co[row, 2]
  list(estimate = est, se = se, p = co[row, 4],
       ci_lo = est - crit * se, ci_hi = est + crit * se)
}

fit_lmm <- function(analysis, adat) {
  if (!requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("lmerTest", quietly = TRUE)) {
    stop("estimator 'lmm_random_intercept' requires the suggested packages ",
         "'lme4' and 'lmerTest'", call. = FALSE)
  }
  if (analysis$inference == "kenward_roger" &&
      !requireNamespace("pbkrtest", quietly = TRUE)) {
    stop("inference = 'kenward_roger' requires the suggested package ",
         "'pbkrtest'", call. = FALSE)
  }
  alpha <- analysis$alpha
  warns <- character(0)
  fit <- withCallingHandlers(
    tryCatch(lmerTest::lmer(analysis$formula, data = adat),
             error = function(e) e),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage")
  )
  if (inherits(fit, "error")) return(list(fatal = TRUE))

  degenerate <- isTRUE(lme4::isSingular(fit, tol = 1e-4))
  # lme4 reports singular/boundary fits through the same message slot as
  # true optimizer failures; those belong to failure class (c) degenerate,
  # not class (b) nonconvergence, and must not be counted as (b) — the
  # taxonomy keeps them separate so that the pre-specified
  # `degenerate_counts` choice governs whether they count.
  conv_msgs <- fit@optinfo$conv$lme4$messages
  conv_msgs <- conv_msgs[!grepl("singular|boundary", conv_msgs,
                                ignore.case = TRUE)]
  nonconv_warns <- warns[grepl("converg", warns, ignore.case = TRUE) &
                           !grepl("singular|boundary", warns,
                                  ignore.case = TRUE)]
  nonconverged <- length(conv_msgs) > 0 || length(nonconv_warns) > 0
  warned <- length(warns) > 0

  if (analysis$inference == "wald_z") {
    co <- tryCatch(lme4::fixef(fit), error = function(e) NULL)
    vc <- tryCatch(as.matrix(stats::vcov(fit)), error = function(e) NULL)
    if (is.null(co) || is.null(vc)) return(list(fatal = TRUE))
    idx <- grep("^treatment", names(co))[1]
    if (is.na(idx)) return(list(fatal = TRUE))
    est <- co[idx]; se <- sqrt(vc[idx, idx])
    crit <- stats::qnorm(1 - alpha / 2)
    p <- 2 * stats::pnorm(abs(est / se), lower.tail = FALSE)
    return(list(estimate = unname(est), se = se, p = unname(p),
                ci_lo = unname(est - crit * se),
                ci_hi = unname(est + crit * se),
                nonconverged = nonconverged, degenerate = degenerate,
                warned = warned))
  }

  ddf <- if (analysis$inference == "kenward_roger") "Kenward-Roger"
         else "Satterthwaite"
  co <- tryCatch(stats::coef(summary(fit, ddf = ddf)),
                 error = function(e) NULL)
  row <- if (is.null(co)) NULL else extract_treatment_row(co)
  if (is.null(co) || is.null(row)) return(list(fatal = TRUE))
  est <- co[row, "Estimate"]
  se <- co[row, "Std. Error"]
  df <- co[row, "df"]
  p <- co[row, "Pr(>|t|)"]
  if (!is.finite(df) || df <= 0) {
    # Approximate df can fail on degenerate fits; fall back to a
    # between-cluster df and recompute the p-value from the t statistic.
    df <- max(length(unique(adat$cluster)) - 2, 1)
    p <- 2 * stats::pt(abs(est / se), df, lower.tail = FALSE)
  }
  crit <- stats::qt(1 - alpha / 2, df)
  list(estimate = est, se = se, p = p,
       ci_lo = est - crit * se, ci_hi = est + crit * se,
       nonconverged = nonconverged, degenerate = degenerate,
       warned = warned)
}

# Two-sample t-test on equal-weighted cluster means, n_clusters - 2 df.
fit_cluster_ttest <- function(analysis, adat) {
  alpha <- analysis$alpha
  res <- tryCatch({
    means <- tapply(adat$y_observed, adat$cluster, mean)
    zc <- tapply(adat$treatment, adat$cluster, function(x) x[1])
    m1 <- means[zc == 1]; m0 <- means[zc == 0]
    n1 <- length(m1); n0 <- length(m0)
    est <- mean(m1) - mean(m0)
    sp2 <- ((n1 - 1) * stats::var(m1) + (n0 - 1) * stats::var(m0)) /
      (n1 + n0 - 2)
    se <- sqrt(sp2 * (1 / n1 + 1 / n0))
    df <- n1 + n0 - 2
    tval <- est / se
    crit <- stats::qt(1 - alpha / 2, df)
    list(estimate = est, se = se,
         p = 2 * stats::pt(abs(tval), df, lower.tail = FALSE),
         ci_lo = est - crit * se, ci_hi = est + crit * se)
  }, error = function(e) NULL)
  if (is.null(res)) return(list(fatal = TRUE))
  res
}

# Multiple imputation of missing outcomes from a normal model drawing on
# the observed data (proper MI: posterior draws of the completer
# regression), baseline-adjusted analysis of each completed dataset, and
# Rubin's rules with Barnard-Rubin degrees of freedom. Identifiable for
# the ITT estimand when dropout is MAR given the observed baseline.
fit_mi <- function(analysis, dat) {
  alpha <- analysis$alpha
  m <- analysis$m_imputations
  res <- tryCatch({
    y <- ifelse(dat$retained, dat$y_observed, NA_real_)
    X <- stats::model.matrix(~ treatment + baseline, data = dat)
    obs <- !is.na(y)
    Xo <- X[obs, , drop = FALSE]
    yo <- y[obs]
    fit0 <- stats::lm.fit(Xo, yo)
    k <- ncol(Xo)
    dfr <- length(yo) - k
    rss <- sum(fit0$residuals^2)
    XtX_inv <- chol2inv(chol(crossprod(Xo)))

    ests <- vars <- numeric(m)
    for (j in seq_len(m)) {
      sigma2 <- rss / stats::rchisq(1, dfr)
      # N(0, V) draw needs the LOWER factor: Var(t(chol(V)) z) = V.
      beta <- fit0$coefficients +
        drop(t(chol(sigma2 * XtX_inv)) %*% rnorm(k))
      y_imp <- y
      y_imp[!obs] <- drop(X[!obs, , drop = FALSE] %*% beta) +
        sqrt(sigma2) * rnorm(sum(!obs))
      fit_j <- stats::lm(y_imp ~ treatment + baseline, data = dat)
      co <- stats::coef(summary(fit_j))
      row <- extract_treatment_row(co)
      ests[j] <- co[row, 1]
      vars[j] <- co[row, 2]^2
    }
    qbar <- mean(ests)
    ubar <- mean(vars)
    b <- stats::var(ests)
    tt <- ubar + (1 + 1 / m) * b
    # Barnard-Rubin small-sample degrees of freedom
    df_com <- nrow(dat) - k
    lam <- (1 + 1 / m) * b / tt
    df_old <- (m - 1) / lam^2
    df_obs <- ((df_com + 1) / (df_com + 3)) * df_com * (1 - lam)
    df_br <- df_old * df_obs / (df_old + df_obs)
    se <- sqrt(tt)
    crit <- stats::qt(1 - alpha / 2, df_br)
    list(estimate = qbar, se = se,
         p = 2 * stats::pt(abs(qbar / se), df_br, lower.tail = FALSE),
         ci_lo = qbar - crit * se, ci_hi = qbar + crit * se)
  }, error = function(e) NULL)
  if (is.null(res)) return(list(fatal = TRUE))
  res
}
