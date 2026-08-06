# Hardcoded absolute path for Rscript compatibility
source("/Users/sigma/SynologyDrive/kuanke/Downloads/jaerevision/mfxtsemipar/R/mfxtsemipar_utils.R")


#' Mixed-frequency semiparametric regression with fixed knots
#'
#' R implementation of Stata's \code{mfxtsemipar}. Estimates a semiparametric
#' regression with user-specified knots (no cross-validation).
#'
#' @param hf data.frame or data.table containing the high-frequency observations.
#' @param lf data.frame or data.table containing one observation per \code{id*tl}.
#' @param y name of the dependent variable (in \code{lf}).
#' @param x character vector of low-frequency covariates (in \code{lf}).
#' @param uvar name of the semiparametric variable (in \code{hf}).
#' @param id name of the panel identifier.
#' @param tl character vector of time-level variables.
#' @param gen name of the generated fitted-values variable.
#' @param hfcov character vector of high-frequency covariates.
#' @param cluster name of the cluster variable.
#' @param bknots numeric vector of length 2 with boundary knots for the spline
#'   basis. If \code{NULL} (default), boundary knots are set to the minimum and
#'   maximum of \code{uvar} (with a small padding). When supplied, the minimum
#'   and maximum of \code{uvar} must lie within \code{bknots}; otherwise an
#'   error is raised. Interior knot locations are still determined from
#'   \code{uvar}'s range or from \code{startp}/\code{endp}, not from
#'   \code{bknots}.
#' @param startp numeric; optional preset minimum internal knot. Together with
#'   \code{endp}, interior knots are placed between \code{startp} and
#'   \code{endp}; otherwise they follow the range of \code{uvar}.
#' @param endp numeric; optional preset maximum internal knot.
#' @param degree polynomial degree.
#' @param knots numeric vector of interior knot locations.
#' @param type spline type: \code{"poly"}, \code{"bs"}, \code{"ms"},
#'   \code{"is"} or \code{"ibs"}; default \code{"poly"}.
#' @param winsor winsorization percentiles (length 2).
#' @param winsor_values logical; see \code{winsor}.
#' @param eqspace use equally spaced knots when generating from \code{nknots}.
#' @param nknots integer; number of knots to generate (alternative to \code{knots}).
#' @param center centering value for the spline basis.
#' @param absorb fixed-effects specification.
#' @param atu alternative evaluation variable (in \code{hf}).
#' @param intercept logical; if \code{TRUE}, include the intercept/base term.
#' @param brep number of wild-bootstrap replications.
#' @param predy name for full LF prediction.
#' @param weights name of a weight variable.
#' @param ucb logical; if \code{TRUE}, compute a uniform confidence band (UCB)
#'   for the estimated curve. Valid under an undersmoothing assumption (large
#'   enough knot count so that bias is asymptotically negligible).
#' @param ucb_level coverage level for the UCB. Default is \code{0.95}.
#' @param ucb_sim_reps number of multivariate-normal draws used to approximate
#'   the sup-t distribution. Default is \code{2000}.
#' @param ucb_grid numeric vector of evaluation points for the UCB. If
#'   \code{NULL} (default), a grid of 200 equally spaced points over the range
#'   of \code{uvar} (or \code{atu} if supplied) is used.
#'
#' @return A list of class \code{mfxtsemipar}, including \code{rmse}, the
#'   in-sample RMSE of the final \code{fixest::feols} fit. When \code{ucb = TRUE},
#'   \code{fitted} contains additional \code{<gen>_lb} and \code{<gen>_ub}
#'   columns, and the component \code{ucb} stores the evaluation grid, sup-t
#'   critical value, and simulation details. The reported standard errors and
#'   UCB are justified under undersmoothing; they do not account for
#'   smoothing bias.
#'
#' @export
mfxtsemipar <- function(hf,
                        lf,
                        y,
                        x = NULL,
                        uvar,
                        id,
                        tl,
                        gen,
                        hfcov = NULL,
                        cluster = NULL,
                        bknots = NULL,
                        startp = NULL,
                        endp = NULL,
                        degree = 1L,
                        knots = NULL,
                        type = "poly",
                        winsor = NULL,
                        winsor_values = FALSE,
                        eqspace = FALSE,
                        nknots = NULL,
                        center = NULL,
                        absorb = NULL,
                        atu = NULL,
                        intercept = TRUE,
                        brep = 0L,
                        predy = NULL,
                        weights = NULL,
                        ucb = FALSE,
                        ucb_level = 0.95,
                        ucb_sim_reps = 2000L,
                        ucb_grid = NULL) {

  # ------------------------------------------------------------------
  # 0. package checks
  # ------------------------------------------------------------------
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("Package 'fixest' is required.")
  }
  if (!type %in% c("poly", "bs", "ms", "is", "ibs")) {
    stop("type must be one of 'poly', 'bs', 'ms', 'is', 'ibs'.")
  }
  if (type %in% c("ms", "is", "ibs") &&
      !requireNamespace("splines2", quietly = TRUE)) {
    stop("Package 'splines2' is required for type 'ms', 'is' and 'ibs'.")
  }

  # ------------------------------------------------------------------
  # 1. input validation
  # ------------------------------------------------------------------
  if (is.null(knots) && is.null(nknots)) {
    stop("Either knots or nknots must be specified.")
  }
  if (!is.null(knots) && !is.null(nknots)) {
    stop("knots and nknots cannot be specified at the same time.")
  }

  hf <- data.table::as.data.table(hf)
  lf <- data.table::as.data.table(lf)

  if (length(tl) == 1L && is.null(names(tl))) tl <- as.character(tl)
  keys <- c(id, tl)

  miss_hf <- setdiff(c(keys, uvar, hfcov, atu), names(hf))
  if (length(miss_hf)) {
    stop("Variables not found in hf: ", paste(miss_hf, collapse = ", "))
  }
  miss_lf <- setdiff(c(keys, y, x, cluster, predy), names(lf))
  if (length(miss_lf)) {
    stop("Variables not found in lf: ", paste(miss_lf, collapse = ", "))
  }

  if (!is.null(weights)) {
    if (!(weights %in% names(hf)) && !(weights %in% names(lf))) {
      stop("Weights variable not found in hf or lf: ", weights)
    }
  }

  if (anyDuplicated(lf, by = keys)) {
    stop("lf must contain exactly one row per id*tl combination.")
  }

  if (!is.numeric(hf[[uvar]])) {
    stop("uvar must be numeric.")
  }

  weights_fml <- if (!is.null(weights)) stats::as.formula(paste0("~ ", weights)) else NULL

  # ------------------------------------------------------------------
  # 2. winsorize
  # ------------------------------------------------------------------
  if (!is.null(winsor)) {
    if (length(winsor) != 2L) {
      stop("winsor must be a numeric vector of length 2.")
    }
    if (winsor_values) {
      wlow <- winsor[1L]
      whigh <- winsor[2L]
    } else {
      quants <- stats::quantile(hf[[uvar]], probs = winsor / 100,
                                na.rm = TRUE, names = FALSE)
      wlow <- quants[1L]
      whigh <- quants[2L]
    }
    hf[, (uvar) := pmin(pmax(get(uvar), wlow), whigh)]
  }

  # ------------------------------------------------------------------
  # 3. knots and boundary knots
  # ------------------------------------------------------------------
  uvals <- hf[[uvar]]
  umin_data <- min(uvals, na.rm = TRUE)
  umax_data <- max(uvals, na.rm = TRUE)
  if (!is.finite(umin_data) || !is.finite(umax_data)) {
    stop("No non-missing observations available in uvar.")
  }

  if (is.null(bknots)) {
    bknots <- c(umin_data - 0.01, umax_data + 0.01)
  } else {
    bknots <- as.numeric(bknots)
    if (length(bknots) != 2L) {
      stop("bknots must be a numeric vector of length 2.")
    }
    if (bknots[1L] >= bknots[2L]) {
      stop("bknots must be strictly ascending.")
    }
    if (umin_data < bknots[1L]) {
      stop("Minimum uvar value (", umin_data,
           ") is below the lower boundary knot (", bknots[1L], ").",
           call. = FALSE)
    }
    if (umax_data > bknots[2L]) {
      stop("Maximum uvar value (", umax_data,
           ") is above the upper boundary knot (", bknots[2L], ").",
           call. = FALSE)
    }
  }

  if (!is.null(startp) && !is.numeric(startp)) {
    stop("startp must be numeric.")
  }
  if (!is.null(endp) && !is.numeric(endp)) {
    stop("endp must be numeric.")
  }
  if (!is.null(startp) && !is.null(endp) && startp >= endp) {
    stop("startp must be strictly less than endp.")
  }

  if (is.null(knots)) {
    knots <- gennknots(hf[[uvar]], nknots = nknots,
                       eqspace = eqspace,
                       startp = startp, endp = endp)
  }

  # ------------------------------------------------------------------
  # 4. generate splines and aggregate
  # ------------------------------------------------------------------
  sp <- make_splines(
    x = hf[[uvar]],
    type = type,
    knots = knots,
    bknots = bknots,
    degree = degree,
    center = center,
    intercept = intercept,
    prefix = ".Spline_"
  )
  spline_vars <- sp$names
  spline_cmd <- sp$call
  add_splines_to_hf(hf, sp$matrix, spline_vars)

  lf_agg <- agg_hf_to_lf(
    hf = hf,
    id = id,
    tl = tl,
    spline_vars = spline_vars,
    hfcov = hfcov,
    weights = weights
  )

  hf[, (spline_vars) := NULL]

  lf_est <- merge(lf, lf_agg, by = keys, all.x = TRUE)

  # ------------------------------------------------------------------
  # 5. estimate
  # ------------------------------------------------------------------
  fml <- build_formula(y = y, varlist = c(spline_vars, x, hfcov),
                       partialout = NULL, absorb = absorb)
  cluster_fml <- if (!is.null(cluster)) stats::as.formula(paste0("~ ", cluster)) else NULL

  est <- fixest::feols(fml, data = lf_est, cluster = cluster_fml,
                       weights = weights_fml, warn = FALSE, notes = FALSE)

  # ------------------------------------------------------------------
  # 6. bootstrap inference
  # ------------------------------------------------------------------
  if (brep > 0L) {
    V <- wildboot_vcov_formula(
      data = lf_est,
      fml = fml,
      y = y,
      cluster = cluster,
      weights_fml = weights_fml,
      brep = brep,
      seed = NULL
    )
    b <- stats::coef(est)
    est$vcov <- V
    est$cov.scaled <- V
  } else {
    b <- stats::coef(est)
    V <- stats::vcov(est)
  }

  info <- fixest::fitstat(est, type = c("ll", "aic", "bic", "n"))
  rmse <- sqrt(mean(stats::residuals(est)^2, na.rm = TRUE))

  # ------------------------------------------------------------------
  # 7. prediction
  # ------------------------------------------------------------------
  eval_var <- if (!is.null(atu)) atu else uvar
  eval_x <- hf[[eval_var]]

  sp_eval <- make_splines(
    x = eval_x,
    type = type,
    knots = knots,
    bknots = bknots,
    degree = degree,
    center = center,
    intercept = intercept,
    prefix = ".Spline_"
  )
  B_eval <- sp_eval$matrix

  retained_spline <- intersect(names(b), spline_vars)
  if (length(retained_spline) == 0L) {
    stop("All spline coefficients were dropped; check collinearity.")
  }
  b_spline <- b[retained_spline]
  B_eval_retained <- B_eval[, retained_spline, drop = FALSE]
  fitted_vals <- as.numeric(B_eval_retained %*% b_spline)

  V_spline <- V[retained_spline, retained_spline, drop = FALSE]
  se_vals <- sqrt(pmax(0, rowSums((B_eval_retained %*% V_spline) * B_eval_retained)))

  out_dt <- data.table::copy(hf[, c(keys, if (!is.null(atu)) atu else uvar), with = FALSE])
  out_dt[, (gen) := fitted_vals]
  out_dt[, (paste0(gen, "_se")) := se_vals]

  # ------------------------------------------------------------------
  # 8. predy
  # ------------------------------------------------------------------
  predy_dt <- NULL
  if (!is.null(predy)) {
    lf_est[, (predy) := stats::predict(est, type = "response")]
    predy_dt <- lf_est[, c(keys, predy), with = FALSE]
    out_dt <- merge(out_dt, predy_dt, by = keys, all.x = TRUE)
  }

  # ------------------------------------------------------------------
  # 9. return
  # ------------------------------------------------------------------
  result <- list(
    nknots = length(knots),
    knots = knots,
    bknots = bknots,
    type = type,
    degree = degree,
    center = if (is.null(center)) 0 else center,
    intercept = intercept,
    uvar = uvar,
    startp = startp,
    endp = endp,
    id = id,
    tl = tl,
    coef = b,
    vcov = V,
    info = info,
    rmse = rmse,
    fitted = out_dt,
    predy = predy_dt,
    estimation = est,
    splinecmd = spline_cmd,
    call = match.call()
  )

  class(result) <- c("mfxtsemipar", "list")

  # ------------------------------------------------------------------
  # 10. uniform confidence band (UCB) under undersmoothing
  # ------------------------------------------------------------------
  if (ucb) {
    if (is.null(ucb_grid)) {
      ucb_grid <- seq(min(eval_x, na.rm = TRUE),
                      max(eval_x, na.rm = TRUE),
                      length.out = 200L)
    }
    ucb_res <- compute_ucb_mfxtsemipar(result,
                                       newdata = ucb_grid,
                                       uvar = eval_var,
                                       level = ucb_level,
                                       sim_reps = ucb_sim_reps)
    crit <- ucb_res$crit

    result$fitted[, (paste0(gen, "_lb")) := fitted_vals - crit * se_vals]
    result$fitted[, (paste0(gen, "_ub")) := fitted_vals + crit * se_vals]

    result$ucb <- list(
      grid = ucb_res$grid,
      crit = crit,
      level = ucb_level,
      sim_reps = ucb_sim_reps
    )
  }

  return(result)
}


# ==============================================================================
# Predict and diagnostic methods
# ==============================================================================
#' @export
predict.mfxtsemipar <- function(object, newdata, uvar = NULL, ucb = FALSE, ...) {
  if (!inherits(object, "mfxtsemipar")) {
    stop("object must be of class 'mfxtsemipar'.")
  }

  if (is.numeric(newdata)) {
    x <- newdata
  } else {
    newdata <- data.table::as.data.table(newdata)
    if (is.null(uvar)) uvar <- object$uvar
    if (!(uvar %in% names(newdata))) {
      stop("Variable '", uvar, "' not found in newdata.")
    }
    x <- newdata[[uvar]]
  }

  if (!is.numeric(x)) {
    stop("Evaluation variable must be numeric.")
  }

  sp <- make_splines(
    x = x,
    type = object$type,
    knots = object$knots,
    bknots = object$bknots,
    degree = object$degree,
    center = object$center,
    intercept = object$intercept,
    prefix = ".Spline_"
  )
  B <- sp$matrix
  spline_vars <- sp$names

  retained <- intersect(names(object$coef), spline_vars)
  if (length(retained) == 0L) {
    stop("No spline coefficients available for prediction.")
  }

  b <- object$coef[retained]
  B_retained <- B[, retained, drop = FALSE]
  fit <- as.numeric(B_retained %*% b)

  V <- object$vcov[retained, retained, drop = FALSE]
  se <- sqrt(pmax(0, rowSums((B_retained %*% V) * B_retained)))

  out <- data.table::data.table(u = x, g = fit, se = se)

  if (ucb) {
    if (is.null(object$ucb)) {
      warning("No UCB information found in object; run mfxtsemipar with ucb = TRUE.",
              call. = FALSE)
    } else {
      crit <- object$ucb$crit
      out[, lb := g - crit * se]
      out[, ub := g + crit * se]
    }
  }

  out
}


#' @export
g_diagnostic.mfxtsemipar <- function(object, newdata, true_g, uvar = NULL) {
  pred <- predict.mfxtsemipar(object, newdata = newdata, uvar = uvar)
  compute_g_diagnostic(pred, true_g)
}


#' @export
print.mfxtsemipar <- function(x, ...) {
  cat("Mixed-frequency semiparametric regression with fixed knots\n")
  cat("  Number of knots: ", x$nknots, "\n", sep = "")
  cat("  Knot locations: ", paste(round(x$knots, 4), collapse = ", "), "\n", sep = "")
  cat("  Boundary knots: ", paste(round(x$bknots, 4), collapse = ", "), "\n", sep = "")
  cat("  Spline type: ", x$type, "\n", sep = "")
  cat("  Final fit RMSE: ", round(x$rmse, 6), "\n", sep = "")
  cat("  Note: reported SE/UCB are valid under undersmoothing (sufficiently large K).\n")
  invisible(x)
}


# ==============================================================================
# Helper: compute uniform confidence band for mfxtsemipar (undersmoothing)
# ==============================================================================
compute_ucb_mfxtsemipar <- function(object,
                                    newdata,
                                    uvar = NULL,
                                    level = 0.95,
                                    sim_reps = 2000L) {
  if (!inherits(object, "mfxtsemipar")) {
    stop("object must be of class 'mfxtsemipar'.")
  }

  if (is.numeric(newdata)) {
    x <- newdata
    if (is.null(uvar)) uvar <- object$uvar
  } else {
    newdata <- data.table::as.data.table(newdata)
    if (is.null(uvar)) uvar <- object$uvar
    if (!(uvar %in% names(newdata))) {
      stop("Variable '", uvar, "' not found in newdata.")
    }
    x <- newdata[[uvar]]
  }

  if (!is.numeric(x)) {
    stop("Evaluation variable must be numeric.")
  }

  sp <- make_splines(
    x = x,
    type = object$type,
    knots = object$knots,
    bknots = object$bknots,
    degree = object$degree,
    center = object$center,
    intercept = object$intercept,
    prefix = ".Spline_"
  )
  B <- sp$matrix
  spline_vars <- sp$names

  retained <- intersect(names(object$coef), spline_vars)
  if (length(retained) == 0L) {
    stop("No spline coefficients available for UCB.")
  }
  B_retained <- B[, retained, drop = FALSE]

  b <- object$coef[retained]
  V <- object$vcov[retained, retained, drop = FALSE]

  y <- as.numeric(B_retained %*% b)
  se <- sqrt(pmax(0, rowSums((B_retained %*% V) * B_retained)))

  draws <- simulate_mvnorm(sim_reps, V)
  deviations <- draws %*% t(B_retained)

  se_safe <- pmax(se, .Machine$double.eps)
  t_stats <- apply(abs(deviations) / matrix(se_safe, nrow = sim_reps,
                                            ncol = length(se_safe), byrow = TRUE),
                   1, max)
  crit <- as.numeric(stats::quantile(t_stats, probs = level, na.rm = TRUE))

  list(
    grid = data.table::data.table(
      u = x,
      g = y,
      se = se,
      lb = y - crit * se,
      ub = y + crit * se
    ),
    crit = crit,
    level = level,
    sim_reps = sim_reps
  )
}
