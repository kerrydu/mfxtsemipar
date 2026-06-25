#' Mixed-frequency cross-validated semiparametric regression: average-then-model
#'
#' Counterpart to \code{mfxtsemipar_cv}. Instead of evaluating the
#' nonparametric component at high frequency and summing the spline basis to the
#' low frequency (\code{sum(g(x))}), this function first averages the
#' high-frequency variable \code{uvar} within each \code{id*tl} cell and then
#' models the semiparametric component directly at the low frequency
#' (\code{g(mean(x))}). It is intended for a comparison/bias check requested by
#' a referee.
#'
#' @inheritParams mfxtsemipar_cv
#'
#' @details
#' All high-frequency variables listed in \code{uvar}, \code{hfcov} and
#' \code{atu} (if supplied) are averaged to the low-frequency level by taking
#' means. The spline basis is then generated on the averaged \code{uvar} and
#' estimation is performed on the low-frequency data set. The returned fitted
#' values live at the LF level (one per \code{id*tl}).
#'
#' Because \code{mfxtsemipar_cv} models the low-frequency dependent variable as
#' the sum of high-frequency contributions, the LF dependent variable \code{y}
#' is divided by the number of high-frequency observations in each
#' \code{id*tl} cell before estimation. This puts \code{y} and the averaged
#' \code{uvar} on the same per-observation scale, so that the estimated
#' semiparametric component can be interpreted as \code{g(mean(u))}.
#'
#' To keep the spline domain comparable to \code{mfxtsemipar_cv}, the boundary
#' knots are taken from the high-frequency \code{uvar} distribution. The interior
#' knot locations are based on the averaged low-frequency \code{uvar}, and the
#' spline basis is evaluated at the averaged \code{uvar}.
#'
#' @return A list of class \code{mfxtsemipar_cv_avg} with the same components as
#'   \code{mfxtsemipar_cv}, including \code{cv_mse}, \code{min_cv_mse}, and
#'   \code{rmse}. The \code{fitted} data.table contains the low-frequency
#'   evaluation point (named \code{uvar}), the fitted semiparametric component,
#'   and \code{.hf_n}, the number of high-frequency observations per
#'   \code{id*tl} cell.
#'
#' @export
mfxtsemipar_cv_avg <- function(hf,
                               lf,
                               y,
                               x = NULL,
                               uvar,
                               id,
                               tl,
                               gen,
                               hfcov = NULL,
                               cluster = NULL,
                               type = "poly",
                               winsor = NULL,
                               winsor_values = FALSE,
                               eqspace = FALSE,
                               maxnk = 5L,
                               minnk = 2L,
                               center = NULL,
                               absorb = NULL,
                               partialout = NULL,
                               cvgroup = NULL,
                               nfold = 10L,
                               seed = NULL,
                               degree = 1L,
                               keepsplines = FALSE,
                               atu = NULL,
                               dropfirstbase = FALSE,
                               sopt = FALSE,
                               brep = 0L,
                               predy = NULL,
                               weights = NULL) {

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
  # 1. basic input validation and conversion
  # ------------------------------------------------------------------
  hf <- data.table::as.data.table(hf)
  lf <- data.table::as.data.table(lf)

  if (length(tl) == 1L && is.null(names(tl))) tl <- as.character(tl)
  keys <- c(id, tl)

  miss_hf <- setdiff(c(keys, uvar, hfcov, atu), names(hf))
  if (length(miss_hf)) {
    stop("Variables not found in hf: ", paste(miss_hf, collapse = ", "))
  }
  miss_lf <- setdiff(c(keys, y, x, cluster, cvgroup), names(lf))
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

  mean_uvar <- ".u_mean"
  mean_atu  <- ".atu_mean"

  # ------------------------------------------------------------------
  # 2. winsorize uvar if requested (on hf, then average)
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
  # 3. boundary knots from high-frequency uvar (keep same domain as mfxtsemipar_cv)
  # ------------------------------------------------------------------
  uvals <- hf[[uvar]]
  umin <- min(uvals, na.rm = TRUE) - 0.01
  umax <- max(uvals, na.rm = TRUE) + 0.01
  bknots <- c(umin, umax)

  # ------------------------------------------------------------------
  # 4. average hf variables to lf
  # ------------------------------------------------------------------
  mean_vars <- unique(c(uvar, hfcov, if (!is.null(atu)) atu))
  lf_mean <- mean_hf_to_lf(hf = hf, id = id, tl = tl,
                           vars = mean_vars, weights = weights)

  # rename averaged uvar / atu to internal names
  if (uvar %in% names(lf_mean)) {
    data.table::setnames(lf_mean, uvar, mean_uvar)
  }
  if (!is.null(atu) && atu != uvar && atu %in% names(lf_mean)) {
    data.table::setnames(lf_mean, atu, mean_atu)
  }

  lf_est <- merge(lf, lf_mean, by = keys, all.x = TRUE)

  if (any(is.na(lf_est[[mean_uvar]]))) {
    stop("Some id*tl combinations in lf have no matching high-frequency observations.")
  }

  # Put LF y on the same per-HF-observation scale as the averaged uvar.
  # mfxtsemipar_cv models y as a sum of HF contributions; the average-then-model
  # counterpart estimates g(mean(u)), so y must be divided by the HF count.
  lf_est[, (y) := .SD[[1L]] / .hf_n, .SDcols = y]

  # ------------------------------------------------------------------
  # 5. generate CV groups
  # ------------------------------------------------------------------
  if (!is.null(cvgroup)) {
    lf_est[, .cv := as.integer(as.factor(get(cvgroup)))]
  } else {
    if (!is.null(seed)) set.seed(seed)
    ids <- unique(lf_est[[id]])
    folds <- sample(rep_len(seq_len(nfold), length(ids)))
    map <- data.table::data.table(.id = ids, .cv = folds)
    data.table::setnames(map, ".id", id)
    lf_est <- merge(lf_est, map, by = id, all.x = TRUE)
  }

  # ------------------------------------------------------------------
  # 6. split controls into main / partialout
  # ------------------------------------------------------------------
  all_controls <- unique(c(x, hfcov))
  if (is.null(partialout)) {
    x_main <- all_controls
    x_partial <- NULL
  } else if (identical(partialout, "all")) {
    x_main <- NULL
    x_partial <- all_controls
  } else {
    partialout <- as.character(partialout)
    x_partial <- intersect(partialout, all_controls)
    if (length(setdiff(partialout, all_controls))) {
      warning("Some partialout variables are not in x or hfcov: ",
              paste(setdiff(partialout, all_controls), collapse = ", "))
    }
    x_main <- setdiff(all_controls, x_partial)
  }

  # ------------------------------------------------------------------
  # 7. cross-validation over nknots
  # ------------------------------------------------------------------
  cv_rmse_vec <- rep(NA_real_, maxnk - minnk + 1L)
  names(cv_rmse_vec) <- paste0("nk=", seq(minnk, maxnk))

  for (i in seq_along(cv_rmse_vec)) {
    nk <- minnk + i - 1L

    knots <- gennknots(lf_est[[mean_uvar]], nknots = nk, eqspace = eqspace,
                       startp = bknots[1L], endp = bknots[2L])

    sp <- make_splines(
      x = lf_est[[mean_uvar]],
      type = type,
      knots = knots,
      bknots = bknots,
      degree = degree,
      center = center,
      intercept = !dropfirstbase,
      prefix = ".Spline_"
    )
    spline_vars <- sp$names
    add_splines_to_hf(lf_est, sp$matrix, spline_vars)

    cv_rmse_vec[i] <- rmse_cv(
      data = lf_est,
      y = y,
      varlist = c(spline_vars, x_main),
      partialout = x_partial,
      absorb = absorb,
      cluster = cluster,
      weights_fml = weights_fml,
      cvvar = ".cv"
    )

    lf_est[, (spline_vars) := NULL]
  }

  # select optimal knots
  minpos <- which.min(cv_rmse_vec)
  soptnk <- minnk + minpos - 1L

  if (sopt) {
    soptnk <- minnk
    for (k in seq(minnk, maxnk)) {
      idx <- k - minnk + 1L
      if (!is.na(cv_rmse_vec[idx]) &&
          (idx == 1L || cv_rmse_vec[idx] < cv_rmse_vec[idx - 1L])) {
        soptnk <- k
      } else {
        break
      }
    }
  }

  nknots <- if (sopt) soptnk else (minnk + minpos - 1L)
  min_cv_mse <- cv_rmse_vec[minpos]

  cat("\nCross-validation RMSE (average-then-model, for knot selection)\n")
  print(data.table::data.table(nk = seq(minnk, maxnk), cv_rmse = cv_rmse_vec))
  cat("\n")

  # ------------------------------------------------------------------
  # 8. final estimation with optimal knots
  # ------------------------------------------------------------------
  knots <- gennknots(lf_est[[mean_uvar]], nknots = nknots, eqspace = eqspace,
                     startp = bknots[1L], endp = bknots[2L])

  sp <- make_splines(
    x = lf_est[[mean_uvar]],
    type = type,
    knots = knots,
    bknots = bknots,
    degree = degree,
    center = center,
    intercept = !dropfirstbase,
    prefix = ".Spline_"
  )
  spline_vars <- sp$names
  spline_cmd <- sp$call
  add_splines_to_hf(lf_est, sp$matrix, spline_vars)

  fml <- build_formula(y = y, varlist = c(spline_vars, x_main),
                       partialout = x_partial, absorb = absorb)
  cluster_fml <- if (!is.null(cluster)) stats::as.formula(paste0("~ ", cluster)) else NULL

  est <- fixest::feols(fml, data = lf_est, cluster = cluster_fml,
                       weights = weights_fml, warn = FALSE, notes = FALSE)

  # ------------------------------------------------------------------
  # 9. bootstrap inference if requested
  # ------------------------------------------------------------------
  if (brep > 0L) {
    V <- wildboot_vcov(
      data = lf_est,
      y = y,
      varlist = c(spline_vars, x_main),
      partialout = x_partial,
      absorb = absorb,
      cluster = cluster,
      weights_fml = weights_fml,
      brep = brep,
      seed = seed
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
  # 10. prediction: semiparametric component at LF
  # ------------------------------------------------------------------
  eval_col <- if (!is.null(atu) && mean_atu %in% names(lf_est)) mean_atu else mean_uvar
  eval_x <- lf_est[[eval_col]]

  sp_eval <- make_splines(
    x = eval_x,
    type = type,
    knots = knots,
    bknots = bknots,
    degree = degree,
    center = center,
    intercept = !dropfirstbase,
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

  out_dt <- data.table::copy(lf_est[, c(keys, eval_col, ".hf_n"), with = FALSE])
  out_dt[, (gen) := fitted_vals]
  out_dt[, (paste0(gen, "_se")) := se_vals]

  # expose the averaged evaluation point under the original uvar name
  data.table::setnames(out_dt, eval_col, uvar)

  if (keepsplines) {
    for (j in seq_along(spline_vars)) {
      out_dt[, (spline_vars[j]) := B_eval[, j]]
    }
  }

  # ------------------------------------------------------------------
  # 11. predy: full prediction at low frequency
  # ------------------------------------------------------------------
  predy_dt <- NULL
  if (!is.null(predy)) {
    lf_est[, (predy) := stats::predict(est, type = "response")]
    predy_dt <- lf_est[, c(keys, predy), with = FALSE]
    out_dt <- merge(out_dt, predy_dt, by = keys, all.x = TRUE)
  }

  # ------------------------------------------------------------------
  # 12. return
  # ------------------------------------------------------------------
  result <- list(
    soptnk = soptnk,
    nknots = nknots,
    knots = knots,
    bknots = bknots,
    min_cv_mse = min_cv_mse,
    cv_mse = data.table::data.table(nk = seq(minnk, maxnk), cv_rmse = cv_rmse_vec),
    rmse = rmse,
    splinecmd = spline_cmd,
    type = type,
    degree = degree,
    center = if (is.null(center)) 0 else center,
    intercept = !dropfirstbase,
    uvar = uvar,
    id = id,
    tl = tl,
    coef = b,
    vcov = V,
    info = info,
    fitted = out_dt,
    predy = predy_dt,
    estimation = est,
    call = match.call()
  )

  class(result) <- c("mfxtsemipar_cv_avg", "list")
  return(result)
}


# ==============================================================================
# Helper: average high-frequency variables to low frequency
# ==============================================================================
mean_hf_to_lf <- function(hf, id, tl, vars, weights = NULL) {
  keys <- c(id, tl)
  out <- hf[, lapply(.SD, mean, na.rm = TRUE), by = keys, .SDcols = vars]

  # number of HF observations per LF cell; needed to put LF y on per-HF scale
  n_dt <- hf[, .(.hf_n = .N), by = keys]
  out <- merge(out, n_dt, by = keys, all.x = TRUE)

  if (!is.null(weights) && weights %in% names(hf)) {
    wmean <- hf[, lapply(.SD, mean, na.rm = TRUE), by = keys, .SDcols = weights]
    out <- merge(out, wmean, by = keys, all.x = TRUE)
  }

  return(out)
}


# ==============================================================================
# Print method
# ==============================================================================
#' @export
print.mfxtsemipar_cv_avg <- function(x, ...) {
  cat("Mixed-frequency semiparametric regression (average-then-model) with cross-validation\n")
  cat("  Selected knots: ", x$nknots, "\n", sep = "")
  cat("  Boundary knots: ", paste(round(x$bknots, 4), collapse = ", "), "\n", sep = "")
  cat("  Knot locations: ", paste(round(x$knots, 4), collapse = ", "), "\n", sep = "")
  cat("  Minimum CV RMSE: ", round(x$min_cv_mse, 6), "\n", sep = "")
  cat("  Final fit RMSE: ", round(x$rmse, 6), "\n", sep = "")
  cat("  Simple optimal knots: ", x$soptnk, "\n", sep = "")
  cat("\nCV RMSE by number of knots:\n")
  print(x$cv_mse, row.names = FALSE)
  invisible(x)
}


# ==============================================================================
# Predict method for g(u)
# ==============================================================================
#' Predict the semiparametric component g(u) from an mfxtsemipar_cv_avg fit
#'
#' @param object an object of class \code{mfxtsemipar_cv_avg}.
#' @param newdata a data.frame/data.table containing the averaged evaluation
#'   variable, or a numeric vector of evaluation points.
#' @param uvar name of the evaluation variable in \code{newdata}. If
#'   \code{NULL}, the \code{uvar} used in estimation is used.
#' @param ... additional arguments (currently ignored).
#'
#' @return A data.table with columns \code{u}, \code{g} and \code{se}.
#'
#' @export
predict.mfxtsemipar_cv_avg <- function(object, newdata, uvar = NULL, ...) {
  if (!inherits(object, "mfxtsemipar_cv_avg")) {
    stop("object must be of class 'mfxtsemipar_cv_avg'.")
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
    stop("No spline coefficients available for prediction.")
  }

  b <- object$coef[retained]
  B_retained <- B[, retained, drop = FALSE]

  fit <- as.numeric(B_retained %*% b)

  V <- object$vcov[retained, retained, drop = FALSE]
  se <- sqrt(pmax(0, rowSums((B_retained %*% V) * B_retained)))

  data.table::data.table(u = x, g = fit, se = se)
}


# ==============================================================================
# Diagnostic: compare estimated g(u) to true g(u)
# ==============================================================================
#' @export
g_diagnostic.mfxtsemipar_cv_avg <- function(object, newdata, true_g, uvar = NULL) {
  pred <- predict.mfxtsemipar_cv_avg(object, newdata = newdata, uvar = uvar)
  compute_g_diagnostic(pred, true_g)
}
