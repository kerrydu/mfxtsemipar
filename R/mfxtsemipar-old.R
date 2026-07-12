# Hardcoded absolute path for Rscript compatibility
source("/Users/sigma/SynologyDrive/kuanke/Downloads/jaerevision/mfxtsemipar/R/mfxtsemipar_cv.R")


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
#'
#' @return A list of class \code{mfxtsemipar}, including \code{rmse}, the
#'   in-sample RMSE of the final \code{fixest::feols} fit.
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
    knots <- gennknots_endpoints(hf[[uvar]], nknots = nknots,
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
    V <- wildboot_vcov(
      data = lf_est,
      y = y,
      varlist = c(spline_vars, x, hfcov),
      partialout = NULL,
      absorb = absorb,
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
  return(result)
}


# ==============================================================================
# Predict and diagnostic methods
# ==============================================================================
#' @export
predict.mfxtsemipar <- function(object, newdata, uvar = NULL, ...) {
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

  data.table::data.table(u = x, g = fit, se = se)
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
  invisible(x)
}


# ==============================================================================
# ==============================================================================
# Helper: parse bc_nknots specification
#
#   NULL     -> same knots as main model
#   6        -> 6 interior BC knots (integer)
#   "6"      -> 6 interior BC knots
#   "1.5k"   -> round(1.5 * CV-selected K) interior BC knots
# ==============================================================================
parse_bc_nknots_spec <- function(bc_nknots) {
  if (is.null(bc_nknots)) {
    return(list(abs = NULL, mult = NULL))
  }
  if (length(bc_nknots) != 1L) {
    stop("bc_nknots must be a single value.")
  }

  if (is.character(bc_nknots)) {
    txt <- trimws(bc_nknots)
    if (grepl("k$", txt, ignore.case = TRUE)) {
      mult <- suppressWarnings(as.numeric(sub("k$", "", txt, ignore.case = TRUE)))
      if (!is.finite(mult) || mult <= 0) {
        stop("bc_nknots multiplier must be positive, e.g. '1.5k'.")
      }
      return(list(abs = NULL, mult = mult))
    }
    val <- suppressWarnings(as.numeric(txt))
    if (!is.finite(val) || val <= 0 || abs(val - round(val)) >= sqrt(.Machine$double.eps)) {
      stop("bc_nknots must be a positive integer, or use '1.5k' for a K multiplier.")
    }
    return(list(abs = as.integer(round(val)), mult = NULL))
  }

  if (is.numeric(bc_nknots)) {
    if (!is.finite(bc_nknots) || bc_nknots <= 0) {
      stop("bc_nknots must be positive.")
    }
    if (abs(bc_nknots - round(bc_nknots)) >= sqrt(.Machine$double.eps)) {
      stop("numeric bc_nknots must be a positive integer; use '1.5k' for a K multiplier.")
    }
    return(list(abs = as.integer(round(bc_nknots)), mult = NULL))
  }

  stop("bc_nknots must be NULL, a positive integer, or a string like '1.5k'.")
}

# ==============================================================================
# Helper: generate knots with optional preset endpoints (local to mfxtsemipar)
#
# Mirrors gennknots() in mfxtsemipar_cv3.R: when startp/endp are supplied they
# are included as the first/last interior knots, with the remaining knots placed
# between them. When both are NULL, behavior falls back to interior quantile
# (or equally spaced) knots over the range of x. Defined under a distinct name
# so it does not override gennknots() from mfxtsemipar_cv.R.
# ==============================================================================
gennknots_endpoints <- function(x, nknots, eqspace = FALSE,
                                startp = NULL, endp = NULL) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    stop("No non-missing observations available to form knots.")
  }

  n_fixed <- (!is.null(startp)) + (!is.null(endp))
  n_middle <- as.integer(nknots - n_fixed)
  if (n_middle < 0L) {
    stop("nknots must be at least the number of specified startp/endp knots.")
  }

  lo <- if (!is.null(startp)) startp else min(x)
  hi <- if (!is.null(endp)) endp else max(x)
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    stop("Invalid knot range: startp must be strictly less than endp.")
  }

  if (n_middle == 0L) {
    cuts <- c(
      if (!is.null(startp)) startp,
      if (!is.null(endp)) endp
    )
    return(cuts)
  }

  if (eqspace) {
    step <- (hi - lo) / (n_middle + 1L)
    middle <- seq(lo + step, hi - step, length.out = n_middle)
  } else {
    if (!is.null(startp) && !is.null(endp)) {
      keep <- x >= startp & x < endp
    } else if (!is.null(endp)) {
      keep <- x < endp
    } else if (!is.null(startp)) {
      keep <- x >= startp
    } else {
      keep <- rep(TRUE, length(x))
    }

    if (sum(keep) < n_middle + 1L) {
      stop("Not enough observations inside [startp, endp] to form knots.")
    }
    probs <- seq_len(n_middle) / (n_middle + 1L)
    middle <- as.numeric(stats::quantile(x[keep], probs = probs,
                                         type = 2, names = FALSE))
  }

  cuts <- c(
    if (!is.null(startp)) startp,
    middle,
    if (!is.null(endp)) endp
  )
  return(cuts)
}
