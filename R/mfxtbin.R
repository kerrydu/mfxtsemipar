.this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.src_dir <- if (!is.null(.this_file)) dirname(.this_file) else "R"
source(file.path(.src_dir, "mfxtbin_cv.R"))


#' Mixed-frequency binned semiparametric regression with fixed bins
#'
#' R implementation of Stata's \code{mfxtbin}. Estimates a binned (step-function)
#' regression with user-specified bins (no cross-validation).
#'
#' @param hf data.frame or data.table containing the high-frequency observations.
#' @param lf data.frame or data.table containing one observation per \code{id*tl}.
#' @param y name of the dependent variable (in \code{lf}).
#' @param x character vector of low-frequency covariates (in \code{lf}).
#' @param uvar name of the binning variable (in \code{hf}).
#' @param id name of the panel identifier.
#' @param tl character vector of time-level variables.
#' @param gen name of the generated fitted-values variable.
#' @param hfcov character vector of high-frequency covariates.
#' @param cluster name of the cluster variable.
#' @param nbin integer; number of bins (alternative to \code{cut}).
#' @param eqspace logical; equally spaced cutpoints when using \code{nbin}.
#' @param cut numeric vector of cutpoints (alternative to \code{nbin} / \code{bw}).
#' @param absorb fixed-effects specification.
#' @param bw numeric; bandwidth for binning (alternative to \code{nbin} / \code{cut}).
#' @param dropbin integer bin number or numeric value whose bin to drop.
#' @param atu alternative evaluation variable (in \code{hf}).
#' @param startp numeric; lower bound for cutpoints.
#' @param endp numeric; upper bound for cutpoints.
#' @param predy name for full LF prediction.
#' @param weights name of a weight variable.
#'
#' @return A list of class \code{mfxtbin}, including \code{rmse}, the
#'   in-sample RMSE of the final \code{fixest::feols} fit.
#'
#' @export
mfxtbin <- function(hf,
                    lf,
                    y,
                    x = NULL,
                    uvar,
                    id,
                    tl,
                    gen,
                    hfcov = NULL,
                    cluster = NULL,
                    nbin = NULL,
                    eqspace = FALSE,
                    cut = NULL,
                    absorb = NULL,
                    bw = NULL,
                    dropbin = NULL,
                    atu = NULL,
                    startp = NULL,
                    endp = NULL,
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

  # ------------------------------------------------------------------
  # 1. mutually exclusive bin specification
  # ------------------------------------------------------------------
  n_spec <- sum(!sapply(list(cut, nbin, bw), is.null))
  if (n_spec == 0L) {
    stop("One of cut, nbin or bw must be specified.")
  }
  if (n_spec > 1L) {
    stop("Only one of cut, nbin and bw can be specified at the same time.")
  }

  # ------------------------------------------------------------------
  # 2. input validation
  # ------------------------------------------------------------------
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
  # 3. generate bins
  # ------------------------------------------------------------------
  gb <- genbins(
    x = hf[[uvar]],
    nbin = nbin,
    cut = cut,
    bw = bw,
    eqspace = eqspace,
    startp = startp,
    endp = endp,
    prefix = ".bin_"
  )
  bin_vars <- gb$names
  cutpoints <- gb$cutpoints

  dropped <- drop_bin(bin_vars, dropbin, cutpoints, hf[[uvar]])
  bin_vars <- dropped$bin_vars
  drop_name <- dropped$drop_name

  add_splines_to_hf(hf, gb$matrix, gb$names)
  if (!is.null(drop_name)) hf[, (drop_name) := NULL]

  # ------------------------------------------------------------------
  # 4. aggregate to LF
  # ------------------------------------------------------------------
  lf_agg <- agg_hf_to_lf(
    hf = hf,
    id = id,
    tl = tl,
    spline_vars = bin_vars,
    hfcov = hfcov,
    weights = weights
  )

  hf[, (gb$names) := NULL]

  lf_est <- merge(lf, lf_agg, by = keys, all.x = TRUE)

  # ------------------------------------------------------------------
  # 5. estimate
  # ------------------------------------------------------------------
  fml <- build_formula(y = y, varlist = c(bin_vars, x, hfcov),
                       partialout = NULL, absorb = absorb)
  cluster_fml <- if (!is.null(cluster)) stats::as.formula(paste0("~ ", cluster)) else NULL

  est <- fixest::feols(fml, data = lf_est, cluster = cluster_fml,
                       weights = weights_fml, warn = FALSE, notes = FALSE)

  b <- stats::coef(est)
  V <- stats::vcov(est)
  info <- fixest::fitstat(est, type = c("ll", "aic", "bic", "n"))
  rmse <- sqrt(mean(stats::residuals(est)^2, na.rm = TRUE))

  # ------------------------------------------------------------------
  # 6. prediction
  # ------------------------------------------------------------------
  eval_var <- if (!is.null(atu)) atu else uvar
  eval_x <- hf[[eval_var]]

  if (!is.null(atu)) {
    atu_min <- min(eval_x, na.rm = TRUE)
    atu_max <- max(eval_x, na.rm = TRUE)
    cp_min <- min(cutpoints)
    cp_max <- max(cutpoints)
    if (atu_min >= cp_min) {
      stop("minimum value of atu should be less than the minimum cutpoint")
    }
    if (atu_max <= cp_max) {
      stop("maximum value of atu should be greater than the maximum cutpoint")
    }
  }

  gb_eval <- genbins(
    x = eval_x,
    cut = cutpoints,
    prefix = ".bin_"
  )
  B_eval <- gb_eval$matrix
  all_bin_vars <- gb_eval$names

  if (!is.null(drop_name)) {
    keep_cols <- setdiff(all_bin_vars, drop_name)
    B_eval <- B_eval[, keep_cols, drop = FALSE]
    all_bin_vars <- keep_cols
  }

  retained_bin <- intersect(names(b), all_bin_vars)
  if (length(retained_bin) == 0L) {
    stop("All bin coefficients were dropped; check collinearity.")
  }
  b_bin <- b[retained_bin]
  B_eval_retained <- B_eval[, retained_bin, drop = FALSE]
  fitted_vals <- as.numeric(B_eval_retained %*% b_bin)

  V_bin <- V[retained_bin, retained_bin, drop = FALSE]
  se_vals <- sqrt(pmax(0, rowSums((B_eval_retained %*% V_bin) * B_eval_retained)))

  out_dt <- data.table::copy(hf[, c(keys, if (!is.null(atu)) atu else uvar), with = FALSE])
  out_dt[, (gen) := fitted_vals]
  out_dt[, (paste0(gen, "_se")) := se_vals]

  # ------------------------------------------------------------------
  # 7. predy
  # ------------------------------------------------------------------
  predy_dt <- NULL
  if (!is.null(predy)) {
    lf_est[, (predy) := stats::predict(est, type = "response")]
    predy_dt <- lf_est[, c(keys, predy), with = FALSE]
    out_dt <- merge(out_dt, predy_dt, by = keys, all.x = TRUE)
  }

  # ------------------------------------------------------------------
  # 8. return
  # ------------------------------------------------------------------
  result <- list(
    nbin = length(retained_bin),
    cutpoints = cutpoints,
    uvar = uvar,
    id = id,
    tl = tl,
    coef = b,
    vcov = V,
    info = info,
    rmse = rmse,
    fitted = out_dt,
    predy = predy_dt,
    estimation = est,
    call = match.call()
  )

  class(result) <- c("mfxtbin", "list")
  return(result)
}


# ==============================================================================
# Predict and diagnostic methods
# ==============================================================================
#' @export
predict.mfxtbin <- function(object, newdata, uvar = NULL, ...) {
  if (!inherits(object, "mfxtbin")) {
    stop("object must be of class 'mfxtbin'.")
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

  gb <- genbins(
    x = x,
    cut = object$cutpoints,
    prefix = ".bin_"
  )
  B <- gb$matrix
  bin_vars <- gb$names

  retained <- intersect(names(object$coef), bin_vars)
  if (length(retained) == 0L) {
    stop("No bin coefficients available for prediction.")
  }

  b <- object$coef[retained]
  B_retained <- B[, retained, drop = FALSE]
  fit <- as.numeric(B_retained %*% b)

  V <- object$vcov[retained, retained, drop = FALSE]
  se <- sqrt(pmax(0, rowSums((B_retained %*% V) * B_retained)))

  data.table::data.table(u = x, g = fit, se = se)
}


#' @export
g_diagnostic.mfxtbin <- function(object, newdata, true_g, uvar = NULL) {
  pred <- predict.mfxtbin(object, newdata = newdata, uvar = uvar)
  compute_g_diagnostic(pred, true_g)
}


#' @export
print.mfxtbin <- function(x, ...) {
  cat("Mixed-frequency binned semiparametric regression with fixed bins\n")
  cat("  Number of bins: ", x$nbin, "\n", sep = "")
  cat("  Cutpoints: ", paste(round(x$cutpoints, 4), collapse = ", "), "\n", sep = "")
  cat("  Final fit RMSE: ", round(x$rmse, 6), "\n", sep = "")
  invisible(x)
}
