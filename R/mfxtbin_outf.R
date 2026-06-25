# Source shared helpers from mfxtbin_cv.R (which sources mfxtsemipar_cv.R).
.this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.src_dir <- if (!is.null(.this_file)) dirname(.this_file) else "R"
source(file.path(.src_dir, "mfxtbin_cv.R"))


#' Out-of-sample prediction for mixed-frequency binned semiparametric regression
#'
#' R implementation of Stata's \code{mfxtbin_outf}. The model is fit on the
#' in-sample observations and predicted values are produced for both in-sample
#' and out-of-sample observations. Results are saved to a file.
#'
#' @param hf data.frame or data.table containing the high-frequency observations.
#' @param lf data.frame or data.table containing one observation per \code{id*tl}.
#' @param y name of the dependent variable (in \code{lf}).
#' @param x character vector of low-frequency covariates (in \code{lf}).
#' @param uvar name of the binning variable (in \code{hf}).
#' @param id name of the panel identifier.
#' @param tl character vector of time-level variables.
#' @param insample name of an indicator variable (in \code{lf}) that equals
#'   \code{1} for in-sample observations and \code{0} for out-of-sample
#'   observations.
#' @param save path to the output file. If the extension is \code{.csv}, a CSV
#'   file is written; otherwise an RDS file is written.
#' @param replace logical; if \code{TRUE}, overwrite an existing file.
#' @param maxnbin maximum number of bins considered in CV.
#' @param minnbin minimum number of bins considered in CV.
#' @param nbin integer; number of bins (alternative to \code{cut} / \code{bw}).
#' @param cut numeric vector of cutpoints (alternative to \code{nbin} / \code{bw}).
#' @param bw numeric; bandwidth for binning (alternative to \code{nbin} /
#'   \code{cut}).
#' @param eqspace logical; equally spaced cutpoints when using \code{nbin}.
#' @param dropbin integer bin number or numeric value whose bin to drop.
#' @param startp numeric; lower bound for cutpoints.
#' @param endp numeric; upper bound for cutpoints.
#' @param hfcov character vector of high-frequency covariates.
#' @param cluster name of the cluster variable for robust standard errors.
#' @param cvgroup name of a variable that defines CV folds (within the
#'   in-sample observations).
#' @param nfold number of CV folds when \code{cvgroup} is \code{NULL}.
#' @param seed integer random seed for fold generation.
#' @param atu name of an alternative evaluation variable. Not used by
#'   \code{mfxtbin_outf} but accepted for compatibility.
#' @param absorb fixed-effects specification passed to \code{fixest::feols}.
#' @param sopt logical; if \code{TRUE}, select the first local minimum of the
#'   CV curve.
#' @param predy name for the full predicted-value column in the output file.
#'   If \code{NULL}, \code{pred_<y>} is used.
#' @param partialout \code{NULL}, \code{"all"}, or a character vector of
#'   variables to partial out.
#' @param weights name of a weight variable (in \code{hf} or \code{lf}).
#'
#' @details
#' If at least one of \code{nbin}, \code{cut} or \code{bw} is supplied, the
#' specified bins are used (no CV). Otherwise the number of bins is selected
#' by cross-validation on the in-sample observations only.
#'
#' The saved dataset contains \code{id}, \code{tl}, \code{insample}, \code{y},
#' \code{predy}, \code{_M_<y>} and \code{_M_ghat}.
#'
#' @return Invisibly returns the saved data.table.
#'
#' @export
mfxtbin_outf <- function(hf,
                         lf,
                         y,
                         x = NULL,
                         uvar,
                         id,
                         tl,
                         insample,
                         save,
                         replace = FALSE,
                         maxnbin = 5,
                         minnbin = 2,
                         nbin = NULL,
                         cut = NULL,
                         bw = NULL,
                         eqspace = FALSE,
                         dropbin = NULL,
                         startp = NULL,
                         endp = NULL,
                         hfcov = NULL,
                         cluster = NULL,
                         cvgroup = NULL,
                         nfold = 10L,
                         seed = NULL,
                         atu = NULL,
                         absorb = NULL,
                         sopt = FALSE,
                         predy = NULL,
                         partialout = NULL,
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
  # 1. basic input validation
  # ------------------------------------------------------------------
  hf <- data.table::as.data.table(hf)
  lf <- data.table::as.data.table(lf)

  if (length(tl) == 1L && is.null(names(tl))) tl <- as.character(tl)
  keys <- c(id, tl)

  miss_hf <- setdiff(c(keys, uvar, hfcov), names(hf))
  if (length(miss_hf)) {
    stop("Variables not found in hf: ", paste(miss_hf, collapse = ", "))
  }
  miss_lf <- setdiff(c(keys, y, insample, x, cluster, cvgroup), names(lf))
  if (length(miss_lf)) {
    stop("Variables not found in lf: ", paste(miss_lf, collapse = ", "))
  }

  if (!is.numeric(lf[[insample]])) {
    stop("insample must be numeric (0/1).")
  }
  if (!all(lf[[insample]] %in% c(0, 1))) {
    stop("insample must contain only 0 and 1.")
  }
  if (!any(lf[[insample]] == 1)) {
    stop("insample must contain at least one in-sample observation.")
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

  if (is.null(predy)) predy <- paste0("pred_", y)

  weights_fml <- if (!is.null(weights)) stats::as.formula(paste0("~ ", weights)) else NULL

  # ------------------------------------------------------------------
  # 2. split controls into main / partialout
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
  # 3. determine bins
  # ------------------------------------------------------------------
  fixed_bins <- !is.null(nbin) || !is.null(cut) || !is.null(bw)

  if (fixed_bins) {
    n_spec <- sum(!sapply(list(cut, nbin, bw), is.null))
    if (n_spec > 1L) {
      stop("Only one of cut, nbin and bw can be specified at the same time.")
    }

    gb <- genbins(
      x = hf[[uvar]], nbin = nbin, cut = cut, bw = bw,
      eqspace = eqspace, startp = startp, endp = endp, prefix = ".bin_"
    )
    cutpoints <- gb$cutpoints
    bin_vars <- gb$names

    dropped <- drop_bin(bin_vars, dropbin, cutpoints, hf[[uvar]])
    bin_vars <- dropped$bin_vars
    drop_name <- dropped$drop_name
  } else {
    # CV selection on in-sample data only
    if (!is.null(cvgroup)) {
      lf[, .cv := as.integer(as.factor(get(cvgroup)))]
    } else {
      if (!is.null(seed)) set.seed(seed)
      in_rows_lf <- which(lf[[insample]] == 1)
      ids_in <- unique(lf[[id]][in_rows_lf])
      folds <- sample(rep_len(seq_len(nfold), length(ids_in)))
      map <- data.table::data.table(.id = ids_in, .cv = folds)
      data.table::setnames(map, ".id", id)
      lf <- merge(lf, map, by = id, all.x = TRUE)
    }

    mse_vec <- rep(NA_real_, maxnbin - minnbin + 1L)
    names(mse_vec) <- paste0("nbin=", seq(minnbin, maxnbin))

    hf <- merge(hf, lf[, c(keys, insample), with = FALSE], by = keys, all.x = TRUE)
    in_rows_hf <- which(hf[[insample]] == 1)
    hf_in <- hf[in_rows_hf]

    for (i in seq_along(mse_vec)) {
      nbin_i <- minnbin + i - 1L
      gb <- genbins(
        x = hf_in[[uvar]], nbin = nbin_i, eqspace = eqspace,
        startp = startp, endp = endp, prefix = ".bin_"
      )
      all_bin_vars <- gb$names
      cutpoints_i <- gb$cutpoints

      dropped <- drop_bin(all_bin_vars, dropbin, cutpoints_i, hf_in[[uvar]])
      bin_vars <- dropped$bin_vars
      drop_name_i <- dropped$drop_name

      add_splines_to_hf(hf_in, gb$matrix, all_bin_vars)
      if (!is.null(drop_name_i)) hf_in[, (drop_name_i) := NULL]

      agg_in <- agg_hf_to_lf(
        hf = hf_in, id = id, tl = tl, spline_vars = bin_vars,
        hfcov = hfcov, weights = weights
      )
      in_rows_lf <- which(lf[[insample]] == 1)
      lf_in_cv <- merge(lf[in_rows_lf], agg_in, by = keys, all.x = TRUE)

      mse_vec[i] <- rmse_cv(
        data = lf_in_cv,
        y = y,
        varlist = c(bin_vars, x_main),
        partialout = x_partial,
        absorb = absorb,
        cluster = cluster,
        weights_fml = weights_fml,
        cvvar = ".cv"
      )

      hf_in[, (all_bin_vars) := NULL]
    }

    minpos <- which.min(mse_vec)
    soptbin <- minnbin + minpos - 1L

    if (sopt) {
      soptbin <- minnbin
      for (k in seq(minnbin, maxnbin)) {
        idx <- k - minnbin + 1L
        if (!is.na(mse_vec[idx]) &&
            (idx == 1L || mse_vec[idx] < mse_vec[idx - 1L])) {
          soptbin <- k
        } else {
          break
        }
      }
    }

    nbin <- if (sopt) soptbin else (minnbin + minpos - 1L)

    gb <- genbins(
      x = hf[[uvar]], nbin = nbin, eqspace = eqspace,
      startp = startp, endp = endp, prefix = ".bin_"
    )
    cutpoints <- gb$cutpoints
    bin_vars <- gb$names

    dropped <- drop_bin(bin_vars, dropbin, cutpoints, hf[[uvar]])
    bin_vars <- dropped$bin_vars
    drop_name <- dropped$drop_name

    cat("\nCross-validation MSE (in-sample)\n")
    print(data.table::data.table(nbin = seq(minnbin, maxnbin), MSE = mse_vec))
    cat("\n")
  }

  # ------------------------------------------------------------------
  # 4. fit and predict in/out of sample
  # ------------------------------------------------------------------
  result <- bin_outf_work(
    hf = hf, lf = lf, y = y, x_main = x_main, x_partial = x_partial,
    uvar = uvar, id = id, tl = tl, insample = insample,
    save = save, replace = replace,
    cutpoints = cutpoints, drop_name = drop_name,
    hfcov = hfcov, cluster = cluster, absorb = absorb,
    weights_fml = weights_fml, weights = weights, predy = predy
  )

  return(invisible(result))
}


# ==============================================================================
# Helper: in/out-of-sample binned prediction and file saving
# ==============================================================================
bin_outf_work <- function(hf, lf, y, x_main, x_partial, uvar, id, tl,
                          insample, save, replace,
                          cutpoints, drop_name,
                          hfcov, cluster, absorb, weights_fml,
                          weights = NULL, predy) {

  keys <- c(id, tl)

  # generate bins on full hf
  gb <- genbins(x = hf[[uvar]], cut = cutpoints, prefix = ".bin_")
  bin_vars <- gb$names
  add_splines_to_hf(hf, gb$matrix, bin_vars)
  if (!is.null(drop_name)) {
    bin_vars <- setdiff(bin_vars, drop_name)
    hf[, (drop_name) := NULL]
  }

  # merge insample indicator to hf
  if (!(insample %in% names(hf))) {
    hf <- merge(hf, lf[, c(keys, insample), with = FALSE], by = keys, all.x = TRUE)
  }

  # aggregate separately
  in_rows_hf <- which(hf[[insample]] == 1)
  out_rows_hf <- which(hf[[insample]] == 0)
  hf_in <- hf[in_rows_hf]
  hf_out <- hf[out_rows_hf]

  agg_in <- agg_hf_to_lf(
    hf = hf_in, id = id, tl = tl, spline_vars = bin_vars,
    hfcov = hfcov, weights = weights
  )
  agg_out <- agg_hf_to_lf(
    hf = hf_out, id = id, tl = tl, spline_vars = bin_vars,
    hfcov = hfcov, weights = weights
  )

  in_rows_lf <- which(lf[[insample]] == 1)
  out_rows_lf <- which(lf[[insample]] == 0)
  lf_in <- merge(lf[in_rows_lf], agg_in, by = keys, all.x = TRUE)
  lf_out <- merge(lf[out_rows_lf], agg_out, by = keys, all.x = TRUE)

  # clean hf
  hf[, (bin_vars) := NULL]
  if (!is.null(drop_name)) hf[, (drop_name) := NULL]

  cluster_fml <- if (!is.null(cluster)) stats::as.formula(paste0("~ ", cluster)) else NULL
  main_varlist <- c(bin_vars, x_main)

  # ------------------------------------------------------------------
  # fit main model on in-sample
  # ------------------------------------------------------------------
  fml_main <- build_formula(
    y = y, varlist = main_varlist, partialout = x_partial, absorb = absorb
  )
  fit_in <- fixest::feols(
    fml_main, data = lf_in, cluster = cluster_fml,
    weights = weights_fml, warn = FALSE, notes = FALSE
  )

  b <- stats::coef(fit_in)
  retained <- intersect(names(b), main_varlist)
  if (length(retained) == 0L) {
    stop("All main-regression variables were dropped; check collinearity.")
  }

  # ------------------------------------------------------------------
  # in-sample predictions
  # ------------------------------------------------------------------
  predy_in <- as.numeric(stats::predict(fit_in, type = "response"))

  X_in <- as.matrix(lf_in[, retained, with = FALSE])
  ghat_in <- as.numeric(X_in %*% b[retained])

  fml_M <- fml_residualize(y = y, x_partial = x_partial, absorb = absorb)
  fml_Mg <- fml_residualize(y = ".ghat", x_partial = x_partial, absorb = absorb)

  fit_My_in <- fixest::feols(
    fml_M, data = lf_in, cluster = cluster_fml,
    weights = weights_fml, warn = FALSE, notes = FALSE
  )
  M_y_in <- as.numeric(stats::residuals(fit_My_in))

  lf_in[, .ghat := ghat_in]
  fit_Mg_in <- fixest::feols(
    fml_Mg, data = lf_in, cluster = cluster_fml,
    weights = weights_fml, warn = FALSE, notes = FALSE
  )
  M_ghat_in <- as.numeric(stats::residuals(fit_Mg_in))

  # ------------------------------------------------------------------
  # out-of-sample predictions
  # ------------------------------------------------------------------
  X_out <- as.matrix(lf_out[, retained, with = FALSE])
  ghat_out <- as.numeric(X_out %*% b[retained])

  lf_out[, .ghat := ghat_out]

  fit_My_out <- fixest::feols(
    fml_M, data = lf_out, cluster = cluster_fml,
    weights = weights_fml, warn = FALSE, notes = FALSE
  )
  M_y_out <- as.numeric(stats::residuals(fit_My_out))
  predy_out <- as.numeric(stats::predict(fit_My_out, type = "response"))

  fit_Mg_out <- fixest::feols(
    fml_Mg, data = lf_out, cluster = cluster_fml,
    weights = weights_fml, warn = FALSE, notes = FALSE
  )
  M_ghat_out <- as.numeric(stats::residuals(fit_Mg_out))

  predy_out <- predy_out + M_ghat_out

  # ------------------------------------------------------------------
  # combine and save
  # ------------------------------------------------------------------
  out_cols <- c(keys, insample, y)
  out_dt <- data.table::rbindlist(
    list(lf_in[, out_cols, with = FALSE], lf_out[, out_cols, with = FALSE]),
    use.names = TRUE
  )

  M_y_name <- paste0("_M_", y)
  M_ghat_name <- "_M_ghat"

  out_dt[, (predy) := c(predy_in, predy_out)]
  out_dt[, (M_y_name) := c(M_y_in, M_y_out)]
  out_dt[, (M_ghat_name) := c(M_ghat_in, M_ghat_out)]

  outf_save(out_dt, save, replace)

  attr(out_dt, "predy") <- predy
  attr(out_dt, "M_y") <- M_y_name
  attr(out_dt, "M_ghat") <- M_ghat_name

  message("Saved out-of-sample predictions to: ", save)
  return(invisible(out_dt))
}
