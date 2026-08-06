# Hardcoded absolute path for Rscript compatibility
source("/Users/sigma/SynologyDrive/kuanke/Downloads/jaerevision/mfxtsemipar/R/mfxtsemipar_cknots.R")


# ==============================================================================
# Shape-constrained (U-shaped) B-spline helpers
# ==============================================================================

#' Greville abscissae for a B-spline basis matching splines2::bSpline
bs_greville_sites <- function(knots, bknots, degree, intercept = TRUE) {
  d <- as.integer(degree)
  ext <- c(rep(bknots[1L], d + 1L), sort(as.numeric(knots)),
           rep(bknots[2L], d + 1L))
  nb <- length(ext) - d - 1L
  sites <- vapply(seq_len(nb), function(i) mean(ext[i + seq_len(d)]), 0)
  if (!intercept) sites <- sites[-1L]
  as.numeric(sites)
}


#' Map free level + nonnegative increments to U-shaped B-spline coefficients
#'
#' Coefficients satisfy
#' \code{beta[1] >= ... >= beta[s] <= ... <= beta[m]}
#' when \code{delta >= 0}, via \code{beta = alpha * 1 + M %*% delta}.
ushape_coef_map <- function(m, s) {
  m <- as.integer(m)
  s <- as.integer(s)
  if (s < 1L || s > m) stop("Valley index s must be in 1..m.")
  nL <- s - 1L
  nR <- m - s
  M <- matrix(0, nrow = m, ncol = nL + nR)
  if (nL > 0L) {
    for (i in seq_len(nL)) {
      M[seq_len(i), i] <- 1
    }
  }
  if (nR > 0L) {
    for (k in seq_len(nR)) {
      M[(s + k):m, nL + k] <- 1
    }
  }
  M
}


#' Parse an absorb formula / character vector into a list of FE factors
absorb_to_flist <- function(absorb, data) {
  if (is.null(absorb)) return(NULL)
  if (inherits(absorb, "formula")) {
    abs_txt <- paste(deparse(absorb), collapse = "")
    abs_txt <- sub("^\\s*~\\s*", "", abs_txt)
    vars <- trimws(strsplit(abs_txt, "\\+", fixed = FALSE)[[1L]])
    vars <- vars[nzchar(vars)]
  } else {
    vars <- as.character(absorb)
  }
  miss <- setdiff(vars, names(data))
  if (length(miss)) {
    stop("Absorb variables not found in data: ", paste(miss, collapse = ", "))
  }
  lapply(vars, function(v) data[[v]])
}


#' Weighted least squares for free columns after subtracting V %*% delta
ushape_profile_eta <- function(y, U, V, delta, w = NULL) {
  r <- as.numeric(y - V %*% delta)
  if (is.null(U) || ncol(U) == 0L) {
    return(numeric(0))
  }
  if (is.null(w)) {
    fit <- stats::lm.fit(U, r)
  } else {
    sw <- sqrt(pmax(as.numeric(w), 0))
    fit <- stats::lm.fit(U * sw, r * sw)
  }
  as.numeric(fit$coefficients)
}


#' Constrained LS for U-shaped B-spline coefficients (and free covariates)
#'
#' Minimizes \code{sum(w * (y - B %*% beta - X %*% gamma)^2)} subject to
#' B-spline coefficients being nonincreasing to the left of a valley index
#' and nondecreasing to the right (sufficient for a U-shaped B-spline).
fit_ushape_bs <- function(y, B, X = NULL, w = NULL, valley = NULL,
                          center = 0, knots = NULL, bknots = NULL,
                          degree = 3L, intercept = TRUE) {
  y <- as.numeric(y)
  B <- as.matrix(B)
  n <- length(y)
  if (nrow(B) != n) stop("B and y have incompatible dimensions.")
  m <- ncol(B)
  if (m < 1L) stop("B must have at least one column.")

  if (is.null(valley)) {
    if (is.null(knots) || is.null(bknots)) {
      stop("valley, or knots+bknots, must be supplied to locate the U-shape trough.")
    }
    sites <- bs_greville_sites(knots, bknots, degree, intercept = intercept)
    if (length(sites) != m) {
      stop("Greville site count (", length(sites),
           ") does not match B columns (", m, ").")
    }
    valley <- which.min(abs(sites - center))
  }
  valley <- as.integer(valley)

  if (!is.null(X)) {
    X <- as.matrix(X)
    if (nrow(X) != n) stop("X and y have incompatible dimensions.")
    colnames_X <- colnames(X)
  } else {
    colnames_X <- character(0)
  }

  M <- ushape_coef_map(m, valley)
  Bone <- as.numeric(B %*% rep(1, m))
  V <- if (ncol(M) > 0L) B %*% M else matrix(0, n, 0L)
  U <- cbind(Bone, X)
  colnames_U <- c(".ushape_alpha", colnames_X)

  # Drop collinear / near-zero free columns
  if (ncol(U) > 0L) {
    keep_U <- rep(TRUE, ncol(U))
    for (j in seq_len(ncol(U))) {
      if (max(abs(U[, j])) < sqrt(.Machine$double.eps)) keep_U[j] <- FALSE
    }
    U <- U[, keep_U, drop = FALSE]
    colnames_U <- colnames_U[keep_U]
    if (ncol(U) > 0L) {
      qru <- qr(if (is.null(w)) U else U * sqrt(pmax(as.numeric(w), 0)))
      keep_rank <- qru$pivot[seq_len(qru$rank)]
      U <- U[, sort(keep_rank), drop = FALSE]
      colnames_U <- colnames_U[sort(keep_rank)]
    }
  }

  n_delta <- ncol(V)
  if (n_delta == 0L) {
    eta <- ushape_profile_eta(y, U, V, numeric(0), w = w)
    names(eta) <- colnames_U
    alpha <- if (".ushape_alpha" %in% names(eta)) eta[[".ushape_alpha"]] else 0
    beta <- rep(alpha, m)
    gamma <- eta[setdiff(names(eta), ".ushape_alpha")]
    fitted <- as.numeric(B %*% beta)
    if (length(gamma)) fitted <- fitted + as.numeric(X[, names(gamma), drop = FALSE] %*% gamma)
    return(list(
      beta = beta,
      gamma = gamma,
      alpha = alpha,
      delta = numeric(0),
      valley = valley,
      fitted = fitted,
      ssr = sum((y - fitted)^2 * if (is.null(w)) 1 else w)
    ))
  }

  obj <- function(delta) {
    eta <- ushape_profile_eta(y, U, V, delta, w = w)
    r <- as.numeric(y - V %*% delta)
    if (length(eta)) r <- r - as.numeric(U %*% eta)
    if (is.null(w)) sum(r * r) else sum(w * r * r)
  }

  opt <- stats::optim(
    par = rep(0, n_delta),
    fn = obj,
    method = "L-BFGS-B",
    lower = rep(0, n_delta),
    control = list(maxit = 2000L)
  )
  delta <- pmax(opt$par, 0)
  eta <- ushape_profile_eta(y, U, V, delta, w = w)
  names(eta) <- colnames_U
  alpha <- if (".ushape_alpha" %in% names(eta)) eta[[".ushape_alpha"]] else 0
  beta <- as.numeric(alpha + M %*% delta)
  gamma <- eta[setdiff(names(eta), ".ushape_alpha")]

  fitted <- as.numeric(B %*% beta)
  if (length(gamma)) {
    fitted <- fitted + as.numeric(X[, names(gamma), drop = FALSE] %*% gamma)
  }

  list(
    beta = beta,
    gamma = gamma,
    alpha = alpha,
    delta = delta,
    valley = valley,
    fitted = fitted,
    ssr = obj(delta),
    convergence = opt$convergence
  )
}


#' Approximate OLS sandwich vcov on demeaned design (ignores inequality)
ushape_ols_vcov <- function(y, Z, w = NULL, fitted = NULL) {
  Z <- as.matrix(Z)
  y <- as.numeric(y)
  n <- length(y)
  if (is.null(fitted)) fitted <- as.numeric(Z %*% qr.coef(qr(Z), y))
  e <- y - fitted
  if (is.null(w)) {
    XtX <- crossprod(Z)
    s2 <- sum(e * e) / max(1, n - ncol(Z))
  } else {
    sw <- sqrt(pmax(as.numeric(w), 0))
    XtX <- crossprod(Z * sw)
    s2 <- sum(w * e * e) / max(1, n - ncol(Z))
  }
  V <- tryCatch(
    s2 * solve(XtX),
    error = function(e) {
      s2 * MASS_ginv(XtX)
    }
  )
  V
}


# Minimal Moore-Penrose fallback without requiring MASS on the search path
MASS_ginv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  pos <- s$d > max(tol * s$d[1L], 0)
  if (!any(pos)) {
    return(matrix(0, nrow(A), ncol(A)))
  }
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}


#' Mixed-frequency semiparametric regression with center-anchored U-shaped B-splines
#'
#' Extension of \code{mfxtsemipar_cknots} that (i) always uses B-splines,
#' (ii) keeps \code{center} as an interior knot with left/right knots allocated
#' by sample mass, and (iii) imposes a U-shape constraint on the B-spline
#' coefficients: nonincreasing to the left of a valley near \code{center} and
#' nondecreasing to the right. This is a standard sufficient condition for the
#' estimated curve to be monotone nonincreasing left of the trough and
#' monotone nondecreasing to the right.
#'
#' Estimation: aggregate the HF B-spline basis to LF, demean fixed effects
#' (if any) with \code{fixest::demean}, then fit the shape-constrained least
#' squares problem by nonnegative reparameterization and \code{L-BFGS-B}.
#'
#' @inheritParams mfxtsemipar_cknots
#' @param degree B-spline degree; default \code{3}.
#'
#' @return A list of class
#'   \code{c("mfxtsemipar_cknots_shape", "mfxtsemipar_cknots", "mfxtsemipar")}.
#'
#' @export
mfxtsemipar_cknots_shape <- function(hf,
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
                                     degree = 3L,
                                     knots = NULL,
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

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("Package 'fixest' is required.")
  }
  if (!requireNamespace("splines2", quietly = TRUE)) {
    stop("Package 'splines2' is required for B-splines.")
  }

  if (is.null(knots) && is.null(nknots)) {
    stop("Either knots or nknots must be specified.")
  }
  if (!is.null(knots) && !is.null(nknots)) {
    stop("knots and nknots cannot be specified at the same time.")
  }

  type <- "bs"
  degree <- as.integer(degree)

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

  # ------------------------------------------------------------------
  # winsorize
  # ------------------------------------------------------------------
  if (!is.null(winsor)) {
    if (length(winsor) != 2L) {
      stop("winsor should be a numeric vector of length 2.")
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
  # knots (center-anchored)
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
    if (length(bknots) != 2L) stop("bknots must be a numeric vector of length 2.")
    if (bknots[1L] >= bknots[2L]) stop("bknots must be strictly ascending.")
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

  if (!is.null(startp) && !is.numeric(startp)) stop("startp must be numeric.")
  if (!is.null(endp) && !is.numeric(endp)) stop("endp must be numeric.")
  if (!is.null(startp) && !is.null(endp) && startp >= endp) {
    stop("startp must be strictly less than endp.")
  }

  if (is.null(center)) center <- 0

  knot_meta <- list(n_left = NA_integer_, n_right = NA_integer_,
                    share_left = NA_real_, share_right = NA_real_)
  if (is.null(knots)) {
    knots <- gennknots_center(hf[[uvar]], nknots = nknots,
                              center = center,
                              eqspace = eqspace,
                              startp = startp, endp = endp)
    knot_meta$n_left <- attr(knots, "n_left")
    knot_meta$n_right <- attr(knots, "n_right")
    knot_meta$share_left <- attr(knots, "share_left")
    knot_meta$share_right <- attr(knots, "share_right")
    attributes(knots) <- NULL
  } else {
    knots <- ensure_center_knot(knots, center)
    knot_meta$n_left <- sum(knots < center)
    knot_meta$n_right <- sum(knots > center)
  }

  # ------------------------------------------------------------------
  # B-spline basis at HF, aggregate to LF
  # ------------------------------------------------------------------
  sp <- make_splines(
    x = hf[[uvar]],
    type = "bs",
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

  other_vars <- c(x, hfcov)
  design_vars <- c(spline_vars, other_vars)
  miss_des <- setdiff(c(y, design_vars), names(lf_est))
  if (length(miss_des)) {
    stop("Variables missing after merge: ", paste(miss_des, collapse = ", "))
  }

  ok <- stats::complete.cases(lf_est[, c(y, design_vars), with = FALSE])
  if (!all(ok)) {
    lf_est <- lf_est[ok]
  }
  n_lf <- nrow(lf_est)
  if (n_lf < 2L) stop("Not enough complete LF observations for estimation.")

  wvec <- NULL
  if (!is.null(weights)) {
    if (!(weights %in% names(lf_est))) {
      stop("Weights variable not found in LF estimation sample: ", weights)
    }
    wvec <- as.numeric(lf_est[[weights]])
  }

  y_raw <- as.numeric(lf_est[[y]])
  B_raw <- as.matrix(lf_est[, spline_vars, with = FALSE])
  X_raw <- if (length(other_vars)) {
    as.matrix(lf_est[, other_vars, with = FALSE])
  } else {
    NULL
  }

  # ------------------------------------------------------------------
  # fixed-effect demeaning
  # ------------------------------------------------------------------
  flist <- absorb_to_flist(absorb, lf_est)
  if (!is.null(flist)) {
    mat_all <- cbind(y = y_raw, B_raw, X_raw)
    dm <- fixest::demean(mat_all, f = flist, weights = wvec, as.matrix = TRUE)
    y_dm <- as.numeric(dm[, 1L])
    B_dm <- dm[, 1L + seq_len(ncol(B_raw)), drop = FALSE]
    colnames(B_dm) <- spline_vars
    if (!is.null(X_raw)) {
      X_dm <- dm[, 1L + ncol(B_raw) + seq_len(ncol(X_raw)), drop = FALSE]
      colnames(X_dm) <- other_vars
    } else {
      X_dm <- NULL
    }
  } else {
    y_dm <- y_raw
    B_dm <- B_raw
    X_dm <- X_raw
  }

  # ------------------------------------------------------------------
  # shape-constrained fit
  # ------------------------------------------------------------------
  fit <- fit_ushape_bs(
    y = y_dm,
    B = B_dm,
    X = X_dm,
    w = wvec,
    center = center,
    knots = knots,
    bknots = bknots,
    degree = degree,
    intercept = intercept
  )

  b_spline <- fit$beta
  names(b_spline) <- spline_vars
  b_other <- fit$gamma
  b <- c(b_spline, b_other)

  fitted_lf_dm <- fit$fitted
  # Rebuild FE-level fitted values for RMSE / bootstrap residuals:
  # residual on demeaned scale; RMSE uses demeaned residual (same as within R^2).
  resid_dm <- y_dm - fitted_lf_dm
  rmse <- sqrt(mean(resid_dm^2, na.rm = TRUE))

  # Approximate within-model vcov (treat active inequalities as free)
  Z_dm <- cbind(B_dm, X_dm)
  colnames(Z_dm) <- names(b)
  V <- ushape_ols_vcov(y_dm, Z_dm, w = wvec, fitted = fitted_lf_dm)
  rownames(V) <- colnames(V) <- names(b)

  # ------------------------------------------------------------------
  # wild bootstrap of the constrained estimator
  # ------------------------------------------------------------------
  if (brep > 0L) {
    bb <- matrix(NA_real_, nrow = brep, ncol = length(b))
    colnames(bb) <- names(b)
    yhat_dm <- fitted_lf_dm
    ehat <- resid_dm
    for (r in seq_len(brep)) {
      if (is.null(cluster)) {
        radw <- sample(c(-1, 1), size = n_lf, replace = TRUE) * ehat
      } else {
        cl <- lf_est[[cluster]]
        ucl <- unique(cl)
        signs <- sample(c(-1, 1), size = length(ucl), replace = TRUE)
        names(signs) <- as.character(ucl)
        radw <- as.numeric(signs[as.character(cl)]) * ehat
      }
      y_star <- yhat_dm + radw
      fit_r <- fit_ushape_bs(
        y = y_star,
        B = B_dm,
        X = X_dm,
        w = wvec,
        valley = fit$valley,
        center = center,
        knots = knots,
        bknots = bknots,
        degree = degree,
        intercept = intercept
      )
      br <- c(fit_r$beta, fit_r$gamma)
      names(br) <- c(spline_vars, names(fit_r$gamma))
      bb[r, names(b)] <- br[names(b)]
    }
    V <- stats::cov(bb)
    rownames(V) <- colnames(V) <- names(b)
  }

  # LF prediction levels: demeaned fit + FE-implied means via y - resid_dm
  # (recovered outcome fitted on original scale)
  yhat_lf <- y_raw - resid_dm
  if (!is.null(predy)) {
    lf_est[, (predy) := yhat_lf]
  }

  info <- list(
    n = n_lf,
    valley = fit$valley,
    shape = "U: nonincreasing left of trough, nondecreasing right",
    convergence = if (!is.null(fit$convergence)) fit$convergence else 0L
  )

  # ------------------------------------------------------------------
  # HF curve prediction
  # ------------------------------------------------------------------
  eval_var <- if (!is.null(atu)) atu else uvar
  eval_x <- hf[[eval_var]]
  sp_eval <- make_splines(
    x = eval_x,
    type = "bs",
    knots = knots,
    bknots = bknots,
    degree = degree,
    center = center,
    intercept = intercept,
    prefix = ".Spline_"
  )
  B_eval <- sp_eval$matrix
  fitted_vals <- as.numeric(B_eval %*% b_spline)
  V_spline <- V[spline_vars, spline_vars, drop = FALSE]
  se_vals <- sqrt(pmax(0, rowSums((B_eval %*% V_spline) * B_eval)))

  out_dt <- data.table::copy(
    hf[, c(keys, if (!is.null(atu)) atu else uvar), with = FALSE]
  )
  out_dt[, (gen) := fitted_vals]
  out_dt[, (paste0(gen, "_se")) := se_vals]

  predy_dt <- NULL
  if (!is.null(predy)) {
    predy_dt <- lf_est[, c(keys, predy), with = FALSE]
    out_dt <- merge(out_dt, predy_dt, by = keys, all.x = TRUE)
  }

  result <- list(
    nknots = length(knots),
    knots = knots,
    bknots = bknots,
    type = type,
    degree = degree,
    center = center,
    n_left = knot_meta$n_left,
    n_right = knot_meta$n_right,
    share_left = knot_meta$share_left,
    share_right = knot_meta$share_right,
    valley = fit$valley,
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
    estimation = fit,
    splinecmd = spline_cmd,
    call = match.call()
  )
  class(result) <- c("mfxtsemipar_cknots_shape", "mfxtsemipar_cknots",
                     "mfxtsemipar", "list")

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

  result
}


#' @export
print.mfxtsemipar_cknots_shape <- function(x, ...) {
  cat("Mixed-frequency semiparametric regression (U-shaped B-splines)\n")
  cat("  Spline type: bs (shape-constrained)\n")
  cat("  Number of knots: ", x$nknots, "\n", sep = "")
  cat("  Center knot: ", x$center, "\n", sep = "")
  cat("  Valley coef index: ", x$valley, "\n", sep = "")
  if (!is.na(x$n_left) || !is.na(x$n_right)) {
    cat("  Knot split (L/C/R): ", x$n_left, " / 1 / ", x$n_right, sep = "")
    if (!is.na(x$share_left)) {
      cat("  (shares ", round(100 * x$share_left, 1), "% / ",
          round(100 * x$share_right, 1), "%)", sep = "")
    }
    cat("\n")
  }
  cat("  Shape: left of trough monotone nonincreasing;",
      " right monotone nondecreasing\n")
  cat("  Knot locations: ", paste(round(x$knots, 4), collapse = ", "), "\n", sep = "")
  cat("  Boundary knots: ", paste(round(x$bknots, 4), collapse = ", "), "\n", sep = "")
  cat("  Final fit RMSE (within): ", round(x$rmse, 6), "\n", sep = "")
  invisible(x)
}
