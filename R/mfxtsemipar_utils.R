# ==============================================================================
# Shared utilities for mfxtsemipar family
#
# This file contains helper functions used by:
#   - mfxtsemipar.R
#   - mfxtsemipar_bc.R
#   - mfxtsemipar2_cv.R
#
# Do not modify mfxtsemipar_cv.R (legacy file).
# ==============================================================================

if (!exists("parse_bc_nknots_spec", mode = "function", inherits = FALSE)) {
  parse_bc_nknots_spec <- function(bc_nknots) {
    if (is.null(bc_nknots)) return(list(abs = NULL, mult = NULL))
    if (length(bc_nknots) != 1L) stop("bc_nknots must be a single value.")
    if (is.character(bc_nknots)) {
      txt <- trimws(bc_nknots)
      if (grepl("k$", txt, ignore.case = TRUE)) {
        mult <- suppressWarnings(as.numeric(sub("k$", "", txt, ignore.case = TRUE)))
        if (!is.finite(mult) || mult <= 0) stop("bc_nknots multiplier must be positive, e.g. '1.5k'.")
        return(list(abs = NULL, mult = mult))
      }
      val <- suppressWarnings(as.numeric(txt))
      if (!is.finite(val) || val <= 0 || abs(val - round(val)) >= sqrt(.Machine$double.eps)) {
        stop("bc_nknots must be a positive integer, or use '1.5k' for a K multiplier.")
      }
      return(list(abs = as.integer(round(val)), mult = NULL))
    }
    if (is.numeric(bc_nknots)) {
      if (!is.finite(bc_nknots) || bc_nknots <= 0) stop("bc_nknots must be positive.")
      if (abs(bc_nknots - round(bc_nknots)) >= sqrt(.Machine$double.eps)) {
        stop("numeric bc_nknots must be a positive integer; use '1.5k' for a K multiplier.")
      }
      return(list(abs = as.integer(round(bc_nknots)), mult = NULL))
    }
    stop("bc_nknots must be NULL, a positive integer, or a string like '1.5k'.")
  }
}

gennknots <- function(x, nknots, eqspace = FALSE,
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


# ==============================================================================
# Helper: build spline basis
# ==============================================================================
make_splines <- function(x, type, knots, bknots, degree = 1L,
                         center = NULL, intercept = TRUE,
                         prefix = ".Spline_") {
  if (is.null(center)) center <- 0

  out <- switch(type,
    poly = polysplines(x, knots, degree, center, intercept, prefix),
    bs   = bs_splines(x, knots, bknots, degree, center, intercept, prefix),
    ms   = ms_splines(x, knots, bknots, degree, intercept, prefix),
    is   = is_splines(x, knots, bknots, degree, intercept, prefix),
    ibs  = ibs_splines(x, knots, bknots, degree, intercept, prefix)
  )

  out$call <- paste0(
    type, " splines: degree=", degree,
    ", knots=c(", paste(round(knots, 4), collapse = ", "),
    "), bknots=c(", paste(round(bknots, 4), collapse = ", "), ")"
  )
  return(out)
}


polysplines <- function(x, knots, degree = 1L, center = 0,
                        intercept = TRUE, prefix = ".Spline_") {
  n <- length(x)
  p <- as.integer(degree)
  j <- 1L
  vars <- character()
  mat <- matrix(0, nrow = n, ncol = 0L)

  for (pow in seq_len(p)) {
    v <- paste0(prefix, j)
    vars <- c(vars, v)
    mat <- cbind(mat, x^pow - center^pow)
    j <- j + 1L
  }

  for (k in knots) {
    v <- paste0(prefix, j)
    vars <- c(vars, v)
    mat <- cbind(mat,
      (x - k)^p * (x > k) - (center - k)^p * (center > k)
    )
    j <- j + 1L
  }

  if (intercept) {
    v0 <- paste0(prefix, "0")
    vars <- c(v0, vars)
    mat <- cbind(1, mat)
  }

  colnames(mat) <- vars
  list(matrix = mat, names = vars)
}


bs_splines <- function(x, knots, bknots, degree = 3L, center = 0,
                       intercept = TRUE, prefix = ".Spline_") {
  # Use splines2::bSpline to match Stata gensplines type(bs)
  B <- splines2::bSpline(x, knots = knots, Boundary.knots = bknots,
                         degree = degree, intercept = intercept)
  B <- as.matrix(B)

  # Apply centering (centerv): B(x) - B(center)
  B_center <- splines2::bSpline(center, knots = knots, Boundary.knots = bknots,
                               degree = degree, intercept = intercept)
  B <- sweep(B, 2, as.numeric(B_center), "-")

  start_idx <- ifelse(intercept, 0L, 1L)
  vars <- paste0(prefix, start_idx:(start_idx + ncol(B) - 1L))
  colnames(B) <- vars
  list(matrix = B, names = vars)
}


ms_splines <- function(x, knots, bknots, degree = 3L,
                       intercept = TRUE, prefix = ".Spline_") {
  B <- splines2::mSpline(x, knots = knots, Boundary.knots = bknots,
                         degree = degree, intercept = intercept)
  vars <- paste0(prefix, seq_len(ncol(B)) - ifelse(intercept, 1L, 0L))
  colnames(B) <- vars
  list(matrix = as.matrix(B), names = vars)
}


is_splines <- function(x, knots, bknots, degree = 2L,
                       intercept = TRUE, prefix = ".Spline_") {
  B <- splines2::iSpline(x, knots = knots, Boundary.knots = bknots,
                         degree = degree, intercept = intercept)
  vars <- paste0(prefix, seq_len(ncol(B)) - ifelse(intercept, 1L, 0L))
  colnames(B) <- vars
  list(matrix = as.matrix(B), names = vars)
}


ibs_splines <- function(x, knots, bknots, degree = 3L,
                        intercept = TRUE, prefix = ".Spline_") {
  B <- splines2::ibs(x, knots = knots, Boundary.knots = bknots,
                     degree = degree, intercept = intercept)
  vars <- paste0(prefix, seq_len(ncol(B)) - ifelse(intercept, 1L, 0L))
  colnames(B) <- vars
  list(matrix = as.matrix(B), names = vars)
}


# ==============================================================================
# Helper: attach spline matrix to hf without copying
# ==============================================================================
add_splines_to_hf <- function(hf, mat, names) {
  for (j in seq_along(names)) {
    data.table::set(hf, j = names[j], value = mat[, j])
  }
}


# ==============================================================================
# Helper: aggregate high-frequency spline basis to low frequency
# ==============================================================================
agg_hf_to_lf <- function(hf, id, tl, spline_vars,
                         hfcov = NULL, weights = NULL) {
  keys <- c(id, tl)

  agg <- hf[, lapply(.SD, sum, na.rm = TRUE), by = keys, .SDcols = spline_vars]

  if (length(hfcov)) {
    means <- hf[, lapply(.SD, mean, na.rm = TRUE), by = keys, .SDcols = hfcov]
    agg <- merge(agg, means, by = keys, all.x = TRUE)
  }

  if (!is.null(weights) && weights %in% names(hf)) {
    wmean <- hf[, lapply(.SD, mean, na.rm = TRUE), by = keys, .SDcols = weights]
    agg <- merge(agg, wmean, by = keys, all.x = TRUE)
  }

  return(agg)
}


# ==============================================================================
# Helper: build fixest formula (no weights here)
# ==============================================================================
build_formula <- function(y, varlist, partialout = NULL, absorb = NULL) {
  rhs <- paste(c(varlist, partialout), collapse = " + ")
  if (nzchar(rhs)) {
    fml <- paste0(y, " ~ ", rhs)
  } else {
    fml <- paste0(y, " ~ 0")
  }

  if (!is.null(absorb)) {
    if (inherits(absorb, "formula")) {
      abs_txt <- paste(deparse(absorb), collapse = "")
      abs_txt <- sub("^\\s*~\\s*", "", abs_txt)
    } else {
      abs_txt <- paste(absorb, collapse = " + ")
    }
    fml <- paste0(fml, " | ", abs_txt)
  }

  stats::as.formula(fml)
}


# ==============================================================================
# Helper: residualization formula (matches reghdfe y depvar [partialout] [, absorb])
# ==============================================================================
fml_residualize <- function(y, x_partial = NULL, absorb = NULL) {
  if (length(x_partial) > 0L || !is.null(absorb)) {
    build_formula(y = y, varlist = x_partial, partialout = NULL, absorb = absorb)
  } else {
    stats::as.formula(paste0(y, " ~ 1"))
  }
}


# ==============================================================================
# Helper: save outf data to file
# ==============================================================================
outf_save <- function(dt, save, replace = FALSE) {
  if (!replace && file.exists(save)) {
    stop("File '", save, "' already exists. Use replace = TRUE to overwrite.")
  }
  ext <- tolower(tools::file_ext(save))
  if (ext == "csv") {
    data.table::fwrite(dt, save)
  } else if (ext %in% c("rdata", "rda")) {
    outf_data <- dt
    save(outf_data, file = save)
  } else {
    saveRDS(dt, save)
  }
}


# ==============================================================================
# Helper: cross-validation RMSE
# ==============================================================================
rmse_cv <- function(data, y, varlist, partialout = NULL,
                    absorb = NULL, cluster = NULL, weights_fml = NULL,
                    cvvar = ".cv") {
  data <- data.table::as.data.table(data)
  folds <- sort(unique(data[[cvvar]]))

  resi2 <- numeric(nrow(data))
  nvars <- length(varlist)

  fml <- build_formula(y = y, varlist = varlist,
                       partialout = partialout, absorb = absorb)
  cluster_fml <- if (!is.null(cluster)) stats::as.formula(paste0("~ ", cluster)) else NULL

  for (f in folds) {
    train <- data[get(cvvar) != f]
    val   <- data[get(cvvar) == f]

    fit <- fixest::feols(fml, data = train, cluster = cluster_fml,
                         weights = weights_fml, warn = FALSE, notes = FALSE)

    b <- stats::coef(fit)
    # allow for collinear drops: keep only coefficients that correspond to varlist
    retained <- intersect(names(b), varlist)
    if (length(retained) == 0L) {
      stop("All varlist variables were dropped; check collinearity.")
    }
    b <- b[retained]

    X_val <- as.matrix(val[, retained, with = FALSE])
    y_val <- val[[y]]
    gjhat <- as.numeric(X_val %*% b)
    resid <- y_val - gjhat

    val_fit <- NULL
    if (length(partialout)) {
      val_dt <- data.table::as.data.table(val)
      val_dt[, .resid := resid]
      val_fml <- build_formula(y = ".resid", varlist = partialout,
                               partialout = NULL, absorb = absorb)
      val_fit <- tryCatch(
        fixest::feols(val_fml, data = val_dt, cluster = cluster_fml,
                      weights = weights_fml, warn = FALSE, notes = FALSE),
        error = function(e) NULL
      )
    }

    if (is.null(val_fit) && !is.null(absorb)) {
      # fall back to fixed-effects only (all partialout vars collinear with FE)
      val_dt <- data.table::as.data.table(val)
      val_dt[, .resid := resid]
      val_fml <- build_formula(y = ".resid", varlist = NULL,
                               partialout = NULL, absorb = absorb)
      val_fit <- fixest::feols(val_fml, data = val_dt, cluster = cluster_fml,
                               weights = weights_fml, warn = FALSE, notes = FALSE)
    }

    if (!is.null(val_fit)) {
      yhat_val <- stats::predict(val_fit, type = "response")
      err <- resid - yhat_val
    } else {
      err <- resid
    }

    resi2[data[[cvvar]] == f] <- err^2
  }

  return(sqrt(mean(resi2, na.rm = TRUE)))
}


# ==============================================================================
# Helper: wild bootstrap variance-covariance matrix from a fitted fixest formula
# ==============================================================================
wildboot_vcov_formula <- function(data, fml, y, cluster = NULL,
                                  weights_fml = NULL, brep = 100L,
                                  seed = NULL) {
  if (brep <= 0L) return(NULL)
  if (!is.null(seed)) set.seed(seed)

  data_dt <- data.table::as.data.table(data)
  cluster_fml <- if (!is.null(cluster)) stats::as.formula(paste0("~ ", cluster)) else NULL

  fit0 <- fixest::feols(fml, data = data_dt, cluster = cluster_fml,
                        weights = weights_fml, warn = FALSE, notes = FALSE)
  b <- stats::coef(fit0)
  k <- length(b)

  yhat <- stats::predict(fit0, type = "response")
  ehat <- data_dt[[y]] - yhat

  bb <- matrix(NA_real_, nrow = brep, ncol = k)

  for (r in seq_len(brep)) {
    if (is.null(cluster)) {
      radw <- sample(c(-1L, 1L), size = nrow(data_dt), replace = TRUE) * ehat
    } else {
      cl <- data_dt[[cluster]]
      ucl <- unique(cl)
      signs <- sample(c(-1L, 1L), size = length(ucl), replace = TRUE)
      names(signs) <- as.character(ucl)
      radw <- signs[as.character(cl)] * ehat
    }

    data_dt[, .ystar := yhat + radw]
    rhs_txt <- as.character(fml)[3L]
    fml_star <- stats::as.formula(paste0(".ystar ~ ", rhs_txt))
    fit_star <- fixest::feols(fml_star, data = data_dt, cluster = cluster_fml,
                              weights = weights_fml, warn = FALSE, notes = FALSE)
    b_star <- stats::coef(fit_star)
    idx <- match(names(b), names(b_star))
    bb[r, !is.na(idx)] <- b_star[idx[!is.na(idx)]]
  }

  data_dt[, .ystar := NULL]

  V <- stats::cov(bb)
  rownames(V) <- colnames(V) <- names(b)
  return(V)
}


# ==============================================================================
# Print method
# ==============================================================================
#' @export
print.mfxtsemipar_cv3 <- function(x, ...) {
  cat("Mixed-frequency semiparametric regression with cross-validation and robust bias correction\n")
  cat("  Selected knots: ", x$nknots, "\n", sep = "")
  cat("  Boundary knots: ", paste(round(x$bknots, 4), collapse = ", "), "\n", sep = "")
  cat("  Knot locations: ", paste(round(x$knots, 4), collapse = ", "), "\n", sep = "")
  cat("  Minimum CV RMSE: ", round(x$min_cv_mse, 6), "\n", sep = "")
  cat("  Final fit RMSE: ", round(x$rmse, 6), "\n", sep = "")
  cat("  Simple optimal knots: ", x$soptnk, "\n", sep = "")
  bc_type <- if (!is.null(x$bc_type)) paste0(x$bc_type, " (degree ", x$bc_degree, ")") else "none"
  cat("  Bias correction: ", bc_type, "\n", sep = "")
  if (!is.null(x$bc_est)) {
    cat("  BC estimation: ", x$bc_est, "\n", sep = "")
  }
  if (!is.null(x$brep) && x$brep > 0L) {
    cat("  Bootstrap replications: ", x$brep, "\n", sep = "")
  }
  cat("\nCV RMSE by number of knots:\n")
  print(x$cv_mse, row.names = FALSE)
  invisible(x)
}


# ==============================================================================
# Predict method for g(u)
# ==============================================================================
#' Predict the semiparametric component g(u) from an mfxtsemipar_cv3 fit
#'
#' Given a new set of values for the semiparametric variable \code{uvar},
#' reconstructs the spline basis using the estimated knots and settings, then
#' returns the fitted semiparametric component and its delta-method standard
#' error.
#'
#' @param object an object of class \code{mfxtsemipar_cv3}.
#' @param newdata a data.frame or data.table containing the evaluation variable,
#'   or a numeric vector of evaluation points.
#' @param uvar name of the evaluation variable in \code{newdata}. If
#'   \code{NULL}, the \code{uvar} used in estimation is used.
#' @param bias_correct logical; if \code{TRUE} (default), use the bias-corrected
#'   estimate when available; otherwise use the raw estimate.
#' @param ucb logical; if \code{TRUE} and \code{object} contains UCB
#'   information, add \code{lb} and \code{ub} columns to the output.
#' @param ... additional arguments (currently ignored).
#'
#' @return A data.table with columns \code{u}, \code{g} (the predicted
#'   g(u)), and \code{se}. When \code{ucb = TRUE}, also \code{lb} and
#'   \code{ub}.
#'
#' @export
predict.mfxtsemipar_cv3 <- function(object, newdata, uvar = NULL,
                                    bias_correct = TRUE, ucb = FALSE, ...) {
  if (!inherits(object, "mfxtsemipar_cv3")) {
    stop("object must be of class 'mfxtsemipar_cv3'.")
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

  B_retained <- B[, retained, drop = FALSE]

  if (bias_correct && !is.null(object$bc_type)) {
    bc_type <- object$bc_type
    retained_bc <- object$bc_retained_bc

    sp_bc <- make_splines(
      x = x,
      type = object$type,
      knots = if (is.null(object$bc_knots)) object$knots else object$bc_knots,
      bknots = object$bknots,
      degree = object$bc_degree,
      center = object$center,
      intercept = object$intercept,
      prefix = ".SplineBC_"
    )
    Btilde <- sp_bc$matrix
    Btilde_retained <- Btilde[, retained_bc, drop = FALSE]

    all_retained <- c(retained, retained_bc)
    b_joint <- object$coef_bc[all_retained]
    V_joint <- object$vcov_bc[all_retained, all_retained, drop = FALSE]

    # Orthogonalize BC evaluation basis using the stored projection.
    # bc_M was computed for all BC basis columns; subset to those retained
    # in the joint regression.
    proj_main <- object$bc_G0 %*% object$bc_M[, retained_bc, drop = FALSE]
    Btilde_orth_retained <- Btilde_retained - B_retained %*% proj_main

    w_bc <- cbind(B_retained, Btilde_orth_retained)

    fit <- as.numeric(w_bc %*% b_joint)
    se <- sqrt(pmax(0, rowSums((w_bc %*% V_joint) * w_bc)))
  } else {
    b <- object$coef[retained]
    fit <- as.numeric(B_retained %*% b)
    V <- object$vcov[retained, retained, drop = FALSE]
    se <- sqrt(pmax(0, rowSums((B_retained %*% V) * B_retained)))
  }

  out <- data.table::data.table(u = x, g = fit, se = se)

  if (ucb) {
    if (!bias_correct) {
      warning("ucb = TRUE requires bias_correct = TRUE in predict(); UCB omitted.",
              call. = FALSE)
    } else if (is.null(object$ucb)) {
      warning("No UCB information found in object; run mfxtsemipar_bc with ucb = TRUE.",
              call. = FALSE)
    } else {
      crit <- object$ucb$crit
      out[, lb := g - crit * se]
      out[, ub := g + crit * se]
    }
  }

  out
}


# ==============================================================================
# Diagnostic: compare estimated g(u) to true g(u)
# ==============================================================================
#' Compare the estimated semiparametric component to a true g(u)
#'
#' Predicts g(u) at the evaluation points provided in \code{newdata} and reports
#' standard accuracy metrics against the supplied true values.
#'
#' @param object an object of class \code{mfxtsemipar_cv3}.
#' @param newdata a data.frame or data.table containing the evaluation variable,
#'   or a numeric vector of evaluation points.
#' @param true_g numeric vector of true g(u) values, aligned with \code{newdata}.
#' @param uvar name of the evaluation variable in \code{newdata}. If
#'   \code{NULL}, the \code{uvar} used in estimation is used.
#' @param bias_correct logical; passed to \code{predict.mfxtsemipar_cv3}.
#'
#' @return A list with components \code{metrics} (a one-row data.table),
#'   \code{pred} (the prediction data.table), and \code{resid} (true - fitted).
#'
#' @export
g_diagnostic <- function(object, ...) {
  UseMethod("g_diagnostic")
}

g_diagnostic.mfxtsemipar_cv3 <- function(object, newdata, true_g,
                                         uvar = NULL, bias_correct = TRUE) {
  pred <- predict.mfxtsemipar_cv3(object, newdata = newdata, uvar = uvar,
                                  bias_correct = bias_correct)
  compute_g_diagnostic(pred, true_g)
}

compute_g_diagnostic <- function(pred, true_g) {
  fit <- pred$g
  se <- pred$se

  if (length(true_g) != length(fit)) {
    stop("length(true_g) must equal the number of evaluation points.")
  }

  resid <- true_g - fit

  metrics <- data.table::data.table(
    n = length(fit),
    mse = mean(resid^2, na.rm = TRUE),
    rmse = sqrt(mean(resid^2, na.rm = TRUE)),
    mae = mean(abs(resid), na.rm = TRUE),
    corr = stats::cor(fit, true_g, use = "pairwise.complete.obs"),
    max_abs_err = max(abs(resid), na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE)
  )

  list(metrics = metrics, pred = pred, resid = resid)
}

# ==============================================================================
# Helper: simulate from a multivariate normal with possibly semi-definite sigma
# ==============================================================================
if (!exists("simulate_mvnorm", mode = "function", inherits = FALSE)) {
  simulate_mvnorm <- function(n_draws, sigma) {
    sigma <- as.matrix((sigma + t(sigma)) / 2)
    eig <- eigen(sigma, symmetric = TRUE)
    values <- pmax(eig$values, 0)
    positive <- values > 0

    if (!any(positive)) {
      return(matrix(0, nrow = n_draws, ncol = ncol(sigma)))
    }

    loadings <- eig$vectors[, positive, drop = FALSE] %*%
      diag(sqrt(values[positive]), nrow = sum(positive))
    random_normals <- matrix(stats::rnorm(n_draws * ncol(loadings)), nrow = n_draws)
    random_normals %*% t(loadings)
  }
}

