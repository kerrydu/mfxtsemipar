# Test mfxtsemipar2_cv (short-/long-run curves with individual means)
.libPaths(c("Rlib", .libPaths()))
source("R/mfxtsemipar2_cv.R")

library(data.table)

set.seed(123456)

n_id <- 50
n_t <- 20
n_hf <- 10

Sigma <- matrix(c(1, 0, 0.42, 0, 1, 0.85, 0.42, 0.85, 1), nrow = 3)
lf_vars <- MASS::mvrnorm(n_id, mu = rep(0, 3), Sigma = Sigma)
lf <- data.table(id = 1:n_id, x2f = lf_vars[, 1], x3f = lf_vars[, 2], d = lf_vars[, 3])
lf <- lf[rep(1:.N, each = n_t)]
lf[, t := rep(1:n_t, length.out = .N)]

N <- n_id * n_t * n_hf
D <- matrix(c(1, 0.2, 0.8, 0.2, 1, 0, 0.8, 0, 1), nrow = 3)
hf_noise <- MASS::mvrnorm(N, mu = rep(0, 3), Sigma = D)
hf <- lf[rep(1:.N, each = n_hf)]
hf[, day := 1:.N, by = .(id, t)]
hf[, x1 := hf_noise[, 1]]
hf[, x2e := hf_noise[, 2]]
hf[, x3e := hf_noise[, 3]]
hf[, x2 := (x2f + x2e) / sqrt(2)]
hf[, x3 := (x3f + x3e) / sqrt(2)]
hf[, gf := 1 * x3 + 2 * x3^2 - 0.25 * x3^3]
hf[, x1_tot := sum(x1), by = .(id, t)]
hf[, x2_tot := sum(x2), by = .(id, t)]
hf[, gf_tot := sum(gf), by = .(id, t)]
hf[, e_tot := sum(rnorm(.N)), by = .(id, t)]

lf <- merge(lf, unique(hf[, .(id, t, x1_tot, x2_tot, gf_tot, e_tot)]), by = c("id", "t"))
lf[, y := x1_tot - x2_tot + gf_tot + d + e_tot]
hf <- hf[, .(id, t, day, x1, x2, x3, gf)]
lf <- lf[, .(id, t, y, x2_tot, d)]
setnames(lf, "x2_tot", "x2")

# Basic run with default bias correction and no iabsorb
res <- mfxtsemipar2_cv(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "fitted",
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 3,
  minnk = 1,
  nfold = 5,
  seed = 123,
  absorb = ~ t,
  partialout = "all",
  iabsorb = FALSE
)
print(res)
stopifnot(all(c("fitted_sr", "fitted_lr", "fitted_sr_se", "fitted_lr_se") %in% names(res$fitted)))
stopifnot(all(c("fitted_raw_sr", "fitted_raw_lr") %in% names(res$fitted)))
cat("Basic run OK.\n")

# Run with UCB
res_ucb <- mfxtsemipar2_cv(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "fitted",
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 3,
  minnk = 1,
  nfold = 5,
  seed = 123,
  absorb = ~ t,
  partialout = "all",
  iabsorb = FALSE,
  ucb = TRUE,
  ucb_level = 0.95,
  ucb_sim_reps = 500L
)
stopifnot(all(c("fitted_sr_lb", "fitted_sr_ub", "fitted_lr_lb", "fitted_lr_ub") %in% names(res_ucb$fitted)))
stopifnot(!is.null(res_ucb$ucb))
stopifnot(res_ucb$ucb$level == 0.95)
stopifnot(res_ucb$ucb$sim_reps == 500L)
grid_ucb <- data.table(x3 = seq(min(hf$x3), max(hf$x3), length.out = 20))
stopifnot(all(c("u", "g_sr", "g_lr", "se_sr", "se_lr", "lb_sr", "ub_sr", "lb_lr", "ub_lr") %in%
  names(predict(res_ucb, newdata = grid_ucb, uvar = "x3", ucb = TRUE))))
cat("UCB run OK.\n")

# UCB should be skipped when bias_correct = FALSE
res_nobc_warn <- mfxtsemipar2_cv(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "fitted_nobc_ucb",
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 3,
  minnk = 1,
  nfold = 5,
  seed = 123,
  absorb = ~ t,
  partialout = "all",
  bias_correct = FALSE,
  iabsorb = FALSE,
  ucb = TRUE,
  ucb_sim_reps = 100L
)
stopifnot(is.null(res_nobc_warn$ucb))
cat("UCB skip without bias_correct OK.\n")

# Run with iabsorb = TRUE
res_ia <- mfxtsemipar2_cv(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "fitted_ia",
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 3,
  minnk = 1,
  nfold = 5,
  seed = 123,
  absorb = ~ t,
  partialout = "all",
  iabsorb = TRUE
)
print(res_ia)
cat("iabsorb run OK.\n")

# Run without bias correction (closer to Stata mfxtsemipar2_cv)
res_nobc <- mfxtsemipar2_cv(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "fitted_nobc",
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 3,
  minnk = 1,
  nfold = 5,
  seed = 123,
  absorb = ~ t,
  partialout = "all",
  bias_correct = FALSE,
  iabsorb = FALSE
)
print(res_nobc)
stopifnot(all(c("fitted_nobc_sr", "fitted_nobc_lr") %in% names(res_nobc$fitted)))
cat("No-BC run OK.\n")

# Test B-splines with iabsorb
res_bs <- mfxtsemipar2_cv(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "fitted_bs",
  type = "bs",
  degree = 3,
  maxnk = 3,
  minnk = 1,
  nfold = 5,
  seed = 123,
  absorb = ~ t,
  partialout = "all",
  iabsorb = TRUE
)
print(res_bs)
cat("B-spline iabsorb run OK.\n")

# Test predict method
grid <- data.table(x3 = seq(min(hf$x3), max(hf$x3), length.out = 20))
pred <- predict(res, newdata = grid, uvar = "x3")
stopifnot(all(c("u", "g_sr", "g_lr", "se_sr", "se_lr") %in% names(pred)))
cat("Predict method OK.\n")

cat("All mfxtsemipar2_cv tests completed.\n")
