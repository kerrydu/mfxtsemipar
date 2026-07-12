# Test mfxtsemipar_bc (single curve with bias correction and UCB)
.libPaths(c("Rlib", .libPaths()))
source("R/mfxtsemipar_bc.R")

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

# Basic run with default bias correction
res <- mfxtsemipar_bc(
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
  absorb = ~ id
)
print(res)
stopifnot(all(c("fitted", "fitted_se") %in% names(res$fitted)))
cat("Basic mfxtsemipar_bc run OK.\n")

# Run with UCB
res_ucb <- mfxtsemipar_bc(
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
  absorb = ~ id,
  ucb = TRUE,
  ucb_level = 0.95,
  ucb_sim_reps = 500L
)
stopifnot(all(c("fitted_lb", "fitted_ub") %in% names(res_ucb$fitted)))
stopifnot(!is.null(res_ucb$ucb))
stopifnot(res_ucb$ucb$level == 0.95)
stopifnot(res_ucb$ucb$sim_reps == 500L)

cat("UCB columns present.\n")
grid <- data.table(x3 = seq(min(hf$x3), max(hf$x3), length.out = 20))
pred_ucb <- predict(res_ucb, newdata = grid, uvar = "x3", ucb = TRUE)
stopifnot(all(c("u", "g", "se", "lb", "ub") %in% names(pred_ucb)))
cat("Predict with UCB OK.\n")

# UCB should be skipped when bias_correct = FALSE
res_nobc <- mfxtsemipar_bc(
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
  absorb = ~ id,
  bias_correct = FALSE,
  ucb = TRUE,
  ucb_sim_reps = 100L
)
stopifnot(is.null(res_nobc$ucb))
cat("UCB skip without bias_correct OK.\n")

cat("All mfxtsemipar_bc tests completed.\n")
