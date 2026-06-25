# Comparison of mfxtsemipar_cv (sum of g(x)) vs mfxtsemipar_cv_avg (g(mean(x)))
.libPaths(c("Rlib", .libPaths()))
source("R/mfxtsemipar_cv.R")
source("R/mfxtsemipar_cv_avg.R")
library(data.table)

set.seed(123456)

# Simulate mixed-frequency panel data (same DGP as test_mfxtsemipar_cv.R)
n_id <- 50
n_t <- 20
n_hf <- 10

Sigma <- matrix(c(1, 0, 0.42,
                  0, 1, 0.85,
                  0.42, 0.85, 1), nrow = 3)
lf_vars <- MASS::mvrnorm(n_id, mu = rep(0, 3), Sigma = Sigma)

lf <- data.table(
  id = 1:n_id,
  x2f = lf_vars[, 1],
  x3f = lf_vars[, 2],
  d = lf_vars[, 3]
)
lf <- lf[rep(1:.N, each = n_t)]
lf[, t := rep(1:n_t, length.out = .N)]

D <- matrix(c(1, 0.2, 0.8,
              0.2, 1, 0,
              0.8, 0, 1), nrow = 3)
N <- n_id * n_t * n_hf
hf_noise <- MASS::mvrnorm(N, mu = rep(0, 3), Sigma = D)

hf <- lf[rep(1:.N, each = n_hf)]
hf[, day := 1:.N, by = .(id, t)]
hf[, x1 := hf_noise[, 1]]
hf[, x2e := hf_noise[, 2]]
hf[, x3e := hf_noise[, 3]]
hf[, x2 := (x2f + x2e) / sqrt(2)]
hf[, x3 := (x3f + x3e) / sqrt(2)]
hf[, gf := 1 * x3 + 2 * x3^2 - 0.25 * x3^3]

# Low-frequency aggregates
hf[, x1_tot := sum(x1), by = .(id, t)]
hf[, x2_tot := sum(x2), by = .(id, t)]
hf[, gf_tot := sum(gf), by = .(id, t)]
hf[, e_tot := sum(rnorm(.N)), by = .(id, t)]
lf <- merge(lf, unique(hf[, .(id, t, x1_tot, x2_tot, gf_tot, e_tot)]),
            by = c("id", "t"))
lf[, y := x1_tot - x2_tot + gf_tot + d + e_tot]

# Keep only necessary columns
hf <- hf[, .(id, t, day, x1, x2, x3, gf)]
lf <- lf[, .(id, t, y, x2_tot, d)]
setnames(lf, "x2_tot", "x2")

# ---------------------------------------------------------------------------
# Original estimator: sum of g(x) at high frequency
# ---------------------------------------------------------------------------
res_sum <- mfxtsemipar_cv(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "g_sum",
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 5,
  minnk = 2,
  nfold = 5,
  seed = 123,
  absorb = ~ id,
  partialout = "all"
)

# ---------------------------------------------------------------------------
# Counterpart estimator: average x3 to LF, then model g(mean(x3))
# ---------------------------------------------------------------------------
res_avg <- mfxtsemipar_cv_avg(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "g_avg",
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 5,
  minnk = 2,
  nfold = 5,
  seed = 123,
  absorb = ~ id,
  partialout = "all"
)

cat("\n=== Original (sum of g(x)) ===\n")
print(res_sum)

cat("\n=== Counterpart (g(mean(x))) ===\n")
print(res_avg)

# ---------------------------------------------------------------------------
# Compare fitted values at the LF level
# ---------------------------------------------------------------------------
# True sum of g(x3) at LF
lf_true <- unique(hf[, .(id, t, x3, gf)])
lf_true[, gf_tot_true := sum(gf), by = .(id, t)]
lf_true[, x3_mean_true := mean(x3), by = .(id, t)]
lf_true[, gf_at_mean_true := 1 * x3_mean_true +
          2 * x3_mean_true^2 - 0.25 * x3_mean_true^3]
lf_true <- unique(lf_true[, .(id, t, gf_tot_true, gf_at_mean_true, x3_mean_true)])

# Aggregate original HF fitted values to LF (sum of estimated g(x))
sum_fit <- res_sum$fitted[, .(id, t, g_sum)]
sum_fit[, g_sum_hat := sum(g_sum), by = .(id, t)]
sum_fit <- unique(sum_fit[, .(id, t, g_sum_hat)])

avg_fit <- res_avg$fitted[, .(id, t, x3 = x3, g_avg = g_avg)]
setnames(avg_fit, "x3", "x3_mean_hat")

comp <- merge(lf_true, sum_fit, by = c("id", "t"))
comp <- merge(comp, avg_fit, by = c("id", "t"))

bias_sum <- comp[, mean(gf_tot_true - g_sum_hat)]
bias_avg <- comp[, mean(gf_at_mean_true - g_avg)]
rmse_sum <- comp[, sqrt(mean((gf_tot_true - g_sum_hat)^2))]
rmse_avg <- comp[, sqrt(mean((gf_at_mean_true - g_avg)^2))]
agg_bias <- comp[, mean(gf_tot_true - gf_at_mean_true)]

cat("\n=== Bias comparison at LF level ===\n")
cat(sprintf("Aggregation bias (true sum(g(x)) - true g(mean(x))): % .4f\n",
            agg_bias))
cat(sprintf("Original  sum(g(x)) vs true sum(g(x)):  bias = % .4f, RMSE = %.4f\n",
            bias_sum, rmse_sum))
cat(sprintf("Counterpart g(mean(x)) vs true g(mean(x)): bias = % .4f, RMSE = %.4f\n",
            bias_avg, rmse_avg))

cat("\nTest completed.\n")
