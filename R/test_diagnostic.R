.libPaths(c("Rlib", .libPaths()))
source("R/mfxtsemipar_cv.R")
library(data.table)
set.seed(123456)

n_id <- 100
n_t <- 30
n_hf <- 10

# LF variables
lf <- data.table(
  id = 1:n_id,
  x2f = rnorm(n_id),
  d = rnorm(n_id)
)
lf <- lf[rep(1:.N, each = n_t)]
lf[, t := rep(1:n_t, length.out = .N)]

# HF variables
N <- n_id * n_t * n_hf
hf <- lf[rep(1:.N, each = n_hf)]
hf[, day := 1:.N, by = .(id, t)]
hf[, x1 := rnorm(.N)]
hf[, x2e := rnorm(.N)]
hf[, x2 := (x2f + x2e) / sqrt(2)]
hf[, x3 := rnorm(.N)]

# True semiparametric component at HF: g(x3) = x3 + 2*x3^2 - 0.25*x3^3
hf[, gf := 1 * x3 + 2 * x3^2 - 0.25 * x3^3]

# Aggregate to LF: totals within id*t
hf[, x1_tot := sum(x1), by = .(id, t)]
hf[, x2_tot := sum(x2), by = .(id, t)]
hf[, gf_tot := sum(gf), by = .(id, t)]
hf[, e_tot := sum(rnorm(.N)), by = .(id, t)]

# LF outcome
lf <- merge(lf, unique(hf[, .(id, t, x1_tot, x2_tot, gf_tot, e_tot)]), by = c("id", "t"))
lf[, y := x1_tot - x2_tot + gf_tot + d + e_tot]

# Keep only needed columns
hf <- hf[, .(id, t, day, x1, x2, x3, gf)]
lf <- lf[, .(id, t, y, x2_tot, d)]
setnames(lf, "x2_tot", "x2")

# Estimate
res <- mfxtsemipar_cv(
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
  maxnk = 8,
  minnk = 2,
  nfold = 5,
  seed = 123,
  absorb = ~ id,
  partialout = "all"
)

# Compare on a grid
grid <- data.table(x3 = seq(quantile(hf$x3, 0.05), quantile(hf$x3, 0.95), length.out = 100))
grid[, gf_true := 1 * x3 + 2 * x3^2 - 0.25 * x3^3]

diag <- g_diagnostic(res, newdata = grid, true_g = grid$gf_true)
print(diag$metrics)

# Compare on original HF observations
diag_hf <- g_diagnostic(res, newdata = hf, true_g = hf$gf)
print(diag_hf$metrics)

# Also compare predicted g to true g at a few specific points
print(head(cbind(
  x3 = hf$x3[1:10],
  true_g = hf$gf[1:10],
  fitted = res$fitted$fitted[1:10]
)))

cat("diagnostic test completed.\n")
