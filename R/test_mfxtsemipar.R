.libPaths(c("Rlib", .libPaths()))
source("R/mfxtsemipar.R")
library(data.table)
set.seed(123456)

n_id <- 100
n_t <- 30
n_hf <- 10

lf <- data.table(id = 1:n_id, x2f = rnorm(n_id), d = rnorm(n_id))
lf <- lf[rep(1:.N, each = n_t)]
lf[, t := rep(1:n_t, length.out = .N)]

N <- n_id * n_t * n_hf
hf <- lf[rep(1:.N, each = n_hf)]
hf[, day := 1:.N, by = .(id, t)]
hf[, x1 := rnorm(.N)]
hf[, x2e := rnorm(.N)]
hf[, x2 := (x2f + x2e) / sqrt(2)]
hf[, x3 := rnorm(.N)]
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

# Estimate with fixed knots
res <- mfxtsemipar(
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
  nknots = 5,
  center = 0,
  absorb = ~ id,
  intercept = TRUE
)

print(res)
head(res$fitted)

# Predict on grid
grid <- data.table(x3 = seq(quantile(hf$x3, 0.05), quantile(hf$x3, 0.95), length.out = 100))
grid[, gf_true := 1 * x3 + 2 * x3^2 - 0.25 * x3^3]
pred <- predict(res, newdata = grid)
head(pred)

diag <- g_diagnostic(res, newdata = grid, true_g = grid$gf_true)
print(diag$metrics)

# Test with explicit knots
res2 <- mfxtsemipar(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "fitted2",
  type = "poly",
  degree = 1,
  knots = c(-1, 0, 1),
  center = 0,
  absorb = ~ id,
  intercept = TRUE
)
print(res2)

cat("mfxtsemipar test completed.\n")
