# Test hfcov, weights, partialout list, and bootstrap
.libPaths(c("Rlib", .libPaths()))
source("R/mfxtsemipar_cv.R")

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
# high-frequency covariate
hf[, zhf := rnorm(.N)]
hf[, gf := 1 * x3 + 2 * x3^2 - 0.25 * x3^3 + 0.5 * zhf]
hf[, x1_tot := sum(x1), by = .(id, t)]
hf[, x2_tot := sum(x2), by = .(id, t)]
hf[, gf_tot := sum(gf), by = .(id, t)]
hf[, e_tot := sum(rnorm(.N)), by = .(id, t)]

lf <- merge(lf, unique(hf[, .(id, t, x1_tot, x2_tot, gf_tot, e_tot)]), by = c("id", "t"))
lf[, y := x1_tot - x2_tot + gf_tot + d + e_tot]

# use hf-level weights
hf[, w := runif(.N, 0.5, 2)]

# low-frequency control x2 and high-frequency covariate zhf
hf <- hf[, .(id, t, day, x1, x2, x3, zhf, gf, w)]
lf <- lf[, .(id, t, y, x2_tot, d)]
setnames(lf, "x2_tot", "x2")

# partialout only zhf (hfcov), keep x2 in main regression
res <- mfxtsemipar_cv(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "fitted",
  hfcov = "zhf",
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 3,
  minnk = 1,
  nfold = 5,
  seed = 123,
  absorb = ~ id,
  partialout = "zhf",
  weights = "w"
)
print(res)
head(res$fitted)

# Test wild bootstrap (small brep)
res_boot <- mfxtsemipar_cv(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  gen = "fitted_boot",
  hfcov = "zhf",
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 3,
  minnk = 1,
  nfold = 5,
  seed = 123,
  absorb = ~ id,
  partialout = "zhf",
  brep = 50,
  cluster = "id"
)
print(res_boot)
cat("Bootstrap vcov dim:", dim(res_boot$vcov), "\n")

cat("All advanced tests completed.\n")
