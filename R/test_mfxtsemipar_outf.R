# Test script for mfxtsemipar_outf
.libPaths(c("Rlib", .libPaths()))
source("R/mfxtsemipar_outf.R")

library(data.table)

set.seed(123456)

# Simulate mixed-frequency panel data
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

hf[, x1_tot := sum(x1), by = .(id, t)]
hf[, x2_tot := sum(x2), by = .(id, t)]
hf[, e_tot := sum(rnorm(.N)), by = .(id, t)]
hf[, gf_tot := sum(gf), by = .(id, t)]

lf <- merge(lf, unique(hf[, .(id, t, x1_tot, x2_tot, gf_tot, e_tot)]), by = c("id", "t"))
lf[, y := x1_tot - x2_tot + gf_tot + d + e_tot]

hf <- hf[, .(id, t, day, x1, x2, x3, gf)]
lf <- lf[, .(id, t, y, x2_tot, d)]
setnames(lf, "x2_tot", "x2")

# Create insample indicator (clustered by id)
set.seed(123)
in_ids <- sample(unique(lf$id), size = floor(n_id * 0.7))
lf[, insample := as.integer(id %in% in_ids)]

cat("In-sample share:", mean(lf$insample), "\n")

# Test 1: fixed knots
save_file <- tempfile(fileext = ".rds")
mfxtsemipar_outf(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  insample = "insample",
  save = save_file,
  nknots = 4,
  type = "poly",
  degree = 1,
  center = 0,
  absorb = ~ id,
  partialout = "all",
  predy = "pred_y"
)

out1 <- readRDS(save_file)
cat("Fixed-knot output dimensions:", nrow(out1), "x", ncol(out1), "\n")
print(names(out1))
print(head(out1))

# Test 2: CV knots
save_file2 <- tempfile(fileext = ".rds")
mfxtsemipar_outf(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  insample = "insample",
  save = save_file2,
  type = "poly",
  degree = 1,
  center = 0,
  maxnk = 4,
  minnk = 2,
  nfold = 5,
  seed = 123,
  absorb = ~ id,
  partialout = "all"
)

out2 <- readRDS(save_file2)
cat("CV-knot output dimensions:", nrow(out2), "x", ncol(out2), "\n")
print(head(out2))

# Test 3: no partialout
save_file3 <- tempfile(fileext = ".csv")
mfxtsemipar_outf(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  insample = "insample",
  save = save_file3,
  replace = TRUE,
  nknots = 3,
  type = "poly",
  degree = 1,
  absorb = ~ id
)

out3 <- data.table::fread(save_file3)
cat("CSV output dimensions:", nrow(out3), "x", ncol(out3), "\n")
print(names(out3))

cat("mfxtsemipar_outf tests completed.\n")
