# Test script for mfxtbin_outf
.libPaths(c("Rlib", .libPaths()))
source("R/mfxtbin_outf.R")

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
hf[, gf := fcase(x3 < 0, 1, x3 < 1, 2, default = 3)]
hf[, x1_tot := sum(x1), by = .(id, t)]
hf[, x2_tot := sum(x2), by = .(id, t)]
hf[, gf_tot := sum(gf), by = .(id, t)]
hf[, e_tot := sum(rnorm(.N)), by = .(id, t)]
lf <- merge(lf, unique(hf[, .(id, t, x1_tot, x2_tot, gf_tot, e_tot)]), by = c("id", "t"))
lf[, y := x1_tot - x2_tot + gf_tot + d + e_tot]
hf <- hf[, .(id, t, day, x1, x2, x3, gf)]
lf <- lf[, .(id, t, y, x2_tot, d)]
setnames(lf, "x2_tot", "x2")

# Create insample indicator
set.seed(123)
in_ids <- sample(unique(lf$id), size = floor(n_id * 0.7))
lf[, insample := as.integer(id %in% in_ids)]

cat("In-sample share:", mean(lf$insample), "\n")

# Test 1: fixed bins
save_file <- tempfile(fileext = ".rds")
mfxtbin_outf(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  insample = "insample",
  save = save_file,
  nbin = 4,
  absorb = ~ id,
  partialout = "all"
)

out1 <- readRDS(save_file)
cat("Fixed-bin output dimensions:", nrow(out1), "x", ncol(out1), "\n")
print(names(out1))
print(head(out1))

# Test 2: CV bins
save_file2 <- tempfile(fileext = ".rds")
mfxtbin_outf(
  hf = hf,
  lf = lf,
  y = "y",
  x = "x2",
  uvar = "x3",
  id = "id",
  tl = "t",
  insample = "insample",
  save = save_file2,
  maxnbin = 5,
  minnbin = 2,
  nfold = 5,
  seed = 123,
  absorb = ~ id,
  partialout = "all"
)

out2 <- readRDS(save_file2)
cat("CV-bin output dimensions:", nrow(out2), "x", ncol(out2), "\n")
print(head(out2))

# Test 3: explicit cut
save_file3 <- tempfile(fileext = ".csv")
mfxtbin_outf(
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
  cut = c(-1, 0, 1),
  absorb = ~ id
)

out3 <- data.table::fread(save_file3)
cat("CSV output dimensions:", nrow(out3), "x", ncol(out3), "\n")
print(names(out3))

cat("mfxtbin_outf tests completed.\n")
