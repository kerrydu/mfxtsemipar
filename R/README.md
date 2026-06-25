# R implementation of the `mfxtsemipar` family

This directory contains R translations of the Stata commands `mfxtsemipar_cv.ado`,
`mfxtbin_cv.ado`, `mfxtsemipar.ado`, `mfxtbin.ado`, `mfxtsemipar_outf.ado` and
`mfxtbin_outf.ado`.

## File

- `mfxtsemipar_cv.R` – main implementation. Source it with `source("R/mfxtsemipar_cv.R")`.

## Dependencies

- `data.table`
- `fixest`
- `splines` (base R)
- `splines2` (only for `type = "bs"`, `"ms"`, `"is"`, `"ibs"`)

Install them, for example, into a project-local library:

```r
install.packages(c("data.table", "fixest", "splines2"), lib = "Rlib")
```

and load with:

```r
.libPaths(c("Rlib", .libPaths()))
source("R/mfxtsemipar_cv.R")
```

## Key design difference from Stata

The high-frequency (`hf`) and low-frequency (`lf`) data are supplied as two
separate `data.frame`/`data.table` objects linked by `id` and `tl`. The spline
basis is generated on `hf`, summed within each `id × tl` cell, and merged into
`lf` for estimation. This avoids storing the expanded spline basis at low
frequency and keeps memory use low.

## Function signature

```r
mfxtsemipar_cv(
  hf, lf, y, x = NULL, uvar, id, tl, gen,
  hfcov = NULL, cluster = NULL, type = "poly",
  winsor = NULL, winsor_values = FALSE, eqspace = FALSE,
  maxnk = 5, minnk = 2, center = NULL, absorb = NULL,
  partialout = NULL, cvgroup = NULL, nfold = 10, seed = NULL,
  degree = 1, keepsplines = FALSE, atu = NULL,
  dropfirstbase = FALSE, sopt = FALSE, brep = 0,
  predy = NULL, weights = NULL
)
```

### Arguments

| Argument | Description |
|----------|-------------|
| `hf` | High-frequency data.table |
| `lf` | Low-frequency data.table (one row per `id × tl`) |
| `y` | Dependent variable name |
| `x` | Low-frequency covariates |
| `uvar` | Variable for the semiparametric component (in `hf`) |
| `id` | Panel identifier |
| `tl` | Time-level variable(s) defining low-frequency units |
| `gen` | Name for generated fitted values |
| `hfcov` | High-frequency covariates (averaged to LF) |
| `cluster` | Cluster variable for robust SE |
| `type` | Spline type: `"poly"`, `"bs"`, `"ms"`, `"is"`, `"ibs"` |
| `winsor` | Winsorization percentiles (length-2 numeric) |
| `winsor_values` | If `TRUE`, `winsor` contains raw cutoffs |
| `eqspace` | Equally spaced knots instead of quantiles |
| `maxnk` / `minnk` | Maximum / minimum number of knots |
| `center` | Centering value for splines |
| `absorb` | Fixed effects: formula (`~ id + year`) or character vector |
| `partialout` | `NULL`, `"all"`, or character vector of variables to partial out |
| `cvgroup` | User-supplied CV fold variable |
| `nfold` / `seed` | Number of folds and seed for automatic CV |
| `degree` | Polynomial degree |
| `keepsplines` | Keep spline columns in the output |
| `atu` | Alternative evaluation variable (in `hf`) |
| `dropfirstbase` | Drop intercept/base from spline basis |
| `sopt` | Use first local minimum of CV curve |
| `brep` | Wild-bootstrap replications |
| `predy` | Name for full LF prediction (`xb + FE`) |
| `weights` | Weight variable (in `hf` or `lf`) |

### Return value

A list of class `mfxtsemipar_cv` containing:

- `nknots`, `knots`, `min_cv_mse`, `soptnk`
- `cv_mse`: data.table of cross-validated RMSE by number of knots (column `cv_rmse`; used to select `nknots`)
- `rmse`: in-sample RMSE of the final `fixest::feols` fit
- `coef`, `vcov`: final model coefficients and variance matrix
- `info`: information criteria from `fixest`
- `fitted`: data.table with `id`, `tl`, `gen`, `gen_se` (HF level)
- `predy`: LF-level predictions if `predy` was requested
- `estimation`: the full `fixest` object

### Predicting `g(u)` on a new grid

After estimation, use the S3 `predict` method to evaluate the semiparametric
component at arbitrary points:

```r
grid <- data.table(x3 = seq(min(hf$x3), max(hf$x3), length.out = 100))
pred <- predict(res, newdata = grid)

# or pass a vector directly
pred <- predict(res, newdata = c(-1, 0, 1))
```

The returned object is a standalone `data.table` with exactly three columns:
`u` (the evaluation point), `g` (the predicted `g(u)`), and `se` (delta-method
standard error).

### Diagnostics: compare estimated `g(u)` to truth

If you know the true `g(u)` (e.g. in simulations), use `g_diagnostic()`:

```r
grid <- data.table(x3 = seq(min(hf$x3), max(hf$x3), length.out = 100))
grid[, gf_true := 1 * x3 + 2 * x3^2 - 0.25 * x3^3]

diag <- g_diagnostic(res, newdata = grid, true_g = grid$gf_true)
print(diag$metrics)
# columns: n, mse, rmse, mae, corr, max_abs_err, mean_se
```

## Example

```r
library(data.table)
source("R/mfxtsemipar_cv.R")

# Assume hf and lf are already built:
# hf: id, t, x3 (high-frequency uvar)
# lf: id, t, y, x2 (low-frequency)

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
  absorb = ~ id,
  partialout = "all",
  maxnk = 5,
  minnk = 2,
  nfold = 5,
  seed = 123
)

print(res)
head(res$fitted)
```

See `test_mfxtsemipar_cv.R` for a complete simulated-data example.

## Differences from the Stata original

- `partialout` unifies Stata's `partialout` and `partialout1(...)`:
  - `partialout = "all"` = partial out all covariates.
  - `partialout = c("x1", "x2")` = partial out only those variables.
- Weights are supplied as a single variable name; if it lives in `hf`, it is
  averaged to the LF level automatically.
- `fixest::feols` replaces `reghdfe`.


---

# R implementation of `mfxtbin_cv`

`R/mfxtbin_cv.R` provides a binned analogue of `mfxtsemipar_cv`. Instead of
splines, the nonparametric component is a step function defined by bins. The
number of bins is chosen by cross-validation.

## File

- `mfxtbin_cv.R` – main implementation. It sources `mfxtsemipar_cv.R` for
  shared helpers, so loading it also makes `mfxtsemipar_cv()` available.

## Function signature

```r
mfxtbin_cv(
  hf, lf, y, x = NULL, uvar, id, tl, gen,
  hfcov = NULL, cluster = NULL,
  maxnbin = 5, minnbin = 2, eqspace = FALSE,
  dropbin = NULL, startp = NULL, endp = NULL,
  cvgroup = NULL, nfold = 10, seed = NULL,
  atu = NULL, absorb = NULL, sopt = FALSE,
  predy = NULL, partialout = NULL, weights = NULL
)
```

Key options:

| Argument | Description |
|----------|-------------|
| `maxnbin` / `minnbin` | Maximum / minimum number of bins |
| `eqspace` | Equally spaced cutpoints instead of quantiles |
| `dropbin` | Integer bin number or numeric value whose bin to drop |
| `startp` / `endp` | Bounds used to form cutpoints |
| `atu` | Alternative evaluation variable |

### Return value

A list of class `mfxtbin_cv` containing:

- `nbin`, `cutpoints`, `min_cv_mse`, `soptbin`
- `cv_mse`: data.table of cross-validated RMSE by number of bins (column `cv_rmse`; used to select `nbin`)
- `rmse`: in-sample RMSE of the final `fixest::feols` fit
- `coef`, `vcov`, `info`, `fitted`, `predy`, `estimation`

### Predict and diagnose

The same S3 methods work for `mfxtbin_cv` objects:

```r
pred <- predict(res_bin, newdata = grid)   # columns u, g, se
diag <- g_diagnostic(res_bin, newdata = grid, true_g = grid$gf_true)
```

## Example

```r
source("R/mfxtbin_cv.R")

res_bin <- mfxtbin_cv(
  hf = hf, lf = lf,
  y = "y", x = "x2", uvar = "x3",
  id = "id", tl = "t", gen = "fitted_bin",
  maxnbin = 6, minnbin = 2,
  absorb = ~ id, partialout = "all",
  nfold = 5, seed = 123
)

print(res_bin)
head(res_bin$fitted)
```

See `test_mfxtbin_cv.R` for a complete simulated-data example.


---

# R implementation of `mfxtsemipar` (fixed knots)

`R/mfxtsemipar.R` is the non-cross-validation version of
`mfxtsemipar_cv`. You supply the knots directly (or the number of knots to
be generated).

## Function signature

```r
mfxtsemipar(
  hf, lf, y, x = NULL, uvar, id, tl, gen,
  hfcov = NULL, cluster = NULL,
  bknots = NULL, degree = 1, knots = NULL, type = "poly",
  winsor = NULL, winsor_values = FALSE, eqspace = FALSE,
  nknots = NULL, center = NULL, absorb = NULL,
  atu = NULL, intercept = TRUE, brep = 0,
  predy = NULL, weights = NULL
)
```

Either `knots` or `nknots` must be supplied (but not both).

## Example

```r
source("R/mfxtsemipar.R")

res <- mfxtsemipar(
  hf = hf, lf = lf,
  y = "y", x = "x2", uvar = "x3",
  id = "id", tl = "t", gen = "fitted",
  type = "poly", degree = 1, nknots = 5, center = 0,
  absorb = ~ id, intercept = TRUE
)

pred <- predict(res, newdata = grid)
```

See `test_mfxtsemipar.R` for a complete example.

---

# R implementation of `mfxtbin` (fixed bins)

`R/mfxtbin.R` is the non-cross-validation version of `mfxtbin_cv`. You
supply the cutpoints or the number of bins.

## Function signature

```r
mfxtbin(
  hf, lf, y, x = NULL, uvar, id, tl, gen,
  hfcov = NULL, cluster = NULL,
  nbin = NULL, eqspace = FALSE, cut = NULL,
  absorb = NULL, bw = NULL, dropbin = NULL,
  atu = NULL, startp = NULL, endp = NULL,
  predy = NULL, weights = NULL
)
```

Exactly one of `cut`, `nbin` or `bw` must be supplied.

## Example

```r
source("R/mfxtbin.R")

res <- mfxtbin(
  hf = hf, lf = lf,
  y = "y", x = "x2", uvar = "x3",
  id = "id", tl = "t", gen = "fitted_bin",
  nbin = 4, absorb = ~ id
)

pred <- predict(res, newdata = grid)
```

See `test_mfxtbin.R` for a complete example.


---

# R implementation of `mfxtsemipar_outf` (out-of-sample spline prediction)

`R/mfxtsemipar_outf.R` fits the spline model on the in-sample observations and
saves predicted values for both in-sample and out-of-sample observations. This
matches the Stata command `mfxtsemipar_outf`.

## Function signature

```r
mfxtsemipar_outf(
  hf, lf, y, x = NULL, uvar, id, tl, insample, save, replace = FALSE,
  hfcov = NULL, cluster = NULL, type = "poly",
  winsor = NULL, winsor_values = FALSE, eqspace = FALSE,
  maxnk = 5, minnk = 2, allknots = NULL, bknots = NULL,
  nknots = NULL, knots = NULL, center = NULL, absorb = NULL,
  cvgroup = NULL, nfold = 10, seed = NULL,
  degree = 1, dropfirstbase = FALSE, sopt = FALSE,
  predy = NULL, partialout = NULL, weights = NULL
)
```

## Key arguments

| Argument | Description |
|----------|-------------|
| `insample` | Name of a 0/1 indicator variable in `lf`; `1` = fit / predict, `0` = predict only |
| `save` | Path to the output file (`.csv` writes CSV, otherwise RDS) |
| `replace` | Overwrite `save` if it already exists |
| `nknots` / `knots` / `allknots` | Fixed knot specification (no CV) |
| `maxnk` / `mink` / `cvgroup` / `nfold` / `seed` / `sopt` | CV knot selection on in-sample data only |

All other arguments are the same as in `mfxtsemipar_cv()` / `mfxtsemipar()`.

If none of `nknots`, `knots` or `allknots` is supplied, the number of knots is
selected by cross-validation on the in-sample observations only.

## Saved output

The output file contains one row per `id × tl` with the columns:

- `id`, `tl`, `insample`, `y`
- `pred_<y>` (or the name supplied in `predy`): full predicted value
- `_M_<y>`: residual of `y` after partialling out `partialout` and fixed effects
- `_M_ghat`: residual of the predicted semiparametric component after partialling
  out `partialout` and fixed effects

## Example

```r
source("R/mfxtsemipar_outf.R")

# assume lf has an insample indicator
mfxtsemipar_outf(
  hf = hf, lf = lf,
  y = "y", x = "x2", uvar = "x3",
  id = "id", tl = "t",
  insample = "insample",
  save = "predictions.rds",
  nknots = 5, type = "poly", degree = 1,
  absorb = ~ id, partialout = "all"
)

pred <- readRDS("predictions.rds")
head(pred)
```

See `test_mfxtsemipar_outf.R` for a complete example.


---

# R implementation of `mfxtbin_outf` (out-of-sample binned prediction)

`R/mfxtbin_outf.R` is the binned analogue of `mfxtsemipar_outf`. It fits the
binned model on the in-sample observations and saves predicted values for both
in-sample and out-of-sample observations, matching Stata's `mfxtbin_outf`.

## Function signature

```r
mfxtbin_outf(
  hf, lf, y, x = NULL, uvar, id, tl, insample, save, replace = FALSE,
  maxnbin = 5, minnbin = 2,
  nbin = NULL, cut = NULL, bw = NULL,
  eqspace = FALSE, dropbin = NULL,
  startp = NULL, endp = NULL,
  hfcov = NULL, cluster = NULL,
  cvgroup = NULL, nfold = 10, seed = NULL,
  atu = NULL, absorb = NULL, sopt = FALSE,
  predy = NULL, partialout = NULL, weights = NULL
)
```

## Key arguments

| Argument | Description |
|----------|-------------|
| `insample` | 0/1 indicator variable in `lf` |
| `save` / `replace` | Output file and overwrite flag |
| `nbin` / `cut` / `bw` | Fixed bin specification (no CV) |
| `maxnbin` / `minnbin` / `cvgroup` / `nfold` / `seed` / `sopt` | CV bin selection on in-sample data only |

Exactly one of `cut`, `nbin` or `bw` may be supplied for fixed bins. If none is
supplied, the number of bins is selected by cross-validation on the in-sample
observations only.

## Saved output

Same as for `mfxtsemipar_outf`: `id`, `tl`, `insample`, `y`, `pred_<y>`,
`_M_<y>`, `_M_ghat`.

## Example

```r
source("R/mfxtbin_outf.R")

mfxtbin_outf(
  hf = hf, lf = lf,
  y = "y", x = "x2", uvar = "x3",
  id = "id", tl = "t",
  insample = "insample",
  save = "predictions_bin.rds",
  nbin = 4, absorb = ~ id, partialout = "all"
)

pred <- readRDS("predictions_bin.rds")
head(pred)
```

See `test_mfxtbin_outf.R` for a complete example.
