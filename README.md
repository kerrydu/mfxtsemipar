## Stata codes for Mixed-Frequency Semiparametric Panel Data Model

### R implementation

R translations are available in the `R/` directory:

- `mfxtsemipar_cv.R` – spline-based mixed-frequency semiparametric regression
  with cross-validated knot selection.
- `mfxtbin_cv.R` – binned (step-function) version with cross-validated bin
  selection.
- `mfxtsemipar.R` – spline-based version with user-specified fixed knots.
- `mfxtbin.R` – binned version with user-specified fixed bins/cutpoints.
- `mfxtsemipar_outf.R` – out-of-sample spline prediction, fitted on in-sample
  observations and saved to a file.
- `mfxtbin_outf.R` – out-of-sample binned prediction, fitted on in-sample
  observations and saved to a file.

All keep high-frequency and low-frequency data in separate `data.table`s
linked by `id` and `tl`, which saves memory compared to the long-format Stata
workflow. See `R/README.md` for usage and examples.
