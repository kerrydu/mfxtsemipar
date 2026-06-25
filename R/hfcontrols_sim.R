#!/usr/bin/env Rscript
# =============================================================================
# Simulation: High-frequency confounders and semiparametric g(z) estimation
# Replicates HFCONTROLS.do using mfxtsemipar_cv.R
# 9 panels x 2 degrees (bs1 + bs2) = each panel shows True, BS(1), BS(2)
# Plotting with ggplot2, following Stata twoway style from DO file
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(splines2)
  library(splines)
  library(ggplot2)
})

source("/Users/sigma/Library/CloudStorage/SynologyDrive-kuanke/Downloads/jaerevision/mfxtsemipar/R/mfxtsemipar_cv.R")

# =============================================================================
# 1. Generate data (replicate gendata4 from Stata)
# =============================================================================
set.seed(123456)

rho    <- 0.5
lambda <- 1

n_id <- 50; n_year <- 2; n_month <- 12; n_day <- 30

dt <- data.table(
  ID    = rep(1:n_id, times = n_year * n_month * n_day),
  year  = rep(rep(1:n_year, each = n_id), times = n_month * n_day),
  month = rep(rep(1:n_month, each = n_id * n_year), times = n_day),
  day   = rep(1:n_day, each = n_id * n_year * n_month)
)

N <- nrow(dt)

dt[, epsilon := rnorm(N)]
dt[, tm := runif(N) * 38]
dt[, z := sqrt(rho) * tm + sqrt(1 - rho^2) * 38 * runif(N)]

b1 <- -0.78536
b2 <-  0.064746
b3 <- -0.003567
b4 <-  0.000079

dt[, f_x := b1*(tm - 20) + b2*(tm^2 - 20^2) + b3*(tm^3 - 20^3) + b4*(tm^4 - 20^4)]
dt[, g_z := lambda * (b1*(z - 20) + b2*(z^2 - 20^2) + b3*(z^3 - 20^3) + b4*(z^4 - 20^4))]
dt[, y_daily := f_x + g_z + rnorm(N)]

dt[, ai := rnorm(N)]
dt[, nt := rnorm(N)]
dt[, ai := ai[1], by = .(ID, year, month)]
dt[, nt := nt[1], by = .(year, month)]
dt[, y_daily := y_daily + ai + nt]

lf <- dt[, .(ysimulate = sum(y_daily),
             tm_m      = mean(tm),
             z_m       = mean(z)),
          by = .(ID, year, month)]

set.seed(123456)
ids <- unique(lf$ID)
folds <- sample(rep_len(1:2, length(ids)))
fold_map <- data.table(ID = ids, group = folds)
lf <- merge(lf, fold_map, by = "ID")
lf[, insample := as.integer(group == 1)]

setorder(lf, ID, year, month)
lf[, tas := -2 + seq_len(.N) / 10]

dt <- merge(dt, lf[, .(ID, year, month, tas)], by = c("ID", "year", "month"))
dt[, z2 := z^2]
dt[, z3 := z^3]
dt[, z4 := z^4]
dt[, ym := year * 100 + month]
lf[, ym := year * 100 + month]

absorb_spec <- c("ID", "ym")

# =============================================================================
# 2. Define the 9 control strategies (matching DO file comments exactly)
# =============================================================================
strategies <- c(
  "1. Mean(z)",
  "2. z, z^2",
  "3. z, z^2, z^3",
  "4. Spline(Z), knots=2",
  "5. Spline(Z), knots=3",
  "6. Spline(Z), knots=4",
  "7. Bins(Z), bins=5",
  "8. Bins(Z), bins=10",
  "9. z, z^2, z^3, z^4"
)

# results[[i]][[d]] where d=1 for degree=1 (BS1/green), d=2 for degree=2 (BS2/red)
results <- vector("list", 9)
for (i in 1:9) results[[i]] <- list(NULL, NULL)

run_strategy <- function(hf, lf, hfcov, type, degree, minnk, maxnk,
                         label, absorb_vars) {
  res <- tryCatch({
    mfxtsemipar_cv(
      hf     = copy(hf),
      lf     = copy(lf),
      y      = "ysimulate",
      x      = NULL,
      uvar   = "tm",
      id     = "ID",
      tl     = c("year", "month"),
      gen    = "gfit",
      hfcov  = hfcov,
      type   = type,
      degree = degree,
      minnk  = minnk,
      maxnk  = maxnk,
      center = 20,
      absorb = absorb_vars,
      atu    = "tas",
      seed   = 42,
      nfold  = 2,
      cvgroup = "group"
    )
  }, error = function(e) {
    message("Error in ", label, " deg=", degree, ": ", e$message)
    return(NULL)
  })

  if (is.null(res)) return(NULL)

  fit_dt <- res$fitted
  fit_unique <- fit_dt[, .(tas, gfit)]
  fit_unique <- unique(fit_unique, by = "tas")
  setorder(fit_unique, tas)

  list(fit = fit_unique, result = res, label = label, degree = degree)
}

create_z_spline <- function(hf, nknots, degree = 3) {
  z_vals <- hf$z
  zmin <- min(z_vals) - 0.01
  zmax <- max(z_vals) + 0.01
  bknots <- c(zmin, zmax)
  probs <- seq_len(nknots) / (nknots + 1)
  knots <- as.numeric(quantile(z_vals, probs = probs, type = 2, names = FALSE))
  B <- splines2::bSpline(z_vals, knots = knots, Boundary.knots = bknots,
                         degree = degree, intercept = TRUE)
  as.matrix(B)
}

create_z_bins <- function(hf, nbins) {
  z_vals <- hf$z
  cuts <- quantile(z_vals, probs = seq(0, 1, length.out = nbins + 1), type = 2, names = FALSE)
  cuts[1] <- cuts[1] - 0.01
  cuts[length(cuts)] <- cuts[length(cuts)] + 0.01
  bin_idx <- cut(z_vals, breaks = cuts, labels = FALSE)
  bin_cols <- paste0(".zbin", 1:nbins)
  for (k in 1:nbins) {
    hf[, (paste0(".zbin", k)) := as.integer(bin_idx == k)]
  }
  bin_cols
}

cat("\n=== Running 9 strategies x 2 degrees ===\n\n")

for (i in 1:9) {
  cat(sprintf("--- Strategy %d: %s ---\n", i, strategies[i]))

  # Prepare hf for this strategy
  if (i == 4) { hf_i <- copy(dt); B <- create_z_spline(hf_i, 2, 3)
    for (j in seq_len(ncol(B))) hf_i[, paste0(".zsp",j-1L) := B[,j]]
    hfcov_i <- paste0(".zsp", 0:(ncol(B)-1L))
  } else if (i == 5) { hf_i <- copy(dt); B <- create_z_spline(hf_i, 3, 3)
    for (j in seq_len(ncol(B))) hf_i[, paste0(".zsp",j-1L) := B[,j]]
    hfcov_i <- paste0(".zsp", 0:(ncol(B)-1L))
  } else if (i == 6) { hf_i <- copy(dt); B <- create_z_spline(hf_i, 4, 3)
    for (j in seq_len(ncol(B))) hf_i[, paste0(".zsp",j-1L) := B[,j]]
    hfcov_i <- paste0(".zsp", 0:(ncol(B)-1L))
  } else if (i == 7) { hf_i <- copy(dt); hfcov_i <- create_z_bins(hf_i, 5)
  } else if (i == 8) { hf_i <- copy(dt); hfcov_i <- create_z_bins(hf_i, 10)
  } else if (i == 1) { hf_i <- copy(dt); hfcov_i <- "z"
  } else if (i == 2) { hf_i <- copy(dt); hfcov_i <- c("z","z2")
  } else if (i == 3) { hf_i <- copy(dt); hfcov_i <- c("z","z2","z3")
  } else if (i == 9) { hf_i <- copy(dt); hfcov_i <- c("z","z2","z3","z4")
  }

  # Run degree=1 (BS1 / green in Stata)
  results[[i]][[1]] <- run_strategy(hf_i, copy(lf), hfcov_i,
    type="bs", degree=1, minnk=2, maxnk=10,
    label=strategies[i], absorb_vars=absorb_spec)

  # Run degree=2 (BS2 / red in Stata)
  results[[i]][[2]] <- run_strategy(hf_i, copy(lf), hfcov_i,
    type="bs", degree=2, minnk=2, maxnk=10,
    label=strategies[i], absorb_vars=absorb_spec)

  cat("\n")
}

# =============================================================================
# 3. Prepare plotting data — 3 lines per panel: True, BS(1), BS(2)
# =============================================================================
true_g <- function(x) lambda * (b1*(x-20) + b2*(x^2-20^2) + b3*(x^3-20^3) + b4*(x^4-20^4))

tmin <- min(dt$tm, na.rm = TRUE)
tmax <- max(dt$tm, na.rm = TRUE)

true_df <- data.frame(x = seq(tmin, tmax, length.out = 500))
true_df$y <- true_g(true_df$x)

extract_fit <- function(res_obj) {
  est_df <- data.frame(x = NA_real_, y = NA_real_)
  if (!is.null(res_obj) && !is.null(res_obj$fit)) {
    fit <- res_obj$fit
    fit <- fit[tas >= tmin & tas <= tmax]
    if (nrow(fit) > 0) {
      est_df <- data.frame(x = fit$tas, y = fit$gfit)
    }
  }
  est_df
}

# Build long-format: each row has x, y, line_type (True/BS1/BS2), panel
plot_list <- vector("list", 9)
for (i in 1:9) {
  panel_name <- strategies[i]
  fit_bs1 <- extract_fit(results[[i]][[1]])  # degree=1
  fit_bs2 <- extract_fit(results[[i]][[2]])  # degree=2

  plot_list[[i]] <- rbind(
    transform(true_df, line_type = "True",  panel = panel_name),
    transform(fit_bs1, line_type = "BS(1)",  panel = panel_name),
    transform(fit_bs2, line_type = "BS(2)",  panel = panel_name)
  )
}
plot_all <- do.call(rbind, plot_list)
plot_all$line_type <- factor(plot_all$line_type, levels = c("True", "BS(1)", "BS(2)"))
plot_all$panel <- factor(plot_all$panel, levels = strategies)

# =============================================================================
# 4. Plot with ggplot2 — matching Stata twoway from DO file exactly
#    Black=True, Green=BS(1), Red=BS(2)
# =============================================================================
line_colors <- c("True" = "black", "BS(1)" = "#2E8B57", "BS(2)" = "#E41A1C")
line_sizes  <- c("True" = 0.7, "BS(1)" = 0.95, "BS(2)" = 0.95)
linetypes   <- c("True" = "solid", "BS(1)" = "solid", "BS(2)" = "solid")

p <- ggplot(plot_all, aes(x = x, y = y,
                          color = line_type, linewidth = line_type,
                          linetype = line_type)) +
  geom_line() +
  facet_wrap(~ panel, ncol = 3, scales = "free") +
  scale_color_manual(values = line_colors) +
  scale_linewidth_manual(values = line_sizes) +
  scale_linetype_manual(values = linetypes) +
  labs(x = "tm", y = "g(tm)") +
  theme_classic(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 9.5, hjust = 0.03),
    strip.background = element_rect(fill = "gray92", color = "gray80"),
    panel.grid.major.y = element_line(color = "gray88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    plot.margin = margin(12, 8, 8, 8)
  ) +
  guides(
    color = guide_legend(override.aes = list(
      linewidth = c(0.7, 0.95, 0.95),
      linetype = c("solid", "solid", "solid")
    ), order = 1, nrow = 1),
    linewidth = "none",
    linetype = "none"
  )

# Save
ggsave("/Users/sigma/.qclaw/workspace-ua58rsb93veqtxl7/hfcontrols_9panel.pdf",
       p, width = 14, height = 13, device = cairo_pdf)
ggsave("/Users/sigma/.qclaw/workspace-ua58rsb93veqtxl7/hfcontrols_9panel.png",
       p, width = 14, height = 13, dpi = 200)

cat("\n=== Done ===\n")
cat("Output:\n")
cat("  PDF: hfcontrols_9panel.pdf\n")
cat("  PNG: hfcontrols_9panel.png\n")

# Summary table
cat("\n=== Strategy Summary (BS(1) vs BS(2)) ===\n")
cat(sprintf("  %-28s | %8s | %8s\n", "Strategy", "BS(1) MSE", "BS(2) MSE"))
cat(sprintf("  %-28s | %8s | %8s\n", "----------------------------", "--------", "--------"))
for (i in 1:9) {
  mse1 <- if (!is.null(results[[i]][[1]]$result)) results[[i]][[1]]$result$minmse else NA
  mse2 <- if (!is.null(results[[i]][[2]]$result)) results[[i]][[2]]$result$minmse else NA
  cat(sprintf("  %-28s | %8.4f | %8.4f\n", strategies[i], mse1, mse2))
}
