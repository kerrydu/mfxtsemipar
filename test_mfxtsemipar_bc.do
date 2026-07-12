* Test file for mfxtsemipar_bc
* July 12, 2026

clear all
set more off
set seed 12345

* add the directory containing mfxtsemipar_bc.ado to the adopath
adopath + "/Users/sigma/SynologyDrive/kuanke/Downloads/jaerevision/mfxtsemipar"

*------------------------------------------------------------------------------
* Generate mixed-frequency panel data
*------------------------------------------------------------------------------
set obs 2000
gen id = ceil(_n / 20)
bys id: gen t = ceil(_n / 2)

* high-frequency uvar
gen uvar = runiform() * 4 - 2

* low-frequency covariate (constant within id-t)
gen x = rnormal()

* true g(u) depends on the LF mean of uvar
bys id t: egen ubar = mean(uvar)
gen g_true = ubar^2

* dependent variable at LF level (constant within id-t)
tempvar noise y0
gen `noise' = rnormal()
bys id t: gen `y0' = g_true[1] + 0.5 * x[1] + `noise'[1]
rename `y0' y

* drop helpers
drop ubar g_true `noise'

*------------------------------------------------------------------------------
* Test 1: basic mfxtsemipar_bc without UCB
*------------------------------------------------------------------------------
di _n "=== Test 1: basic mfxtsemipar_bc ==="
mfxtsemipar_bc y x, uvar(uvar) id(id) gen(g) tl(t) type(poly) maxnk(4) minnk(2)

sum g_raw g g_se g_bc_se
di "e(nknots) = " e(nknots) ", e(min_cv_mse) = " e(min_cv_mse) ", e(rmse) = " e(rmse)

*------------------------------------------------------------------------------
* Test 2: mfxtsemipar_bc with UCB
*------------------------------------------------------------------------------
di _n "=== Test 2: mfxtsemipar_bc with UCB ==="
mfxtsemipar_bc y x, uvar(uvar) id(id) gen(g2) tl(t) type(poly) maxnk(4) minnk(2) ///
    ucb ucbgrid(-2(0.1)2) ucbsim(500)

sum g2_raw g2 g2_se g2_bc_se g2_lb g2_ub
di "e(ucb_crit) = " e(ucb_crit) ", e(ucb_level) = " e(ucb_level)

*------------------------------------------------------------------------------
* Test 3: mfxtsemipar_bc with bootstrap SE
*------------------------------------------------------------------------------
di _n "=== Test 3: mfxtsemipar_bc with bootstrap SE ==="
mfxtsemipar_bc y x, uvar(uvar) id(id) gen(g3) tl(t) type(poly) maxnk(4) minnk(2) ///
    brep(50)

sum g3_raw g3 g3_se g3_bc_se g3_bc_boot_se

*------------------------------------------------------------------------------
* Test 4: mfxtsemipar_bc with atu evaluation and bcnknots multiplier
*------------------------------------------------------------------------------
di _n "=== Test 4: mfxtsemipar_bc with atu and bcnknots(1.5k) ==="
gen atu = uvar + 0.1
mfxtsemipar_bc y x, uvar(uvar) id(id) gen(g4) tl(t) type(poly) maxnk(4) minnk(2) ///
    atu(atu) bcnknots(1.5k) bcdegree(3)

sum g4_raw g4 g4_se g4_bc_se

di _n "All tests completed successfully."
