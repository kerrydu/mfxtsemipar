*! mfxtsemipar_bc - Mixed-frequency semiparametric regression with robust bias correction
*! version 1.0 - July 12, 2026
*!
*! Syntax:
*!   mfxtsemipar_bc depvar [indepvars] [if] [in] [weight],
*!                    uvar(varname) id(varname) gen(string) tl(varlist)
*!                    [, options]
*!
*! Description:
*!   Mixed-frequency cross-validated semiparametric regression with a
*!   joint-orthogonalized robust bias correction (RBC) for the nonparametric
*!   component. The command follows the implementation style of Cattaneo,
*!   Farrell and Feng (2020): the bias-correction basis is projected off the
*!   main basis, and the main and orthogonalized bias-correction bases are
*!   estimated jointly. Uniform confidence bands (UCB) for the bias-corrected
*!   curve are optional.
*!
*! Required options:
*!   uvar(varname)     - high-frequency semiparametric variable
*!   id(varname)       - panel identifier
*!   gen(string)       - name stub for generated fitted-value/SE variables
*!   tl(varlist)       - time-level variables for collapse to low frequency
*!
*! Key options:
*!   cluster(string)   - cluster variable for robust SEs
*!   bknots(numlist)   - boundary knots (must cover range of uvar)
*!   type(string)      - spline type: poly, bs, ms, is, ibs (default: poly)
*!   winsor(string)    - winsorization percentiles for uvar
*!   WINSORValues      - winsor() contains raw values, not percentiles
*!   EQSPACE           - equally spaced interior knots
*!   MAXNK(integer 5)  - maximum number of knots in CV
*!   MINNK(integer 2)  - minimum number of knots in CV
*!   center(numlist)   - centering value for spline basis
*!   Absorb(string)    - fixed-effects specification for reghdfe
*!   cvgroup(varname)  - variable defining CV folds
*!   nfold(integer 10) - number of CV folds when cvgroup() absent
*!   seed(numlist)     - seed for fold generation
*!   degree(integer 1) - polynomial/spline degree
*!   keepsplines       - keep generated spline variables
*!   atu(varname)      - alternative evaluation variable
*!   hfcov(varlist)    - high-frequency covariates
*!   dropfirstbase     - drop intercept/base term from spline basis
*!   sopt              - select first local minimum of CV curve
*!   brep(real 0)      - wild-bootstrap replications for BC inference
*!   predy(name)       - low-frequency full prediction variable
*!   PARTIALOUT        - partial out all covariates
*!   PARTIALOUT1(varlist) - partial out selected covariates
*!   ucb               - compute uniform confidence band
*!   ucblevel(real 95) - UCB coverage level
*!   ucbgrid(numlist)  - grid points for UCB
*!   ucbsim(integer 2000) - simulations for sup-t critical value
*!   bcdegree(integer) - degree of BC basis (default: degree+1)
*!   bcnknots(string)  - BC knots: integer, or e.g. "1.5k" for multiplier
*!
*! Generated variables:
*!   <gen>_raw        - raw spline prediction (main model)
*!   <gen>            - bias-corrected prediction
*!   <gen>_se         - raw standard error
*!   <gen>_bc_se      - bias-corrected standard error
*!   <gen>_bc_boot_se - bootstrap BC SE (if brep>0)
*!   <gen>_lb, <gen>_ub - UCB bounds (if ucb)
*!
*! Returned results:
*!   e(soptnk)        - first-local-minimum knot count
*!   e(knots)         - selected interior knots
*!   e(nknots)        - selected number of interior knots
*!   e(min_cv_mse)    - minimum CV RMSE
*!   e(cv_mse)        - CV RMSE matrix
*!   e(rmse)          - in-sample RMSE of main model
*!   e(splinecmd)     - command used to generate main splines
*!   e(info)          - model information criteria
*!   e(ucb_crit)      - sup-t critical value (if ucb)
*!   e(ucb_level)     - UCB level (if ucb)

*! Sunday, July 12, 2026

cap program drop mfxtsemipar_bc
program define mfxtsemipar_bc, eclass

version 16

* bcdegree() conflicts with degree() in Stata's syntax parser (one option
* name is a suffix of the other). Parse it manually from the full command
* line and remove it before invoking syntax.
local bcdegree ""
local fullcmd `"`0'"'
if regexm(`"`fullcmd'"', "bcdegree\s*\(\s*([0-9]+)\s*\)") {
    local bcdegree = regexs(1)
    local fullcmd = regexr(`"`fullcmd'"', "bcdegree\s*\(\s*[0-9]+\s*\)", "")
    local 0 `"`fullcmd'"'
}

syntax varlist(min=1) [if] [in] [fw aw pw/], ///
                            uvar(varname) ///
                            id(varname) ///
                            gen(string) ///
                            tl(varlist) ///
                            [cluster(string) ///
                             bknots(numlist ascending min=2 max=2) ///
                             type(string) ///
                             winsor(string) ///
                             WINSORValues ///
                             EQSPACE ///
                             MAXNK(integer 5) ///
                             MINNK(integer 2) ///
                             center(numlist max=1) ///
                             Absorb(string) ///
                             cvgroup(varname) ///
                             nfold(integer 10) ///
                             seed(numlist) ///
                             degree(integer 1) ///
                             keepsplines ///
                             atu(varname) ///
                             hfcov(varlist) ///
                             dropfirstbase ///
                             sopt ///
                             brep(real 0) ///
                             predy(name) ///
                             PARTIALOUT ///
                             PARTIALOUT1(varlist) ///
                             ucb ///
                             ucblevel(real 95) ///
                             ucbgrid(numlist) ///
                             ucbsim(integer 2000) ///
                             bcnknots(string) ///
                             startp(numlist max=1) ///
                             endp(numlist max=1)]

* dependencies
cap which gensplines
if _rc ssc install gensplines
cap which savesome
if _rc ssc install savesome
cap which reghdfe
if _rc ssc install reghdfe
cap which ftools
if _rc ssc install ftools
cap which splitsample
if _rc ssc install splitsample

* predy validation
if `"`predy'"' != "" checkpredy `predy'
local xbd `r(xbd)'
local predy `r(predy)'

* weight variable
if (`"`weight'"' != "") {
    tempvar weightvar
    qui gen `weightvar' = `exp'
    local weightexp [`weight'=`weightvar']
}

* intercept / dropfirstbase
if "`dropfirstbase'" != "" {
    local intercept
}
else local intercept intercept

* parse varlist and partialout
local varlist0 `varlist'
gettoken ydep varlist0: varlist0
local varlist0 `varlist0' `hfcov'
if "`partialout'" != "" & `"`partialout1'"' == "" {
    local partialout1 `varlist0'
    local varlist `ydep'
    local varlist0
}
if `"`partialout1'"' != "" {
    local varlist0 : list varlist0 - partialout1
    local varlist `ydep' `varlist0'
}
local partialout `partialout1'

* spline type and command
local gensplines gensplines
if "`type'" == "" | "`type'" == "poly" {
    local type poly
    local gensplines polysplines
}
if inlist("`type'", "bs", "ms", "is", "ibs", "poly") == 0 {
    di as err "type() should be one of bs, ms, is, ibs and poly"
    exit
}

* absorb variables that need to be aggregated (exclude id and tl)
if (`"`absorb'"' != "") {
    local absorb0 `absorb'
    gettoken absorb1 absorb0 : absorb0, p(#)
    local absorbvars `absorb1'
    while "`absorb0'" != "" {
        gettoken absorb1 absorb0 : absorb0, p(#)
        if "`absorb1'" != "#" local absorbvars `absorbvars' `absorb1'
    }
    local absorbvars : list uniq absorbvars
    local absorbvars : list absorbvars - tl
    local absorbvars : list absorbvars - id
}

* sample marker
marksample touse
markout `touse' `uvar' `id' `tl' `hfcov', strok

local xvar `uvar'

* winsorization of uvar
if "`winsor'" != "" {
    if strpos("`winsor'", ",") == 0 local winsor `winsor',
    if "`winsorvalues'" != "" local valuesopt values
    GetWinsorOpts `winsor' xvar(`xvar') touse(`touse') `valuesopt'
    if "`xvar'" != "" {
        tempvar xvar_winsor
        local vartype: type `xvar'
        gen `vartype' `xvar_winsor' = cond(`xvar' <= `winsor_low', ///
                                        `winsor_low', ///
                                        cond(`xvar' >= `winsor_high', `winsor_high', `xvar'))
        local xvar `xvar_winsor'
    }
}

* boundary knots (with validation)
qui su `xvar' if `touse'
local rawmin = r(min)
local rawmax = r(max)
local min = `rawmin' - 0.01
local max = `rawmax' + 0.01

if "`bknots'" != "" {
    local bkmn : word 1 of `bknots'
    local bkmx : word 2 of `bknots'
    if `rawmin' < `bkmn' {
        di as err "Minimum uvar value (" `rawmin' ") is below the lower boundary knot (" `bkmn' ")."
        exit 498
    }
    if `rawmax' > `bkmx' {
        di as err "Maximum uvar value (" `rawmax' ") is above the upper boundary knot (" `bkmx' ")."
        exit 498
    }
}
if "`bknots'" == "" {
    local bknots `min' `max'
}

if `"`cluster'"' != "" {
    local setype cluster(`cluster')
}

* CV groups
if `"`cvgroup'"' != "" {
    tempvar cv
    qui egen `cv' = group(`cvgroup')
}
else {
    tempvar cv
    if "`seed'" != "" qui set seed `seed'
    qui splitsample, cluster(`id') nsplit(`nfold') gen(`cv')
}

* startp / endp options for gennknots
local startpop = cond("`startp'" == "", "", "startp(`startp')")
local endpop   = cond("`endp'" == "", "", "endp(`endp')")

* clone uvar for safe restores
tempvar uvar2
clonevar `uvar2' = `uvar'

*===============================================================================
* 1. Cross-validation over nk
*===============================================================================
mat mse = J(`minnk' - 1, 1, .)

forv i = `minnk' / `maxnk' {
    preserve
    gennknots `xvar' if `touse', nknots(`i') `eqspace' `startpop' `endpop'
    local knots `r(knots)'
    `gensplines' `xvar', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type')
    local allbins `r(splinevarlist)'

    qui collapse (mean) `ydep' `varlist0' `partialout1' `cv' `absorbvars' `weightvar' (sum) `allbins' if `touse', by(`id' `tl')

    qui rmse_cv `ydep' `allbins' `varlist0' `weightexp', cv(`cv') partialout(`partialout1') absorb(`absorb')
    local rmse = r(rmse)
    mat mse = mse \ `rmse'
    local mserowname "`mserowname' nk=`i'"
    restore
}

minpos mse
local nknots = r(pos)
local min = r(min)

* simple-optimal (first local minimum)
local soptnk = `minnk'
local mse0 1e9
forv j = `minnk' / `maxnk' {
    local binj = mse[`j', 1]
    if `binj' != . & `binj' < `mse0' {
        local mse0 = `binj'
        local soptnk = `j'
    }
    else continue, break
}

mat mse = mse[`minnk'..`maxnk', 1]
mat rownames mse = `mserowname'
mat colnames mse = "cv_rmse"
di _n "Cross-validation RMSE (for knot selection)"
matlist mse
di _n

if "`sopt'" != "" local nknots `soptnk'

*===============================================================================
* 2. Final main regression with optimal knots
*===============================================================================
gennknots `xvar' if `touse', nknots(`nknots') `eqspace' `startpop' `endpop'
local knots `r(knots)'

`gensplines' `xvar', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type')
local allbins `r(splinevarlist)'
local splinecmd `gensplines' `xvar', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type')

preserve
qui collapse (mean) `ydep' `varlist0' `partialout1' `absorbvars' `weightvar' (sum) `allbins' if `touse', by(`id' `tl')

tempvar res0
qui reghdfe `ydep' `varlist0' `partialout1' `allbins' `weightexp', ///
            absorb(`absorb') `setype' residuals(`res0')

tempvar __rmse_sq
qui gen double `__rmse_sq' = `res0'^2
qui summarize `__rmse_sq'
local rmse = sqrt(r(mean))
cap drop `__rmse_sq'

if `"`predy'"' != "" {
    qui predict `predy', `xbd'
    tempfile predy_file
    qui savesome `id' `tl' `predy' using `predy_file', replace
}

qui estimate store mfxtsemipar_bc_main
mat b_main = e(b)
mat V_main = e(V)
qui estat ic, all
tempname info_main
mat `info_main' = r(S)

* retained main spline coefficients
local allcoefs : colfullnames e(b)
local retained_main : list allcoefs & allbins
restore

if "`retained_main'" == "" {
    di as err "All main spline coefficients were dropped; check collinearity."
    exit 498
}

*===============================================================================
* 3. Robust bias correction basis and joint regression
*===============================================================================
* BC degree
if "`bcdegree'" == "" local bcdegree = `degree' + 1

* BC knots
if "`bcnknots'" == "" {
    local bc_nk `nknots'
}
else {
    local lastchr = substr("`bcnknots'", length("`bcnknots'"), 1)
    if inlist("`lastchr'", "k", "K") {
        local mult = substr("`bcnknots'", 1, length("`bcnknots'") - 1)
        local bc_nk = round(`mult' * `nknots')
        if `bc_nk' < 1 local bc_nk 1
    }
    else {
        cap confirm integer number `bcnknots'
        if _rc {
            di as err "bcnknots() must be a positive integer or a multiplier ending in k (e.g. 1.5k)"
            exit 198
        }
        local bc_nk `bcnknots'
    }
}

gennknots `xvar' if `touse', nknots(`bc_nk') `eqspace' `startpop' `endpop'
local bc_knots `r(knots)'

`gensplines' `xvar', gen(__SplineBC_) knots(`bc_knots') bknots(`bknots') degree(`bcdegree') centerv(`center') `intercept' type(`type')
local bc_spline_vars `r(splinevarlist)'

preserve
qui collapse (mean) `ydep' `varlist0' `partialout1' `absorbvars' `weightvar' (sum) `allbins' `bc_spline_vars' if `touse', by(`id' `tl')

* compute projection P = (B'WB)^{-1} B'W Btilde on complete cases
local bc_orth_vars
foreach v of local bc_spline_vars {
    local vorth `v'_orth
    qui gen double `vorth' = .
    local bc_orth_vars `bc_orth_vars' `vorth'
}

mata: mfxt_bc_proj("`retained_main'", "`bc_spline_vars'", "`weightvar'", "__mfxt_bc_G0", "__mfxt_bc_M", "__mfxt_bc_P")

* orthogonalize BC basis in LF data
mata: mfxt_bc_orth("`retained_main'", "`bc_spline_vars'", "__mfxt_bc_P", "`bc_orth_vars'")

* joint regression: main + orthogonalized BC
tempvar res_joint
qui reghdfe `ydep' `varlist0' `partialout1' `retained_main' `bc_orth_vars' `weightexp', ///
            absorb(`absorb') `setype' residuals(`res_joint')

qui estimate store mfxtsemipar_bc_joint
mat b_joint = e(b)
mat V_joint = e(V)
qui estat ic, all
tempname info_joint
mat `info_joint' = r(S)

* retained BC orthogonal variables in joint regression
local allcoefs : colfullnames e(b)
local retained_bc_orth : list allcoefs & bc_orth_vars

* corresponding raw BC variables
local retained_bc
foreach v of local retained_bc_orth {
    local vraw = subinstr("`v'", "_orth", "", 1)
    local retained_bc `retained_bc' `vraw'
}

if "`retained_bc_orth'" == "" {
    restore
    di as err "All bias-correction spline coefficients were dropped; check collinearity."
    exit 498
}

*--------------------------------------------------------------------------
* wild bootstrap VCE for joint regression (if requested)
*--------------------------------------------------------------------------
if `brep' > 0 {
    tempvar yhat ehat ystar
    qui predict double `yhat', xbd
    qui gen double `ehat' = `ydep' - `yhat'
    qui gen double `ystar' = .
    local k_joint = colsof(b_joint)
    mata: bb = J(`brep', `k_joint', .)
    forv b = 1 / `brep' {
        qui wildboot `ystar' `yhat' `ehat', cluster(`cluster')
        qui reghdfe `ystar' `varlist0' `partialout1' `retained_main' `bc_orth_vars' `weightexp', ///
                    absorb(`absorb') `setype'
        mat bi = e(b)
        mata: bb[`b', .] = st_matrix("bi")
    }
    mata: bb = quadvariance(bb)
    mata: st_matrix("V_bc_boot", bb)
    local bjoint_names : colnames b_joint
    mat rownames V_bc_boot = `bjoint_names'
    mat colnames V_bc_boot = `bjoint_names'
}

restore

*===============================================================================
* 4. Prediction at high frequency
*===============================================================================
* generate evaluation splines (atu if requested, otherwise xvar/uvar)
if "`atu'" != "" {
    cap drop `allbins'
    cap drop `bc_spline_vars'
    `gensplines' `atu', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type')
    local allbins `r(splinevarlist)'
    `gensplines' `atu', gen(__SplineBC_) knots(`bc_knots') bknots(`bknots') degree(`bcdegree') centerv(`center') `intercept' type(`type')
    local bc_spline_vars `r(splinevarlist)'
}

* confirm / drop generated variables
cap drop `gen'_raw
cap drop `gen'
cap drop `gen'_se
cap drop `gen'_bc_se
cap drop `gen'_bc_boot_se
cap drop `gen'_lb
cap drop `gen'_ub
qui gen double `gen'_raw = .
qui gen double `gen' = .
qui gen double `gen'_se = .
qui gen double `gen'_bc_se = .
qui gen double `gen'_lb = .
qui gen double `gen'_ub = .

mata: mfxt_bc_predict("`gen'", "`retained_main'", "`retained_bc_orth'", ///
                      "__mfxt_bc_P", "b_main", "V_main", "b_joint", "V_joint", "_orth")

if `brep' > 0 {
    cap drop `gen'_bc_boot_se
    qui gen double `gen'_bc_boot_se = .
    mata: mfxt_bc_boot_se("`gen'", "`retained_main'", "`retained_bc_orth'", ///
                          "__mfxt_bc_P", "V_bc_boot", "_orth")
}

*===============================================================================
* 5. Uniform confidence band
*===============================================================================
if "`ucb'" != "" {
    local evalvar `xvar'
    if "`atu'" != "" local evalvar `atu'
    qui su `evalvar' if `touse'
    local umin = r(min)
    local umax = r(max)
    if "`ucbgrid'" != "" {
        local ucbgridpts `ucbgrid'
    }
    else {
        local ngrid = 200
        local step = (`umax' - `umin') / (`ngrid' - 1)
        qui numlist "`umin'(`step')`umax'", sort
        local ucbgridpts `r(numlist)'
    }

    preserve
    clear
    local ngrid : word count `ucbgridpts'
    qui set obs `ngrid'
    gen double __ucb_u = .
    tokenize `ucbgridpts'
    forv i = 1 / `ngrid' {
        qui replace __ucb_u = ``i'' in `i'
    }
    `gensplines' __ucb_u, gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type')
    `gensplines' __ucb_u, gen(__SplineBC_) knots(`bc_knots') bknots(`bknots') degree(`bcdegree') centerv(`center') `intercept' type(`type')

    local bc_orth_gridvars
    foreach v of local retained_bc_orth {
        local vgrid `v'
        qui gen double `vgrid' = .
        local bc_orth_gridvars `bc_orth_gridvars' `vgrid'
    }
    mata: mfxt_bc_orth_grid("`retained_main'", "`retained_bc'", "__mfxt_bc_P", "`bc_orth_gridvars'")

    * drop unused spline columns to avoid passing them to UCB routine
    local droplist : list allbins - retained_main
    foreach v of local droplist {
        cap drop `v'
    }
    local droplist2 : list bc_spline_vars - retained_bc
    foreach v of local droplist2 {
        cap drop `v'
    }

    mata: mfxt_ucb_crit("b_joint", "V_joint", "`retained_main' `bc_orth_gridvars'", `ucblevel' / 100, `ucbsim', "ucb_crit")
    restore
    local crit = ucb_crit
    qui replace `gen'_lb = `gen' - `crit' * `gen'_bc_se
    qui replace `gen'_ub = `gen' + `crit' * `gen'_bc_se
}

*===============================================================================
* 6. Cleanup and return
*===============================================================================
if "`keepsplines'" == "" {
    cap drop `allbins'
    cap drop `bc_spline_vars'
}

* merge predy if requested
if `"`predy'"' != "" {
    qui merge m:1 `id' `tl' using `predy_file', nogen
}

* restore joint estimation results for ereturn
qui estimate restore mfxtsemipar_bc_joint

ereturn scalar soptnk = `soptnk'
ereturn local knots `knots'
ereturn scalar nknots = `nknots'
ereturn scalar min_cv_mse = `min'
ereturn matrix cv_mse = mse
ereturn scalar rmse = `rmse'
ereturn local splinecmd `splinecmd'
ereturn matrix info = `info_joint'

* extra BC information
ereturn local bc_knots `bc_knots'
ereturn scalar bc_nknots = `bc_nk'
ereturn scalar bc_degree = `bcdegree'
ereturn matrix bc_P = __mfxt_bc_P

if "`ucb'" != "" {
    ereturn scalar ucb_crit = `crit'
    ereturn scalar ucb_level = `ucblevel'
}

end


*===============================================================================
* Helper programs
*===============================================================================

cap program drop polysplines
program define polysplines, rclass
version 14
syntax varlist(min=1 max=1), knots(numlist) gen(string) [degree(integer 1) centerv(integer 0) *]

local power `degree'
local xvar `varlist'
local j = 1
forv p = 1 / `power' {
    qui gen double `gen'`j' = `xvar'^`p' - `centerv'^`p'
    local splinevarlist `splinevarlist' `gen'`j'
    local j = `j' + 1
}

foreach num of numlist `knots' {
    qui gen double `gen'`j' = (`xvar' - `num')^`power' * (`xvar' > `num') - (`centerv' - `num')^`power' * (`centerv' > `num')
    local splinevarlist `splinevarlist' `gen'`j'
    local j = `j' + 1
}

return local knots `knots'
return local splinevarlist `splinevarlist'
end


cap program drop gennknots
program define gennknots, rclass
version 16
syntax varlist(min=1 max=1) [if] [in], nknots(integer) [eqspace startp(numlist) endp(numlist)]

marksample touse
local uvar `varlist'
qui su `varlist' if `touse'
local nbin = `nknots' + 1
if "`eqspace'" != "" {
    if "`startp'" != "" local min = `startp'
    else local min = r(min)
    if "`endp'" != "" local max = `endp'
    else local max = r(max)
    local step = (`max' - `min') / `nbin'
    qui numlist "`=`min'+`step''(`step')`=`max'-`step''", sort
    local cutpoints `r(numlist)'
}
else {
    if "`startp'" != "" & "`endp'" != "" local ifcond if `uvar' >= `startp' & `uvar' < `endp'
    else if "`endp'" != "" {
        local ifcond if `uvar' < `endp'
    }
    else if "`startp'" != "" {
        local ifcond if `uvar' >= `startp'
    }
    else {
        local ifcond
    }

    qui _pctile `varlist' `ifcond', nq(`nbin')

    forv i = 1 / `=`nbin'-1' {
        local cutpoints `cutpoints' `=r(r`i')'
    }
}

return local knots `cutpoints'
end


program define GetWinsorOpts
syntax anything [if][in], xvar(string) [values] touse(varname)
marksample touse
numlist "`anything'", ascending min(2) max(2)

capture confirm number `xvar'
if !_rc & "`values'" == "" {
    di as error "If using the winsor() option with a scalar you need to use the values option."
    exit 198
}

if "`values'" == "" {
    _pctile `xvar' if `touse' `fw', percentiles(`r(numlist)')
    c_local winsor_low  `r(r1)'
    c_local winsor_high `r(r2)'
}
else {
    tokenize `r(numlist)'
    c_local winsor_low  `1'
    c_local winsor_high `2'
}
end


cap program drop rmse_cv
program define rmse_cv, rclass
version 16.0
syntax varlist [fw aw pw/], [CV(varname) n(integer 10) SEED(integer 1234) opt(string) cluster(varname) LOGLik absorb(string) PARTIALOUT(varlist)]

if ("`weight'" != "") local weightexp [`weight'=`exp']

tempvar gid resi2 gjhat valy
qui gen double `resi2' = .
qui gen double `gjhat' = .
if "`cv'" != "" {
    qui egen int `gid' = group(`cv')
}
else {
    set seed `seed'
    if "`cluster'" != "" {
        qui splitsample, cluster(`cluster') nsplit(`n') gen(`gid')
    }
    else {
        qui splitsample, nsplit(`n') gen(`gid')
    }
}
qui su `gid'
local numgid = r(max)

gettoken depvar varlist: varlist
local nvars : word count `varlist'

forv i = 1 / `numgid' {
    qui reghdfe `depvar' `varlist' `partialout' `weightexp' if `gid' != `i', a(`absorb')
    mat b = e(b)
    mat b = b[1, 1..`nvars']

    qui computeghat `gjhat' `varlist', bmat(b)
    qui replace `gjhat' = `depvar' - `gjhat' if `gid' == `i'
    qui reghdfe `gjhat' `partialout' `weightexp' if `gid' == `i', a(`absorb') residuals(`valy')
    qui replace `resi2' = (`valy')^2 if `gid' == `i'
    cap drop `valy'
}

qui su `resi2', meanonly
return scalar rmse = sqrt(r(mean))
end


cap program drop computeghat
program define computeghat
version 16
syntax varlist, bmat(name)
gettoken ydep varlist: varlist
mata: b = st_matrix("`bmat'")
mata: st_view(sdata = ., ., "`varlist'")
mata: st_view(ydep = ., ., "`ydep'")
mata: ydep[., .] = sdata * (b')
end


cap program drop minpos
program define minpos, rclass
version 14
args m

confirm matrix `m'

mata: __mm__ = st_matrix("`m'")
mata: __mm__ = (range(1, rows(__mm__), 1), __mm__)
mata: __mm__ = sort(__mm__, 2)
mata: st_numscalar("r(pos)", __mm__[1, 1])
mata: st_numscalar("r(min)", __mm__[1, 2])
return scalar pos = r(pos)
return scalar min = r(min)
end


cap program drop wildboot
program define wildboot
version 16
syntax varlist(min=3 max=3), [cluster(varlist)]

local ystar : word 1 of `varlist'
local yhat : word 2 of `varlist'
local ehat : word 3 of `varlist'

tempvar radw rn
qui gen double `rn' = runiform()
if "`cluster'" == "" {
    qui gen double `radw' = cond(`rn' <= 0.5, 1, -1) * `ehat'
}
else {
    qui bys `cluster': gen double `radw' = cond(`rn'[1] <= 0.5, 1, -1) * `ehat'
}
qui replace `ystar' = `yhat' + `radw'
end


program define checkpredy, rclass
syntax varlist, [xbd xb]
confirm new var `varlist'
return local xbd `xbd'
return local predy `varlist'
end


*===============================================================================
* Mata functions
*===============================================================================
version 16
mata:

real rowvector name2idx(string scalar mname, string rowvector names)
{
    string matrix stripe
    real scalar i, j, k
    real vector idx
    stripe = st_matrixcolstripe(mname)
    k = cols(names)
    idx = J(1, k, .)
    for (i = 1; i <= k; i++) {
        for (j = 1; j <= rows(stripe); j++) {
            if (stripe[j, 2] == names[i]) {
                idx[i] = j
                break
            }
        }
    }
    return(idx)
}

void mfxt_bc_proj(string scalar mainvars_s, string scalar bcvars_s,
                  string scalar wvar_s, string scalar G0name,
                  string scalar Mname, string scalar Pname)
{
    real matrix B, Bt, Bw, Btw
    real colvector w
    real scalar sw
    string rowvector mainvars, bcvars

    mainvars = tokens(mainvars_s)
    bcvars   = tokens(bcvars_s)

    B  = st_data(., mainvars)
    Bt = st_data(., bcvars)
    if (wvar_s != "") w = st_data(., wvar_s)
    else w = J(rows(B), 1, 1)

    real colvector cc
    cc = rowmissing(B) :== 0 :& rowmissing(Bt) :== 0 :& w :< . :& w :> 0
    B  = select(B, cc)
    Bt = select(Bt, cc)
    w  = select(w, cc)

    sw = sum(w)
    if (sw <= 0) {
        printf("No valid observations for BC projection\n")
        exit(2000)
    }

    Bw  = B  :* sqrt(w)
    Btw = Bt :* sqrt(w)

    real matrix G0, M, P
    G0 = invsym(cross(Bw, Bw) / sw)
    M  = cross(Bw, Btw) / sw
    P  = G0 * M

    st_matrix(G0name, G0)
    st_matrix(Mname, M)
    st_matrix(Pname, P)
    st_matrixrowstripe(Pname, (J(length(mainvars), 1, ""), mainvars'))
    st_matrixcolstripe(Pname, (J(length(bcvars), 1, ""), bcvars'))
}

void mfxt_bc_orth(string scalar mainvars_s, string scalar bcvars_s,
                  string scalar Pname, string scalar outvars_s)
{
    real matrix B, Bt, P, Btorth
    string rowvector outvars

    B  = st_data(., tokens(mainvars_s))
    Bt = st_data(., tokens(bcvars_s))
    P  = st_matrix(Pname)
    outvars = tokens(outvars_s)

    Btorth = Bt - B * P
    st_store(., outvars, Btorth)
}

void mfxt_bc_orth_grid(string scalar mainvars_s, string scalar bcvars_s,
                       string scalar Pname, string scalar outvars_s)
{
    real matrix B, Bt, P, Btorth
    real rowvector idx
    string rowvector bcvars, outvars

    B  = st_data(., tokens(mainvars_s))
    bcvars = tokens(bcvars_s)
    Bt = st_data(., bcvars)
    P  = st_matrix(Pname)
    idx = name2idx(Pname, bcvars)
    outvars = tokens(outvars_s)

    Btorth = Bt - B * P[., idx]
    st_store(., outvars, Btorth)
}

void mfxt_bc_predict(string scalar gen, string scalar retained_main_s,
                     string scalar retained_bc_orth_s, string scalar Pname,
                     string scalar bmainname, string scalar Vmainname,
                     string scalar bjointname, string scalar Vjointname,
                     string scalar orth_suffix)
{
    string rowvector rmain, rbc_orth, rbc
    real matrix B, Bt, P, Psub, Btorth, wbc
    real matrix Vmain, Vjoint, Vmain_r, Vjoint_r
    real rowvector bmain, bjoint, bmain_r, bjoint_r
    real colvector fit_raw, se_raw, fit_bc, se_bc
    real rowvector idx_main_b, idx_main_V, idx_bc_b, idx_joint_main_b, idx_joint_all_V, idx_Pcols
    real scalar eps
    string scalar gen_raw, gen_se, gen_bc_se

    rmain     = tokens(retained_main_s)
    rbc_orth  = tokens(retained_bc_orth_s)
    rbc       = subinstr(rbc_orth, orth_suffix, "", 1)

    B  = st_data(., rmain)
    Bt = st_data(., rbc)

    bmain  = st_matrix(bmainname)
    Vmain  = st_matrix(Vmainname)
    bjoint = st_matrix(bjointname)
    Vjoint = st_matrix(Vjointname)
    P      = st_matrix(Pname)

    idx_main_b       = name2idx(bmainname, rmain)
    idx_main_V       = name2idx(Vmainname, rmain)
    idx_bc_b         = name2idx(bjointname, rbc_orth)
    idx_joint_main_b = name2idx(bjointname, rmain)
    idx_joint_all_V  = idx_joint_main_b, idx_bc_b
    idx_Pcols        = name2idx(Pname, rbc)

    if (missing(idx_main_b) > 0 | missing(idx_main_V) > 0 | missing(idx_bc_b) > 0 | missing(idx_joint_main_b) > 0 | missing(idx_Pcols) > 0) {
        printf("Name matching failed in mfxt_bc_predict\n")
        exit(3498)
    }

    bmain_r  = bmain[1, idx_main_b]
    Vmain_r  = Vmain[idx_main_V, idx_main_V]
    bjoint_r = bjoint[1, idx_joint_main_b], bjoint[1, idx_bc_b]
    Vjoint_r = Vjoint[idx_joint_all_V, idx_joint_all_V]
    Psub     = P[., idx_Pcols]

    Btorth = Bt - B * Psub
    wbc    = B, Btorth

    fit_raw = B * bmain_r'
    se_raw  = sqrt(rowsum((B * Vmain_r) :* B))
    fit_bc  = wbc * bjoint_r'
    se_bc   = sqrt(rowsum((wbc * Vjoint_r) :* wbc))

    eps = epsilon(1)
    se_raw = se_raw :* (se_raw :> eps) :+ eps :* (se_raw :<= eps)
    se_bc  = se_bc  :* (se_bc  :> eps) :+ eps :* (se_bc  :<= eps)

    gen_raw   = gen + "_raw"
    gen_se    = gen + "_se"
    gen_bc_se = gen + "_bc_se"

    st_store(., gen_raw, fit_raw)
    st_store(., gen, fit_bc)
    st_store(., gen_se, se_raw)
    st_store(., gen_bc_se, se_bc)
}

void mfxt_bc_boot_se(string scalar gen, string scalar retained_main_s,
                     string scalar retained_bc_orth_s, string scalar Pname,
                     string scalar Vbootname, string scalar orth_suffix)
{
    string rowvector rmain, rbc_orth, rbc
    real matrix B, Bt, P, Psub, Btorth, wbc, Vboot, Vboot_r
    real rowvector idx_boot_main, idx_boot_bc, idx_all, idx_Pcols
    real colvector se_boot
    string scalar gen_bc_boot_se

    rmain     = tokens(retained_main_s)
    rbc_orth  = tokens(retained_bc_orth_s)
    rbc       = subinstr(rbc_orth, orth_suffix, "", 1)

    B  = st_data(., rmain)
    Bt = st_data(., rbc)
    P  = st_matrix(Pname)
    Vboot = st_matrix(Vbootname)

    idx_boot_main = name2idx(Vbootname, rmain)
    idx_boot_bc   = name2idx(Vbootname, rbc_orth)
    idx_all       = idx_boot_main, idx_boot_bc
    idx_Pcols     = name2idx(Pname, rbc)

    if (missing(idx_all) > 0 | missing(idx_Pcols) > 0) {
        printf("Name matching failed in mfxt_bc_boot_se\n")
        exit(3498)
    }

    Vboot_r = Vboot[idx_all, idx_all]
    Psub    = P[., idx_Pcols]
    Btorth  = Bt - B * Psub
    wbc     = B, Btorth

    se_boot = sqrt(rowsum((wbc * Vboot_r) :* wbc))
    gen_bc_boot_se = gen + "_bc_boot_se"
    st_store(., gen_bc_boot_se, se_boot)
}

void mfxt_ucb_crit(string scalar bname, string scalar vname,
                   string scalar varlist, real scalar level,
                   real scalar sim_reps, string scalar critname)
{
    string matrix vars
    real matrix V, B, Vsym, Vv, loadings, draws, dev
    real vector b, se, lambda, positive, sorted, tmax
    real scalar k, n, i, crit, eps
    real rowvector colidx

    vars = tokens(varlist)
    k = cols(vars)
    colidx = name2idx(bname, vars)
    if (missing(colidx) > 0) {
        printf("Variable name not found in matrix columns\n")
        exit(3498)
    }
    b = st_matrix(bname)[1, colidx]

    V = st_matrix(vname)[colidx, colidx]
    B = st_data(., vars)
    n = rows(B)

    se = sqrt(rowsum((B * V) :* B))
    eps = epsilon(1)
    se = se :* (se :> eps) :+ eps :* (se :<= eps)

    Vsym = (V + V') / 2
    symeigensystem(Vsym, Vv, lambda)
    lambda = lambda :* (lambda :> 0)
    positive = lambda :> 0
    if (sum(positive) > 0) {
        loadings = Vv[., positive] * diag(sqrt(lambda[positive]))
        draws = rnormal(sim_reps, cols(loadings), 0, 1) * loadings'
        dev = draws * B'
        tmax = rowmax(abs(dev) :/ se')
        sorted = sort(tmax, 1)
        crit = sorted[ceil(level * sim_reps)]
    }
    else {
        crit = invttail(1e6, (1 - level) / 2)
    }

    st_numscalar(critname, crit)
}

end
