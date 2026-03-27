*! mfxtsemipar_cv - Mixed-frequency cross-validated semiparametric regression
*! version 1.0 - July 28, 2024
*! 
*! Syntax:
*!   mfxtsemipar_cv depvar [indepvars] [if] [in] [weight], 
*!                  uvar(varname) id(varname) gen(string) tl(varlist)
*!                  [options]
*!
*! Description:
*!   This command performs mixed-frequency semiparametric regression with 
*!   cross-validation for optimal knot selection. It fits spline or polynomial
*!   functions to a univariate variable while controlling for other covariates
*!   and fixed effects in a panel data setting.
*!
*! Required options:
*!   uvar(varname)     - Variable for semiparametric transformation
*!   id(varname)       - Panel identifier variable
*!   gen(string)       - Name for generated fitted values variable
*!   tl(varlist)       - Time-level variables for collapse operation
*!
*! Options:
*!   cluster(string)   - Cluster variable for robust standard errors
*!   type(string)      - Spline type: bs, ms, is, ibs, or poly (default: poly)
*!   winsor(string)    - Winsorization percentiles for uvar
*!   eqspace          - Use equally spaced knots instead of quantile-based
*!   maxnk(integer)   - Maximum number of knots to consider (default: 5)
*!   minnk(integer)   - Minimum number of knots to consider (default: 2)
*!   center(numlist)  - Centering value for splines
*!   absorb(string)   - Fixed effects specification
*!   partialout       - Use partial-out estimation method
*!   cvgroup(varname) - Variable defining CV groups (alternative to automatic)
*!   nfold(integer)   - Number of CV folds (default: 10)
*!   seed(numlist)    - Random seed for CV group generation
*!   degree(integer)  - Polynomial degree (default: 1)
*!   keepsplines      - Keep generated spline variables
*!   atu(varname)     - Alternative variable for spline evaluation
*!   hfcov(varlist)   - High-frequency covariates
*!   dropfirstbase    - Drop intercept from spline basis
*!   sopt             - Use simple optimal knot selection (first local minimum)
*!   brep(real)       - Number of bootstrap replications for inference (default: 0)
*!
*! Returns:
*!   e(soptnk)        - Simple optimal number of knots
*!   e(knots)         - Selected knot locations
*!   e(minmse)        - Minimum cross-validation MSE
*!   e(splinecmd)     - Command used to generate optimal splines
*!   e(info)          - Model information criteria
*!
*! Examples:
*!   . mfxtsemipar_cv y x1 x2, uvar(z) id(firm) gen(fitted) tl(year quarter)
*!   . mfxtsemipar_cv y, uvar(x) id(id) gen(f) tl(t) type(bs) maxnk(8) absorb(year)
*!
*! Author: Kerui Du
*! Support: kerrydu@xmu.edu.cn

*! Sunday, July 28, 2024 at 17:48:54

cap program drop mfxtsemipar_cv
program define mfxtsemipar_cv, eclass

version 16
syntax varlist(min=1) [if] [in] [fw aw pw/], ///
                            uvar(varname) ///
							id(varname) ///
							gen(string)     ///
							tl(varlist)  ///
							[cluster(string)   ///
                             bknots(numlist ascending min=2 max=2) ///
                             type(string) ///
                             winsor(string) ///
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
                              atu(varname)  ///
                              hfcov(varlist) ///
                              dropfirstbase sopt ///
                              brep(real 0) ///
                              predy(name) ///
                              PARTIALOUT ///
                              PARTIALOUT1(varlist)]

cap which gensplines
if _rc ssc install gensplines
cap which savesome 
if _rc ssc install savesome
cap which reghdfe 
if _rc ssc install reghdfe
cap which ftools
if _rc ssc install ftools

if "`predy'"!="" confirm new var `predy'

if (`"`weight'"'!= "" ){
    tempvar weightvar
    qui gen `weightvar' = `exp'
    local weightexp [`weight'=`weightvar']
} 
 if "`dropfirstbase'"!=""{
    local intercept 
 }
 else local intercept intercept         


 local varlist0 `varlist'
 gettoken ydep varlist0: varlist0
 local varlist0 `varlist0' `hfcov'
 if "`partialout'"!="" & `"`partialout1'"'==""{
    local partialout1 `varlist0'
    local varlist `ydep'
    local varlist0 
 }
 if `"`partialout1'"'!=""{
    local varlist0 : list varlist0 - partialout1
    local varlist `ydep' `varlist0'
 }
 local partialout `partialout1'


local gensplines gensplines
if "`type'" == "" | "`type'"=="poly" {
   local type poly
   local gensplines polysplines
}

// inlist("`type'","bs","ms","is","ibs","poly")==0 error
if inlist("`type'","bs","ms","is","ibs","poly")==0 {
    di as err "type() should be one of bs, ms, is, ibs and poly"
    exit
}

if (`"`absorb'"'!= "" ){
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
      //local absorbvars : list uniq absorbvars
  }

marksample touse
markout `touse' `uvar' `id' `tl' `hfcov', strok

local xvar `uvar'
// winsor options
if "`winsor'" != "" {
    if strpos("`winsor'",",") == 0 local winsor `winsor',
    GetWinsorOpts `winsor' xvar(`xvar')  touse(`touse')
    if "`xvar'" != "" {
      tempvar xvar_winsor
      local vartype: type `xvar'
      gen `vartype' `xvar_winsor' = cond(`xvar'<=`winsor_low', ///
                                      `winsor_low',        ///
                                      cond(`xvar'>=`winsor_high',`winsor_high',`xvar'))
      local xvar `xvar_winsor'
    }
  }  

qui su `xvar' if `touse'
local min = r(min) - 0.01
local max = r(max) + 0.01

if "`bknots'" == "" {
    local bknots `min' `max'
}


if `"`cluster'"'!="" {
    local setype cluster(`cluster')
}

if `"`cvgroup'"'!=""{
	tempvar cv
	qui egen `cv' = group(`cvgroup')
}
else {
	tempvar cv
	if "`seed'"!="" qui set seed `seed'
	qui splitsample, cluster(`id') nsplit(`nfold') gen(`cv')
}

tempvar uvar2 
clonevar `uvar2' = `uvar'

mat mse = J(`minnk'-1,1,.)


////////

forv i=`minnk'/`maxnk'{
	preserve
    gennknots `xvar' if `touse', nknots(`i') `eqspace'
    local knots `r(knots)'
    `gensplines' `xvar', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type') 
    local allbins  `r(splinevarlist)'
    

	qui collapse (mean) `ydep' `varlist0' `partialout1' `cv' `absorbvars' `weightvar' (sum) `allbins'  if `touse', by(`id' `tl')

	qui rmse_cv `ydep' `allbins' `varlist0' `weightexp',  cv(`cv') partialout(`partialout1') absorb(`absorb')
	local rmse = r(rmse)
	mat mse = mse \ `rmse'
    local mserowname "`mserowname' nk=`i'"
	restore

}


minpos mse 
//mat list mse 
// plot the trend of mse which stored in matrix mse

local nknots = r(pos)
//di "nknots: " `nknots'
local min = r(min)

local soptnk = `minnk'
local mse0 1e9
forv j=`minnk'/`maxnk'{
    local binj = mse[`j',1]
    if `binj' != . & `binj'<`mse0' {
        local mse0 = `binj'
        local soptnk = `j'
    }
    else continue, break
}


mat mse = mse[`minnk'..`maxnk',1]
mat rownames mse = `mserowname'
mat colnames mse = "MSE"
matlist mse 
di _n 
//////////
if "`sopt'"!="" local nknots `soptnk'


gennknots `xvar' if `touse', nknots(`nknots') `eqspace'
local knots `r(knots)'

`gensplines' `xvar', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type') 
local allbins  `r(splinevarlist)'
local splinecmd `gensplines' `xvar', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type') 

 preserve
 qui collapse (mean)  `ydep' `varlist0' `partialout1'  `absorbvars' `weightvar' (sum) `allbins'  if `touse', by(`id' `tl')
 if `brep'==0{
    tempvar res0
    reghdfe `ydep' `varlist0' `partialout1' `allbins' `weightexp', absorb(`absorb') `setype' residuals(`res0')
    if `"`predy'"'!=""{
        qui predict `predy', xbd
        tempfile predy_file
        qui savesome `id' `tl' `predy' using `predy_file', replace
    }

    // for safe
    qui estimate store mfxtsemipar_cv_restore
    mat b = e(b)
    qui estat ic,all
    tempname info
    mat `info'= r(S)

 }
 else{
    tempvar res0
    qui reghdfe `ydep' `varlist0' `partialout1' `allbins' `weightexp', absorb(`absorb') residuals(`res0')
    mat b = e(b)
    mata: k = length(st_matrix("b"))
    // estat ic,all
    // tempname info
    // mat `info'= r(S)
    tempvar yhat ehat ystar 
    qui predict double `yhat', xbd
    //local depvar: word 1 of `varlist'
    qui gen double `ehat' = `ydep' - `yhat'
    qui gen double `ystar' = . 
    mata: bb = J(`brep',k,.)
    forv b=1/`brep'{
        qui wildboot `ystar' `yhat' `ehat', cluster(`cluster')
        qui reghdfe `ystar' `varlist0' `partialout1' `allbins' `weightexp', absorb(`absorb')
        mat bi = e(b)
        mata: bb[`b',.] = st_matrix("bi")
    }
    mata: bb= quadvariance(bb)
    mata: st_matrix("V",bb)
    // 将V ereturn 为 e(V)
    tempvar res_yhat
    qui reghdfe `ydep' `varlist0' `partialout1' `allbins' `weightexp', absorb(`absorb') residuals(`res_yhat')
    mat b = e(b)
    if `"`predy'"'!=""{
        qui predict `predy', xbd
        tempfile predy_file
        qui savesome `id' `tl' `predy' using `predy_file', replace
    }
    ereturn repost b V
    ereturn display
    qui estat ic,all
    tempname info
    mat `info'= r(S)
    qui estimate store mfxtsemipar_cv_restore
 }

 //2025-03-23

 restore
 qui estimate restore mfxtsemipar_cv_restore

 if "`atu'"!=""{
     cap drop `allbins'
    `gensplines' `atu', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type') 
    local allbins  `r(splinevarlist)'
 }


 local gf 0
 foreach b in `allbins'{
	 local gf  `gf'+_b[`b']*`b'
 }
 qui predictnl `gen' = `gf', se(`gen'_se)
 if "`keepsplines'"=="" drop `allbins'

ereturn scalar soptnk = `soptnk'
ereturn local knots `knots'
ereturn scalar minmse =`min'
ereturn local splinecmd `splinecmd'
ereturn matrix info = `info'
ereturn matrix bmat = b

if `"`predy'"'!=""{
   qui merge m:1 `id' `tl' using `predy_file', nogen
}
end



cap program drop polysplines
program define polysplines, rclass
version 14
syntax varlist(min=1 max=1), knots(numlist) gen(string) [degree(integer 1) centerv(integer 0) *]

local power `degree' 
local xvar `varlist'
local j =1
forv p =1/`power'{
	qui gen double `gen'`j' = `xvar'^`p' - `centerv'^`p'
    local splinevarlist `splinevarlist' `gen'`j'
	local j= `j'+1
}


foreach num of numlist `knots'{
	qui gen double `gen'`j' = (`xvar' - `num')^`power'*(`xvar'>`num') - (`centerv' - `num')^`power'*(`centerv'>`num')
	local splinevarlist `splinevarlist' `gen'`j'
    local j=`j'+1
}

return local knots `knots'
return local splinevarlist `splinevarlist'

end


cap program drop gennknots
program define gennknots , rclass
version 16
syntax varlist(min=1 max=1) [if] [in], nknots(integer) [eqspace startp(numlist) endp(numlist)]

marksample touse
local uvar `varlist'
qui su `varlist' if `touse'
local nbin = `nknots' + 1
if "`eqspace'"!=""{
    if "`startp'"!="" local min = `startp'
    else local min = r(min)
	if "`endp'"!="" local max = `endp'
    else local max = r(max)
	  local step = (`max'-`min')/(`nbin'-1)
	  qui numlist "`=`min'+`step''(`step')`max'", sort
	  local cutpoints  `r(numlist)'
  }
  else{
      if "`startp'"!="" & "`endp'"!="" local ifcond  if `uvar'>=`startp' & `uvar'<`endp'
      else if "`endp'"!="" {
            local ifcond if `uvar'<`endp'
      }
      else if "`startp'"!="" {
            local ifcond if `uvar'>=`startp'
      }
      else {
            local ifcond
      }  

      qui _pctile `varlist' `ifcond', nq(`nbin')

	  forv i=1/`=`nbin'-1'{
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
    di as error "If usig the winsor() option with a scalar you need to use the values option."
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


///////

* 2024-05-10
/* 10 fold cross validation with full fixed effects 
gen randnum = runiform()
bys year month date: replace randnum = randnum[1]
bys id year month (randnum): gen cv = mod(_n,10)+1
*/
cap program drop rmse_cv
program define rmse_cv, rclass
    version 16.0
    syntax varlist [fw aw pw/],  [CV(varname) n(integer 10) SEED(integer 1234)  opt(string) cluster(varname) LOGLik  absorb(string)  PARTIALOUT(varlist)]

    if ("`weight'"!="") local weightexp [`weight'=`exp']

    tempvar gid resi2 resi mae
    qui gen double `resi2' = .
    // qui gen double `mae' = .
    if "`loglik'"!=""{
        tempvar logL
        qui gen `logL'=.
    }
    if(`"`cv'"'!=""){
        qui egen int  `gid' = group(`cv')
    }
    else{
        set seed `seed'
        if "`cluster'"!="" {
            qui splitsample, cluster(`cluster') nsplit(`n') gen(`gid')
        }
        else{
            qui splitsample, nsplit(`n') gen(`gid')
        }
    }
    qui su `gid'
    local numgid = r(max)
    // 特别注意这个varlist 只是y和splines，不包括控制变量
    gettoken depvar varlist: varlist 
    local nvars : word count `varlist'

    // forv j=1/`nvars'{
    //     tempvar ms`j'
    // }
    tempvar valy  eterm  gjhat
    // qui gen double `eterm' = .
   qui gen double `gjhat' = .

    forv i=1/`numgid'{
        qui reghdfe `depvar' `varlist' `partialout' `weightexp' if `gid'!=`i', a(`absorb')
        mat b = e(b)
        mat b = b[1,1..`nvars']

        ////
        // qui reghdfe `depvar'  `controls' `weightexp' if `gid'==`i', a(`absorb') residuals(`valy')
        // qui replace `eterm' = `valy'  if `gid'==`i'
        * compute `varlist'*b for numgid == `i'
        // forv j=1/`nvars'{
        //     local Mvars
        //     local varj : word `j' of `varlist'
        //     qui reghdfe `varj'  `controls' `weightexp' if `gid'==`i', a(`absorb') residuals(`ms`j'')
        //     local Mvars `Mvars' `ms`j''
        // }
    //    // qui replace `mae' = abs(`valy'/`eterm') if `numgid' == `i'
    //    qui replace `resi2' = (`eterm' - `valy')^2 if `gid'==`i'
    //    cap drop `valy'
    //    cap drop `Mvars'


        // qui computeghat `valy' `Mvars', bmat(b)
        //// tricky (MY-Mg(X))'(MY-Mg(X))
        // reg Y-g(X) on M
        qui computeghat `gjhat' `varlist', bmat(b)
        qui replace `gjhat' = `depvar' - `gjhat' if `gid'==`i'
        qui reghdfe `gjhat' `partialout1' `weightexp' if `gid'==`i', a(`absorb') residuals(`valy')
        qui replace `resi2' = (`valy')^2 if `gid'==`i'
        cap drop `valy'
        //cap drop `gjhat'
    }

    qui su `resi2',meanonly
    //su `resi2'
    return scalar rmse = sqrt(r(mean))
    
end

////////////////////

cap program drop computeghat
program define computeghat
    version 16
    syntax varlist , bmat(name)
    gettoken ydep varlist: varlist
    mata: b =  st_matrix("`bmat'")
    mata: st_view(sdata=.,.,"`varlist'")
    mata: st_view(ydep=.,.,"`ydep'")
    mata: ydep[.,.] = sdata*(b')

end


/////
//////////////////
cap program drop minpos
program define minpos,rclass
version 14
args m

confirm matrix `m'

mata: __mm__ = st_matrix("`m'")
mata: __mm__ = (range(1,rows(__mm__),1), __mm__)
mata: __mm__=sort(__mm__,2)
mata: st_numscalar("r(pos)",__mm__[1,1])
mata: st_numscalar("r(min)",__mm__[1,2])
return scalar pos =r(pos)
return scalar min =r(min)
end

/////
program define wildboot 
version 16
syntax varlist(min=3 max=3) , [cluster(varlist)]

local ystar : word 1 of `varlist'
local yhat : word 2 of `varlist'
local ehat : word 3 of `varlist'

tempvar radw rn
qui gen double `rn' = runiform()
if "`cluster'"=="" {
    qui gen double `radw' = cond(`rn'<=0.5,1,-1) * `ehat'
}
else {
    qui bys `cluster': gen double `radw' = cond(`rn'[1]<=0.5,1,-1) * `ehat'
}
qui replace `ystar' = `yhat' + `radw'


end