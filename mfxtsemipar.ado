/*
  Program Name: mfxtsemipar
  Description: This program estimates a mixed-frequency fixed effects model with semiparametric splines.
  Author: Unknown
  Date: Sunday, July 28, 2024 at 17:48:54
  Version: 16

  Syntax:
    mfxtsemipar varlist(min=1) [if] [in] [fw aw pw/], ///
                 uvar(varname) ///
                 id(varname) ///
                 gen(string)     ///
                 tl(varlist)  ///
                 [cluster(string)   ///
                 BKnots(numlist ascending) ///
                 DEGree(numlist max=1) ///
                 KNots(numlist ascending) ///
                 type(string) ///
                 winsor(string) ///
                 EQSPACE ///
                 NKnots(numlist integer max=1 >0) ///
                 center(numlist max=1) ///
                 Absorb(string) ///
                 atu(varname)  INTERcept ///
                 brep(real)]

  Remarks:
  - This program estimates a mixed-frequency fixed effects model using semiparametric splines.
  - The program requires Stata version 16 or higher.
  - The program supports various options for specifying the model and estimation settings.
  - The program checks for syntax errors and provides error messages if necessary.
  - The program generates spline variables based on the specified knots and degree.
  - The program supports the winsor option for winsorizing the variables.
  - The program estimates the model using the reghdfe command.
  - The program supports the absorb option for absorbing fixed effects.
  - The program supports the atu option for estimating the model for a specific set of units.
  - The program supports the brep option for wild bootstrap inference.

  Example:
    mfxtsemipar y x, uvar(id) id(id) gen(model) tl(time) knots(0.5 1 2) degree(3) type(bs) winsor(0.01, 0.99) Absorb(id) atu(unit) INTERcept brep(500)

  References:
  - [Reference 1]
  - [Reference 2]
*/
*! Sunday, July 28, 2024 at 17:48:54

cap program drop mfxtsemipar
program define mfxtsemipar, eclass

version 16
syntax varlist(min=1) [if] [in] [fw aw pw/], ///
                             uvar(varname) ///
							               id(varname) ///
							               gen(string)     ///
							               tl(varlist)  ///
							              [cluster(string)   ///
                            BKnots(numlist ascending) ///
                            DEGree(numlist max=1) ///
                            KNots(numlist ascending) ///
                            ALLKNots(numlist ascending) ///
                            type(string) ///
                            winsor(string) ///
							              EQSPACE  predy(name) ///
                            NKnots(numlist integer max=1 >0) ///
                            center(numlist max=1) ///
                            Absorb(string) ///
                            atu(varname)  INTERcept ///
                            hfcov(varlist) ///
                            brep(real 0)]

if (`"`weight'"'!= "" ){
    tempvar weightvar
    qui gen `weightvar' = `exp'
    local weightexp [`weight'=`weightvar']
} 

if "`predy'"!="" confirm new var `predy'


// 检查语法：allknots、knots和nknots至少一个需要指定
if "`allknots'" == "" & "`knots'" == "" & "`nknots'" == "" {
    di as err "At least one of allknots(), knots() and nknots() should be specified"
    exit
}

//检查语法：knots和nknots不能同时指定
if "`knots'" != "" & "`nknots'" != "" {
    di as err "knots() and nknots() cannot be specified at the same time"
    exit
}
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

if "`knots'"==""{
    gennknots `xvar' if `touse', nknots(`nknots') `eqspace'
    local knots `r(knots)'
}

`gensplines' `xvar', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type') 
 local splinecmd `gensplines' `xvar', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type') 
 local allbins  `r(splinevarlist)'

 preserve
 qui collapse (mean) `varlist' `weightvar' `absorbvars' `hfcov' (sum) `allbins'  if `touse', by(`id' `tl')
 if `brep'==0{
    tempvar res0
    reghdfe `varlist' `hfcov' `allbins' `weightexp', absorb(`absorb') `setype' residuals(`res0')
    if `"`predy'"'!=""{    if `"`predy'"'!=""{
      qui predict `predy', xbd
      tempfile predy_file
      qui savesome `id' `tl' `predy' using `predy_file', replace
    }
    mat b = e(b)
    estimate store mfxtsemipar_restore
    estat ic,all
    tempname info
    mat `info' = r(S)
 }
 else{
    tempvar res0
    qui reghdfe `varlist' `hfcov' `allbins' `weightexp', absorb(`absorb') residuals(`res0')
    if `"`predy'"'!=""{
      qui predict `predy', xbd
      tempfile predy_file
      qui savesome `id' `tl' `predy' using `predy_file', replace
    }

    mat b = e(b)
    mata: k = length(st_matrix("b"))
    tempvar yhat ehat ystar
    qui predict double `yhat', xdb
    local depvar: word 1 of `varlist'
    qui gen double `ehat' = `depvar' - `yhat'
    qui gen double `ystar' = . 
    mata: bb = J(`brep',k,.)
    forv b=1/`brep'{
        qui wildboot `ystar' `yhat' `ehat', cluster(`cluster')
        qui reghdfe `ystar' `hfcov' `allbins' `weightexp', absorb(`absorb')
        mat bi = e(b)
        mata: bb[`b',.] = st_matrix("bi")
    }
    mata: bb= quadvariance(bb)
    mata: st_matrix("V",bb)
    // 将V ereturn 为 e(V)
    cap drop `res0'
    qui reghdfe `varlist' `hfcov' `allbins' `weightexp', absorb(`absorb') residuals(`res0')
    if `"`predy'"'!=""{
      qui predict `predy', xbd
      tempfile predy_file
      qui savesome `id' `tl' `predy' using `predy_file', replace
    }
    mat b = e(b)
    ereturn repost b V
    ereturn display
    estimate store mfxtsemipar_restore
    estat ic,all
    tempname info
    mat `info' = r(S)
 }
 restore
 qui estimate restore mfxtsemipar_restore


 if "`atu'"!=""{
     drop `allbins'
    `gensplines' `atu', gen(__Spline_) knots(`knots') bknots(`bknots') degree(`degree') centerv(`center') `intercept' type(`type') 
    local allbins  `r(splinevarlist)'
 }


 local gf 0
 foreach b in `allbins'{
	 local gf  `gf'+_b[`b']*`b'
 }
 qui predictnl `gen' = `gf', se(`gen'_se)
 drop `allbins'

 ereturn local knots `knots'
 ereturn local splinecmd `splinecmd'
 ereturn mat info = `info'

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

/////
program define wildboot 
version 16
syntax varlist(min=3 max=3) , [cluster(varlist)]

local ystar : word 1 of `varlist'
local yhat : word 2 of `varlist'
local ehat : word 3 of `varlist'

tempvar radw rn
qui gen double `rn' = uniform()
if "`cluster'"=="" {
    qui gen double `radw' = cond(`rn'<=0.5,1,-1) * `ehat'
}
else {
    qui bys `cluster': gen double `radw' = cond(`rn'[1]<=0.5,1,-1) * `ehat'
}
qui replace `ystar' = `yhat' + `radw'

end
