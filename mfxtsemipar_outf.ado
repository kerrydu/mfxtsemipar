*! 2026-03-27
cap program drop mfxtsemipar_outf
program define mfxtsemipar_outf

version 16
syntax varlist(min=1) [if] [in] [fw aw pw/], ///
                            uvar(varname) ///
							id(varname) ///
							save(name)     ///
							tl(varlist)  ///
                            insample(varname) ///
							[cluster(string)   ///
                             type(string) ///
                             winsor(string) ///
							 EQSPACE ///
                             MAXNK(integer 5) ///
                             MINNK(integer 2) ///
                             ALLKNots(numlist ascending) ///
                             bknots(numlist ascending min=2 max=2) ///
                             nknots(numlist integer max=1 >0) ///
                             knots(numlist ascending) ///
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
                              PARTIALOUT1(varlist) ///
                              replace ]

if "`nknots'"!=""  | "`allknots'"!="" | "`knots'"!=""{
    outsample_pred0 `0' 
}
else{
    outsample_predk `0' 
}



end
/////////////////////////////////////////////////
program define outsample_predk

version 16
syntax varlist(min=1) [if] [in] [fw aw pw/], ///
                            uvar(varname) ///
							id(varname) ///
							save(name)     ///
							tl(varlist)  ///
                            insample(varname) ///
							[cluster(string)   ///
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
                              PARTIALOUT1(varlist) ///
                              replace ]

cap which gensplines
if _rc ssc install gensplines
cap which savesome 
if _rc ssc install savesome
cap which reghdfe 
if _rc ssc install reghdfe
cap which ftools
if _rc ssc install ftools

if "`save'"!=""{
    cap confirm new file `save'
    if _rc & "`replace'"=="" {
        di as error "file `save' already exists, please use replace option to overwrite"
        exit
    }

}

local predy_file `save'


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


 if "`predy'"!="" confirm new var `predy'
 else local predy pred_`ydep'
// 检查语法：allknots、knots和nknots至少一个需要指定

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
    

	qui collapse (mean) `ydep' `varlist0' `partialout1' `cv' `absorbvars' `weightvar' (sum) `allbins'  if `touse' & `insample' == 1, by(`id' `tl' `insample')

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

local predeq 0
foreach v in `allbins' `varlist0' {
    local predeq `predeq' + _b[`v']*`v'
}

 preserve
 qui collapse (mean)  `ydep' `varlist0' `partialout1'  `absorbvars' `weightvar' (sum) `allbins'  if `touse' & `insample' == 1, by(`id' `tl' `insample')
 tempvar res0
 reghdfe `ydep' `varlist0' `partialout1' `allbins' `weightexp', absorb(`absorb') `setype' residuals(`res0')
 qui predict `predy' if `insample' == 1, xbd
 qui estimate store mfxtsemipar_outf_insample_eststore
 tempvar ghat 
 qui predictnl double `ghat' = `predeq' if `insample' == 1
 qui reghdfe `ydep' `partialout1' `weightexp' if `insample'==1, absorb(`absorb') `setype' residuals(_M_`ydep')
 qui reghdfe `ghat' `partialout1' `weightexp' if `insample'==1, absorb(`absorb') `setype' residuals(_M_ghat)
 qui replace `predy' = `predy' + _M_ghat if `insample' == 1

 qui savesome `id' `tl' `insample' `ydep' `predy' _M_`ydep' _M_ghat using `predy_file', `replace'

 restore 

 preserve 
 qui collapse (mean)  `ydep' `varlist0' `partialout1'  `absorbvars' `weightvar' (sum) `allbins'  if `touse' & `insample' == 0, by(`id' `tl' `insample')

 qui estimate restore mfxtsemipar_outf_outsample_eststore
 tempvar ghat 
 qui predictnl double `ghat' = `predeq' if `insample' == 0
 qui reghdfe `ydep' `partialout1' `weightexp' if `insample'==0, absorb(`absorb') `setype' residuals(_M_`ydep')
 qui predict `predy' if `insample' == 0, xbd
 qui reghdfe `ghat' `partialout1' `weightexp' if `insample'==0, absorb(`absorb') `setype' residuals(_M_ghat)
 qui replace `predy' = `predy' + _M_ghat if `insample' == 0

 qui keep `id' `tl' `insample' `ydep' `predy' _M_`ydep' _M_ghat 
 append  using `predy_file'
 label var `predy' "Predicted `ydep'"
 label var _M_`ydep' "MY=Y - partialoutvar * \beta - FEs"
 label var _M_ghat   "Predicted MY based on M[g(x) + controls * \alpha]"
 save `predy_file', `replace'

 restore 

end

//////////////////////////////////////////
program define outsample_pred0

version 16
syntax varlist(min=1) [if] [in] [fw aw pw/], ///
                            uvar(varname) ///
							id(varname)  ///
							insample(varname)   ///
							tl(varlist)  ///
                            save(name)     ///
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
                            brep(real 0) ///
                            replace    ///
                            PARTIALOUT PARTIALOUT1(varlist)  ]

if "`save'"!="" & "`replace'"=="" {
    cap confirm new file `save'
    if _rc  {
        di as error "file `save' already exists, please use replace option to overwrite"
        exit
    }
}

local pred_file `save'

if (`"`weight'"'!= "" ){
    tempvar weightvar
    qui gen `weightvar' = `exp'
    local weightexp [`weight'=`weightvar']
} 

local varlist0 `varlist'
gettoken ydep varlist0: varlist0
local varlist0 `varlist0' `hfcov'

if "`partialout'"!="" & `"`partialout1'"'==""{
    local partialout1 `varlist0'
}

local  controls: list varlist0 - partialout1


if "`predy'"!="" confirm new var `predy'
else local predy pred_`ydep'


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

  local predeq 0
  foreach v in `allbins'  `controls' {
    local predeq `predeq' + _b[`v']*`v'
  }

    preserve
    qui collapse (mean) `ydep' `controls' `partialout1' `weightvar'  ///
       `absorbvars'  (sum) `allbins'  if `touse' & `insample'==1, ///
                                        by(`id' `tl' `insample')

    tempvar res0
    reghdfe `ydep' `controls' `partialout1' `allbins' `weightexp', absorb(`absorb') `setype' residuals(`res0')
    qui predict `predy' if `insample' == 1, xbd
    qui estimate store mfxtsemipar_restore
    tempvar ghat
    qui predictnl double `ghat' = `predeq' if `insample' == 1
    cap drop `res0'
    qui reghdfe `ydep' `partialout1' `weightexp' if `insample'==1, absorb(`absorb') `setype' residuals(_M_`ydep')
    qui reghdfe `ghat' `partialout1' `weightexp' if `insample'==1, absorb(`absorb') `setype' residuals(_M_ghat)
    qui savesome `id' `tl' `insample' `ydep' `predy' _M_`ydep' _M_ghat using `predy_file', `replace'
    restore

    collapse (mean) `ydep' `controls' `partialout1' `weightvar'  ///
       `absorbvars'  (sum) `allbins'  if `touse' & `insample'==0, ///
                                        by(`id' `tl' `insample')
    qui estimate store mfxtsemipar_outsample_eststore
    tempvar ghat
    qui predictnl double `ghat' = `predeq' if `insample' == 0
    cap drop `res0'
    qui reghdfe `ydep' `partialout1' `weightexp' if `insample'==0, absorb(`absorb') `setype' residuals(_M_`ydep')
    qui predict `predy' if `insample' == 0, xbd
    qui reghdfe `ghat' `partialout1' `weightexp' if `insample'==0, absorb(`absorb') `setype' residuals(_M_ghat)
    qui replace `predy' = `predy' + _M_ghat if `insample' == 0
    qui keep `id' `tl' `insample' `ydep' `predy' _M_`ydep' _M_ghat
    append using `predy_file'
    save `predy_file', `replace'

 
end


////////////////////////////////////////////


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
        qui reghdfe `gjhat' `partialout' `weightexp' if `gid'==`i', a(`absorb') residuals(`valy')
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