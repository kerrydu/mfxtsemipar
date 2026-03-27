*2026-03-24

cap program drop mfxtbin_cv
program define mfxtbin_cv,eclass
version 16

cap which gensplines
if _rc ssc install gensplines
cap which savesome 
if _rc ssc install savesome
cap which reghdfe 
if _rc ssc install reghdfe
cap which ftools
if _rc ssc install ftools

syntax  varlist(min=1) [fw aw pw/], ///
                            uvar(varname) ///
							id(varname)  ///
							tl(varlist)  ///
                            gen(string)    ///
							[MAXnbin(real  5) ///
                             MINnbin(real  2) ///
							 EQSPACE  ///
							 dropbin(string) ///
							 cluster(string)  ///
                             hfcov(varlist) ///
							 startp(numlist  max=1) ///
							endp(numlist max=1)  ///
							cvgroup(varname) ///
                            nfold(integer 10) ///
							seed(numlist max=1) ///
                            atu(varname) ///
                            absorb(string) ///
                            sopt ///
                            predy(name) ///
                            PARTIALOUT(varlist) ///
                            PARTIALOUT1(varlist)]


marksample touse
markout `touse' `uvar' `id' `tl', strok

if "`predy'"!="" confirm new var `predy'

if (`"`weight'"'!= "" ){
    tempvar weightvar
    qui gen `weightvar' = `exp'
    local weightexp [`weight'=`weightvar']
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


  local varlist0 `varlist'
  gettoken ydep varlist0: varlist0
  local varlist0  `varlist0' `hfcov'
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

local dropbincopy `dropbin'
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

mat mse = J(`minnbin'-1,1,.)
forv i=`minnbin'/`maxnbin'{
    //su `uvar2'  `uvar'
	preserve
	genbins `uvar2' if `touse', gen(__bin_) nbin(`i') `eqspace' startp(`startp') endp(`endp')
	local allbins `r(varlist)'
    local cutpoints `r(cutpoints)'
    local dropbin  `dropbincopy'
	if `"`dropbin'"'!=""{
		cap confirm integer `dropbin'
		if _rc==0{
            local dropbin __bin_`dropbin'
			local allbins: list allbins - dropbin 
		}
		else{
			numinbin `cutpoints', `dropbin'
			local dropbin = `r(inbin)'
            local dropbin __bin_`dropbin'
			local allbins: list allbins -  dropbin 
		}
		
	}

	qui collapse (mean) `ydep' `varlist0' `partialout1' `cv' `absorbvars' `weightvar' (sum) `allbins'  if `touse', by(`id' `tl')
	qui rmse_cv `ydep' `allbins' `varlist0' `weightexp',  cv(`cv') partialout(`partialout1') absorb(`absorb')
	local rmse = r(rmse)
	mat mse = mse \ `rmse'
	local mserowname `mserowname' "#_of_bins=`i'"
	restore

}


minpos mse 
//mat list mse 

local nstar = r(pos)
local min = r(min)
local soptbin `minnbin'
local mse0 1e9
forv j=`minnbin'/`maxnbin'{
    local binj = mse[`j',1]
    if `binj' != . & `binj'<`mse0' {
        local mse0 = `binj'
        local soptbin = `j'
    }
    else continue, break
}

mat mse = mse[`minnbin'..`maxnbin',1]

mat rownames mse = `mserowname'
mat colnames mse = "MSE"
matlist mse 
di _n 
if "`sopt'"!=""{
    local nstar = `soptbin'
}
genbins `uvar2' if `touse', gen(__bin_) nbin(`nstar') `eqspace' startp(`startp') endp(`endp')
local allbins `r(varlist)'
local cutpoints `r(cutpoints)'
local genbinscmd genbins `uvar2' if `touse', gen(__bin_) nbin(`nstar') `eqspace' startp(`startp') endp(`endp')
//di "`cutpoints'"
local dropbin  `dropbincopy'
if `"`dropbin'"'!=""{
    cap confirm integer `dropbin'
    if _rc==0{
        local dropbin __bin_`dropbin'
        local allbins: list allbins - dropbin 
    }
    else{
        numinbin `cutpoints', `dropbin'
        local dropbin = `r(inbin)'
        local dropbin __bin_`dropbin'
        local allbins: list allbins -  dropbin 
    }
    
}
tempvar res_yhat
preserve

qui collapse (mean)  `ydep' `varlist0' `partialout1'  `absorbvars' `weightvar' (sum) `allbins'  if `touse', by(`id' `tl')
reghdfe `ydep' `varlist0' `partialout1' `allbins' `weightexp', absorb(`absorb') `setype' residuals(`res_yhat')
if `"`predy'"'!=""{
    qui predict `predy', xbd
    tempfile predy_file
    qui savesome `id' `tl' `predy' using `predy_file', replace
}

// for safe
qui estimate store mfxtbin_cv_eststore
mat b = e(b)
estat ic,all
tempname info
mat `info' = r(S)
restore 
qui estimate restore mfxtbin_cv_eststore

if "`atu'"!=""{
    local mincuts : word 1 of `cutpoints'
    local ncuts: word count `cutpoints'
    local maxcuts : word `ncuts' of `cutpoints'
    qui su `atu' if `touse'
    local rmin = r(min)
    local rmax = r(max)
    if (`rmin'>=`mincuts'){
        di as err "mimum value of `atu' should be less than the minimum cutpoints"
        exit
    }
    if (`rmax'<=`maxcuts'){
        di as err "maximum value of `atu' should be greater than the maximum cutpoints"
        exit
    }

    qui drop __bin_*
    genbins `atu', gen(__bin_)  cut(`cutpoints')  
 }

local gf 0
foreach b in `allbins'{
	local gf  `gf'+_b[`b']*`b'
}

qui predictnl  `gen' = `gf', se(`gen'_se)

ereturn local cutpoints `cutpoints'
ereturn scalar soptbin = `soptbin'
ereturn scalar minmse = `min'
eret local genbinscmd `genbinscmd'
ereturn matrix info = `info'
ereturn matrix bmat = b
if `"`predy'"'!=""{
    qui merge m:1 `id' `tl' using `predy_file', nogen
 }
end



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


///////////////

cap program drop genbins
program define genbins,rclass
version 16
syntax  varlist(min=1 max=1) [if] [in], ///
							gen(string)     ///
							[nbin(numlist  max=1 integer >0) ///
							cut(numlist) ///
							bw(numlist >0  max=1) ///
							eqspace ///
							startp(numlist max=1) ///
							endp(numlist max=1) ]
marksample touse
local uvar `varlist'
qui su `uvar' if `touse'
if "`bw'" != "" {
	if (`=r(min)+`bw''>r(max)){
	  di as error "bw() is too large"
	  exit
	}
    if "`startp'"!="" local min = `startp'
    else local min = r(min)
	if "`endp'"!="" local max = `endp'
    else local max = r(max)
	qui numlist "`=`min'+`bw''(`bw')`max'", sort
	local cutpoints  `r(numlist)'
}
else if "`nbin'" != "" {
  if "`eqspace'"!=""{
    if "`startp'"!="" local min = `startp'
    else local min = r(min)
	if "`endp'"!="" local max = `endp'
    else local max = r(max)
	  local step = (`max'-`min')/(`nbin'-1)
	  if "`endp'"!="" qui numlist "`=`min'+`step''(`step')`max'", sort
      else qui numlist "`=`min'+`step''(`step')`=`max'-0.98*`step''", sort
	  local cutpoints  `r(numlist)'
  }
  else{
      if "`startp'"!="" & "`endp'"!="" local ifcond  if `uvar'>=`startp' & `uvar'<`endp'
      else if "`endp'"!="" {
            local ifcond if `uvar'<`endp' & `touse'
      }
      else if "`startp'"!="" {
            local ifcond if `uvar'>=`startp' & `touse'
      }
      else {
            local ifcond if `touse'
      }  

      qui _pctile `uvar' `ifcond', nq(`nbin')

	  forv i=1/`=`nbin'-1'{
		  local cutpoints `cutpoints' `=r(r`i')'
	  }
  }

}
else local cutpoints `cut'

if ("`startp'"!="") {
	local cutpoints `startp' `cutpoints'
}
if ("`endp'"!="") {
	local cutpoints `cutpoints' `endp'
}


tempvar uvarbin
classify_intervals `uvar', cutpoints(`cutpoints') gen(`uvarbin')

local nbin: word count `cutpoints'
local nbin =`nbin'+1
fvrevar i.`uvarbin' if `touse'
local allbins  `r(varlist)'
local i = 1
foreach b in `allbins' {
    if `i' == 1 {
        local p1: word 1 of `cutpoints'
        qui gen int `gen'`i' = (`uvar' < `p1') if `uvar' != .
    }
    else {
        qui gen int `gen'`i' = `b' if `uvar' != .
    }
	local genvars `genvars' `gen'`i'
	local ++i
}

return scalar nbin = `nbin'
return local cutpoints  `cutpoints'
return local varlist `genvars'

end


cap program drop numinbin
program define numinbin,rclass
version 16
gettoken cutpoints 0:0,p(,)
syntax, base(numlist min=1 max=1)

local bin =1
foreach cut in `cutpoints'{

	if `base'>`cut'{
		local bin = `bin'+1
	}
	else continue, break
}

return scalar inbin = `bin'

end 


/////////


* 定义命令
cap program drop classify_intervals
program define classify_intervals
    syntax varlist(fv max=1 min=1) [if] [in], CUTpoints(numlist) gen(name)
    * 创建区间标记变量
	confirm new var `gen'
    qui gen int `gen' = 0

	qui numlist "`cutpoints'", sort
	local cutpoints `r(numlist)'
    
    local n = 1
    foreach cut in `cutpoints' {
        if `n' == 1 {
            qui replace `gen' = `n' if `varlist' < `cut' & `varlist' != .
        }
        else {
            local prevCut = `: word `=`n'-1' of `cutpoints''
            qui replace `gen' = `n' if `varlist' >= `prevCut' & `varlist' < `cut' & `varlist' != .
        }
        local ++n
    }
    * 处理最后一个区间
    local lastCut = `: word `=`n'-1' of `cutpoints''
    qui replace `gen' = `n' if `varlist' >= `lastCut' & `varlist' != .
    
    * 标记完成
    label var `gen' "Interval for `varlist'"
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