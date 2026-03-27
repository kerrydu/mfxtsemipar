/*
    PROGRAM: mfxtbin
    AUTHOR: Unknown
    DATE: Unknown

    PURPOSE:
    This program implements a mixed-frequency fixed effects panel regression model with binning. It allows for the estimation of fixed effects models with variables at different frequencies and handles missing values and clustering.

    SYNTAX:
    mfxtbin varlist [fw aw pw/] , uvar(varname) id(varname) gen(string) tl(varlist) [cluster(string) nbin(numlist max=1 integer >0) EQSPACE cut(numlist) Absorb(varlist) bw(numlist >0 max=1) dropbin(string) atu(varname) startp(numlist max=1) endp(numlist max=1)]

    ARGUMENTS:
    - varlist: A list of variables to include in the regression model.
    - fw: The frequency weight variable.
    - aw: The annual weight variable.
    - pw: The panel weight variable.
    - uvar(varname): The variable used to define the units in the panel data.
    - id(varname): The variable used to identify the units in the panel data.
    - gen(string): The prefix for the generated variables.
    - tl(varlist): The time-varying variables to include in the model.
    - cluster(string): The variable used for clustering.
    - nbin(numlist max=1 integer >0): The number of bins to use for binning.
    - EQSPACE: Option to use equal spacing for binning.
    - cut(numlist): The cut points for binning.
    - Absorb(varlist): The variables to absorb in the fixed effects estimation.
    - bw(numlist >0 max=1): The bandwidth for the binning.
    - dropbin(string): The bin to drop from the model.
    - atu(varname): The variable used to define the time units in the panel data.
    - startp(numlist max=1): The start point for binning.
    - endp(numlist max=1): The end point for binning.

    EXAMPLES:
    mfxtbin y x1 x2, uvar(unit) id(id) gen(model) tl(time) cluster(cluster) nbin(5) Absorb(fe) bw(0.5) dropbin(bin1) atu(time) startp(2000) endp(2020)

    NOTES:
    - This program requires Stata version 16 or higher.
    - The program handles missing values and clustering.
    - The program supports fixed effects estimation with variables at different frequencies.
    - The program allows for binning of variables.
*/

cap program drop mfxtbin
program define mfxtbin, eclass
version 16
syntax varlist(min=1) [fw aw pw/], ///
                            uvar(varname) ///
							id(varname) ///
							gen(string)     ///
							tl(varlist)  ///
							[ cluster(string)   ///
							 nbin(numlist max=1 integer >0) ///
							 EQSPACE ///
							 cut(numlist) ///
                             hfcov(varlist) ///
							Absorb(varlist ts fv) ///
							bw(numlist >0  max=1) ///
                            dropbin(string) ///
                            atu(varname) ///
                            startp(numlist max=1) ///
                            endp(numlist max=1) ///
                            predy(name)]

marksample touse
markout `touse' `uvar' `id' `tl', strok

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

// if cut() is specified, display nbin() and bw() is ignored
 if "`cut'" != ""  & "`nbin'"!="" {
 	di as err "cut() and nbin() cannot be specified at the same time"
 	exit
 }
 if "`cut'" != ""  & "`bw'"!="" {
 	di as err "cut() and bw() cannot be specified at the same time"
 	exit
 }
 if "`bw'" != ""  & "`nbin'"!="" {
 	di as err "bw() and nbin() cannot be specified at the same time"
 	exit
 }

 if `"`cluster'"'!="" {
	local setype cluster(`cluster')
 }

 genbins `uvar' if `touse', gen(__bin_) nbin(`nbin') cut(`cut') bw(`bw') `eqspace' startp(`startp') endp(`endp')

 local binvars `r(varlist)'
 local cutpoints `r(cutpoints)'
 local genbinscmd  genbins `uvar' if `touse', gen(__bin_) nbin(`nbin') cut(`cut') bw(`bw') `eqspace' startp(`startp') endp(`endp')

 if "`dropbin'"!="" {
    cap confirm numeric `dropbin'
    if _rc==0 local dropbin `dropbin'
    else{
        numinbin `cutpoints', `dropbin'
        local dropbin `r(inbin)'
    }
    local dropbin __bin_`dropbin'
    local binvars: list binvars - dropbin
 
 }

 local allbin `binvars'

 preserve


 qui collapse (mean) `varlist' `weightvar' `hfcov' (sum) `allbin' if `touse', by(`id' `tl')
 tempvar res0
 reghdfe `varlist' `hfcov' `allbin' `weightexp', absorb(`absorb') `setype' residuals(`res0')
 if `"`predy'"'!=""{
    qui predict `predy', xbd
    tempfile predy_file
    qui savesome `id' `tl' `predy' using `predy_file', replace
 }


 mat b = e(b)
 estat ic,all
 /* ereturn mat info = r(S) */
 tempname info
 mat `info' = r(S)
 restore
 if "`atu'"!=""{
    local mincuts : word 1 of `cutpoints'
    local ncuts: word count `cutpoints'
    local maxcuts : word `ncuts' of `cutpoints'
    qui su `atu' if `touse'
    local rmin = r(min)
    local rmax = r(max)
    if (`rmin'>=`mincuts'){
        di as err "minimum value of `atu' should be less than the minimum cutpoints"
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
 foreach b in `allbin'{
	 local gf  `gf'+_b[`b']*`b'
 }
 qui predictnl `gen' = `gf', se(`gen'_se)

 eret local cutpoints `cutpoints'
 eret local genbinscmd `genbinscmd'
 eret matrix info = `info'
 
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
fvrevar i.`uvarbin' if `touse'  // mising is a category
local allbin  `r(varlist)'
local i = 1
foreach b in `allbin' {
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




* 2024-05-10
/* 10 fold cross validation with full fixed effects 
gen randnum = runiform()
bys year month date: replace randnum = randnum[1]
bys id year month (randnum): gen cv = mod(_n,10)+1
*/
cap program drop rmse_cv
program define rmse_cv, rclass
    version 16.0
    syntax ,Cmd(string) [CV(varname) n(integer 10) SEED(integer 1234)  opt(string) cluster(varname) LOGLik]
    tempvar gid resi2 resi mae
    qui gen double `resi2' = .
    qui gen double `mae' = .
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
    forv i=1/`numgid'{
        cap drop  __hdfe*
        qui `cmd' if `gid'!=`i', `opt'
        //di "`cmd'"
        //`cmd' if `gid'!=`i', `opt'  
        local sigma2 = e(rss)/e(df_r) 
        local absorb = e(extended_absvars) 
        local y = e(depvar)
        cap ds __hdfe*
        local felist `r(varlist)'
        if "`felist'"!=""{
            tempvar feall xb 
            qui gen double `feall' = 0
           foreach fe of local felist {
                tempvar gfe
                gettoken fe1 absorb:absorb
                local fe1 = subinstr("`fe1'","#"," ",.)
                qui egen int `gfe' = group(`fe1')
                qui bys `gfe': fillmissing `fe', with(any)
                cap drop `gfe'
                qui replace `feall' = `feall' + `fe' 
           }
           qui predict `xb' if `gid'==`i', xb 
           qui gen double `resi'   = `y' - `xb' - `feall' if `gid'==`i'
           cap drop `xb'
        }
        else{
            qui predict `resi' if `gid'==`i',r 
        }
        qui replace `resi2' = `resi'^2 if `gid'==`i'
        qui replace `mae' = abs(`resi')/abs(`y') if `gid'==`i'
        if "`loglik'"!=""  qui replace `logL' = lnnormalden(`resi',0,sqrt(`sigma2')) if `gid'==`i'
        cap drop `resi'
    }
    qui su `resi2',meanonly
    //su `resi2'
    return scalar rmse = sqrt(r(mean))
    qui su `mae',meanonly
    //su `mae',d
    return scalar rmae = r(mean)
    if "`loglik'"!=""{
        qui su `logL'
        return scalar loglik = r(sum)/`numgid'
    }
    
end