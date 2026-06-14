{smcl}
{* *! version 1.0  27mar2026}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{vieweralsosee "[R] reghdfe" "help reghdfe"}{...}
{viewerjumpto "Syntax" "mfxtbin_cv##syntax"}{...}
{viewerjumpto "Description" "mfxtbin_cv##description"}{...}
{viewerjumpto "Options" "mfxtbin_cv##options"}{...}
{viewerjumpto "Remarks" "mfxtbin_cv##remarks"}{...}
{viewerjumpto "Examples" "mfxtbin_cv##examples"}{...}
{viewerjumpto "Stored results" "mfxtbin_cv##results"}{...}
{viewerjumpto "Author" "mfxtbin_cv##author"}{...}
{title:Title}

{phang}
{bf:mfxtbin_cv} {hline 2} Mixed-frequency cross-validated binscatter-type regression


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:mfxtbin_cv}
{depvar} [{indepvars}]
{ifin}
{weight}{cmd:,}
{opt uvar(varname)}
{opt id(varname)}
{opt gen(string)}
{opt tl(varlist)}
[{it:options}]

{synoptset 23 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt uvar(varname)}}variable entering the binned (step) nonlinear component{p_end}
{synopt:{opt id(varname)}}panel identifier{p_end}
{synopt:{opt gen(string)}}base name for fitted nonlinear part; standard errors stored as {it:name}{cmd:_se}{p_end}
{synopt:{opt tl(varlist)}}time-level variables used with {cmd:id} in the collapse to mixed frequency{p_end}

{syntab:Bins and CV}
{synopt:{opt maxnbin(real)}}maximum number of bins tried in cross-validation; default is {cmd:5}{p_end}
{synopt:{opt minnbin(real)}}minimum number of bins tried; default is {cmd:2}{p_end}
{synopt:{opt eqspace}}equally spaced cutpoints on [{cmd:startp},{cmd:endp}] or sample range; default is quantile-based bins{p_end}
{synopt:{opt startp(numlist)}}optional lower endpoint for bin range{p_end}
{synopt:{opt endp(numlist)}}optional upper endpoint for bin range{p_end}
{synopt:{opt dropbin(string)}}reference bin to exclude (integer bin index or a {it:base} value passed to internal cut logic){p_end}
{synopt:{opt cvgroup(varname)}}user-supplied CV fold variable{p_end}
{synopt:{opt nfold(integer)}}number of folds if {cmd:cvgroup()} omitted; default is {cmd:10}{p_end}
{synopt:{opt seed(numlist)}}seed for {helpb splitsample} when forming folds{p_end}
{synopt:{opt sopt}}select the first local minimum along the CV curve instead of the global minimum{p_end}

{syntab:Model}
{synopt:{opt cluster(string)}}cluster-robust SEs in final {cmd:reghdfe}{p_end}
{synopt:{opt absorb(string)}}fixed effects ({cmd:reghdfe} absorb syntax){p_end}
{synopt:{opt hfcov(varlist)}}high-frequency covariates collapsed with the rest of the model{p_end}
{synopt:{opt partialout(varlist)}}partial-out specification; see {help mfxtbin_cv##remarks:Remarks}{p_end}
{synopt:{opt partialout1(varlist)}}explicit list of variables to partial out{p_end}
{synopt:{opt predy(name)}}store linear-predictor part ({cmd:xbd}) merged on {cmd:id} {cmd:tl}{p_end}
{synopt:{opt atu(varname)}}evaluate bins on another variable using estimated cutpoints{p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}
{cmd:fweight}s, {cmd:aweight}s, and {cmd:pweight}s are allowed; see {help weight}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:mfxtbin_cv} fits a mixed-frequency panel regression in which the effect of
{cmd:uvar()} is modeled as a step function (bin indicators), similar in spirit to
binscatter with covariate adjustment. High-frequency observations are collapsed to
the frequency implied by {cmd:id()} and {cmd:tl()} before estimation.

{pstd}
For each candidate number of bins between {cmd:minnbin()} and {cmd:maxnbin()},
the command runs {it:K}-fold cross-validation: within each fold it estimates
coefficients on the bin dummies (and other regressors) on the training folds
using {helpb reghdfe}, then scores prediction error on the held-out fold. The
number of bins minimizing the cross-validated RMSE is selected, and the final
model is re-estimated at that resolution.

{pstd}
The fitted value of the nonlinear part (linear combination of bin dummies with
estimated coefficients) and its standard error are computed with {cmd:predictnl}
and stored under {cmd:gen()} and {cmd:gen()}{cmd:_se}.

{pstd}
The first run may install dependencies from SSC: {helpb gensplines}, {helpb savesome},
{helpb reghdfe}, and {helpb ftools}.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt uvar(varname)}, {opt id(varname)}, {opt gen(string)}, and {opt tl(varlist)}
are required. {cmd:varlist} must contain the dependent variable followed by
low-frequency (or already aggregated) regressors; see {opt hfcov()} for extra
high-frequency controls.

{dlgtab:Bins and CV}

{phang}
{opt maxnbin()} and {opt minnbin()} bound the grid of bin counts. Internally,
bin indicators are constructed with the same {cmd:genbins} logic as in {helpb mfxtbin}.

{phang}
{opt eqspace} places interior cutpoints at equal spacing; otherwise interior
cutpoints are empirical quantiles of {cmd:uvar()} (subject to {cmd:if}/{cmd:in}
and {cmd:startp()}/{cmd:endp()}).

{phang}
{opt dropbin()} omits one bin from the design so the others are identified
relative to that reference. Supply an integer bin index, or a numeric {it:base}
interpreted against the cutpoint list (see the helper {cmd:numinbin} in the ado).

{phang}
{opt cvgroup()} supplies integer fold IDs. If omitted, folds are created by
{cmd:splitsample, cluster(}{cmd:id()}{cmd:) nsplit(}{cmd:nfold()}{cmd:)}.

{phang}
{opt sopt} replaces the global CV minimum by the {it:first} bin count at which
the CV error stops decreasing when moving from {cmd:minnbin()} upward (useful when
the CV curve is flat or noisy).

{dlgtab:Model}

{phang}
{opt absorb()} and {opt cluster()} are passed to {cmd:reghdfe} for the final fit
and for CV training folds (clustering in CV follows {cmd:rmse_cv} as implemented
in the ado).

{phang}
{opt hfcov()} lists variables that enter the collapse and the regression at high
frequency together with {cmd:uvar()} and the dependent variable.

{phang}
{opt partialout()} without {cmd:partialout1()} sets the partial-out set to all
regressors other than the dependent variable. {opt partialout1()} overrides that
list explicitly. These options mirror the partial-out path in {cmd:mfxtbin_cv}'s
{cmd:reghdfe} / CV steps.

{phang}
{opt predy()} saves the {cmd:xbd} linear predictor from the final {cmd:reghdfe}
on the collapsed data and merges it back to the original data by {cmd:id()}
and {cmd:tl()}.

{phang}
{opt atu()} rebuilds bin dummies from the {it:selected} cutpoints on a different
evaluation variable. The command requires values of {cmd:atu()} to extend below
the smallest and above the largest cutpoint so every bin can occur.


{marker remarks}{...}
{title:Remarks}

{pstd}
Cross-validation RMSE by bin count is printed as a matrix before the final fit.
Stored quantities include the chosen cutpoint list and the minimizing bin count;
see {help mfxtbin_cv##results:Stored results}.

{pstd}
Requires Stata 16+. The estimation stack relies on {cmd:reghdfe}, Mata, and
{cmd:predictnl}.


{marker examples}{...}
{title:Examples}

{pstd}Setup: mixed-frequency collapse with CV over bin count{p_end}
{phang2}{cmd:. mfxtbin_cv y x1 x2, uvar(x3) id(panelid) tl(date) gen(gfit) absorb(panelid)}{p_end}

{pstd}Equally spaced bins and user fold variable{p_end}
{phang2}{cmd:. mfxtbin_cv y x1, uvar(x3) id(id) tl(week) gen(g1) eqspace cvgroup(mycv)}{p_end}

{pstd}Narrow the CV grid and drop a reference bin{p_end}
{phang2}{cmd:. mfxtbin_cv y x1 x2, uvar(x3) id(id) tl(t) gen(g2) minnbin(3) maxnbin(8) dropbin(1)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:mfxtbin_cv} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(soptbin)}}selected number of bins after CV (or under {cmd:sopt}, the simple-optimum count){p_end}
{synopt:{cmd:e(minmse)}}minimum cross-validated RMSE reported by the internal CV routine{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cutpoints)}}full cutpoint list used to generate final bins{p_end}
{synopt:{cmd:e(genbinscmd)}}{cmd:genbins} command line for the chosen bin count{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(bmat)}}coefficient vector from the final {cmd:reghdfe}{p_end}
{synopt:{cmd:e(info)}}information matrix from {cmd:estat ic, all}{p_end}


{marker author}{...}
{title:Author}

{pstd}Kerui Du{p_end}
{pstd}Xiamen University{p_end}
{pstd}kerrydu@xmu.edu.cn{p_end}


{title:Also see}

{p 4 14 2}
Online:  {helpb mfxtsemipar_cv}, {helpb reghdfe}, {helpb splitsample}, {helpb predictnl}

{p 4 14 2}
Related command:  {hi:mfxtbin} (companion ado; install with the package)
