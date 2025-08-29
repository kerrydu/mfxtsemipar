{smcl}
{* *! version 1.0  28jul2024}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{vieweralsosee "[R] reghdfe" "help reghdfe"}{...}
{viewerjumpto "Syntax" "mfxtsemipar_cv##syntax"}{...}
{viewerjumpto "Description" "mfxtsemipar_cv##description"}{...}
{viewerjumpto "Options" "mfxtsemipar_cv##options"}{...}
{viewerjumpto "Examples" "mfxtsemipar_cv##examples"}{...}
{viewerjumpto "Stored results" "mfxtsemipar_cv##results"}{...}
{viewerjumpto "References" "mfxtsemipar_cv##references"}{...}
{viewerjumpto "Author" "mfxtsemipar_cv##author"}{...}
{title:Title}

{phang}
{bf:mfxtsemipar_cv} {hline 2} Mixed-frequency cross-validated semiparametric regression


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:mfxtsemipar_cv}
{depvar} [{indepvars}]
{ifin}
{weight}{cmd:,}
{opt uvar(varname)}
{opt id(varname)}
{opt gen(string)}
{opt tl(varlist)}
[{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt uvar(varname)}}variable for semiparametric transformation{p_end}
{synopt:{opt id(varname)}}panel identifier variable{p_end}
{synopt:{opt gen(string)}}name for generated fitted values variable{p_end}
{synopt:{opt tl(varlist)}}time-level variables (low-frequency identifiers) for collapse operation{p_end}

{syntab:Model}
{synopt:{opt cluster(string)}}cluster variable for robust standard errors{p_end}
{synopt:{opt type(string)}}spline type: {cmd:bs}, {cmd:ms}, {cmd:is}, {cmd:ibs}, or {cmd:poly}; default is {cmd:poly}{p_end}
{synopt:{opt winsor(string)}}winsorization percentiles for uvar{p_end}
{synopt:{opt eqspace}}use equally spaced knots instead of quantile-based{p_end}
{synopt:{opt maxnk(integer)}}maximum number of knots to consider; default is {cmd:5}{p_end}
{synopt:{opt minnk(integer)}}minimum number of knots to consider; default is {cmd:2}{p_end}
{synopt:{opt center(numlist)}}centering value for splines{p_end}
{synopt:{opt absorb(string)}}fixed effects specification{p_end}
{synopt:{opt partialout}}use partial-out estimation method{p_end}
{synopt:{opt degree(integer)}}polynomial degree; default is {cmd:1}{p_end}
{synopt:{opt hfcov(varlist)}}high-frequency covariates{p_end}
{synopt:{opt atu(varname)}}alternative variable for spline evaluation{p_end}

{syntab:Cross-validation}
{synopt:{opt cvgroup(varname)}}variable defining CV groups (alternative to automatic){p_end}
{synopt:{opt nfold(integer)}}number of CV folds; default is {cmd:10}{p_end}
{synopt:{opt seed(numlist)}}random seed for CV group generation{p_end}

{syntab:Advanced}
{synopt:{opt keepsplines}}keep generated spline variables{p_end}
{synopt:{opt dropfirstbase}}drop intercept from spline basis{p_end}
{synopt:{opt sopt}}use simple optimal knot selection (first local minimum){p_end}
{synopt:{opt brep(real)}}number of bootstrap replications for inference; default is {cmd:0}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{cmd:fweight}s, {cmd:aweight}s, and {cmd:pweight}s are allowed; see {help weight}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:mfxtsemipar_cv} performs mixed-frequency semiparametric regression with 
cross-validation for optimal knot selection. It fits spline or polynomial
functions to a univariate variable while controlling for other covariates
and fixed effects in a panel data setting.

{pstd}
The command uses cross-validation to automatically select the optimal number
of knots for the semiparametric component, minimizing out-of-sample prediction
error. The method is particularly useful for modeling nonlinear relationships
in panel data with mixed-frequency variables.

{pstd}
The estimation proceeds in two steps: (1) cross-validation is used to determine
the optimal number of knots by comparing mean squared errors across different
knot specifications, and (2) the final model is estimated using the selected
number of knots.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt uvar(varname)} specifies the variable for semiparametric transformation.
This is the variable for which the nonlinear relationship with the dependent
variable will be modeled using splines or polynomials. This option is required.

{phang}
{opt id(varname)} specifies the panel identifier variable. This variable
defines the cross-sectional units in the panel data. This option is required.

{phang}
{opt gen(string)} specifies the name for the generated fitted values variable.
The command will create a new variable with this name containing the fitted
nonlinear component and its standard error (with suffix "_se"). This option is required.

{phang}
{opt tl(varlist)} specifies the time-level variables for the collapse operation.
These variables define the time dimension and are used when aggregating
high-frequency data to lower frequencies. This option is required.

{dlgtab:Model}

{phang}
{opt cluster(string)} specifies the cluster variable for computing robust
standard errors clustered at the specified level.

{phang}
{opt type(string)} specifies the type of spline to use. Options are:
{cmd:bs} (B-splines), {cmd:ms} (M-splines), {cmd:is} (I-splines), 
{cmd:ibs} (integrated B-splines), or {cmd:poly} (polynomial splines).
The default is {cmd:poly}.

{phang}
{opt winsor(string)} specifies winsorization percentiles for the {cmd:uvar}.
For example, {cmd:winsor(1 99)} will winsorize the variable at the 1st and
99th percentiles to reduce the influence of outliers.

{phang}
{opt eqspace} requests that knots be placed at equally spaced intervals
rather than at quantiles of the {cmd:uvar} distribution.

{phang}
{opt maxnk(integer)} specifies the maximum number of knots to consider
in the cross-validation procedure. The default is 5.

{phang}
{opt minnk(integer)} specifies the minimum number of knots to consider
in the cross-validation procedure. The default is 2.

{phang}
{opt center(numlist)} specifies the centering value for splines. This can
help with numerical stability and interpretation.

{phang}
{opt absorb(string)} specifies the fixed effects to absorb. The syntax
follows that of {helpb reghdfe}.

{phang}
{opt partialout} requests the use of a partial-out estimation method,
which can be computationally more efficient for models with many fixed effects.

{phang}
{opt degree(integer)} specifies the polynomial degree for polynomial splines.
The default is 1 (linear splines).

{phang}
{opt hfcov(varlist)} specifies high-frequency covariates that should be
included in the model.

{phang}
{opt atu(varname)} specifies an alternative variable for spline evaluation.
This allows the spline to be evaluated at different points than those used
for estimation.

{dlgtab:Cross-validation}

{phang}
{opt cvgroup(varname)} specifies a variable that defines cross-validation
groups. If not specified, the command will automatically generate CV groups.

{phang}
{opt nfold(integer)} specifies the number of cross-validation folds.
The default is 10.

{phang}
{opt seed(numlist)} specifies the random seed for cross-validation group
generation, ensuring reproducible results.

{dlgtab:Advanced}

{phang}
{opt keepsplines} requests that the generated spline variables be kept
in the dataset after estimation.

{phang}
{opt dropfirstbase} requests that the intercept be dropped from the
spline basis functions.

{phang}
{opt sopt} requests the use of simple optimal knot selection, which
selects the first local minimum in the cross-validation curve rather
than the global minimum.

{phang}
{opt brep(real)} specifies the number of bootstrap replications for
wild bootstrap inference. The default is 0 (no bootstrap).


{marker examples}{...}
{title:Examples}

    {hline}
{pstd}Simulate mixed-frequency panel data with U-shaped nonlinear function{p_end}
{phang2}{cmd:. clear all}{p_end}
{phang2}{cmd:. set seed 123456}{p_end}
{phang2}{cmd:. set obs 50}{p_end}
{phang2}{cmd:. matrix C = (1,0,0.42\0,1,0.85\0.42,0.85,1)}{p_end}
{phang2}{cmd:. drawnorm x2f x3f d, corr(C)}{p_end}
{phang2}{cmd:. generate id = _n}{p_end}
{phang2}{cmd:. expand 20}{p_end}
{phang2}{cmd:. bysort id: generate t = _n}{p_end}
{phang2}{cmd:. xtset id t}{p_end}
{phang2}{cmd:. matrix D = (1,0.2,0.8\0.2,1,0\0.8,0,1)}{p_end}
{phang2}{cmd:. expand 10}{p_end}
{phang2}{cmd:. drawnorm x1 x2e x3e, corr(D)}{p_end}
{phang2}{cmd:. generate x2 = (x2f+x2e)/sqrt(2)}{p_end}
{phang2}{cmd:. bys id t: gen day = _n}{p_end}
{phang2}{cmd:. generate x3 = (x3f+x3e)/sqrt(2)}{p_end}
{phang2}{cmd:. generate gf = 1*x3 + 2*x3^2 - 0.25*(x3)^3}{p_end}
{phang2}{cmd:. bys id t: egen gftot = total(gf)}{p_end}
{phang2}{cmd:. bys id t: ereplace x1 = total(x1)}{p_end}
{phang2}{cmd:. bys id t: ereplace x2 = total(x2)}{p_end}
{phang2}{cmd:. drawnorm e}{p_end}
{phang2}{cmd:. bys id t: ereplace e = total(e)}{p_end}
{phang2}{cmd:. bys id t: replace y = x1 - x2 + gftot + d + e}{p_end}

    {hline}
{pstd}Basic mixed-frequency semiparametric regression with cross-validation{p_end}
{phang2}{cmd:. mfxtsemipar_cv y x1 x2, uvar(x3) id(id) tl(t) gen(fitted) center(0)}{p_end}

    {hline}
{pstd}With fixed effects, clustering, and polynomial degree 1{p_end}
{phang2}{cmd:. mfxtsemipar_cv y x1 x2, uvar(x3) id(id) tl(t) gen(poly1_fit) center(0) absorb(id) cluster(id) type(poly) degree(1) partialout}{p_end}

    {hline}
{pstd}Using higher-order polynomial splines{p_end}
{phang2}{cmd:. mfxtsemipar_cv y x1 x2, uvar(x3) id(id) tl(t) gen(poly2_fit) center(0) type(poly) degree(2) maxnk(10) minnk(2)}{p_end}

    {hline}
{pstd}With evaluation at specific grid points{p_end}
{phang2}{cmd:. gen x3grid = -3 + _n/10 + 0.1 if _n <= 61}{p_end}
{phang2}{cmd:. mfxtsemipar_cv y x1 x2, uvar(x3) id(id) tl(t) gen(grid_fit) center(0) atu(x3grid) absorb(id) partialout}{p_end}

    {hline}
{pstd}Simulate data with sinusoidal nonlinear function{p_end}
{phang2}{cmd:. clear all}{p_end}
{phang2}{cmd:. set seed 123456}{p_end}
{phang2}{cmd:. set obs 50}{p_end}
{phang2}{cmd:. matrix C = (1,0,0.42\0,1,0.85\0.42,0.85,1)}{p_end}
{phang2}{cmd:. drawnorm x2f x3f d, corr(C)}{p_end}
{phang2}{cmd:. generate id = _n}{p_end}
{phang2}{cmd:. expand 20}{p_end}
{phang2}{cmd:. bysort id: generate t = _n}{p_end}
{phang2}{cmd:. matrix D = (1,0.2,0.8\0.2,1,0\0.8,0,1)}{p_end}
{phang2}{cmd:. expand 10}{p_end}
{phang2}{cmd:. drawnorm x1 x2e x3e, corr(D)}{p_end}
{phang2}{cmd:. generate x2 = (x2f+x2e)/sqrt(2)}{p_end}
{phang2}{cmd:. generate x3 = (x3f+x3e)/sqrt(2)}{p_end}
{phang2}{cmd:. generate gf = 3*sin(_pi/3*x3)}{p_end}
{phang2}{cmd:. bys id t: egen gftot = total(gf)}{p_end}
{phang2}{cmd:. bys id t: ereplace x1 = mean(x1)}{p_end}
{phang2}{cmd:. bys id t: ereplace x2 = mean(x2)}{p_end}
{phang2}{cmd:. drawnorm e}{p_end}
{phang2}{cmd:. bys id t: ereplace e = total(e)}{p_end}
{phang2}{cmd:. bys id t: replace y = x1 - x2 + gftot + d + e}{p_end}

    {hline}
{pstd}Estimate sinusoidal function with cross-validation{p_end}
{phang2}{cmd:. mfxtsemipar_cv y x1 x2, uvar(x3) id(id) tl(t) gen(sin_fit) center(0) absorb(id) type(poly) degree(2) partialout}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:mfxtsemipar_cv} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(soptnk)}}simple optimal number of knots{p_end}
{synopt:{cmd:e(minmse)}}minimum cross-validation MSE{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(knots)}}selected knot locations{p_end}
{synopt:{cmd:e(splinecmd)}}command used to generate optimal splines{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(info)}}model information criteria{p_end}


{marker references}{...}
{title:References}

{phang}
Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2024. 
On binscatter. {it:American Economic Review} 114: 1488-1514.

{phang}
Henderson, D. J., and C. F. Parmeter. 2015. 
{it:Applied Nonparametric Econometrics}. 
Cambridge: Cambridge University Press.

{phang}
Robinson, P. M. 1988. Root-N-consistent semiparametric regression. 
{it:Econometrica} 56: 931-954.


{marker author}{...}
{title:Author}

{pstd}Kerui Du{p_end}
{pstd}Xiamen University{p_end}
{pstd}kerrydu@xmu.edu.cn{p_end}


{title:Also see}

{p 4 14 2}
Manual:  {hi:[R] regress}, {hi:[XT] xtreg}

{p 4 14 2}
Online:  {helpb regress}, {helpb reghdfe}, {helpb xtreg}, {helpb spline}
