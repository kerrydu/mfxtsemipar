{smcl}
{* *! version 1.0  28jul2024}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{vieweralsosee "[R] reghdfe" "help reghdfe"}{...}
{viewerjumpto "Syntax" "mfxtsemipar##syntax"}{...}
{viewerjumpto "Description" "mfxtsemipar##description"}{...}
{viewerjumpto "Options" "mfxtsemipar##options"}{...}
{viewerjumpto "Examples" "mfxtsemipar##examples"}{...}
{viewerjumpto "Stored results" "mfxtsemipar##results"}{...}
{viewerjumpto "References" "mfxtsemipar##references"}{...}
{viewerjumpto "Author" "mfxtsemipar##author"}{...}
{title:Title}

{phang}
{bf:mfxtsemipar} {hline 2} Mixed-frequency semiparametric regression with fixed knots


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:mfxtsemipar}
{depvar} [{indepvars}]
{ifin}
{weight}{cmd:,}
{opt uvar(varname)}
{opt id(varname)}
{opt gen(string)}
{opt tl(varlist)}
{cmd:{c -(}}{opt knots(numlist)} {c |} {opt nknots(integer)}{cmd:{c )-}}
[{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt uvar(varname)}}variable for semiparametric transformation{p_end}
{synopt:{opt id(varname)}}panel identifier variable{p_end}
{synopt:{opt gen(string)}}name for generated fitted values variable{p_end}
{synopt:{opt tl(varlist)}}time-level variables (low-frequency identifiers) for collapse operation{p_end}

{syntab:Knot specification (choose one)}
{synopt:{opt knots(numlist)}}manually specify knot locations{p_end}
{synopt:{opt nknots(integer)}}specify number of knots (automatic placement){p_end}

{syntab:Model}
{synopt:{opt cluster(string)}}cluster variable for robust standard errors{p_end}
{synopt:{opt type(string)}}spline type: {cmd:bs}, {cmd:ms}, {cmd:is}, {cmd:ibs}, or {cmd:poly}; default is {cmd:poly}{p_end}
{synopt:{opt bknots(numlist)}}boundary knots for spline basis{p_end}
{synopt:{opt degree(numlist)}}polynomial degree for splines{p_end}
{synopt:{opt winsor(string)}}winsorization percentiles for uvar{p_end}
{synopt:{opt eqspace}}use equally spaced knots instead of quantile-based{p_end}
{synopt:{opt center(numlist)}}centering value for splines{p_end}
{synopt:{opt absorb(string)}}fixed effects specification{p_end}
{synopt:{opt hfcov(varlist)}}high-frequency covariates{p_end}
{synopt:{opt atu(varname)}}alternative variable for spline evaluation{p_end}
{synopt:{opt intercept}}include intercept in spline basis{p_end}

{syntab:Inference}
{synopt:{opt brep(real)}}number of bootstrap replications for inference; default is {cmd:0}{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{cmd:fweight}s, {cmd:aweight}s, and {cmd:pweight}s are allowed; see {help weight}.{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:mfxtsemipar} estimates a mixed-frequency fixed effects model with 
semiparametric splines. Unlike {helpb mfxtsemipar_cv}, this command requires 
the user to specify the knot locations manually or specify the number of 
knots without cross-validation.

{pstd}
The command fits spline or polynomial functions to a univariate variable 
while controlling for other covariates and fixed effects in a panel data 
setting. This is useful when you have prior knowledge about the appropriate 
number or location of knots, or when computational efficiency is more 
important than optimal knot selection.

{pstd}
The estimation aggregates high-frequency data to the specified time level 
using the collapse operation, then estimates the semiparametric model using 
{helpb reghdfe} for efficient handling of fixed effects.


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

{dlgtab:Knot specification}

{phang}
{opt knots(numlist)} manually specifies the knot locations for the spline
basis. The knots should be specified in ascending order and within the range
of the {cmd:uvar}. This option cannot be used together with {cmd:nknots()}.

{phang}
{opt nknots(integer)} specifies the number of knots to be automatically
placed. By default, knots are placed at quantiles of the {cmd:uvar} distribution.
This option cannot be used together with {cmd:knots()}.

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
{opt bknots(numlist)} specifies the boundary knots for the spline basis.
If not specified, boundary knots are automatically set to slightly outside
the range of the {cmd:uvar}.

{phang}
{opt degree(numlist)} specifies the polynomial degree for splines.
Higher degrees allow for more flexible functional forms but may lead to
overfitting.

{phang}
{opt winsor(string)} specifies winsorization percentiles for the {cmd:uvar}.
For example, {cmd:winsor(1 99)} will winsorize the variable at the 1st and
99th percentiles to reduce the influence of outliers.

{phang}
{opt eqspace} requests that knots be placed at equally spaced intervals
rather than at quantiles of the {cmd:uvar} distribution. This option only
applies when using {cmd:nknots()}.

{phang}
{opt center(numlist)} specifies the centering value for splines. This can
help with numerical stability and interpretation of the coefficients.

{phang}
{opt absorb(string)} specifies the fixed effects to absorb. The syntax
follows that of {helpb reghdfe}, allowing for multiple levels of fixed
effects separated by spaces or hash symbols.

{phang}
{opt hfcov(varlist)} specifies high-frequency covariates that should be
summed during the collapse operation and included in the model.

{phang}
{opt atu(varname)} specifies an alternative variable for spline evaluation.
This allows the spline to be evaluated at different points than those used
for estimation, which can be useful for prediction or policy analysis.

{phang}
{opt intercept} includes an intercept term in the spline basis functions.
By default, the intercept may be dropped to avoid collinearity issues.

{dlgtab:Inference}

{phang}
{opt brep(real)} specifies the number of bootstrap replications for
wild bootstrap inference. When {cmd:brep()} is greater than 0, the command
performs wild bootstrap to obtain robust standard errors that account for
potential heteroskedasticity and clustering. The default is 0 (no bootstrap).


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
{pstd}Basic mixed-frequency semiparametric regression with manual knots{p_end}
{phang2}{cmd:. mfxtsemipar y x1 x2, uvar(x3) id(id) gen(fitted) tl(t) knots(-1 0 1) center(0)}{p_end}

    {hline}
{pstd}Automatic knot placement with specified number{p_end}
{phang2}{cmd:. mfxtsemipar y x1 x2, uvar(x3) id(id) gen(auto_fit) tl(t) nknots(4) center(0)}{p_end}

    {hline}
{pstd}With fixed effects and clustering{p_end}
{phang2}{cmd:. mfxtsemipar y x1 x2, uvar(x3) id(id) gen(fe_fit) tl(t) nknots(3) center(0) absorb(id) cluster(id)}{p_end}

    {hline}
{pstd}Using polynomial splines with higher degree{p_end}
{phang2}{cmd:. mfxtsemipar y x1 x2, uvar(x3) id(id) gen(poly_fit) tl(t) knots(-1.5 -0.5 0.5 1.5) type(poly) degree(2) center(0)}{p_end}

    {hline}
{pstd}With evaluation at specific grid points{p_end}
{phang2}{cmd:. gen x3grid = -3 + _n/10 + 0.1 if _n <= 61}{p_end}
{phang2}{cmd:. mfxtsemipar y x1 x2, uvar(x3) id(id) gen(grid_fit) tl(t) nknots(4) center(0) atu(x3grid)}{p_end}

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
{pstd}Estimate sinusoidal function with manual knots{p_end}
{phang2}{cmd:. mfxtsemipar y x1 x2, uvar(x3) id(id) gen(sin_fit) tl(t) knots(-2 -1 0 1 2) center(0) type(poly) degree(2)}{p_end}

    {hline}
{pstd}With bootstrap inference{p_end}
{phang2}{cmd:. mfxtsemipar y x1 x2, uvar(x3) id(id) gen(boot_fit) tl(t) nknots(4) center(0) brep(100) cluster(id)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:mfxtsemipar} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(knots)}}knot locations used in estimation{p_end}
{synopt:{cmd:e(splinecmd)}}command used to generate splines{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(info)}}model information criteria{p_end}

{pstd}
Additionally, {cmd:mfxtsemipar} returns all standard results from {helpb reghdfe},
including coefficient estimates, standard errors, and fit statistics.


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

{phang}
Wood, S. N. 2017. {it:Generalized Additive Models: An Introduction with R}. 
2nd ed. Boca Raton, FL: CRC Press.


{marker author}{...}
{title:Author}

{pstd}Kerui Du{p_end}
{pstd}Xiamen University{p_end}
{pstd}kerrydu@xmu.edu.cn{p_end}


{title:Also see}

{p 4 14 2}
Manual:  {hi:[R] regress}, {hi:[XT] xtreg}

{p 4 14 2}
Online:  {helpb mfxtsemipar_cv}, {helpb regress}, {helpb reghdfe}, {helpb xtreg}, {helpb spline}
