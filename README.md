---
title: "Case I interval censored data with the Addams family"
output: pdf_document
date: "2024-05-02"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## The Addamfamily package

The Addamfamily packge contains routines to fit univariate and mutlitvariate time-to-event models. 
Its main focus are multivariate frailty models. 
Among the available distributions are the $\mathcal{AF}$ and the power variance family (PVF) as well as their special cases. 
Routines are available for current status data (also known as case I interval censored data) and right censored data. 
Clusters do not need to be balanced. A model on the frailty distribution parameters can be imposed along multilevel factors.
The package also contains routines to estimate shifted $\mathcal{NB}$, $\mathcal{B}$, $\mathcal{P}$ 
models for current status data. 
The baseline hazard can be Breslow, piecewise constant or the generalized gamma ($\text{gen-}\mathcal{G}$) distribution as well as 
its special cases. Breslow hazards are not available for current status data, but piecewise hazards might 
be adapted to mimic the non-parametric approach. 
We do not made good experience with the $\text{gen-}\mathcal{G}$ for current status data so far. 
Covariates might also be included in common proportional hazard manner or by stratifying baseline hazards. 
Predict routines are available but not yet put in the predict command framework.

Please note that this package is under development. 
We have full confidence in the numeric results. 
However, it was a project that grew rather unexpectedly and hence, the coding style is certainly a bit inefficient. 
In the moment the convergence criterion is an absolute change in the likelihood. 
Hence one reason for a very slow optimization process is an exaggregated convergence criterion. 
Please adapt where necessary.
```{r conv, echo= TRUe, eval = TRUE}
converge <- 1e-1
```
Basically all optimization routines that are available in R are available for optimization.
Fitting procedures will issue many warnings which is on purpose for the development phase. 

This overview might become a future vignette. 
For the time being, we fit some models in order to illustrate the functionality of the Addamfamily package. 
For questions, bug reports, hints, assistance etc please contact me via maximilian.bardo\@proton.me.

# Installation

The package can be installed from git directly 
```{r install_git, echo=TRUE, eval = FALSE}
devtools::install_github("time-to-MaBo/Addamsfamily")
```
or from local machine
```{r install_loc, echo=TRUE, eval = FALSE}
devtools::install.packages('your_path/Addamsfamily', repos = NULL, type="source")
```
where your path has to be added. 

## Data

The Addamsfamily package contains a simulated dataset that is supposed to mimic the HPV data of the XXAX PAper. 
It can be loaded by
```{r HPV, echo=TRUE, eval = TRUE}
library(Addamsfamily)
data(HPV)
head(data)
#HPV_sim <- HPV_sim[!HPV_sim$HPV %in% c('HPV31', 'HPV58', 'HPV33'),] #reduce sample size as optimization takes a while
#HPV_sim$HPV <- droplevels(HPV_sim$HPV)
```

## Models
We will fit some models in order to explain functionality. 
Note that we do not claim that the fitted models make any sense from a content-related perspective. 
There is no formual option in the function so far. 
Hence the arguments have to passed. 
Note that your dataset should only consist of numeric, logical or factor data types. 
For factors only dummies should be passed as we handle factors with as.numeric()-1 at that point.
Consider building more dummies.

### Univariate Models current status models
Unfortunately the current status data models are particular. 
This might be caused by the Likelihood where many combinations of sums of cumulative hazard 
rates need to be computed which is currently solved via a loop. 

#### piecewise constant hazards
A model without frailty and piecewise constant hazards (stratified by HPV type) can be estimated by
```{r HPV, echo=TRUE, eval = TRUE}
uni.p10 <- ph.surv(data = HPV, X = c('sex'), y = 'age', d = 'status', type = 'nofrail', 
        strata = 'HPV', ID = NULL, baseline = list('piecewise','Intervall', by = 10, NULL ), 
        current.status = TRUE, converge = converge  )
```
Note that current.status has to be set to TRUE otherwise right censored models are fitted. 
This is something we might change in future as this is certainly error prone. 
The baseline option as specified above says that there will be a unique parameter for the time from $[0,10), [10,20),...$.
This might be changed by replacing "Intervall" with "Breslow" where the intervals are set between the distinct event times (what does not make so much sense for current status data as the event times are not known, but technically this is no problem). 
The breaks for the distinct hazard parameters might be set manually as follows 

```{r HPV, echo=TRUE, eval = TRUE}
breaks <- list(c(seq(from = 0, to = 10, by = 5), seq(from = 20, to = 40, by = 10), 80))
breaks <- replicate(breaks1, n=7)
#breaks <- breaks[1:4]
uni.pbreaks <- ph.surv(data = HPV, X = c('sex'), y = 'age', d = 'status', type = 'nofrail', 
        strata = 'HPV', ID = NULL, baseline = list('piecewise','Intervall', by = 10, NULL ), 
        current.status = TRUE, converge = converge  )
```
When a list of length one is provided all hazard models will have the same break points. 
For this data set the survival function reaches a plateau at around 30 years. 
Hence, hazard parameters will often be zero and the Hessian will be singular.
Can be handled by H.singular.
If the baseline hazard should not be stratified the argument strata can simply be 
removed or set to NULL (its default value).

The output of the function is currently just a print command. 
We are working on it to build class-specific appropriate summary, predict .. commands. 

#### Generalized Gamma hazards
A $\text{gen-}\mathcal{G}$ model might be fitted by 
```{r GG, echo=TRUE, eval = TRUE}
uni.GG <- ph.surv(data = HPV, X = c('sex'), y = 'age', d = 'status', type = 'nofrail', 
        strata = 'HPV', ID = NULL, baseline = list('GG',NULL,NULL, NULL ), 
        current.status = TRUE, converge = converge  )
```
We did not make particular good experience with the $\text{gen-}\mathcal{G}$ for current status data so far. 
It works better for right censored data. 
To see which distribution has been chosen through the fitting process we can do the following
```{r GGtree, echo=TRUE, eval = TRUE}
uni.GG$GG.tree
```

Note that the the special cases of the GG are also available.
You can find the options in the following vector:
```{recho=TRUE}
GG_options <- ( 'GG', 'gamma', 'gamma^-1', 'Weibull', 'Weibull^-1', 
                             'exponential', 'exponential^-1', 'ammag', 'ammag^-1', 'lognormal', 'location', 'N/2', 'GGpos', 'GGinv' )
```
Also note that you can specify different distribution for different strata. 
Here we have 7 strata, lets assign each of them a different distribution:
```{r GGmix, echo=TRUE, eval = TRUE}
uni.GGmix <- ph.surv(data = HPV, X = c('sex'), y = 'age', d = 'status', type = 'nofrail', 
        strata = 'HPV', ID = NULL, baseline = list(c('Weibull', 'gamma', 'gamma^-1', 'exponential', 'exponential^-1', 'ammag', 'ammag^-1'),NULL,NULL, NULL ), current.status = TRUE  )
```
The options GGpos and GGinv restirct the GG to positive or negative values of the bla parameter.
See the tutorial of Cox on the $\text{gen-}\mathcal{G}$ distribution. We utilized the same parameterization here. 
The option "location" shares all parameters among the strata ecept for the location parameter of the $\text{gen-}\mathcal{G}$ where each stratum has its own parameter. The strata among which the parameters are shared can also be put in group such that not all have the same bla and bla parameter.

Now let us estimate some frailty models
### Frailty models
Note that initial values for the parameters can be set. See the help. 
Further optimization options for current status can be set by
```{r echo = TRUE, eval = FALSE}
pp.options.cs <- list(EM.cs=FALSE, univariate=TRUE,uni.model=NA,uni.initial=TRUE)
```
where the fist element of  the list specifies whether the EM algorithm should be used which does not 
work well in the case of current status data. The second option specifies whether a 
univariate model should also be fitted, the third option can be the aformentioned univariate model if this 
has already been computed and the last option specifies whether PH and baseline paramers should be intialized with 
the estimates from the univariate model. 

When fitting frailty models it is often important to handle the threshold until 
one of its special cases (including no frailty) is assumed to be the choice of 
the current iteration to handle that case specifically. This can be done by the option 
frail.thresh which is 1e-6 by default.
If an error message like "missing value where TRUE/FALSE needed" is issued, 
frail.thresh will probably do the job (e.g. setting to 1e-5 or 1e-4). 
#### $\mathcal{AF}$
The option type specifies the frailty distribution. 
The option type specifies the frailty distribution. It can be of length(unique(frailty)) when 
different distributional choices and not just the frailty parameters should be stratified. 
A $\mathcal{AF}$ model can be fitted as follows:
```{r Addams, echo=TRUE, eval = TRUE}
ag.p <- ph.surv(data = HPV_sim, X = 'sex', y = 'age', d = 'status', ID = 'ID',strata = 'HPV',  
                current.status = TRUE, converge = converge, H.singular = 'gcholesky',
                baseline = list('piecewise','Intervall',by,NULL), type = 'alpha_gamma', 
                centerX = FALSE, prnt = TRUE,breaks=breaks,cs.options = pp.options.cs, 
                frail.thresh = 1e-5)
pp.options.cs$uni.model <- ag.p$cond
pp.options.cs$univariate <- FALSE
```
To look at the parameters and to see which distribution was chosen the following can be done:
```{r Addamspars, echo=TRUE, eval = TRUE}
ag.p$theta
ag.p$Z.dist
ag.p$lam_0
ag.p$T.dist
ag.p$beta
```
## scaled Poisson model
```{r Poisson, echo=TRUE, eval = TRUE}
a.p <- ph.surv(data = HPV_sim, X = 'sex', y = 'age', d = 'status', ID = 'ID',strata = 'HPV',  
                current.status = TRUE, converge = converge, H.singular = 'gcholesky',
                baseline = list('piecewise','Intervall',by,NULL), type = 'alpha', 
                centerX = FALSE, prnt = TRUE,breaks=breaks,cs.options = pp.options.cs, 
               frail.thresh = 1e-5 )

```
## gamma model
```{r gamma, echo=TRUE, eval = TRUE}
a.p <- ph.surv(data = HPV_sim, X = 'sex', y = 'age', d = 'status', ID = 'ID',strata = 'HPV',  
                current.status = TRUE, converge = converge, H.singular = 'gcholesky',
                baseline = list('piecewise','Intervall',by,NULL), type = 'gamma', 
                centerX = FALSE, prnt = TRUE,breaks=breaks,cs.options = pp.options.cs, 
               frail.thresh = 1e-5 )

```
## PVF model
```{r PVF, echo=TRUE, eval = TRUE}
PVF.p <- ph.surv(data = HPV_sim, X = 'sex', y = 'age', d = 'status', ID = 'ID',strata = 'HPV',  
                current.status = TRUE, converge = converge, H.singular = 'gcholesky',
                baseline = list('piecewise','Intervall',by,NULL), type = 'PVF', 
                centerX = FALSE, prnt = TRUE,breaks=breaks,cs.options = pp.options.cs, 
                frail.thresh = 1e-5 )

```
## IG model
```{r IG, echo=TRUE, eval = TRUE}
IG.p <- ph.surv(data = HPV_sim, X = 'sex', y = 'age', d = 'status', ID = 'ID',strata = 'HPV',  
                current.status = TRUE, converge = converge, H.singular = 'gcholesky',
                baseline = list('piecewise','Intervall',by,NULL), type = 'IG', 
                centerX = FALSE, prnt = TRUE,breaks=breaks,cs.options = pp.options.cs, 
                frail.thresh = 1e-5 )

```
## stratified frailty model
```{r stratifiedfrialty, echo=TRUE, eval = TRUE}
agPVF.p <- ph.surv(data = HPV_sim, X = 'sex', y = 'age', d = 'status', ID = 'ID',strata = 'HPV',  
                current.status = TRUE, converge = converge, H.singular = 'gcholesky',
                baseline = list('piecewise','Intervall',by,NULL), 
               frailty = 'sex', type = c('alpha_gamma','IG'), 
                centerX = FALSE, prnt = TRUE,breaks=breaks,
               cs.options = pp.options.cs, frail.thresh = 1e-5 )

```
Note that the distribution itself has to be the same for all strata (type).
For right censored data the distribution (and not jsut the distribution parameters) 
is allowed to differ across the frailty strata.

#### Overdispersion models
Overdispersion parameters for current status data moels can be computed by setting overdispersion = TRUE.

There is more functionality which supposed to be documented soon.

## Right censored data
This is in preparation. Additionally to the above hazard Breslow hazards can be chosen. 
Synthax is basically the same as for current status data except for setting current.status to FALSE (default).
EM algorithm is chosen for optimization purposes. 
