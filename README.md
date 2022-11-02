## Loading the model

This repository contains the R code (1) and associated data for the
zoonotic hookworm transmission dynamics model described in Walker et al.
(2). The model requires installation of the deSolve (3) package and can
is loaded, with default parameter values, as follows:

    library("deSolve")
    source("hookworm.R")
    source("pars.R")

## Running the model to endemic equilibrium

The model is best used by running simulations to equilibrium and running
interventions using the equilibrium (baseline) value for the mean number
of worms per host as the initial value.

Running the model to equilibrium requires only that parameter
`pars["equib"]<-1` (the default setting):

    pars["equib"]<-1
    base <- funcs$runmod(pars)

The first two coloums of the output `base` contain the mean number of
worms (worm burden) in host 1 (dogs) and host 2 (humans):

    head(base[,1:2])

    ##    time        W1
    ## 1     0  5.000000
    ## 13    1  8.046955
    ## 25    2 11.304249
    ## 37    3 14.524048
    ## 49    4 17.535181
    ## 61    5 20.241719

## Simulating mass drug administration

Mass drug adminsitration (MDA) is simulated by using the final rows of
`base[nrow(base),1:2]` (the equilibrium worm burdens) as initial values
for a second run. The start time, stop time, frequency, duration,
coverage and treatment efficacy are set using parameters:
`pars["start.tx"]`; `pars["stop.t"]`; `pars["freq.tx1"]` and
`pars["freq.tx2"]`; `pars["n.tx1"]` and `pars["n.tx2"]`; `pars["c1"]`
and `pars["c2"]`, `pars["epsilon1"]` and `pars["epsilon2"]`. Hence, to
run 8 years of human-only annual MDA at 75% coverage using the default
efficacy of 90% we can set the parameters as:

    pars["equib"] <- 0
    pars["start.tx"] <- 1
    pars["stop.t"] <- 10
    pars["freq.tx2"] <- 1
    pars["n.tx2"] <- 8
    pars["c2"] <- 0.75
    pars["c1"] <- 0

We also need to reduce the integration step size to avoid numerical
errors during the MDA phase:

    pars["dt"] <- 0.005
    pars["dtout"] <- pars["dt"]*20

We can now run the model:

    HOMDA <- funcs$runmod(pars,inits=as.numeric(base[nrow(base),2:3]))

The prevalence and effective reproduction numbers and stored in
`HOMDA[,"Wp1"]`, `HOMDA[,"Wp2"]` and `HOMDA[,"RE"]`:

![](README_files/figure-markdown_strict/plotHO-1.png)

Running a One Health intervention with 50% coverage in dogs requires
setting up the MDA parameters for dogs (using the same endemic
equilibrium starting values):

    pars["freq.tx1"] <- 1
    pars["n.tx1"] <- 8
    pars["c1"] <- 0.5
    OHMDA <- funcs$runmod(pars,inits=as.numeric(base[nrow(base),2:3]))

![](README_files/figure-markdown_strict/plotOH-1.png)

## References

<span class="csl-left-margin">1. </span><span class="csl-right-inline">R
Core Team. *R: A language and environment for statistical computing*.
Vienna, Austria: R Foundation for Statistical Computing (2021).
<https://www.R-project.org/></span>

<span class="csl-left-margin">2. </span><span
class="csl-right-inline">Walker M, Lambert S, Grillo MB, Neve MI,
Worsley A, Traub R, Colella V. Modelling the effectiveness of one health
interventions against zoonotic hookworm. *Front Med* (2022) under
review:</span>

<span class="csl-left-margin">3. </span><span
class="csl-right-inline">Soetaert K, Petzoldt T, Setzer RW. Solving
differential equations in R: Package deSolve. *J Stat Softw* (2010)
33:1–25. doi:
[10.18637/jss.v033.i09](https://doi.org/10.18637/jss.v033.i09)</span>
