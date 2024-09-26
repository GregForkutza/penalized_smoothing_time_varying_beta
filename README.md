
# Forkutza (2024)

This repository contains the code for the MSc. Thesis:

Greg Forkutza. *Inferring the time-varying transmission rate and
effective reproduction number by fitting semi-mechanistic compartmental
models to incidence data*.

The examples in this thesis can be reproduced by executing any of the
files in the R/Examples directory.

## Abstract

This thesis presents a novel approach to ecological dynamic modeling
using non-stochastic compartmental models. Estimating the transmission
rate ($\beta$) and the effective reproduction number ($R_t$) is
essential for understanding disease spread and guiding public health
interventions. We extend this method to infectious disease models, where
the transmission rate varies dynamically due to external factors. Using
Simon Wood’s partially specified modeling framework, we introduce
penalized smoothing to estimate time-varying latent variables within the
`R` package `macpan2`. This integration provides an accessible tool for
complex estimation problems. The efficacy of our approach is first
validated via a simulation study and then demonstrated with real-world
datasets on Scarlet Fever, COVID-19, and Measles. We infer the effective
reproduction number ($R_t$) using the estimated $\beta$ values,
providing further insights into the dynamics of disease transmission.
Model fit is compared using the Akaike Information Criterion (AIC), and
we evaluate the performance of different smoothing bases derived using
the `mgcv` package. Our findings indicate that this methodology can be
extended to various ecological and epidemiological contexts, offering a
versatile and robust approach to parameter estimation in dynamic models.

# Resources

- [macpan2](https://github.com/canmod/macpan2)

- [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html)
