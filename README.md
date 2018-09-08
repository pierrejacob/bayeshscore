# bayeshscore
**Bayesian model comparison using the Hyvärinen score**.  
*\[more details in S. Shao, P.E. Jacob, J. Ding, and V. Tarokh (2017), available at [https://arxiv.org/abs/1711.00136](https://arxiv.org/abs/1711.00136)\]*

This package provides functions that compute the **Hyvärinen score** (and the log-evidence as an aside) of a Bayesian model. This is achieved by using either **SMC** \[e.g. [Chopin (2002)](https://academic.oup.com/biomet/article-abstract/89/3/539/251804) and [Del Moral, Doucet, Jasra (2006)](http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2006.00553.x/abstract)\] or **SMC<sup>2</sup>** \[cf. [Chopin, Jacob, Papaspiliopoulos (2013)](http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2012.01046.x/abstract)\] depending on whether the likelihood can be evaluated.

## Table of contents
1. [Computing the Hyvärinen score](#howto_hscore)
2. [Defining a model](#howto_model)
3. [Setting algorithmic parameters](#howto_algoparam)
4. [Output description](#output)
5. [Particle filter](#particlefilter)
6. [Examples](#examples)

### <a name="howto_hscore"></a> Computing the Hyvärinen score

After installing and loading the package, computing the Hyvärinen score (abbrv. H-score) of a model is done by calling
```R
hscore(observations, model, algorithmic_parameters)
```
where the inputs consist of
* `observations`: sequence of dimY-dimensional observations Y<sub>1</sub> , ... , Y<sub>T</sub> encoded as a dimY by T `matrix`
* `model`: Bayesian model encoded as a `list` (see [Defining a model](#howto_model))
* `algorithmic_parameters` : algorithmic parameters provided as a `list` (see [Setting algorithmic parameters](#howto_algoparam))

The function `hscore` is a wrapper that either calls `smc(observations, model, algorithmic_parameters)` or `smc2(observations, model, algorithmic_parameters)` depending on whether the `model` specifies the likelihood or not.  
The output is detailed below (see [Output description](#output)).

### <a name="howto_model"></a> Defining a model
The `model` needs to be provided as a `list`. A description of the required fields is provided in `inst/_template_model.R`. Missing fields may be automatically filled-in via the function `set_default_model` defined in `R/util_default.R`.

### <a name="howto_algoparam"></a> Setting algorithmic parameters
The `algorithmic_parameters` need to be provided as a `list`. A complete description of the required fields along with their default values can be found in the function `set_default_algorithmic_parameters` defined in `R/util_default.R`.

### <a name="output"></a> Output description
The output of `hscore`, `smc`, or `smc2` is a `list`. Depending on the specified `algorithmic_parameters`, some of its fields may be set to `NULL`. In its most exhaustive form, the output will contain the following objects:
* `thetas` : last set of particles thetas (dimtheta by Ntheta `matrix`)
* `normw` : corresponding normalized weights (`vector` of length Ntheta)
* `byproducts` or `PFs` : corresponding `list` of byproducts (e.g. particle filters for `smc2`, see [Particle filter](#particlefilter))
* `logtargetdensities` : corresponding evaluations of target log-densities (`vector` of length Ntheta)
* `thetas_history` : `list` of successive sets of particles thetas (one `matrix` per timestep, starting with the prior, so the length of the `list` is T+1)
* `normw_history` : `list` of successive normalized weights (i.e. `list` of `vector`)
* `logtargetdensities_history` : `list` of successive target log-densities evaluations (i.e. `list` of `vector`)
* `byproducts_history` of `PF_history` : `list` of successive byproducts or particle filters (i.e. `list` of `list`)
* `logevidence` : cumulative logevidence (`vector` of length T)
* `hscore` : cumulative H-score using Fisher/Louis type identities (`vector` of length T)
* `hscoreDDE` : cumulative H-score using kernel density estimation (`vector` of length T)
* `ESS` : successive ESS (`vector` of length T)
* `rejuvenation_times` : successive times when rejuvenation occured (`vector` of random length less than T)
* `rejuvenation_rate` : associated acceptance rates (`vector` of random length less than T)
* `method` : method called (`string` equal to either `"SMC"` or `"SMC2"`)
* `algorithmic_parameters` : `list` of algorithmic parameters used

### <a name="particlefilter"></a> Particle filter
The particle filter is implemented in `R/conditional_particle_filter.R`. The function `conditional_particle_filter` uses the bootstrap particle filter. It takes an optional argument `path` when conditioning on a particular path is needed. The output is a `list` containing the number of x-particles (`Nx`), the last set of x-particles (`X`), their respective normalized weights (`xnormW`), an estimator of the log-likelihood (`log_p_y_hat`), the last incremental log-likelihood (`incremental_ll`), a tree encoding the paths and ancestral lineages (`tree`), and a realized path sampled from the final sets of paths (`path`).

### <a name="examples"></a> Examples
Some examples are provided in the folder `inst`. These scripts reproduce the figures presented in [S. Shao, P.E. Jacob, J. Ding, and V. Tarokh (2017)](https://doi.org/10.1080/01621459.2018.1518237) and its [supplementary material](https://www.tandfonline.com/doi/suppl/10.1080/01621459.2018.1518237). Open access to the manuscript and its supplement is also available [here](https://stephaneshao.github.io/publication/hyvarinen/).
