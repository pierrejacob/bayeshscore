# bayeshscore
**Bayesian model comparison using the Hyvärinen score**.  
*\[more details in Shao, Jacob, Ding, Tarokh (2017) at [https://arxiv.org/abs/1711.00136](https://arxiv.org/abs/1711.00136)\]*

This package provides functions that compute the **Hyvärinen score** (and the log-evidence as an aside) of a Bayesian model. This is achieved by using either **SMC** \[e.g. [Chopin (2002)](https://academic.oup.com/biomet/article-abstract/89/3/539/251804) and [Del Moral, Doucet, Jasra (2006)](http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2006.00553.x/abstract)\] or **SMC<sup>2</sup>** \[cf. [Chopin, Jacob, Papaspiliopoulos (2013)](http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2012.01046.x/abstract)\], depending on whether the likelihood can be evaluated.

### <a name="howto_hscore"></a> Computing the Hyvärinen score

After installing and loading the package, computing the Hyvärinen score of a model is done by calling
```R
hscore(observations, model, algorithmic_parameters)
```
where the inputs consist of
* `observations`: sequence of d<sub>Y</sub>-dimensional observations Y<sub>1</sub> , ... , Y<sub>T</sub> encoded as a d<sub>Y</sub> by T `matrix`
* `model`: Bayesian model encoded as a `list` ([learn more](#howto_model))
* `algorithmic_parameters` : algorithmic parameters provided as a `list` ([learn more](#howto_algoparam))

The function `hscore` is essentially a wrapper that calls either 
```R
smc(observations, model, algorithmic_parameters)
```


### <a name="howto_model"></a> Defining a model

### <a name="howto_algoparam"></a> Setting algorithmic parameters
