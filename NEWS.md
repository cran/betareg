# betareg 3.2-1

* New working paper "Extended-Support Beta Regression for [0, 1] Responses"
  by Ioannis Kosmidis and Achim Zeileis in the _arXiv.org E-Print Archive_,
  [doi:10.48550/arXiv.2409.07233](https://doi.org/10.48550/arXiv.2409.07233).

* New package web page (via `altdoc`/`quarto`) at
  <https://topmodels.R-Forge.R-project.org/betareg/>

* Extended functionality of `predict()` method for `betareg` objects and
  enhanced the corresponding documentation, see `?predict.betareg`.

* Turned `vignette("betareg", package = "betareg")` and
  `vignette("betareg-ext", package = "betareg")` from Sweave into Quarto
  vignettes. Some improvements/updates in the text.


# betareg 3.2-0

* Major extension in `betareg()`: In addition to classic beta regression for
  responses in the open interval (0, 1), extended-support beta regression is
  added which can model responses in the closed interval [0, 1] (i.e., including
  boundary observations at 0 and/or 1). This is accomplished by adding two new
  response distributions:  The extended-support beta distribution (`"xbeta"`)
  leverages an underlying symmetric four-parameter beta distribution with
  exceedence parameter `nu` to obtain support [-nu, 1 + nu] that is subsequently
  censored to [0, 1] in order to obtain point masses at the boundary values 0
  and 1. The extended-support beta mixture distribution  (`"xbetax"`) is a
  continuous mixture of extended-support beta distributions where the exceedence
  parameter follows an exponential distribution with mean `nu` (rather than a
  fixed value of `nu`). The latter `"xbetax"` specification is used by default
  in case of boundary observations at 0 and/or 1. The `"xbeta"` specification
  with fixed `nu` is mostly for testing and debugging purposes.

* Quantile residuals are added to the `residuals()` method for `betareg` objects.
  They are easy to compute and have good distributional properties. Hence,
  they are the new default residuals.

* Bug fix in `pseudo.r.squared` computation for weighted models where previously
  the weights were erroneously ignored (reported by Ray Tayek).

* Bug fixes in `betatree()`: Split points were computed incorrectly due to wrong
  sign of the log-likelihood (reported by Se-Wan Jeong). And trees with only
  intercepts for both `mu` and `phi` could not be fitted (reported by Ludwig
  Hothorn).


# betareg 3.1-4

* In `betatree()` the `"xlevels"` attribute from `partykit::mob` is now correctly
  stored in `$levels` (rather than `$xlevels`) of the returned object.


# betareg 3.1-3

* Added `IGNORE_RDIFF` flags in some examples in order to avoid showing
  diffs due to small numeric deviations in some checks (especially on CRAN).


# betareg 3.1-2

* Added `suppressWarnings(RNGversion("3.5.0"))` in those places where `set.seed()`
  was used to assure exactly reproducible results from R 3.6.0 onwards.


# betareg 3.1-1

* Conditional registration of `sctest()` method for `betatree` objects when
  `strucchange` package is loaded.


# betareg 3.1-0

* The `betatree()` function now uses the new `mob()` implementation from the
  `partykit` package (instead of the old `party` package). The user interface
  essentially remained the same but now many more options are available through
  the new `mob()` function. The returned model object is now inheriting from
  `modelparty`/`party`.

* Included `grDevices` in Imports.

* Fixed `model.frame()` method for `betareg` objects which do not store the
  model frame in `$model`.

* `betamix()` gained arguments `weights` (case weights for observations) and
  `offset` (for the mean linear predictor).


# betareg 3.0-5

* The `Formula` package is now only in Imports but not Depends (see below).

* Method `FLXgetModelmatrix` for `FLXMRbeta` objects modified due to
  changes in `flexmix` 2.3.12. 


# betareg 3.0-4

* For some datasets `betareg()` would just "hang" because `dbeta()` "hangs"
  for certain extreme parameter combinations (in current R versions).
  `betareg()` now tries to catch these cases in order to avoid the problem.
  
* Depends/Imports/Suggests have been rearranged to conform with current
  CRAN check policies. This is the last version of `betareg` to have the
  `Formula` package in Depends - from the next version onwards it will
  only be in Imports.


# betareg 3.0-3

* The `predict()` method gained support for `type = "quantile"`, so that
  quantiles of the response distribution can be predicted.

* The `Formula` package is now not only in the list of dependencies
  but is also imported in the `NAMESPACE`, in order to facilitate
  importing `betareg` in other packages.


# betareg 3.0-2

* Avoid `.Call()`-ing logit link functions directly, instead
  use elements of `make.link("logit")`.


# betareg 3.0-1

* Small consistency updates in labeling coefficients for
  current R-devel.


# betareg 3.0-0

* New release accompanying the second JSS paper: "Extended Beta
  Regression in R: Shaken, Stirred, Mixed, and Partitioned" by Gruen,
  Kosmidis, and Zeileis which appears as Journal of Statistical
  Software 48(11). See also `citation("betareg")`. The paper presents
  the recently introduced features: bias correction/reduction in
  `betareg()`, recursive partitioning via `betatree()`, and finite
  mixture modeling via `betamix()`. See also `vignette("betareg-ext",
  package = "betareg")` for the vignette version within the package.


# betareg 2.4-1

* Formula interface for `betamix()` changed to allow for three parts
  in the right hand side where the third part relates to the
  concomitant variables.

* Modified the internal structure of vignettes/tests. The original
  vignettes are now moved to the vignettes directory, containing also
  .Rout.save files. Similarly, an .Rout.save for the examples is
  added in the tests directory.


# betareg 2.4-0

* Support bias-corrected (BC) and bias-reduced (BR) maximum likelihood
  estimation of beta regressions. See the `type` argument of `betareg()`.
  To enable BC/BR, an additional Fisher scoring iteration was added
  that (by default) also enhances the usual ML results.
  
* New `vignette("betareg-ext", package = "betareg")` introducing BC/BR
  estimation along with the recent additions beta regression trees and latent
  class beta regression (aka finite mixture beta regression models).

* Enabled fitting of beta regression models without coefficients in the
  mean equation.

* Enabled usage of offsets in both parts of the model, i.e., one can use
  `betareg(y ~ x + offset(o1) | z + offset(o2))` which is also equivalent to
  `betareg(y ~ x | z + offset(o2), offset = o1)`, i.e., the `offset`
  argument of betareg is employed for the mean equation only. Consequently,
  `betareg_object$offset` is now a list with two elements (`mean`/`precision`).

* Added warning and ad-hoc workaround in the starting value selection
  of `betareg.fit()` for the precision model. If no valid starting value can be
  obtained, a warning is issued and `c(1, 0, ..., 0)` is employed.

* Added `betareg_object$nobs` in the return object containing the number
  of observations with non-zero weights. Then `nobs()` can be used to extract
  this and consequently `BIC()` can be used to compute the BIC.


# betareg 2.3-0

* New `betatree()` function for beta regression trees based
  on model-based recursive partitioning. `betatree()` leverages
  the `mob()` function from the `party` package. For enabling this
  plug-in, a `StatModel` constructor `betaReg()` is provided
  based on the `modeltools` package.

* New `betamix()` function for latent class beta regression, or
  finite mixture beta regression models. `betamix()` leverages the
  `flexmix()` function from the `flexmix` package. For enabling this
  plug-in, the driver `FLXMRbeta()` is provided.

* Added tests/vignette-betareg.R based on the models fitted
  in `vignette("betareg", package = "betareg")`.


# betareg 2.2-3

* The `"levels"` element of a `betareg` object is now a list
  with components `"mean"`, `"precision"`, and `"full"` to match
  the `"terms"` of the object.
  
* Improved data handling bug in `predict()` method.


# betareg 2.2-2

* Documentation updates for `?gleverage`.


# betareg 2.2-1

* Package now published in Journal of Statistical Software,
  see https://www.jstatsoft.org/v34/i02/
  and `citation("betareg")` within R.

* Bug fix and improvements in `gleverage()` method for `betareg`
  objects: Analytic second derivatives are now used and
  variable dispersion models are handled correctly.


# betareg 2.2-0

* `dbeta(..., log = TRUE)` is now used for computing the
  log-likelihood which is numerically more stable
  than the previous hand-crafted version.
  
* The starting values in the dispersion regression are
  now chosen differently, resulting in a somewhat more
  robust specification of starting values. The intercept
  is computed as described in Ferrari & Cribari-Neto
  (2004), plus a link transformation (if any). All further
  parameters (if any) are initially set to zero. See also
  the vignette for details.
  
* Various documentation improvements, especially in the
  vignette.


# betareg 2.1-2

* New vignette (written by Francisco Cribari-Neto and Z)  
  introducing the package and replicating a range of
  publications related to beta regression:
  `vignette("betareg", package = "betareg")`
  provides some theoretical background, a discussion of the
  implementation and several hands-on examples.

* Implemented an optional precision model, yielding
  variable dispersion. The precision parameter `phi` may
  depend on a linear predictor, as suggested by
  Simas, Barreto-Souza, and Rocha (2010). In single part
  formulas of type `y ~ x1 + x2`, `phi` is by default assumed to
  be constant, i.e., an intercept plus identity link. But
  it can be extended to `y ~ x1 + x2 | z1 + z2` where `phi`
  depends on `z1 + z2`, by default through a log link.    

* Allowed all link functions (in mean model) that are
  available in `make.link()` for binary responses, and added
  log-log link.

* Added data and replication code for Smithson & Verkuilen
  (2006, Psychological Methods). See `?ReadingSkills`,
  `?MockJurors`, `?StressAnxiety` as well as the complete
  replication code in `demo("SmithsonVerkuilen2006")`.

* Default in `residuals()` (as well as in the related `plot()`
  and `summary()` components) is now to use standardized
  weighted residuals 2 (`type = "sweighted2"`).


# betareg 2.0-0

* Package `betareg` was orphaned on CRAN, Z took over
  as maintainer, ended up re-writing the whole package.
  The package still provides all functionality as before
  but the interface is not fully backward-compatible.
  
* `betareg()`: More standard formula-interface arguments;
  `betareg` objects do _not_ inherit from `lm` anymore.

* `betareg.fit()`: Renamed from `br.fit()`, enhanced interface
  with more arguments and returned information. Untested
  support of weighted regressions is enabled.
  
* `betareg.control()`: New function encapsulating control
  of `optim()`, slightly modified default values.

* `anova()` method was removed, use `lrtest()` from `lmtest`
  package instead.
  
* `gen.lev.betareg()` was changed to `gleverage()` method
  (with new generic) and a bug in the method was fixed.
  
* `envelope.beta()` was removed and is now included in
  `plot()` method for `betareg` objects.
  
* Datasets `prater` and `pratergrouped` were incorporated
  into a single `GasolineYield` dataset.

* New data set `FoodExpenditure` from Griffiths et al. (1993),
  replicating second application from Ferrari and Cribari-Neto
  (2004).

* Added `NAMESPACE`.

* The `residuals()` method now has three further types of
  residuals suggested by Espinheira et al. (2008) who recommend
  to use "standardized weighted residuals 2" (`type = "sweighted2"`).
  The default are Pearson (aka standardized) residuals.
