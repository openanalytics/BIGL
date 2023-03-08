# BIGL 1.7.0 (2023-03-08)

- switch from `rgl` to `plotly` for the 3D plots
- move NEWS to NEWS.md

# BIGL 1.6.7 (2023-01-19)

- replace `rgl.*` functions due to upcoming deprecation in the `rgl` package

# BIGL 1.6.6 (2022-05-04)

- skip further tests on CRAN that may fail randomly

# BIGL 1.6.5 (2021-10-26)

- skip some tests on CRAN that may fail randomly

# BIGL 1.6.4 (2021-09-01)

- Fix vignettes to avoid CRAN warnings on systems without pandoc or webshot2

- Make sure NAs in bootConfint are handled

# BIGL 1.6.3 (2021-06-29)

- Fix code for compatibility with new ggplot2 version

# BIGL 1.6.2 (2021-03-10)

- Fix issue with testing when variance at off-axis points equals zero

# BIGL 1.6.0 (2020-12-16)

- Residual resampling according to mean-variance trend

- Use bootstrapped confidence intervals for point-by-point and overall effect sizes

- Internal redesign for speed-up, using "d1d2" identifiers

# BIGL 1.5.3 (2020-09-10)

- fix issue with `generalizedLoewe` and `hsa` functions when data is a tibble

# BIGL 1.5.2 (2020-07-24)

- fix issue with smooth monotherapy plot when dotted line was not shown due to numeric inaccuracies

# BIGL 1.5.1 (2020-07-02)

- show smooth prediction curve on the marginal fit plot

- show warning for agonist/antagonist combination when using generalized Loewe model

# BIGL 1.4.3 (2020-01-30)

- fixes for special case of 2 flat monotherapy profiles

# BIGL 1.4.2 (2020-01-17)

- even better support for flat monotherapies, consistent calls/colors for all representations (contour plot, summary table and 3d plot)

- added clarification to description of Bliss method in the "Methodology" vignette

- fixed typos

# BIGL 1.4.1 (2019-06-20)

- improve tests

# BIGL 1.4.0 (2019-06-18)

- added Alternative Loewe independence null model

# BIGL 1.3.0 (2019-02-28)

- added Bliss independence null model

# BIGL 1.2.3 (2018-11-21)

- better support for flat monotherapies

# BIGL 1.2.2 (2018-09-07)

- make sure if no replicates are there model/unequal method is not allowed

# BIGL 1.2.1 (2018-08-24)

- fix regression bug to allow custom axis labels for the contour plot

# BIGL 1.2.0 (2018-06-13)

- support heterogeneous variance with new argument 'method' for 'fitSurface'

# BIGL 1.1.3 (2018-05-14)

- show EC50 in log10 terms for monotherapy summary

- allow custom compound names for plotting

# BIGL 1.1.2 (2018-05-07)

- fix monotherapy formula in the vignette

# BIGL 1.1.1 (2018-01-02)

- added reference to the original article

- fix DESCRIPTION file according to new R CMD check policy

# BIGL 1.1.0 (2017-09-22)

- isobologram() plot is now available for ResponseSurface objects It is also now included in the shiny app.

- summary.ResponseSurface does not fail if marginal fits have constraints

- Change of package maintainer

# BIGL 1.0.3 (2017-08-18)

- Contour plot for maxR now has correct p-value labels

- Summary for marginal fits does not fail if Hessian computation failed

# BIGL 1.0.2 (2017-07-31)

- extra arguments in fitMarginals are now passed to the optimizer functions

# BIGL 1.0.1 (2017-06-29)

- `predict` method added for the MarginalFit objects

# BIGL 1.0 (2017-06-28)

- Initial public release

