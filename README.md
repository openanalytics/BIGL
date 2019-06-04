[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-ago/BIGL)](https://cran.r-project.org/package=BIGL)

# Biochemically Intuitive Generalized Loewe Model

`BIGL` package for R can be used in studying synergistic effects between two drugs or compounds. It uses response surface approach to study deviations of the observed effects from a predicted response surface.

To install the package from CRAN, run:

```r
install.packages("BIGL")
```

To install the latest development version from Github, we suggest using `devtools`:

```r
devtools::install_github("OpenAnalytics/BIGL", build_vignettes = TRUE)
```

`BIGL` methodology currently allows generalized Loewe, classical Loewe and Highest Single Agent models for predicted response surface generation. Generalized Loewe approach allows in particular to treat compounds with differing maximal response and is detailed in the accompanying BIGL paper (pending submission).

Scientific workflow is briefly explained in the `methodology` vignette provided with the package. `analysis` vignette will guide the user through the main functionality of the package. If you have built the vignettes for the package during the installation, these can be accessed by

```r
## User guide
vignette("analysis", package = "BIGL")

## Methodology overview
vignette("methodology", package = "BIGL")
```

If you have `shiny` and `DT` packages installed, you might be interested in trying out the demo for the type of response surfaces that the package works with.

```r
library(BIGL)
runBIGL()
```

Alternatively, the demo can be accessed [here](https://bigl.openanalytics.eu/app/bigl).
