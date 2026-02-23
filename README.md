
<!-- This brief doc/description is only for the github page. README.md is generated from README.Rmd. Please edit the latter file - rmarkdown::render('README.Rmd', output_format = 'github_document', output_file = 'README.md') -->
<!-- SetUp process create_package("~/Documents/Synched/Uni/Work/2020/SpatialInference"), use_readme_rmd(), use_git() -->
<!--  -->
<!-- Get in Rcpp http://adv-r.had.co.nz/Rcpp.html#rcpp-package -->
<!-- Get in conley https://github.com/potterzot/vcovConley -->

# SpatialInference

<!-- badges: start -->
<!-- badges: end -->

Package containes helper functions to do inference with spatial data in
R. Most importantly, the fast C++ implementation for Conley (spatial
HAC) standard errors by Darin Christensen and a function to estimate the
correlation range of regression residuals.

## Some note on pkg development

- for this readme on Github, edit the .Rmd and then use
  `rmarkdown::render()`
- after you edited or included a function, just use the cmd+shift+B
  short cut to install an rebuild the package (rinse and repeat)

### Functions (actual and sheduled)

<!-- * grid_identifier() [for the FE] -->
<!-- * put lon/lat (or the current CRS coord) as columns in df [to control for "smooth position in space"] -->
<!-- * -->
<!-- * decluster() -->
<!-- * -->
<!-- * spatbag() -->
<!-- * spat_did() [boils down to spec of the term with W] -->
<!-- * SpatFD from Martin73 [again just figuring out the Wlagged vectors and substract the lag from the y] -->
<!-- * -->
<!-- * conley() -->
<!-- * conleySANDWICHfiona? -->
<!-- * conleyGMM? -->
<!-- * kelly() -->
<!-- * acreg() -->
<!-- * sphac_strap()? -->
<!-- * kimsun11 (not implemented anywhere acc to KP17) -->
<!-- * Bester et al 16 fixed-b assymptot? -->
<!-- *  -->
<!-- * block bootstraps? (conley19,...) -->
<!-- * -->
<!-- * robust Morans I -->
<!-- * -->
<!-- * center of gravity (https://twitter.com/undertheraedar/status/1285904073825366016/photo/1) [from the quant geo textbook (arlinghaus?), using the ArcGIS formula] -->

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
#install.packages("devtools")
#devtools::install_github("axlehner/SpatialInference")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SpatialInference)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

# Known Errors

Fortran compiler on MaxOS:
<https://github.com/RubD/Giotto_site/issues/11>

worked without the last one: the pointer to the direction
