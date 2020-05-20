# xtreg2way: An Algorithm to Estimate the Two-Way Fixed Effects Model 

The `xtreg2way` package is an algorithm to efficiently estimate a two-way fixed effects model.  This algorithm is adapted from the [Matlab code](https://soma.people.stanford.edu/research) written for [this paper](https://soma.people.stanford.edu/sites/g/files/sbiybj4171/f/jem-2014-0008.pdf) by Paulo Somaini and Frank A Wolak.  

## Abstract
We present an algorithm to estimate the two-way fixed effect linear model. The algorithm relies on
the Frisch-Waugh-Lovell theorem and applies to ordinary least squares (OLS), two-stage least squares (TSLS)
and generalized method of moments (GMM) estimators. The coefficients of interest are computed using the
residuals from the projection of all variables on the two sets of fixed effects. Our algorithm has three desirable
features. First, it manages memory and computational resources efficiently which speeds up the computation of the estimates. Second, it allows the researcher to estimate multiple specifications using the same set
of fixed effects at a very low computational cost. Third, the asymptotic variance of the parameters of interest
can be consistently estimated using standard routines on the residualized data. 

## How do I use it?

There are three ways you can interface with the algorithm.

In the most basic way, you can input all variables, including your `y` independent variable, matrix of `X` dependent variables, and `iid`, `tid` as your corresponding group identifiers:

`output <- xtreg2way(y,X,iid,tid,w)`

After running the algorithm once, the output will contain a list `output$struc`.  The list `struc` contains matrices that do not depend on the columns passed.  Thus, you can save time by running the algorithm with `struc` defined on a new set of columns.

`output2 <- xtreg2way(y,newColumnsX, struc=struc)`

Finally, instead of passing your variables explicitly as `y` and `X`, you can pass a formula and a dataframe.

`output3 <- xtreg2way(y~x1+x2, myDataFrame, iid, tid)`

## Whats in this repo?

* R/
    * This folder contains all R code for the `xtreg2way` function, as well as helper functions
* man/
    * This folder contains function documentation.
    * This documentation is automatically built using `roxygen2`. To edit the documentation:
        * Edit the documentation in each function's R/... file
        * Run `devtools::document()` from an R console
        * Save and push this updated documentation
* tests/
     * This folder currently contains a single test, under tests/testthat/test_xtreg2way.R
* vignettes/
     * This folder contains a single vignette, which is an RMarkdown file containing example code.
