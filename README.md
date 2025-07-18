
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DualProbitCDF Repository

<!--**DUA**l probit marginal **L**ikelihood **E**stimation via **EP** **A**pproximation, -->

[![Last
Commit](https://img.shields.io/github/last-commit/fradenti/DualProbitCDF)](https://github.com/fradenti/DualProbitCDF)

The repository **DualProbitCDF** contains the R package and code
implementing the estimation strategy proposed in the article:

> *“Multivariate Gaussian cumulative distribution functions as the
> marginal likelihood of their dual Bayesian probit model”* by A. Fasano
> and F. Denti, accepted for publication in **Biometrika (2025, In
> press)**.

In addition to the package itself, the repository includes all the
scripts needed to reproduce the simulation studies and the figures
presented in the paper.

## Getting Started

1.  **Open the R Project**  
    Open the R project file `DualProbitCDF.Rproj` and set it as your
    working directory.

2.  **Install the R Package**  
    Install the included R package `EPmvnCDF`, which stands for  
    *Expectation Propagation for the estimation of Multivariate Normal
    Cumulative Distribution Functions*, via

    ``` r
    install.packages("EPmvnCDF_0.2.0.tar.gz", repos = NULL, type = "source")
    ```

    Otherwise open the `R` project contained in the `EPmvnCDF/` folder
    and build the package from there.

Then you can run the code.

In particular:

- Scripts starting with the letter **`A`** contain utility functions
  that are reused throughout the project.

- Scripts starting with the letter **`B`** contain the code to run the
  four estimation methods described in the main paper:

  - our **Expectation Propagation (EP)** proposal, implemented using
    both Cholesky and Eigenvalue decompositions;
  - the method in the **`TruncatedNormal`** package;
  - the methods in the **`tlrmvnmvt`** package;
  - the **Orthant method** by Ridgway (2016), implemented in **C++**
    (see the `cpp_source/` directory).

  The output of each method is stored in the corresponding subfolder
  within the `RDS/` directory.

- The `.Rmd` file starting with the letter **`C`** extracts and
  processes the simulation results stored in the `.RDS` files. It
  creates summary data frames (saved in `RDS/`) and saves the figures
  used in the paper (in `NewFigures/`).

- Scripts starting with the letter **`D`** reproduce an additional
  simulation study included in the Supplementary Material.
