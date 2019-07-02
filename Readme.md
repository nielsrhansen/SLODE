## Statistical Learning of Ordinary Differential Equation Systems

This repository contains the code for the simulations from the paper

*Statistical Learning of Ordinary Differential Equation Systems*

by Frederik Vissing Mikkelsen and Niels Richard Hansen. 

The main package developed for the paper is [`episode`](https://cran.r-project.org/package=episode).

## Installation

The code contained in this respository is separated into two parts. One part
consists of the R package `tsars`, that implements an alternative to the 
main learning algorithms in `episode`, and the other part consists of a 
collection of R scripts that run all simulations. 

Install the `tsars` package by running the command below. This will 
also install all packages, including `episode`, that the scripts depend 
upon. 

``` r
# install.packages("devtools") # If `devtools` is not installed
devtools::install_github("nielsrhansen/SLODE/tsars", dependencies = TRUE)
```

## Usage

To run the R scripts, first install `tsars` as descriped above, then 
clone the repository and follow the instructions 
from the [Readme](SimStudies/README.md) file in the SimStudies directory. 

