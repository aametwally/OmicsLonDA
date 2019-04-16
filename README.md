# OmicsLonDA

OmicsLonDA (Omics Longitudinal Differential Analysis) is a statistical framework that provides robust identification of time intervals where omics features are significantly different between groups. OmicsLonDA is based on 5 main steps: (a) adjust measurements based on each subject's specific baseline, (b) global testing using linear mixed-effect model to select candidate features and covariates for time intervals analysis, (c) fitting smoothing spline regression model, (d) Monte Carlo simulation to generate the empirical distribution of the test statistic, and (e) inference of significant time intervals of omics features. 


<br>

# Getting Started


## Prerequisites

* R(>= 3.6)


## Installation

Download the latest development code of OmicsLonDA from GitHub using devtools
```
library(devtools)
install_github("aametwally/OmicsLonDA", ref = "master")
```




## Example:


### Bugs and Suggestions
OmicsLonDA is under active research development. Please report any bugs/suggestions to Ahmed Metwally (ametwall@stanford.edu).
