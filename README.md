# OmicsLonDA

OmicsLonDA (Omics Longitudinal Differential Analysis) is a statistical framework that provides robust identification of time intervals where omics features are significantly different between groups. OmicsLonDA is based on 5 main steps: (a) adjust measurements based on each subject's specific baseline, (b) global testing using linear mixed-effect model to select candidate features and covariates for time intervals analysis, (c) fitting smoothing spline regression model, (d) Monte Carlo simulation to generate the empirical distribution of the test statistic, and (e) inference of significant time intervals of omics features. 


<br>

# Getting Started


## Prerequisites

* R(>= 3.3.2)

<br>

## Installation


#### Install OmicsLonDA from GitHub:
Download the latest development code of OmicsLonDA from GitHub using devtools
```
library(devtools)
install_github("aametwally/OmicsLonDA", ref = "master")
```


<br>



## Example:
```
library(OmicsLonDA)

## Load 1000 simulated features
load("data/simulatedDataset_diff_OmicsLonDA_normalized_1000.RData", envir = parent.frame(), verbose = FALSE)
head(diff_simulatedDataset_norm[[1]])

```


```
## Define the prediction timepoints 
points = seq(1, 10, length.out = 100)
```

### Apply OmicsLonDA on feature #1 

```
output.omicslonda_diff_1 = omicslonda(formula = normalizedCount ~ Time, df = diff_simulatedDataset_norm[[1]], n.perm = 1000, 
                                      fit.method = "ssgaussian", points = points,
                                      text = "sim_f1", parall = FALSE, pvalue.threshold = 0.05,
                                      adjust.method = "BH", col = c("blue", "green"),
                                      prefix = "OmicsLonDA_clr_f1-2", ylabel = "CLR-NormalizedCount",
                                      DrawTestStatDist = FALSE, time.unit = "days")
```


<br>

### Bugs and Suggestions
OmicsLonDA is under active research development. Please report any bugs/suggestions to Ahmed Metwally (ametwall@stanford.edu).
