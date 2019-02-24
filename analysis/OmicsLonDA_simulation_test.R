# install.packages("devtools", repos = c("http://rstudio.org/_packages", "http://cran.rstudio.com"))

library(ggplot2)
library(gss)
library(plyr)
library(pracma)
library(parallel)
library(doParallel)
library(zoo)
library(pheatmap)
library(viridis)


setwd("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/OmicsLonDA_github/")
#setwd("C:/Users/ametwall/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")


source("R/OmicsLonDA.R")
source("R/CurveFitting.R")
source("R/Visualization.R")
source("R/Permutation.R")
source("R/Normalization.R")
#source("R/OmicsLonDA_Evaluation.R")



### load dataset
#load("data/simulatedDataset_diff_OmicsLonDA_1000.RData", envir = parent.frame(), verbose = FALSE)
load("data/simulatedDataset_diff_OmicsLonDA_normalized_1000.RData", envir = parent.frame(), verbose = FALSE)
head(diff_simulatedDataset_norm[[1]])


## Test OmicsLonDA
points = seq(1, 200, length.out = 200)
output.omicslonda_diff_1 = omicslonda(formula = normalizedCount ~ Time, df = diff_simulatedDataset_norm[[1]], n.perm = 50, 
                                      fit.method = "ssgaussian", points = points,
                                      text = "sim_f1", parall = FALSE, pvalue.threshold = 0.05,
                                      adjust.method = "BH", col = c("blue", "green"),
                                      prefix = "OmicsLonDA_clr_f1-2", ylabel = "CLR-NormalizedCount",
                                      DrawTestStatDist = FALSE, time.unit = "days")
