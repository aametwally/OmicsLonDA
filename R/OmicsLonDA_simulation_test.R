## OmicsLonDA
library(ggplot2)
library(gss)
library(plyr)
library(pracma)
library(parallel)
library(doParallel)
library(zoo)
library(pheatmap)
library(viridis)


setwd("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")
#setwd("C:/Users/ametwall/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")


source("OmicsLonDA/R/OmicsLonDA.R")
source("OmicsLonDA/R/CurveFitting.R")
source("OmicsLonDA/R/Visualization.R")
source("OmicsLonDA/R/Permutation.R")
source("OmicsLonDA/R/Normalization.R")
source("OmicsLonDA/R/OmicsLonDA_Evaluation.R")



### Normalized
load("simulatedDataset_diff_OmicsLonDA_optimized_normalized.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_nondiff_OmicsLonDA_optimized_normalized.RData", envir = parent.frame(), verbose = FALSE)
points = seq(1, 200, length.out = 200)
evaluation_summary_diff = automateEval_OmicsLonDA(data = diff_simulatedDataset_norm, n.perm = 100, points = points,
                                                  pvalue_threshold = 0.05, prefix = "Thu_SU_diff_noSubjRandom_normalized_all")

evaluation_summary_nondiff = automateEval_OmicsLonDA(data = nondiff_simulatedDataset_norm, n.perm = 100, points = points,
                                                     pvalue_threshold = 0.05, prefix = "Thu_SU_nondiff_noSubjRandom_normalized_all")


### 1000 permutations
points = seq(1, 200, length.out = 200)
evaluation_summary_diff = automateEval_OmicsLonDA(data = diff_simulatedDataset_norm, n.perm = 1000, points = points,
                                                  pvalue_threshold = 0.05, prefix = "Thu_SU_diff_noSubjRandom_normalized_all_1000p")

evaluation_summary_nondiff = automateEval_OmicsLonDA(data = nondiff_simulatedDataset_norm, n.perm = 1000, points = points,
                                                     pvalue_threshold = 0.05, prefix = "Thu_SU_nondiff_noSubjRandom_normalized_all_1000p")






# ### Workout the evaluation dimension
points = seq(1, 200, length.out = 200)
evaluation_summary_nondiff_f1_f2 = automateEval_OmicsLonDA(data = diff_simulatedDataset_norm[1:2], n.perm = 100, points = points,
                                                     pvalue_threshold = 0.05, prefix = "Thu_SU_diff_noSubjRandom_normalized_1_2")





## Test omicslonda alone
points = seq(1, 200, length.out = 200)
output.omicslonda_diff_1 = omicslonda(formula = Count ~ Time, df = diff_simulatedDataset_norm[[1]], n.perm = 50, fit.method = "ssgaussian", points = points,
                                      text = "simdf_1", parall = FALSE, pvalue.threshold = 0.05,
                                      adjust.method = "BH", col = c("blue", "green"),
                                      prefix = "TestOmicsLonDA_ssgaussian_norm_stepbystep_f1", ylabel = "NormalizedCounts",
                                      DrawTestStatDist = FALSE, time.unit = "days")
