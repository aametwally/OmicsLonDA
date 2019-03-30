context("OmicsLonDA Test")
library(OmicsLonDA)

test_that("Element test", {
  data(diff_simulatedDataset_norm)
  head(diff_simulatedDataset_norm[[1]])
  points = seq(1, 200, length.out = 200)
  output.omicslonda_diff_1 = omicslonda(formula = normalizedCount ~ Time, df = diff_simulatedDataset_norm[[1]], n.perm = 10, 
                                        fit.method = "ssgaussian", points = points,
                                        text = "sim_f1", parall = FALSE, pvalue.threshold = 0.05,
                                        adjust.method = "BH", col = c("blue", "green"),
                                        prefix = "OmicsLonDA_clr_f1", ylabel = "CLR-NormalizedCount",
                                        DrawTestStatDist = FALSE, time.unit = "days")
  #correct <- c("sim_f1")
  #expect_match(output.omicslonda_diff_1[[1]][[1]][1], correct)
  expect_true(identical(output.omicslonda_diff_1[[1]][[1]][1], "sim_f1"))
})
