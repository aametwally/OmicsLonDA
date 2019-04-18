context("OmicsLonDA Test")
library(OmicsLonDA)

test_that("Element test", {
  set.seed(34443)
  data("omicslonda_data_example")
  omicslonda_se_object_adjusted = adjustBaseline(se_object = omicslonda_data_example$omicslonda_se_object)
  omicslonda_test_object = omicslonda_se_object_adjusted[1,]
  points = seq(1, 500, length.out = 500)
  res = omicslonda(se_object = omicslonda_test_object, n.perm = 10,
                   fit.method = "ssgaussian", points = points, text = "Feature_1",
                   parall = FALSE, pvalue.threshold = 0.05, 
                   adjust.method = "BH", time.unit = "days",
                   ylabel = "Normalized Count",
                   col = c("blue", "firebrick"), prefix = "OmicsLonDA_example")
  #correct <- c("sim_f1")
  #expect_match(output.omicslonda_diff_1[[1]][[1]][1], correct)
  expect_true(identical(res$details$feature[1], "Feature_1"))
})
