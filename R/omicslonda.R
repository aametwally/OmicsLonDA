#' Omics Longitudinal Differential Analysis for one feature
#'
#' Find significant time intervals of omic feature
#' 
#' @param formula formula to be passed to the regression model
#' @param df dataframe that contains (subject, time, group, count, normalizedCount, or any any other covariates) for each feature
#' @param n.perm number of permutations.
#' @param fit.method fitting method (ssguassian).
#' @param points points at which the prediction should happen.
#' @param text Feature's name.
#' @param parall boolean to indicate whether to use multicore.
#' @param pvalue.threshold p-value threshold cutoff for identifing significant time intervals.
#' @param adjust.method multiple testing correction method.
#' @param time.unit time unit used in the Time vector (hours, days, weeks, months, etc.)
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures (default: "Normalized Count")
#' @param prefix prefix to be used to create directory for the analysis results
#' @param DrawTestStatDist boolean to indicate if the histogram of the testStat needs to be plotted or not. Default is "FALSE" since the histogram is usually a big file
#' @return returns a list of the significant time intervals for the tested feature.
#' @import parallel
#' @import doParallel
#' @import stats
#' @import zoo 
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
omicslonda = function(formula = Count ~ Time, df, n.perm = 500, fit.method = "ssnbinomial", 
                      points, text = 0, parall = FALSE, pvalue.threshold = 0.05, 
                      adjust.method = "BH", time.unit = "days", ylabel = "Normalized Count", col = c("blue", "firebrick"),
                      prefix = "Test", DrawTestStatDist = FALSE)
{
  cat("Start OmicsLonDA \n")
  # all.vars(formula)[1]
  
  if (!dir.exists(prefix)){
    dir.create(file.path(prefix))
  }
  
  # Extract groups
  Group = as.character(df$Group)
  group.levels = sort(unique(Group))
  #print("Group levels = \n")
  print(group.levels)
  #print("\n")
  
  
  
  if(length(group.levels) > 2){
    stop("You have more than two phenotypes.")
  } else if(length(group.levels) < 2){
    stop("You have less than two phenotypes.")
  }
  
  
  gr.1 = as.character(group.levels[1])
  gr.2 = as.character(group.levels[2])

  
  
  levels(df$Group) = c(levels(df$Group), "0", "1")
  df$Group[which(df$Group == gr.1)] = 0
  df$Group[which(df$Group == gr.2)] = 1
  
  
  
  ## Preprocessing: add pseudo counts
  ## TODO: Make sure this is applied to the corrsponding column (Count, rawCount, normalized Count)
  #df$Count = df$Count + 1e-8
  
  

  
  ## Visualize feature's abundance accross different time points  
  visualizeFeature(formula = formula, df, text, group.levels, unit = time.unit, ylabel = ylabel, col = col, prefix = prefix)
  
  
  group.0 = df[df$Group == 0, ]
  group.1 = df[df$Group == 1, ]
  points.min = max(sort(group.0$Time)[1], sort(group.1$Time)[1])
  points.max = min(sort(group.0$Time)[length(group.0$Time)], sort(group.1$Time)[length(group.1$Time)])
  points = points[which(points >= points.min & points <= points.max)]
  
  
  cat("points.min = ", points.min, "\n")
  cat("points.max = ", points.max, "\n")
  cat("points = ", points, "\n")
  
  cat("Start Curve Fitting \n") 
  
  if (fit.method == "ssgaussian")
  {
    cat("Fitting: Smoothing Spline Gaussian Regression \n")
    model = tryCatch({
      curveFitting(formula = formula, df, method= "ssgaussian", points)
    },  error = function(err) {
      print(paste("ERROR in gss = ", err, sep="")); 
      return("ERROR")
    })
  } else {
    cat("You have entered unsupported fitting method\n")
    quit()
  }
  
  ## Visualize feature's trajectories spline
  visualizeFeatureSpline2(formula = formula, df, model, fit.method, text, group.levels, unit = time.unit, ylabel = ylabel, 
                         col = col, prefix = prefix)

  
  ### Test Statistic
  stat = testStat(model)$testStat
  
  
  
  ### TODO: Test Bootstrapping
  # library(boot)
  # #bt = boot(formula = formula, bs.df = df, method = fit.method, points = points, parall = parall, prefix = prefix)
  # set.seed(27262)
  # #index = sample(1:726, 726)
  # #data = df
  # bt = boot(df, bootstrapOmicslonda, 10)
  # 
  # set.seed(27262)
  # bootstrapOmicslonda(df,1:726)
  
  
  
  
  
  ## TODO: fix the occasional warning  # 1: <anonymous>: ... may be used in an incorrect context: ‘.fun(piece, ...)’
  ## Permutation 
  perm  = permutationMC2(formula = formula, perm.dat = df, n.perm = n.perm, method = fit.method, points = points, parall = parall, prefix = prefix)

  test.stat.prem = testStatPermutation(perm)
  t1 = do.call(rbind, test.stat.prem)
  

  t2 = unlist(t1[,1])
  t3 = as.vector(t2)
  length(t3)
  pvalue.test.stat = vapply(1:(length(points)-1), function(i){
    if(stat[i]>=0)
    {
      sum(t3 > stat[i])/length(t3)
    }
    else if(stat[i]<0)
    {
      sum(t3 < stat[i])/length(t3)
    }
  } )

  ## Adjust p-values
  if(adjust.method == "qvalue"){
    #adjusted.pvalue = qvalue(pvalue.test.stat, pfdr = TRUE)
    adjusted.pvalue = .data$qvalue_truncp(pvalue.test.stat)$qvalues
  }else{
    adjusted.pvalue = p.adjust(pvalue.test.stat, method = adjust.method)
  }
    
  
  
  if(DrawTestStatDist)
  {
    ##  Visualize testStat empirical distribution for the null distribution
    visualizeTestStatHistogram(t3, text, fit.method, prefix = prefix, modelStat = stat)
  }


  interval = findSigInterval2(adjusted.pvalue, threshold = pvalue.threshold, sign = sign(stat))
  
  st = points[interval$start]
  en = points[interval$end + 1]
  
  
  if(length(st) > 0)
  {
    ## Visualize sigificant area
    visualizeArea(formula = formula, model, fit.method, st, en, text, group.levels, unit = time.unit, ylabel = ylabel,
                  col = col, prefix = prefix)
  }
  
  
  ## Calculate start, end, dominant for each interval
  interval.start = points[-length(points)]
  interval.end = points[-1]
  dominant = sign(stat)
  dominant[which(dominant == 1)] = gr.1
  dominant[which(dominant == -1)] = gr.2
  
  
  ## Calculate FoldChange
  avg.mod0.count = rollapply(model$dd.0$Count, 2, mean)
  avg.mod1.count = rollapply(model$dd.1$Count, 2, mean)
  foldChange = avg.mod0.count/avg.mod1.count
  
  
  output.details = list(feature = rep(text, length(interval.start)), 
                        interval.start = interval.start, interval.end = interval.end,
                        avg.mod0.count = avg.mod0.count, avg.mod1.count = avg.mod1.count, 
                        foldChange = foldChange,
                        testStat = stat, testStat.abs = abs(stat), testStat.sign = sign(stat), dominant = dominant,
                        intervals.pvalue = pvalue.test.stat, adjusted.pvalue = adjusted.pvalue, points = points)
  output.summary = data.frame(feature = rep(text, length(interval$start)), start = st, end = en,
                              dominant = interval$dominant, pvalue = interval$pvalue)
  
  
  ## Output table that summarize time intervals statistics
  feature.summary = as.data.frame(do.call(cbind, output.details), stringsAsFactors = FALSE)
  write.csv(feature.summary, file = sprintf("%s/Feature_%s_Summary_%s.csv", prefix, text, fit.method), row.names = FALSE)
  cat("\n\n")
  
  return(list(detailed = output.details, summary = output.summary))
}
