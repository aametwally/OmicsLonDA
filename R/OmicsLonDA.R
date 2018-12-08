#' Omics Longitudinal Differential Abundance Analysis for one feature
#'
#' Find significant time intervals of the one feature
#' 
#' @param Count matrix has the number of reads that mapped to each feature in each sample.
#' @param Time vector of the time label of each sample.
#' @param Group vector of the group label of each sample.
#' @param Subject vector of the subject Subject label of each sample.
#' @param n.perm number of permutations.
#' @param fit.method fitting method (ssnbinomial, ssguassian, lowess).
#' @param points points at which the prediction should happen.
#' @param text Feature's name.
#' @param parall boolean to indicate whether to use multicore.
#' @param pvalue.threshold p-value threshold cutoff for identifing significant time intervals.
#' @param adjust.method multiple testing correction method.
#' @param time.unit time unit used in the Time vector (hours, days, weeks, months, etc.)
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures (default: "Normalized Count")
#' @return returns a list of the significant time intervals for the tested feature.
#' @import parallel
#' @import doParallel
#' @import stats
#' @import zoo 
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' data(omicslonda_test_data)
#' n.sample = 5
#' n.timepoints = 10
#' n.group = 2
#' Group = factor(c(rep(0, n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' Subject = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' \dontrun{
#' output.nbinomial = omicslonda(Count = omicslonda_test_data[1,], Time = Time, Group = Group,
#' Subject = Subject, fit.method =  "ssnbinomial", n.perm = 10, points = points,
#' text = rownames(omicslonda_test_data)[1], parall = FALSE, pvalue.threshold = 0.05, 
#' adjust.method = "BH", time.unit = "hours", ylabel = "Normalized Count", col = c("black", "green"))
#' }
#' @export
omicslonda = function(formula = Count ~ Time, df, n.perm = 500, fit.method = "ssnbinomial", 
                      points, text = 0, parall = FALSE, pvalue.threshold = 0.05, 
                      adjust.method = "BH", time.unit = "days", ylabel = "Normalized Count", col = c("blue", "firebrick"),
                      prefix = "Test")
{
  cat("Start OmicsLonDA \n")
  
  if (!dir.exists(prefix)){
    dir.create(file.path(prefix))
  }
  
  # Extract groups
  Group = as.character(df$Group)
  group.levels = sort(unique(Group))
  
  
  
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
  
  
  
  ## Preprocessing: add pseudo counts for zero abudnance features
  df$Count = df$Count + 1e-8
  
  
  ## Form OmicsLonDA dataframe
  #aggregate.df = data.frame(Count = Count, Time = Time, Group = Group, Subject = Subject)
  
  
  
  print(df)
  
  ## Visualize feature's abundance accross different time points  
  visualizeFeature(df, text, group.levels, unit = time.unit, ylabel = ylabel, col = col, prefix = prefix)
  
  
  group.0 = df[df$Group == 0, ]
  group.1 = df[df$Group == 1, ]
  points.min = max(sort(group.0$Time)[1], sort(group.1$Time)[1])
  points.max = min(sort(group.0$Time)[length(group.0$Time)], sort(group.1$Time)[length(group.1$Time)])
  points = points[which(points >= points.min & points <= points.max)]
  
  
  cat("points.min = ", points.min, "\n")
  cat("points.max = ", points.max, "\n")
  cat("points = ", points, "\n")
  
  cat("Start Curve Fitting \n") 
  
  if (fit.method == "ssnbinomial")
  {
    cat("Fitting: Smoothing Spline Negative-Binomial Regression  \n")
    model = tryCatch({
      curveFitting(formula = formula, df, method= "ssnbinomial", points)
    },  error = function(err) {
      print(paste("ERROR in gss = ", err, sep="")); 
      return("ERROR")
    })
  } else if (fit.method == "ssgaussian")
  {
    cat("Fitting: Smoothing Spline Gaussian Regression \n")
    model = tryCatch({
      curveFitting(formula = formula, df, method= "ssgaussian", points)
    },  error = function(err) {
      print(paste("ERROR in gss = ", err, sep="")); 
      return("ERROR")
    })
  } else if (fit.method == "lowess")
  {
    cat("Fitting: LOWESS \n")
    model = curveFitting(formula = formula, df = df, method= "lowess", points)
  } else {
    cat("You have entered unsupported fitting method\n")
    quit()
  }
  
  ## Visualize feature's trajectories spline
  visualizeFeatureSpline(df, model, fit.method, text, group.levels, unit = time.unit, ylabel = ylabel, 
                         col = col, prefix = prefix)
  
  
  
  
  ## Calculate area under the fitted curve for each time interval
  # cat("Calculate Area Under the Fitted Curves \n")
  # area = intervalArea(model)
  # modelAR.abs = area$ar.abs
  
  
  
  ### DEVELOPMENT: TEST the new test Statistic
  stat = testStat(model)$testStat
  
  
  
  
  ## TODO: remove this warning
  # 1: <anonymous>: ... may be used in an incorrect context: ‘.fun(piece, ...)’
  
  ## Permutation 
  perm  = permutation(formula = formula, df, n.perm, fit.method, points, parall = parall)
  
  
  # ## Area p-value per unit interval
  # area.perm = areaPermutation(perm)
  # a1 = do.call(rbind, area.perm)
  # a2 = do.call(rbind, a1[,2])
  # 
  # ##  Visualize AR empirical distribution
  # visualizeARHistogram(a2, text, fit.method, prefix = prefix, modelAR = modelAR.abs)
  # 
  # ## Calculate AR p-value 
  # pvalue.area = sapply(1:(length(points)-1), function(i){
  #   sum(a2[,i] >= area$ar.abs[i])/length(a2[,i])
  # } )
  # 
  # 
  # ## Identify significant time inetrval based on the adjusted p-value 
  # cat("p-value Adjustment Method = ", adjust.method, "\n")
  # adjusted.pvalue = p.adjust(pvalue.area, method = adjust.method)
  
  
  ### DEVELOPEMENT: Test for having the null distribution from all timepoints
  test.stat.prem = testStatPermutation(perm)
  t1 = do.call(rbind, test.stat.prem)
  t2 = do.call(rbind, t1[,1])
  t3 = as.vector(t2)
  length(t3)
  pvalue.test.stat = sapply(1:(length(points)-1), function(i){
    if(stat[i]>=0)
    {
      sum(t3 > stat[i])/length(t3)
    }
    else if(stat[i]<0)
    {
      sum(t3 < stat[i])/length(t3)
    }
  } )
  adjusted.pvalue2 = p.adjust(pvalue.test.stat, method = adjust.method)
  adjusted.pvalue = adjusted.pvalue2
  ##  Visualize AR empirical distribution for the null distribution
  visualizeARHistogram2(t3, text, fit.method, prefix = prefix, modelStat = stat)
  
  
  
  #interval = findSigInterval(adjusted.pvalue, threshold = pvalue.threshold, sign = area$ar.sign)
  interval = findSigInterval2(adjusted.pvalue, threshold = pvalue.threshold, sign = sign(stat))
  
  st = points[interval$start]
  en = points[interval$end + 1]
  
  
  if(length(st) > 0)
  {
    ## Visualize sigificant area
    visualizeArea(aggregate.df, model, fit.method, st, en, text, group.levels, unit = time.unit, ylabel = ylabel,
                  col = col, prefix = prefix)
  }
  
  
  ## Calculate start, end, dominant for each interval
  interval.start = points[-length(points)]
  interval.end = points[-1]
  dominant = sign(stat)
  dominant[which(dominant == 1)] = gr.1
  dominant[which(dominant == -1)] = gr.2
  
  
  ## Prepare Log2FoldChange
  avg.mod0.count = rollapply(model$dd.0$Count, 2, mean)
  avg.mod1.count = rollapply(model$dd.1$Count, 2, mean)
  foldChange = avg.mod0.count/avg.mod1.count
  log2FoldChange = log2(foldChange)
  
  
  output.details = list(feature = rep(text, length(interval.start)), 
                        interval.start = interval.start, interval.end = interval.end,
                        avg.mod0.count = avg.mod0.count, avg.mod1.count = avg.mod1.count, 
                        foldChange = foldChange, log2FoldChange = log2FoldChange, 
                        testStat = stat, testStat.abs = abs(stat), testStat.sign = sign(stat), dominant = dominant,
                        intervals.pvalue = pvalue.test.stat, adjusted.pvalue = adjusted.pvalue)
  output.summary = data.frame(feature = rep(text, length(interval$start)), start = st, end = en,
                              dominant = interval$dominant, pvalue = interval$pvalue)
  
  
  ## Output table and volcano plot that summarize time intervals statistics
  feature.summary = as.data.frame(do.call(cbind, output.details), stringsAsFactors = FALSE)
  write.csv(feature.summary, file = sprintf("%s/Feature_%s_Summary_%s.csv", prefix, text, fit.method), row.names = FALSE)
  x = as.data.frame(sapply(feature.summary[, c("foldChange","log2FoldChange", "intervals.pvalue", "adjusted.pvalue")], as.numeric))
  visualizeVolcanoPlot(df = x, text, prefix = prefix, fit.method = fit.method)
  cat("\n\n")
  
  return(list(detailed = output.details, summary = output.summary))
}


#' Omics Longitudinal Differential Abundance Analysis for all Features
#'
#' Identify significant features and their significant time interval
#' 
#' @param Count Count matrix of all features
#' @param Time Time label of all samples
#' @param Group Group label of all samples
#' @param Subject individual Subject label for samples
#' @param n.perm number of permutations
#' @param fit.method The fitting method (ssnbinomial, ssgaussian, lowess)
#' @param num.intervals The number of time intervals at which omicslonda test differential abundance 
#' @param parall logic to indicate whether to use multicore
#' @param pvalue.threshold p-value threshold cutoff
#' @param adjust.method Multiple testing correction methods
#' @param time.unit time unit used in the Time vector (hours, days, weeks, months, etc.)
#' @param norm.method normalization method to be used to normalize count matrix (css, tmm, ra, log10, median_ratio) 
#' @param prefix prefix for the output figure
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures (default: "Normalized Count")
#' @return Returns a list of the significant features a long with their significant time intervals
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' \dontrun{
#' data(omicslonda_test_data)
#' n.sample = 5
#' n.timepoints = 10
#' n.group = 2
#' Group = factor(c(rep(0, n.sample*n.timepoints), rep(1,n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' Subject = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' output.nbinomial = omicslondaAll(Count = omicslonda_test_data, Time = Time, Group = Group,
#' Subject = Subject, n.perm = 10, fit.method =  "ssnbinomial", num.intervals = 100, 
#' parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", 
#' time.unit = "hours", norm.method = "none", prefix = "Test",  time.unit = "hours", 
#' ylabel = "Normalized Count", col = c("black", "green"))
#' }
#' @export
omicslondaAll = function(formula = Count ~ Time, countTable, metadata, n.perm = 500, fit.method = "ssnbinomial", 
                        num.intervals = 100, parall = FALSE, pvalue.threshold = 0.05, 
                        adjust.method = "BH", time.unit = "days", norm.method = "none", prefix = "TestOmicsLonDA",
                        ylabel = "Normalized Count", col = c("blue", "firebrick"))
{
  ## Check the dimentions of the annotation vectors and count matrix
  if(ncol(countTable) == nrow(metadata))
  {
      cat("Dimensionality check passed\n")
  } else
  {
    stop("Number of samples in count table does not match the number of samples in metadata table.")
  }
  
  
  #### Create Prefix folder
  dir.create(file.path(prefix), showWarnings = FALSE)
  
  
  ## Normalization
  if(norm.method != "none")
  {
    if(norm.method == "css" | norm.method == "tmm" | norm.method == "ra" | norm.method == "log10" | norm.method == "median_ratio")
    {
      cat("Normalizaton method = ", norm.method, "\n")
      countTable = normalize(countTable, method = norm.method)
    } else{
      stop("You have entered a wrong normalization method")
    }
  }
  
  ## Specify the test/prediction timepoints for omicslonda
  if(num.intervals == "none")
    points = seq(min(metadata$Time), max(metadata$Time))
  else
    points = seq(min(metadata$Time), max(metadata$Time), length.out = num.intervals + 1)
  
  cat("Prediction Points = ")
  print(points)
  cat("\n")
  
  ## TODO: Filter out the taxa that always have zero of one/both group
  metadata$Group = as.character(metadata$Group)
  group.levels = sort(unique(metadata$Group))
  if(length(group.levels) > 2){
    stop("You have more than two phenotypes.")
  }
  gr.1 = group.levels[1]
  gr.2 = group.levels[2]
  
  data.count.filt = as.matrix(countTable)
  
  ## Apply omicslonda for each feature
  n.features = nrow(data.count.filt)
  detailed = list()
  summary = list()
  for (i in 1:n.features)
  {
    cat ("Feature  = ", rownames(data.count.filt)[i], "\n")
    x = cbind(Count = data.count.filt[i,], metadata)
    x = na.omit(x)
    print(x)
    out = omicslonda(formula = formula, df = x, n.perm = n.perm, fit.method = fit.method, points = points,
               text = rownames(data.count.filt)[i], parall = parall, pvalue.threshold = pvalue.threshold,     
               adjust.method = adjust.method, time.unit = time.unit, ylabel = ylabel, col = col, prefix = prefix)
    detailed[[i]] = out$detailed
    summary[[i]] = out$summary
  }
  
  summary.tmp = do.call(rbind, summary)
  summary.tmp$dominant[which(summary.tmp$dominant == 1)] = gr.1
  summary.tmp$dominant[which(summary.tmp$dominant == -1)] = gr.2
  
  ## Output table and figure that summarize the significant time intervals
  write.csv(summary.tmp, file = sprintf("%s/OmicsLonDA_TimeIntervals_%s_%s.csv", prefix, fit.method, prefix), row.names = FALSE)
  visualizeTimeIntervals(interval.details = summary.tmp, prefix, unit = time.unit, col = col, fit.method = fit.method)
  
  
  aggregateData = list(output.detail = detailed, output.summary = summary.tmp)
  save(aggregateData, file = sprintf("%s/OmicsLonDA_Summary_%s_%s.RData", prefix, fit.method, prefix))
  return(aggregateData)
}










omicslondaCorr = function (data1, data2,
                  num.intervals = 100, parall = FALSE, pvalue.threshold = 0.05, 
                  time.unit = "days", prefix = "TestOmicsLonDA_Correlation",
                 col = c("blue", "firebrick"))
{
  cat("Start omicslondaCorr \n")
  
  ## TODO: Export in MetaLonDA number of rows and features
  ## TODO: Separate the positive and negative association, and mixed association
  ## TODO: HOW to deal with unequal length of intervals? This might happen when measurements of one variable is not taken at the same time
  ### TODO: Differentiate between each dominant (1,-1)
  ## TODO: limit the X, Y, Q matrices to the significant features
  ## TODO: Think of visualization technique to visualize time series
  ## TODO: Check from the timepoints for each file as they may not be the same due to filteration step of omicslonda
  
  
  ## Put 1 in every cell that has significant
  nfeatures1 = length(data1$output.detail)
  ninterval1 = length(data1$output.detail[[1]]$adjusted.pvalue)
  X = matrix(0, nrow = nfeatures1, ncol = ninterval1)
  for(i in 1:nfeatures1)
  {
    X[i, which(data1$output.detail[[i]]$adjusted.pvalue <= 0.05)] = 1
    
  }
  
  nfeatures2 = length(data2$output.detail)
  ninterval2 = length(data2$output.detail[[1]]$adjusted.pvalue)
  Y = matrix(0, nrow = nfeatures2, ncol = ninterval2)
  for(i in 1:nfeatures2)
  {
    Y[i, which(data2$output.detail[[i]]$adjusted.pvalue <= 0.05)] = 1
    
  }
  
  cat("v0.1 \n")


  ### Design Matrix 
  Z = array(0, c(nfeatures1, nfeatures2, ninterval1))  
  
  
  ## TODO: How to do this in linear algebra form
  for(i in 1:nfeatures1)
  {
    for(j in 1:nfeatures2)
    {
      Z[i,j,]=X[i,]*Y[j,]
    }
  }
  
  cat("v0.2 \n")
  ## Create 2d matrix between data1xdata2 where each element represent the number of intervals that are significant between both
  Q = matrix(0, nrow = nfeatures1, ncol = nfeatures2)  
  for(i in 1:nfeatures1)
  {
    for(j in 1:nfeatures2)
    {
      Q[i,j]= sum (Z[i,j,])
    }
  }
  cat("v0.3 \n")
  
  
  ### Get names of Matrices X, Y Features
  Xfeatures = c()
  for (i in 1:nfeatures1)
  {
    tmp = unique(data1$output.detail[[i]]$feature)
    Xfeatures = c(Xfeatures, tmp)
  }
  
  Yfeatures = c()
  for (i in 1:nfeatures2)
  {
    tmp = unique(data2$output.detail[[i]]$feature)
    Yfeatures = c(Yfeatures, tmp)
  }
  
  Q = as.data.frame(Q)
  rownames(Q) = Xfeatures
  colnames(Q) = Yfeatures
  cat("v0.4 \n")
  
  
  ### Visualization
  # Remove the insignificant rows and columns
  Q = Q[-which(rowSums(Q)==0), -which(colSums(Q)==0)]
  
  
  cat("v0.5 \n")
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(as.matrix(xs), probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  
  custom_breaks <- quantile_breaks(Q, n = 10)
  
  cat("v0.6 \n")
  xx = paste(prefix, "_Heatmap.jpg", sep = "")
  jpeg(filename = xx, res = 1200, height = 20, width = 20, units = 'cm')
  pheatmap(
    mat               = Q,
    color             = inferno(15),
    #breaks            = custom_breaks,
    border_color      = "black",
    # cluster_cols      = daily_cluster_cols,
    # cluster_rows      = daily_cluster_rows,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    # annotation_row    = daily_col,
    # annotation_colors = daily_colors,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = prefix
  )
  dev.off()
  cat("v0.7 \n")
}
