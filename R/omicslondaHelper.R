#' Fit longitudinal data
#'
#' Fits longitudinal samples from the same group using negative binomial smoothing splines or LOWESS
#' 
#' @param formula formula to be passed to the regression model
#' @param df dataframe has the Count, Group, Subject, Time
#' @param method fitting method (ssgaussian)
#' @param points points at which the prediction should happen
#' @return returns the fitted model
#' @import gss
#' @import stats
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
curveFitting = function(formula = Count ~ Time, df, method = "ssnbinomial", points){
  
  ## Seprate the two groups
  group.null = df
  group.0 = df[df$Group==0, ]
  group.1 = df[df$Group==1, ]
  
  ## Fitting 
  if(method == "ssgaussian"){
    
    ### Difference between random = mkran(~1|Subject, group.null) and random = ~1|Subject
    ran.null = mkran(~1|Subject, group.null)
    ran.0 = mkran(~1|Subject, group.0)
    ran.1 = mkran(~1|Subject, group.1)
    
    
    ### With Subjet Random Effect
    ## null model
    # mod.null = ssanova(formula, data = group.null, skip.iter=TRUE, random = ran.null, na.action = na.omit)
    # 
    # ## full model
    # mod.0 = ssanova(formula, data = group.0, skip.iter=TRUE, random = ran.0, na.action = na.omit)
    # mod.1 = ssanova(formula, data = group.1, skip.iter=TRUE, random = ran.1, na.action = na.omit)
    
    # ## null model
    mod.null = ssanova(formula, data = group.null, skip.iter=TRUE)
    
    ## full model
    mod.0 = ssanova(formula, data = group.0, skip.iter=TRUE)
    mod.1 = ssanova(formula, data = group.1, skip.iter=TRUE)
    
    # ### Troubleshoot ssanova
    # mod.0.test = ssanova(formula, data = group.0, skip.iter=TRUE, nbasis=10)
    # mod.1.test = ssanova(formula, data = group.1, skip.iter=TRUE, nbasis=10)
    
    ## Calculate goodness of fit F-statistic the fitted model
    rss.null = summary(mod.null)$rss
    rss.full = summary(mod.0)$rss+summary(mod.1)$rss
    f.stat = (rss.null - rss.full)/rss.null
  }
  
  
  ## Estimate values at the provided time points
  est.null = predict(mod.null, data.frame(Time = points), include = c("1", "Time"), se = TRUE)
  est.0 = predict(mod.0, data.frame(Time = points), include = c("1", "Time"), se = TRUE)
  est.1 = predict(mod.1, data.frame(Time = points), include = c("1", "Time"), se = TRUE)
  
  
  
  ## prepare dataframe for plotting
  if (method == "ssgaussian")
  {
    ## Curve dataframe
    dd.null = data.frame(Time = points, Count = est.null$fit, Group = "NULL", Subject = "NULL", SE = est.null$se)
    dd.0 = data.frame(Time = points, Count = est.0$fit, Group = "fit.0", Subject = "fit.0", SE = est.0$se)
    dd.1 = data.frame(Time = points, Count = est.1$fit, Group = "fit.1", Subject = "fit.1", SE = est.1$se)
    
    ## Confidence interval dataframe
    dd.null.u95 = data.frame(Time = points, Count = (est.null$fit + 1.96*est.null$se), Group = "null.u", Subject = "null.u")
    dd.null.l95 = data.frame(Time = points, Count = (est.null$fit - 1.96*est.null$se), Group = "null.l", Subject = "null.l")
    dd.0.u95 = data.frame(Time = points, Count = (est.0$fit + 1.96*est.0$se), Group = "fit.0.u", Subject = "fit.0.u")
    dd.0.l95 = data.frame(Time = points, Count = (est.0$fit - 1.96*est.0$se), Group = "fit.0.l", Subject = "fit.0.l")
    dd.1.u95 = data.frame(Time = points, Count = (est.1$fit + 1.96*est.1$se), Group = "fit.1.u", Subject = "fit.1.u")
    dd.1.l95 = data.frame(Time = points, Count = (est.1$fit - 1.96*est.1$se), Group = "fit.1.l", Subject = "fit.1.l")
  } 
  
  
  
  ## Return the results
  if(method == "ssgaussian")
  {
    output = list(f.stat = f.stat, rss.null = rss.null, rss.full = rss.full, 
                  dd.null = dd.null, dd.0 = dd.0, dd.1 = dd.1, mod.null = mod.null, 
                  mod.0 = mod.0, mod.1 = mod.1, dd.null.u95 = dd.null.u95, dd.null.l95 = dd.null.l95,
                  dd.0.u95 = dd.0.u95, dd.0.l95 = dd.0.l95, dd.1.u95 = dd.1.u95, dd.1.l95= dd.1.l95)
  }
  
  
  return(output)
}


#' Calculate Test-Statistic of each feature's time interval
#'
#' Calculate Test-Statistic of each feature's time interval
#' 
#' @param curve.fit.df gss data object of the fitted spline
#' @return returns testStat for all time intervals
#' @import pracma
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
testStat = function(curve.fit.df){
  size = length(curve.fit.df$dd.null$Time)
  testStat = numeric(size - 1)
  
  
  for(i in 1:(size - 1)){
    testStat.null = trapz(curve.fit.df$dd.null$Time[i:(i+1)], curve.fit.df$dd.null$Count[i:(i+1)])
    testStat.0 = trapz(curve.fit.df$dd.0$Time[i:(i+1)], curve.fit.df$dd.0$Count[i:(i+1)])
    testStat.1 = trapz(curve.fit.df$dd.1$Time[i:(i+1)], curve.fit.df$dd.1$Count[i:(i+1)])
    
    testStat[i] = (testStat.0 - testStat.1) / sqrt((mean(curve.fit.df$dd.0$SE[i:(i+1)]))^2 + (mean(curve.fit.df$dd.1$SE[i:(i+1)]))^2)
    if(is.na(testStat[i])){
      testStat[i] = 0
    }
  }
  
  return(list(testStat = testStat))
}

# 
# 
# #' Find Significant Interval based on testStat
# #'
# #' Find Significant Interval based on testStat
# #' 
# #' @param adjusted.pvalue vector of the adjusted p-value
# #' @param threshold p-value cut off
# #' @param sign vector hold area sign of each time interval 
# #' @return returns a list of the start and end points of all significant time intervals
# #' @references
# #' Ahmed Metwally (ametwall@stanford.edu)
# #' @export
findSigInterval2 = function(adjusted.pvalue, threshold = 0.05, sign)
{
  sig = which(adjusted.pvalue < threshold/2)
  sign = sign[sig]
  padj = adjusted.pvalue[sig]
  start = numeric()
  end = numeric()
  p = numeric()
  dom = numeric()
  
  if(length(sig) == 0)
  {
    cat("No Significant Intevals Found \n")
  }
  else if(length(sig) == 1)
  {
    start = sig[1]
    end = sig [1]
    p = padj[1]
    dom = sign[1]
  }
  else
  {
    start = sig[1]
    
    if((sig[2] - sig[1]) != 1 | sign[2] != sign[1])
    {
      end = c(end, sig[1])
      dom = c(dom, sign[1])
      p = c(p, padj[1])
    }
    
    for(i in 2:length(sig))
    {
      if(i != length(sig))
      {
        if((sig[i] - sig[i-1]) > 1 | sign[i] != sign[i-1])
        {
          start= c(start, sig[i])
        }
        
        if((sig[i+1] - sig[i]) != 1 | sign[i+1] != sign[i])
        {
          end = c(end, sig[i])
          dom = c(dom, sign[i])
          p = c(p, mean(adjusted.pvalue[start[length(start)] : end[length(end)]]))
        }
      }
      else
      {
        if((sig[i]-sig[i-1]) > 1 | sign[i] != sign[i-1])
        {
          start= c(start, sig[i])
        }
        end= c(end, sig[i])
        dom = c(dom, sign[i])
        p = c(p, mean(adjusted.pvalue[start[length(start)] : end[length(end)]]))
      }
    }
  }
  
  return(list(start = start, end = end, pvalue = p, dominant = dom))
}




#' Calculate testStat of each feature's time interval for all permutations
#' 
#' @param perm list has all the permutated models
#' @return returns a list of all permutation area ratio
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
testStatPermutation = function(perm)
{
  testStat.list = list()
  list.len = length(perm)
  for (j in 1:list.len)
  {
    testStat.list[[j]] = testStat(perm[[j]])
  }
  
  return(testStat.list)
}





#' Normalize count matrix 
#'
#' Normalize count matrix
#'
#' @param count count matrix
#' @param method normalization method
#' @return the normalized count matrix
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @export
normalize = function(count, method = "css"){
  # col.data=0 ## this line is for CRAN package
  if(method == "clr")
  {
    cat("Normalization using CLR method \n")
  }
  else if(method == "css")
  {
    cat("Normalization using CSS method \n")
    otu = metagenomeSeq::newMRexperiment(count)
    p.1 = metagenomeSeq::cumNormStatFast(otu, pFlag = TRUE)
    otu.2 = metagenomeSeq::cumNorm(otu, p = p.1)
    count.normalized = metagenomeSeq::MRcounts(otu.2, norm = TRUE)
  }
  else if(method == "tmm")
  {
    cat("Normalization using TMM method \n")
    factors = edgeR::calcNormFactors(count, method="TMM")
    eff.lib.size = colSums(count) * factors
    ref.lib.size = mean(eff.lib.size) #Use the mean of the effective library sizes as a reference library size
    count.normalized = sweep(count, MARGIN = 2, eff.lib.size, "/") * ref.lib.size 
  }
  else if(method == "ra")
  {
    cat("Normalization using Relative Abundance (RA) method \n")
    count.normalized  = apply(count, 2, function(x) (x/sum(x)))
  }
  else if(method == "log10")
  {
    cat("Normalization using log10 of the RA method \n")
    count.normalized  = apply(count, 2, function(x) log10(x/sum(x) + 1)) 
  }
  else if(method == "median_ratio")
  {
    cat("Normalization using Median-Ratio method \n")
    col.data = as.data.frame(cbind(colnames(count), rep(1:2, length.out= ncol(count)), rep(1:2, length.out= ncol(count))))
    rownames(col.data) = col.data[,1]
    col.data =  col.data[,-1]
    colnames(col.data) = c("test","condition")
    data.deseq = DESeq2::DESeqDataSetFromMatrix(countData = count, colData = col.data, ~ condition)
    cds = DESeq2::estimateSizeFactors( data.deseq )
    count.normalized = t( t(DESeq2::counts(cds)) / DESeq2::sizeFactors(cds) )
  }
  return(count.normalized)
}
