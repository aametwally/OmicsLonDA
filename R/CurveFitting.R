#' Fit longitudinal data
#'
#' Fits longitudinal samples from the same group using negative binomial smoothing splines or LOWESS
#' 
#' @param df dataframe has the Count, Group, Subject, Time
#' @param method fitting method (ssnbinomial, ssgaussian, lowess)
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
  ## TODO Check the projection
  if(method == "ssnbinomial"){
    
    ran.null = mkran(~1|Subject, group.null)
    ran.0 = mkran(~1|Subject, group.0)
    ran.1 = mkran(~1|Subject, group.1)
    
    ## null model
    mod.null = gssanova(formula, data = group.null, family = "nbinomial", skip.iter=TRUE, random = ran.null, na.action = na.omit)
    mod.null.nbinomial.project = project(mod.null, c("Time"))
    
    ## full model
    mod.0 = gssanova(formula, data = group.0, family = "nbinomial", skip.iter=TRUE, random = ran.0, na.action = na.omit)
    mod.1 = gssanova(formula, data = group.1, family = "nbinomial", skip.iter=TRUE, random = ran.1, na.action = na.omit)
    mod.0.nbinomial.project = project(mod.0, c("Time"))
    mod.1.nbinomial.project = project(mod.1, c("Time"))
  }
  else if(method == "ssgaussian"){
    
    ### Different between random = mkran(~1|Subject, group.null) and random = ~1|Subject
    
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
    
    
    
    # ##### TEST
    # ## null model
    # mod.null = ssanova(formula, data = group.null, skip.iter=TRUE)
    # 
    # ## full model
    # mod.0 = ssanova(formula, data = group.0, skip.iter=TRUE)
    # mod.1 = ssanova(formula, data = group.1, skip.iter=TRUE)
    
    ## Calculate goodness of fit F-statistic the fitted model
    rss.null = summary(mod.null)$rss
    rss.full = summary(mod.0)$rss+summary(mod.1)$rss
    f.stat = (rss.null - rss.full)/rss.null
    
  }
  else if(method == "lowess"){
    ## TODO: LOWESS need regularization
    ## TODO: LOWESS can only take up to four covariates
    ## loess is diffefrent from lowess that it considers more covariates versus only one for the lowess
    
    ## null model
    mod.null = loess(formula, data = group.null, model = TRUE)
    
    ## full model
    mod.0 = loess(formula, data = group.0, model = TRUE)
    mod.1 = loess(formula, data = group.1, model = TRUE)
  }
  
  
  
  
  
  ## TODO: make sure this inc works for lowess
  ## Estimate values at the provided time points
  #cat("in prediction \n")
  est.null = predict(mod.null, data.frame(Time = points), include = c("1", "Time"), se = TRUE)
  est.0 = predict(mod.0, data.frame(Time = points), include = c("1", "Time"), se = TRUE)
  est.1 = predict(mod.1, data.frame(Time = points), include = c("1", "Time"), se = TRUE)
  #cat("After prediction \n")
  
  
  ## prepare dataframe for plotting
  if (method == "ssnbinomial"){
    ## Curve dataframe
    dd.null = data.frame(Time = points, Count = mod.null$nu/exp(est.null$fit), Group = "null", Subject = "null", SE = est.null$se)
    dd.0 = data.frame(Time = points, Count = mod.0$nu/exp(est.0$fit), Group = "fit.0", Subject = "fit.0", SE = est.0$se)
    dd.1 = data.frame(Time = points, Count = mod.1$nu/exp(est.1$fit), Group = "fit.1", Subject = "fit.1", SE = est.1$se)
    
    ## Confidence interval dataframe
    dd.null.u95 = data.frame(Time = points, Count = mod.null$nu/exp(est.null$fit +1.96*est.null$se), Group = "null.u", Subject = "null.u")
    dd.null.l95 = data.frame(Time = points, Count = mod.null$nu/exp(est.null$fit -1.96*est.null$se), Group = "null.l", Subject = "null.l")
    dd.0.u95 = data.frame(Time = points, Count = mod.0$nu/exp(est.0$fit + 1.96*est.0$se), Group = "fit.0.u", Subject = "fit.0.u")
    dd.0.l95 = data.frame(Time = points, Count = mod.0$nu/exp(est.0$fit - 1.96*est.0$se), Group = "fit.0.l", Subject = "fit.0.l")
    dd.1.u95 = data.frame(Time = points, Count = mod.1$nu/exp(est.1$fit + 1.96*est.1$se), Group = "fit.1.u", Subject = "fit.1.u")
    dd.1.l95 = data.frame(Time = points, Count = mod.1$nu/exp(est.1$fit - 1.96*est.1$se), Group = "fit.1.l", Subject = "fit.1.l")
  }
  else if (method == "ssgaussian" | method == "lowess")
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
  if (method == "ssnbinomial"){
    output = list(dd.null = dd.null, dd.0 = dd.0, dd.1 = dd.1, mod.null = mod.null, mod.0 = mod.0, 
                  mod.1 = mod.1, mod.0.nbinomial.project = mod.0.nbinomial.project, 
                  mod.1.nbinomial.project = mod.1.nbinomial.project, 
                  mod.null.nbinomial.project = mod.null.nbinomial.project,
                  dd.null.u95 = dd.null.u95, dd.null.l95= dd.null.l95,
                  dd.0.u95 = dd.0.u95, dd.0.l95= dd.0.l95, dd.1.u95 = dd.1.u95, dd.1.l95= dd.1.l95)
  }
  else if(method == "ssgaussian")
  {
    output = list(f.stat = f.stat, rss.null = rss.null, rss.full = rss.full, 
                  dd.null = dd.null, dd.0 = dd.0, dd.1 = dd.1, mod.null = mod.null, 
                  mod.0 = mod.0, mod.1 = mod.1, dd.null.u95 = dd.null.u95, dd.null.l95 = dd.null.l95,
                  dd.0.u95 = dd.0.u95, dd.0.l95 = dd.0.l95, dd.1.u95 = dd.1.u95, dd.1.l95= dd.1.l95)
  }
  else if (method == "lowess"){
    output = list(dd.null = dd.null, dd.0 = dd.0, dd.1 = dd.1, mod.null = mod.null, 
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
#' @examples 
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
# #' @examples 
# #' p = c(0.04, 0.01, 0.02, 0.04, 0.06, 0.2, 0.06, 0.04)
# #' sign = c(1, 1, 1, 1, -1, -1, 1, 1)
# #' findSigInterval2(p, threshold = 0.05, sign)
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
#' Fits longitudinal samples from the same group using SSANOVA
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

