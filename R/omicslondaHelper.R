#' Adjust for baseline differences using Centered-Log Ratio
#'
#' Adjust for baseline differences using Centered-Log Ratio
#' 
##' @param se_object SummarizedExperiment object contains omics count/level matrix 
#' and  metadata 
#' @return a SummarizedExperiment object with the adjusted baseline assay
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment assays
#' @importFrom methods is
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' data(omicslonda_data_example)
#' adjusted_se = adjustBaseline(omicslonda_data_example$omicslonda_se_object)
#' @export
adjustBaseline = function(se_object = NULL){
  
    ### validate se_object
    stopifnot(is(se_object, "SummarizedExperiment"))
    stopifnot(all(c("Subject", "Time", "Group") %in% colnames(colData(se_object))))

    subjects = unique(colData(se_object)$Subject)
    
    df = assay(se_object)
    metadata = colData(se_object)
    
    ### Add psuedo count to prevent NA in CLR
    df = df + 0.000000001
    
    for(i in seq_len(nrow(df)))
    {
        for(subj in subjects)
        {
            m = min(metadata[metadata$Subject==subj, "Time"])
            subj_samples = rownames(metadata)[metadata$Subject==subj]
            df[i, subj_samples] = log(df[i, subj_samples]/df[i, metadata$Subject == subj & metadata$Time == m])
        }
    }
    
    se_object_normalized = se_object
    assay(se_object_normalized) = as.matrix(df)
    return(se_object_normalized)
}




#' Fit longitudinal data
#'
#' Fits longitudinal samples from the same group using negative binomial
#' smoothing splines or LOWESS
#' 
#' @param formula formula to be passed to the regression model
#' @param df dataframe has the Count, Group, Subject, Time
#' @param fit.method fitting method (ssgaussian)
#' @param points points at which the prediction should happen
#' @return a list that contains fitted smoothing spline for each group 
#' along with 95% confidence intervals
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment assays
#' @importFrom stats predict
#' @import gss
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' library(SummarizedExperiment)
#' data("omicslonda_data_example")
#' omicslonda_se_object_adjusted = adjustBaseline(
#'                  se_object = omicslonda_data_example$omicslonda_se_object)
#' se_object = omicslonda_se_object_adjusted[1,]
#' dt = data.frame(colData(se_object))
#' dt$Count = as.vector(assay(se_object))
#' Group = as.character(dt$Group)
#' group.levels = sort(unique(Group))
#' gr.1 = as.character(group.levels[1])
#' gr.2 = as.character(group.levels[2])
#' df = dt
#' levels(df$Group) = c(levels(df$Group), "0", "1")
#' df$Group[which(df$Group == gr.1)] = 0
#' df$Group[which(df$Group == gr.2)] = 1
#' group.0 = df[df$Group == 0, ]
#' group.1 = df[df$Group == 1, ]
#' points = seq(1, 500)
#' model = curveFitting(formula = Count ~ Time, df = df, fit.method = "ssgaussian", points = points)
#' @export
curveFitting = function(formula = Count ~ Time, df = "NULL", fit.method = "ssgaussian",
                        points = NULL){
    
    ## Seprate the two groups
    group.null = df
    group.0 = df[df$Group==0, ]
    group.1 = df[df$Group==1, ]
    
    ## Fitting 
    if(fit.method == "ssgaussian"){
        
        ### Difference between random = mkran(~1|Subject, group.null) and
        ### random = ~1|Subject
        # ran.null = mkran(~1|Subject, group.null)
        # ran.0 = mkran(~1|Subject, group.0)
        # ran.1 = mkran(~1|Subject, group.1)
        
        
        # ## null model
        mod.null = ssanova(formula, data = group.null, skip.iter=TRUE)
        
        ## full model
        mod.0 = ssanova(formula, data = group.0, skip.iter=TRUE)
        mod.1 = ssanova(formula, data = group.1, skip.iter=TRUE)
        
        
        ## Calculate goodness of fit F-statistic the fitted model
        rss.null = summary(mod.null)$rss
        rss.full = summary(mod.0)$rss+summary(mod.1)$rss
        f.stat = (rss.null - rss.full)/rss.null
    }
    
    
    ## Estimate values at the provided time points
    est.null = predict(mod.null, data.frame(Time = points), include =
                        c("1", "Time"), se = TRUE)
    est.0 = predict(mod.0, data.frame(Time = points), include = c("1", "Time"),
                    se = TRUE)
    est.1 = predict(mod.1, data.frame(Time = points), include = c("1", "Time"),
                    se = TRUE)
    
    
    
    ## prepare dataframe for plotting
    if (fit.method == "ssgaussian")
    {
        ## Curve dataframe
        dd.null = data.frame(Time = points, Count = est.null$fit,
                            Group = "NULL", Subject = "NULL", SE = est.null$se)
        dd.0 = data.frame(Time = points, Count = est.0$fit, Group = "fit.0",
                        Subject = "fit.0", SE = est.0$se)
        dd.1 = data.frame(Time = points, Count = est.1$fit, Group = "fit.1",
                        Subject = "fit.1", SE = est.1$se)
        
        ## Confidence interval dataframe
        dd.null.u95 = data.frame(Time = points,
                                Count = (est.null$fit + 1.96*est.null$se),
                                Group = "null.u", Subject = "null.u")
        dd.null.l95 = data.frame(Time = points,
                                Count = (est.null$fit - 1.96*est.null$se),
                                Group = "null.l", Subject = "null.l")
        dd.0.u95 = data.frame(Time = points,
                            Count = (est.0$fit + 1.96*est.0$se),
                            Group = "fit.0.u", Subject = "fit.0.u")
        dd.0.l95 = data.frame(Time = points,
                            Count = (est.0$fit - 1.96*est.0$se),
                            Group = "fit.0.l", Subject = "fit.0.l")
        dd.1.u95 = data.frame(Time = points,
                            Count = (est.1$fit + 1.96*est.1$se),
                            Group = "fit.1.u", Subject = "fit.1.u")
        dd.1.l95 = data.frame(Time = points,
                            Count = (est.1$fit - 1.96*est.1$se),
                            Group = "fit.1.l", Subject = "fit.1.l")
    } 
    
    if(fit.method == "ssgaussian")
    {
        output = list(f.stat = f.stat, rss.null = rss.null, rss.full = rss.full,
                    dd.null = dd.null, dd.0 = dd.0, dd.1 = dd.1,
                    mod.null = mod.null, mod.0 = mod.0, mod.1 = mod.1,
                    dd.null.u95 = dd.null.u95, dd.null.l95 = dd.null.l95,
                    dd.0.u95 = dd.0.u95, dd.0.l95 = dd.0.l95,
                    dd.1.u95 = dd.1.u95, dd.1.l95= dd.1.l95)
    }
    
    
    return(output)
}


#' Calculate Test-Statistic of each feature's time interval
#'
#' Calculate Test-Statistic of each feature's time interval
#' 
#' @param curve.fit.df gss data object of the fitted spline
#' @return a list that has the test statistic for each time interval
#' @import pracma
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' library(SummarizedExperiment)
#' data("omicslonda_data_example")
#' model = omicslonda_data_example$omicslonda_results$model
#' stat = testStat(model)$testStat
#' @export
testStat = function(curve.fit.df){
    size = length(curve.fit.df$dd.null$Time)
    testStat = numeric(size - 1)
    
    testStat = sapply(seq_len(size-1), function(i){
      testStat.null = trapz(curve.fit.df$dd.null$Time[i:(i+1)],
                            curve.fit.df$dd.null$Count[i:(i+1)])
      testStat.0 = trapz(curve.fit.df$dd.0$Time[i:(i+1)],
                         curve.fit.df$dd.0$Count[i:(i+1)])
      testStat.1 = trapz(curve.fit.df$dd.1$Time[i:(i+1)],
                         curve.fit.df$dd.1$Count[i:(i+1)])
      
      testStat = (testStat.0 - testStat.1) /
        sqrt((mean(curve.fit.df$dd.0$SE[i:(i+1)]))^2 +
               (mean(curve.fit.df$dd.1$SE[i:(i+1)]))^2)
      if(is.na(testStat)){
        testStat = 0
      }
      
      testStat
    })

    return(list(testStat = testStat))
}


#' Calculate testStat of each feature's time interval for all permutations
#' 
#' @param perm list has all the permutated models
#' @return a list of test statistic for each time interval for all all permutations
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' library(SummarizedExperiment)
#' data("omicslonda_data_example")
#' omicslonda_se_object_adjusted = adjustBaseline(
#'     se_object = omicslonda_data_example$omicslonda_se_object)
#' omicslonda_test_object = omicslonda_se_object_adjusted[1,]
#' se_object = omicslonda_test_object
#' dt = data.frame(colData(se_object))
#' dt$Count = as.vector(assay(se_object))
#' Group = as.character(dt$Group)
#' group.levels = sort(unique(Group))
#' gr.1 = as.character(group.levels[1])
#' gr.2 = as.character(group.levels[2])
#' df = dt
#' levels(df$Group) = c(levels(df$Group), "0", "1")
#' df$Group[which(df$Group == gr.1)] = 0
#' df$Group[which(df$Group == gr.2)] = 1
#' group.0 = df[df$Group == 0, ]
#' group.1 = df[df$Group == 1, ]
#' points = seq(100, 130)
#' perm  = permutationMC(formula = Count ~ Time, perm.dat = df, n.perm = 10,
#'                       fit.method = "ssgaussian", points = points,
#'                       parall = FALSE, prefix = tempfile())
#' test.stat.prem = testStatPermutation(perm)
#' @export
testStatPermutation = function(perm)
{
  testStat.list <- lapply(perm, testStat)
  return(testStat.list)
}


#' Find Significant Interval based on testStat
#'
#' Find Significant Interval based on testStat
#' 
#' @param adjusted.pvalue vector of the adjusted p-value
#' @param threshold p-value cut off
#' @param sign vector hold area sign of each time interval 
#' @return returns a list of the start and end points of all significant
#' time intervals
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' library(SummarizedExperiment)
#' padjusted = abs(rnorm(10, mean = 0.05, sd = 0.04))
#' sign = sample(x = c(1,-1), 10, replace = TRUE)
#' intervals = findSigInterval(adjusted.pvalue = padjusted, 
#'                         threshold = 0.05, sign = sign)
#' @export
findSigInterval = function(adjusted.pvalue, threshold = 0.05, sign)
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
        message("No Significant Intevals Found")
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
                p = c(p, mean(adjusted.pvalue[start[length(start)] :
                                                end[length(end)]]))
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