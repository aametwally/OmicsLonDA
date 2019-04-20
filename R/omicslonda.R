#' Omics Longitudinal Differential Analysis for one feature
#'
#' Find significant time intervals of omic feature
#' 
#' @param se_object SummarizedExperiment object contains omics count/level matrix 
#' and  metadata contains (subject, time, group, and any any other covariates)
#' @param n.perm number of permutations.
#' @param fit.method fitting method (ssguassian).
#' @param points points at which the prediction should happen.
#' @param text Feature's name.
#' @param parall boolean to indicate whether to use multicore.
#' @param pvalue.threshold p-value threshold cutoff for identifing significant
#' time intervals.
#' @param adjust.method multiple testing correction method.
#' @param time.unit time unit used in the Time vector (hours, days, weeks,
#' months, etc.)
#' @param col two color to be used for the two groups (eg., c("red", "blue")).
#' @param ylabel text to be shown on the y-axis of all generated figures
#' (default: "Normalized Count")
#' @param prefix prefix to be used to create directory for the analysis results
#' @return a list of the significant time intervals for the tested
#' feature, fitted model for each group, null distribution of the test statistic 
#' of the tested feature, and the original input data.
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment
#' @importFrom methods is
#' @importFrom BiocGenerics rowSums colSums rowMeans colMeans
#' @import BiocParallel
#' @importFrom stats p.adjust
#' @import zoo
#' @import SummarizedExperiment
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples
#' library(SummarizedExperiment)
#' data(omicslonda_data_example)
#' omicslonda_se_object_adjusted = adjustBaseline(
#'                  se_object = omicslonda_data_example$omicslonda_se_object)
#' omicslonda_test_object = omicslonda_se_object_adjusted[1,]
#' points = seq(1, 500, length.out = 500)
#' res = omicslonda(se_object = omicslonda_test_object, n.perm = 10,
#'                  fit.method = "ssgaussian", points = points, text = "Feature_1",
#'                  parall = FALSE, pvalue.threshold = 0.05, 
#'                  adjust.method = "BH", time.unit = "days",
#'                  ylabel = "Normalized Count",
#'                  col = c("blue", "firebrick"), prefix = tempfile())
#' @export
omicslonda = function(se_object = NULL, n.perm = 500,
                        fit.method = "ssgaussian", points = NULL, text = "FeatureName",
                        parall = FALSE, pvalue.threshold = 0.05, 
                        adjust.method = "BH", time.unit = "days",
                        ylabel = "Normalized Count",
                        col = c("blue", "firebrick"), prefix = "Test")
{
    message("Start OmicsLonDA")
    
    rowSums <- BiocGenerics::rowSums
    colSums <- BiocGenerics::colSums
    rowMeans <- BiocGenerics::rowMeans
    colMeans <- BiocGenerics::colMeans
    
    ### validate se_object
    stopifnot(is(se_object, "SummarizedExperiment"))
    stopifnot(all(c("Subject", "Time", "Group") %in% colnames(colData(se_object))))
    ## validate fit.method
    stopifnot(fit.method %in% c("ssgaussian"))
    ## validate pvalue.threshold
    stopifnot(pvalue.threshold >= 0, pvalue.threshold <= 1)
    ## validate pvalue adjustment method
    stopifnot(adjust.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
    ## validate col
    stopifnot(length(col) == 2)
    ## validate points
    stopifnot(length(points) >= 2, is(points, "numeric"))
   
  
    if (!dir.exists(prefix)){
        dir.create(file.path(prefix))
    }
    
  
    dt = data.frame(colData(se_object))
    dt$Count = as.vector(assay(se_object))
  
    # Extract groups
    Group = as.character(dt$Group)
    group.levels = sort(unique(Group))
    print(group.levels)
    
    
    if(length(group.levels) > 2){
        stop("You have more than two phenotypes.")
    } else if(length(group.levels) < 2){
        stop("You have less than two phenotypes.")
    }
    
    
    gr.1 = as.character(group.levels[1])
    gr.2 = as.character(group.levels[2])

    
    df = dt
    levels(df$Group) = c(levels(df$Group), "0", "1")
    df$Group[which(df$Group == gr.1)] = 0
    df$Group[which(df$Group == gr.2)] = 1
    
    
    group.0 = df[df$Group == 0, ]
    group.1 = df[df$Group == 1, ]
    points.min = max(sort(group.0$Time)[1], sort(group.1$Time)[1])
    points.max = min(sort(group.0$Time)[length(group.0$Time)],
                    sort(group.1$Time)[length(group.1$Time)])
    points = points[which(points >= points.min & points <= points.max)]
    
    
    message("points.min = ", points.min)
    message("points.max = ", points.max)
    
    message("Start Curve Fitting") 
    
    if (fit.method == "ssgaussian")
    {
        message("Fitting: Smoothing Spline Gaussian Regression")
        model = tryCatch({
            curveFitting(formula = Count ~ Time, df, fit.method = "ssgaussian", points)
            },  error = function(err) {
              stop("ERROR in gss = ", err)
              #print(paste("ERROR in gss = ", err, sep="")); 
              #return("ERROR")
            })
    } else {
        stop("You have entered unsupported fitting method")
    }
    
    
    ### Test Statistic
    stat = testStat(model)$testStat
    
    ## Permutation 
    perm  = permutationMC(formula = Count ~ Time, perm.dat = df, n.perm = n.perm,
                            fit.method = fit.method, points = points,
                            parall = parall, prefix = prefix)

    test.stat.prem = testStatPermutation(perm)
    t1 = do.call(rbind, test.stat.prem)
    

    t2 = unlist(t1[,1])
    t3 = as.vector(t2)
    length(t3)
    pvalue.test.stat = vapply(head(seq_along(points), -1), function(i){
        if(stat[i]>=0)
        {
        sum(t3 > stat[i])/length(t3)
        }
        else if(stat[i]<0)
        {
        sum(t3 < stat[i])/length(t3)
        }
    }, 1)

    ## Adjust p-values
    if(adjust.method == "qvalue"){
        #adjusted.pvalue = qvalue(pvalue.test.stat, pfdr = TRUE)
        adjusted.pvalue = .data$qvalue_truncp(pvalue.test.stat)$qvalues
    }else{
        adjusted.pvalue = p.adjust(pvalue.test.stat, method = adjust.method)
    }
        
    



    interval = findSigInterval(adjusted.pvalue, threshold = pvalue.threshold,
                                sign = sign(stat))
    
    st = points[interval$start]
    en = points[interval$end + 1]
  
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
                        interval.start = interval.start,
                        interval.end = interval.end,
                        avg.mod0.count = avg.mod0.count,
                        avg.mod1.count = avg.mod1.count, 
                        foldChange = foldChange,
                        testStat = stat, testStat.abs = abs(stat),
                        testStat.sign = sign(stat), dominant = dominant,
                        intervals.pvalue = pvalue.test.stat,
                        adjusted.pvalue = adjusted.pvalue, points = points)
    output.details$points = output.details$points[-length(output.details$points)]
    
    output.summary = data.frame(feature = rep(text, length(interval$start)),
                                start = st, end = en,
                                dominant=interval$dominant,
                                pvalue = interval$pvalue)
    
    omicslonda_output = list(df = dt, details = output.details, summary = output.summary, 
                             model = model, distribution = t3, start = st, end = en)

    return(omicslonda_output)
}
