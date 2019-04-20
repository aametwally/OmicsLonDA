#' Permute group labels
#'
#' Permutes the group label of the samples in order to construct the
#' testStatistics empirical distibution
#'
#' @param formula formula to be passed to the regression model
#' @param perm.dat dataframe has the Count, Group, Subject, Time
#' @param n.perm number of permutations
#' @param fit.method The fitting method (ssgaussian)
#' @param points The points at which the prediction should happen
#' @param parall boolean to indicate whether to use multicore.
#' @param prefix prefix to be used to create directory for the analysis results
#' @return a list of the fitted model for each group for all the permutations
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment
#' @import plyr
#' @import utils
#' @import gss
#' @importFrom parallel detectCores
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
#' @export
permutationMC = function(formula = Count ~ Time, perm.dat = NULL, n.perm = 500,
                        fit.method = "ssgaussian", points, parall = FALSE, prefix = "Test"){
    
    
    message("Start Permutation")
    message("Number of permutation = ", n.perm)
    
    pp = list() 
    perm = 0 # to be able to store the value
    n.subjects = length(unique(perm.dat$Subject))
    message("# of Subjects = ", n.subjects)
    
    
    ## Parallelization
    if(parall == TRUE) {
        max.cores = detectCores()
        desired.cores = max.cores - 1
        message("# of used cores = ", desired.cores)
       
        param = SnowParam(workers = desired.cores, progressbar = TRUE, type = "SOCK")
        #param = MulticoreParam(workers = desired.cores)
    } else{
        param = SerialParam()
    }
    
    
    sample_group = unique(perm.dat[,c("Subject","Group")])
    permute = function(j, curveFitting){
        sample_group$Group = sample(sample_group$Group, 
                                    length(sample_group$Group), replace = FALSE)
        for (i in seq_len(nrow(perm.dat)))
        {
            perm.dat[i, "Group"] = sample_group[which(sample_group$Subject ==
                                                          perm.dat[i, "Subject"]),]$Group
        }
        
        g.0 = perm.dat[perm.dat$Group == 0, ]
        g.1 = perm.dat[perm.dat$Group == 1, ]
        g.min = max(sort(g.0$Time)[1], sort(g.1$Time)[1])
        g.max = min(sort(g.0$Time)[length(g.0$Time)],
                    sort(g.1$Time)[length(g.1$Time)])
        
        ### This part to handle the situation when min/max timepoint of each
        ### group from the permuted subjects, lies outside the range of the
        ### points vector
        if(g.min > min(points) | g.max < max(points))
        {
            points = points[which(points>=g.min & points<=g.max)]
            perm = curveFitting(formula, df = perm.dat, fit.method = fit.method, points)
            assign(paste("Model", j, sep = "_"), perm)
        }
        else
        {
            perm = curveFitting(formula, df = perm.dat, fit.method = fit.method, points)
            assign(paste("Model", j, sep = "_"), perm)
        }
    }
    
    pp = bplapply(seq_len(n.perm), permute, BPPARAM = param, curveFitting = curveFitting)
    pp[vapply(pp, is.null, logical(1))] = NULL
    return(pp)
}  