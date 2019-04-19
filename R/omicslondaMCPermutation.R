#' Permute group labels
#'
#' Permutes the group label of the samples in order to construct the
#' testStatistics empirical distibution
#'
#' @param formula formula to be passed to the regression model
#' @param perm.dat dataframe has the Count, Group, Subject, Time
#' @param n.perm number of permutations
#' @param method The fitting method (ssgaussian)
#' @param points The points at which the prediction should happen
#' @param parall boolean to indicate whether to use multicore.
#' @param prefix prefix to be used to create directory for the analysis results
#' @return a list of the fitted model for each group for all the permutations
#' @importFrom SummarizedExperiment colData assay SummarizedExperiment
#' @import plyr
#' @import utils
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
#'                       method = "ssgaussian", points = points,
#'                       parall = FALSE, prefix = tempfile())
#' @export
permutationMC = function(formula = Count ~ Time, perm.dat, n.perm = 500,
                        method="ssgaussian", points, parall = FALSE, prefix){
    
    
    message("Start Permutation \n")
    
    message("Number of permutation = ", n.perm, "\n")
    
    pp = list() 
    perm = 0 # to be able to store the value
    n.subjects = length(unique(perm.dat$Subject))
    message("# of Subjects = ", n.subjects, "\n")
    
    
    ## Run in Parallel
    if(parall == TRUE) {max.cores = detectCores()
        message("# cores = ", max.cores, "\n")
        desired.cores = max.cores - 1
        cl = makeCluster(desired.cores)
        registerDoParallel(cl)
    } 
    
    
    sample_group = unique(perm.dat[,c("Subject","Group")])
    pp = llply(seq_len(n.perm), function(j){
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
        #message("\n")
        #message("old: g.min = ", g.min, "   min(points) = ", min(points), "\n")
        #message("old: g.max = ", g.max, "   max(points) = ", max(points), "\n")
        points = points[which(points>=g.min & points<=g.max)]
        #message("new: g.min = ", g.min, "   min(points) = ", min(points), "\n")
        #message("new: g.max = ", g.max, "   max(points) = ", max(points), "\n")
        
        perm = curveFitting(formula, df = perm.dat, method = method, points)
        assign(paste("Model", j, sep = "_"), perm)
        
        #message("Special Case: generated permutation is out of range \n")
        # assign(paste("Model", j, sep = "_"), NULL)
        }
        # else if (length(which(sum(g.0$Count) == 0 | sum(g.1$Count)==0)))
        # {
        #   message("Special Case: zero for all variable of one group \n")
        #   assign(paste("Model", j, sep = "_"), NULL)
        # }
        else
        {
            perm = curveFitting(formula, df = perm.dat, method = method, points)
            assign(paste("Model", j, sep = "_"), perm)
        }
    }, .parallel = parall, .progress = "text", .inform = TRUE,
    .paropts = list(.export=ls(.GlobalEnv),
                    .packages=.packages(all.available=TRUE)))
    
    
    if(parall == TRUE) {
        stopCluster(cl)
    }
    
    pp[vapply(pp, is.null, logical(1))] = NULL
    return(pp)
}  