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
#' @return returns the fitted model for all the permutations
#' @import plyr
#' @import utils
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' data(diff_simulatedDataset_norm)
#' df = diff_simulatedDataset_norm[[1]]
#' gr.1 = as.character(group.levels[1])
#' gr.2 = as.character(group.levels[2])
#' levels(df$Group) = c(levels(df$Group), "0", "1")
#' df$Group[which(df$Group == gr.1)] = 0
#' df$Group[which(df$Group == gr.2)] = 1
#' group.0 = df[df$Group == 0, ]
#' group.1 = df[df$Group == 1, ]
#' points.min = max(sort(group.0$Time)[1], sort(group.1$Time)[1])
#' points.max = min(sort(group.0$Time)[length(group.0$Time)],
#'                  sort(group.1$Time)[length(group.1$Time)])
#' points = points[which(points >= points.min & points <= points.max)]
#' model = curveFitting(formula = formula, df, method= "ssgaussian", points)
#' perm  = permutationMC2(formula = Count ~ Time, perm.dat = df, n.perm = 10,
#' method = "ssgaussian", points = points, parall = "FALSE",
#' prefix = "Test")
#' @export
permutationMC2 = function(formula = Count ~ Time, perm.dat, n.perm = 500,
                        method="ssgaussian", points, parall = FALSE, prefix){
    
    cat("Number of permutation = ", n.perm, "\n")
    
    ## Start permutation
    cat("Start Permutation \n")

    
    pp = list() 
    perm = 0 # to be able to store the value
    n.subjects = length(unique(perm.dat$Subject))
    cat("# of Subjects = ", n.subjects, "\n")
    
    
    ## Run in Parallel
    if(parall == TRUE) {max.cores = detectCores()
        cat("# cores = ", max.cores, "\n")
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
        #cat("Group after = ", perm.dat$Group , "\n")
        
        
        
        g.0 = perm.dat[perm.dat$Group == 0, ]
        g.1 = perm.dat[perm.dat$Group == 1, ]
        g.min = max(sort(g.0$Time)[1], sort(g.1$Time)[1])
        g.max = min(sort(g.0$Time)[length(g.0$Time)],
                sort(g.1$Time)[length(g.1$Time)])
        
        
        
        #cat("\n", "2nd time points = ", points, "\n")
        #cat("g.min = ", g.min, "\n")
        #cat("g.max = ", g.max, "\n")
        
        #cat("min(points) = ", min(points), "\n")
        
        
        
        ### This part to handle the situation when min/max timepoint of each
        ### group from the permuted subjects, lies outside the range of the
        ### points vector
        if(g.min > min(points) | g.max < max(points))
        {
        #cat("\n")
        #cat("old: g.min = ", g.min, "   min(points) = ", min(points), "\n")
        #cat("old: g.max = ", g.max, "   max(points) = ", max(points), "\n")
        points = points[which(points>=g.min & points<=g.max)]
        #cat("new: g.min = ", g.min, "   min(points) = ", min(points), "\n")
        #cat("new: g.max = ", g.max, "   max(points) = ", max(points), "\n")
        
        perm = curveFitting(formula, df = perm.dat, method = method, points)
        assign(paste("Model", j, sep = "_"), perm)
        
        #cat("Special Case: generated permutation is out of range \n")
        # assign(paste("Model", j, sep = "_"), NULL)
        }
        # else if (length(which(sum(g.0$Count) == 0 | sum(g.1$Count)==0)))
        # {
        #   cat("Special Case: zero for all variable of one group \n")
        #   assign(paste("Model", j, sep = "_"), NULL)
        # }
        else
        {
        #cat("In else", "\n")
        perm = curveFitting(formula, df = perm.dat, method = method, points)
        # visualizeFeatureSpline_permute(df = perm.dat, model = perm,
        #                                method = method, text = "ssgaussian",
        #                                group.levels = c("A", "B"),
        #                                prefix = paste(prefix, "/", "perm_", j,
        #                                "_", sep = ""))
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






# bootstrapOmicslonda = function(data, index){
#   #print("here\n")
#   bs.df = data[index, ]
#   ## TODO: This needs to be changed to reflect the nromalized formula
#   #formula = Count ~ Time
#   method = "ssgaussian"
#   
#   
#   ## Run in Parallel
#   if(.data$parall == TRUE) {
#     max.cores = detectCores()
#     cat("# cores = ", max.cores, "\n")
#     desired.cores = max.cores - 1
#     cl = makeCluster(desired.cores)
#     registerDoParallel(cl)
#   } 
#   
#   
#   g.0 = bs.df[bs.df$Group == 0, ]
#   g.1 = bs.df[bs.df$Group == 1, ]
#   g.min = max(sort(g.0$Time)[1], sort(g.1$Time)[1])
#   g.max = min(sort(g.0$Time)[length(g.0$Time)],sort(g.1$Time)[length(
#                                                               g.1$Time)])
#   
#   
#   
#   ### This part to handle the situation when min/max timepoint of eac group
#   ### from the permuted subjects, lies outside the range of the points vector 
#   if(g.min > min(points) | g.max < max(points))
#   {
#     #cat("\n")
#     #cat("old: g.min = ", g.min, "   min(points) = ", min(points), "\n")
#     #cat("old: g.max = ", g.max, "   max(points) = ", max(points), "\n")
#     points = points[which(points>=g.min & points<=g.max)]
#     #cat("new: g.min = ", g.min, "   min(points) = ", min(points), "\n")
#     #cat("new: g.max = ", g.max, "   max(points) = ", max(points), "\n")
#     
#     bootstrapped = curveFitting(formula, df = bs.df, method = method, points)
#     #assign(paste("Model", j, sep = "_"), bootstrapped)
#     
#     #cat("Special Case: generated permutation is out of range \n")
#     # assign(paste("Model", j, sep = "_"), NULL)
#   }  else
#   {
#     #cat("In else", "\n")
#     bootstrapped = curveFitting(formula, df = bs.df, method = method, points)
#     # visualizeFeatureSpline_permute(df = perm.dat, model = perm,
#     #                        method = method, text = "ssgaussian",
#     #                        group.levels = c("A", "B"),
#     #                        prefix = paste(prefix, "/", "perm_", j, "_",
#     #                                       sep = ""))
#     #assign(paste("Model", j, sep = "_"), bootstrapped)
#   }
#   
#   
#   bs.stat = testStat(bootstrapped)$testStat
#   
#   if(.data$parall == TRUE) {
#     stopCluster(cl)
#   }
#   
#   print(bs.stat[1])
#   return(bs.stat[1])
# }