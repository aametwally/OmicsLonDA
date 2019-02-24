#' Permute group labels
#'
#' Permutes the group label of the samples in order to construct the testStatistics empirical distibution
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
#' @export
permutationMC2 = function(formula = Count ~ Time, perm.dat, n.perm = 500, method = "ssgaussian", points, parall = FALSE, prefix){
  
  cat("Number of permutation = ", n.perm, "\n")
  
  ## Start permutation
  cat("Start Permutation \n")

  
  pp = list() 
  perm = 0 # to be able to store the value
  n.subjects = length(unique(perm.dat$Subject))
  cat("# of Subjects = ", n.subjects, "\n")
  
  
  ## Run in Parallel
  if(parall == TRUE) {
    max.cores = detectCores()
    cat("# cores = ", max.cores, "\n")
    desired.cores = max.cores - 1		
    cl = makeCluster(desired.cores)
    registerDoParallel(cl)
  } 
  
  
  sample_group = unique(perm.dat[,c("Subject","Group")])
  pp = llply(1:n.perm, function(j){
    sample_group$Group = sample(sample_group$Group, length(sample_group$Group), replace = F)
    for (i in 1:nrow(perm.dat))
    {
      perm.dat[i, "Group"] = sample_group[which(sample_group$Subject == perm.dat[i, "Subject"]),]$Group
    }
    #cat("Group after = ", perm.dat$Group , "\n")
    
    
    
    g.0 = perm.dat[perm.dat$Group == 0, ]
    g.1 = perm.dat[perm.dat$Group == 1, ]
    g.min = max(sort(g.0$Time)[1], sort(g.1$Time)[1])
    g.max = min(sort(g.0$Time)[length(g.0$Time)], sort(g.1$Time)[length(g.1$Time)])
    
    
    
    #cat("\n", "2nd time points = ", points, "\n")
    #cat("g.min = ", g.min, "\n")
    #cat("g.max = ", g.max, "\n")
    
    #cat("min(points) = ", min(points), "\n")
    
    
    
    ### This part to handle the situation when min/max timepoint of eac group from the permuted subjects, lies outside the range of the points vector 
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
      # visualizeFeatureSpline_permute(df = perm.dat, model = perm, method = method, text = "ssgaussian", group.levels = c("A", "B"),
      #                        prefix = paste(prefix, "/", "perm_", j, "_", sep = ""))
      assign(paste("Model", j, sep = "_"), perm)
    }
  }, .parallel = parall, .progress = "text", .inform = TRUE,
  .paropts = list(.export=ls(.GlobalEnv),
                  .packages=.packages(all.available=T)))
  
  
  if(parall == TRUE) {
    stopCluster(cl)
  }
  
  pp[sapply(pp, is.null)] = NULL
  return(pp)
}  




# permutationMC = function(formula = Count ~ Time, perm.dat, n.perm = 500, method = "ssnbinomial", points, parall = FALSE, prefix){
#   
#   cat(" number of permutation = ", n.perm, "\n")
#   #cat("beging group=  ", perm.dat$Group, "\n")
#   
#   ## Start permutation
#   cat("Start Permutation \n")
#   
#   #cat("points = ", points, "\n")
#   
#   pp = list() 
#   perm = 0 # to be able to store the value
#   n.subjects = length(unique(perm.dat$Subject))
#   cat("# of Subjects = ", n.subjects, "\n")
#   
#   
#   ## Run in Parallel
#   if(parall == TRUE) {
#     max.cores = detectCores()
#     cat("# cores = ", max.cores, "\n")
#     desired.cores = max.cores - 1		
#     cl = makeCluster(desired.cores)
#     registerDoParallel(cl)
#   } 
#   
#   
#   sample_group = unique(perm.dat[,c("Subject","Group")])
#   pp = llply(1:n.perm, function(j){
#     
#     # for (i in levels(perm.dat$Subject)){
#     #   perm.uniq.len = 1
#     #   while(perm.uniq.len == 1){
#     #     perm.dat[which(perm.dat$Subject == i),]$Group = rep(sample(c(0,1),1), sum(perm.dat$Subject == i))# ,  replace = TRUE) #rep(sample(1:2), each = time.point)
#     #     perm.uniq.len = length(unique(perm.dat$Group))
#     #   }
#     # }
#     
#     # cat("Group before = ", perm.dat$Group , "\n")
#     # perm.dat$Group = sample(perm.dat$Group, length(perm.dat$Group), replace = F)
#     # cat("Group after = ", perm.dat$Group , "\n")
#     
#     ## Perumute group labels
#     #cat("Group before = ", perm.dat$Group , "\n")
#     sample_group$Group = sample(sample_group$Group, length(sample_group$Group), replace = F)
#     for (i in 1:nrow(perm.dat))
#     {
#       perm.dat[i, "Group"] = sample_group[which(sample_group$Subject == perm.dat[i, "Subject"]),]$Group
#     }
#     #cat("Group after = ", perm.dat$Group , "\n")
#     
#     
#     
#     g.0 = perm.dat[perm.dat$Group == 0, ]
#     g.1 = perm.dat[perm.dat$Group == 1, ]
#     g.min = max(sort(g.0$Time)[1], sort(g.1$Time)[1])
#     g.max = min(sort(g.0$Time)[length(g.0$Time)], sort(g.1$Time)[length(g.1$Time)])
#     
#     
#     
#     #cat("\n", "2nd time points = ", points, "\n")
#     #cat("g.min = ", g.min, "\n")
#     #cat("g.max = ", g.max, "\n")
#     
#     #cat("min(points) = ", min(points), "\n")
#     
#     
#     
#     ### TODO: Fix this issue. The method should be able to handle this situation 
#     if(g.min > min(points) | g.max < max(points))
#     {
#       cat("\n")
#       cat("g.min = ", g.min, "   min(points) = ", min(points), "\n")
#       cat("g.max = ", g.max, "   max(points) = ", max(points), "\n")
# 
#       cat("Special Case: generated permutation is out of range \n")
#       assign(paste("Model", j, sep = "_"), NULL)
#     }
#     # else if (length(which(sum(g.0$Count) == 0 | sum(g.1$Count)==0)))
#     # {
#     #   cat("Special Case: zero for all variable of one group \n")
#     #   assign(paste("Model", j, sep = "_"), NULL)
#     # }
#     else
#     {
#       #cat("In else", "\n")
#       perm = curveFitting(formula, df = perm.dat, method = method, points)
#       # visualizeFeatureSpline_permute(df = perm.dat, model = perm, method = method, text = "permutation", group.levels = c("A", "B"),
#       #                        prefix = paste(prefix, "/", "perm_", j, "_", sep = ""))
#       assign(paste("Model", j, sep = "_"), perm)
#     }
#   }, .parallel = parall, .progress = "text", .inform = TRUE,
#   .paropts = list(.export=ls(.GlobalEnv),
#                   .packages=.packages(all.available=T)))
#   
#   
#   if(parall == TRUE) {
#     stopCluster(cl)
#   }
#   
#   pp[sapply(pp, is.null)] = NULL
#   return(pp)
# }  
