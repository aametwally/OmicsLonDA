#' Permute group labels
#'
#' Permutes the group label of the samples in order to construct the AR empirical distibution
#'
#' @param perm.dat dataframe has the Count, Group, Subject, Time
#' @param n.perm number of permutations
#' @param method The fitting method (negative binomial, LOWESS)
#' @param points The points at which the prediction should happen
#' @param parall boolean to indicate whether to use multicore.
#' @return returns the fitted model for all the permutations
#' @import plyr
#' @import utils
#' @references
#' Ahmed Metwally (ametwall@stanford.edu)
#' @examples 
#' data(metalonda_test_data)
#' n.sample = 5
#' n.timepoints = 10
#' n.perm = 3
#' n.group = 2
#' Group = factor(c(rep(0, n.sample*n.timepoints), rep(1, n.sample*n.timepoints)))
#' Time = rep(rep(1:n.timepoints, times = n.sample), 2)
#' Subject = factor(rep(1:(2*n.sample), each = n.timepoints))
#' points = seq(1, 10, length.out = 10)
#' aggregate.df = data.frame(Count = metalonda_test_data[1,], Time = Time, Group = Group, Subject = Subject)
#' prm = permutation(aggregate.df, n.perm = 3, method = "ssnbinomial", points)
#' @export
permutation = function(formula = Count ~ Time, perm.dat, n.perm = 500, method = "ssnbinomial", points, parall = FALSE){

  ## Start permutation
  cat("Start Permutation \n")
  
  #cat("points = ", points, "\n")
  
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
  
  
  pp = llply(1:n.perm, function(j){
    for (i in levels(perm.dat$Subject)){
      perm.uniq.len = 1
      
      while(perm.uniq.len == 1){
        perm.dat[which(perm.dat$Subject == i),]$Group = rep(sample(c(0,1),1), sum(perm.dat$Subject == i))# ,  replace = TRUE) #rep(sample(1:2), each = time.point)
        perm.uniq.len = length(unique(perm.dat$Group))
      }
    }
    
    g.0 = perm.dat[perm.dat$Group == 0, ]
    g.1 = perm.dat[perm.dat$Group == 1, ]
    g.min = max(sort(g.0$Time)[1], sort(g.1$Time)[1])
    g.max = min(sort(g.0$Time)[length(g.0$Time)], sort(g.1$Time)[length(g.1$Time)])
    
    
    
    #cat("\n", "2nd time points = ", points, "\n")
    #cat("g.min = ", g.min, "\n")
    #cat("g.max = ", g.max, "\n")
    
    #cat("min(points) = ", min(points), "\n")
    
    
    if(g.min > min(points) | g.max < max(points))
    {
      cat("Special Case: generated permutation is out of range \n")
      assign(paste("Model", j, sep = "_"), NULL)
    } 
    else if (length(which(sum(g.0$Count) == 0 | sum(g.1$Count)==0)))
    {
      cat("Special Case: zero for all variable of one group \n")
      assign(paste("Model", j, sep = "_"), NULL)
    }
    else
    {
      #cat("In else", "\n")
      perm = curveFitting(formula, df = perm.dat, method = method, points)
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






# 
# perm.dat = df
# n.perm = 5
# method = "ssnbinomial"
# 
# ### Extract Group vector for each subject
# 
# groups_subjects = unique(df[, c("Group", "Subject")])$Group
# 
# 
# combn(c("A", "B", "A"), 2)
# 
# 
# require(gtools)
# permutations(n = length(groups_subjects), r = length(groups_subjects), v = groups_subjects)
# permutations(n = length(groups_subjects), r = length(groups_subjects), v = groups_subjects, set = FALSE)
# permutations(n = 3, r = 3, v = c("B", "B", "C"), set = FALSE)
# 
# 
# 
# ss = table(interaction(df$Group,df$Subject))
# mm = ss[ss>0]
# ## generate differen permutations
# 
# 
# ### Sample n.perm from the generated cominations
# 
# 
# ## for each perm, assign the new groups to the subjects
# 
# 
# ## Do curve fitting and proceed
# 
# 
# permutationMCMC = function(formula = Count ~ Time, perm.dat, n.perm = 500, method = "ssnbinomial", points, lev, parall = FALSE){
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
#   pp = llply(1:n.perm, function(j){
#     for (i in levels(perm.dat$Subject)){
#       perm.uniq.len = 1
#       
#       while(perm.uniq.len == 1){
#         perm.dat[which(perm.dat$Subject == i),]$Group = rep(sample(c(0,1),1), sum(perm.dat$Subject == i))# ,  replace = TRUE) #rep(sample(1:2), each = time.point)
#         perm.uniq.len = length(unique(perm.dat$Group))
#       }
#     }
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
#     if(g.min > min(points) | g.max < max(points))
#     {
#       cat("Special Case: generated permutation is out of range \n")
#       assign(paste("Model", j, sep = "_"), NULL)
#     } 
#     else if (length(which(sum(g.0$Count) == 0 | sum(g.1$Count)==0)))
#     {
#       cat("Special Case: zero for all variable of one group \n")
#       assign(paste("Model", j, sep = "_"), NULL)
#     }
#     else
#     {
#       #cat("In else", "\n")
#       perm = curveFitting(formula, df = perm.dat, method = method, points)
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