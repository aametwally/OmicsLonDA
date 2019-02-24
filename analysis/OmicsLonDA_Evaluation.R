automateEval_OmicsLonDA = function(data = data, n.perm = 100, points = points, pvalue_threshold = 0.05, prefix = "MMM",
                                   pattern = "nondifferential")
{
  
  ##############################
  ###### Test all features  ####
  ##############################
  
  omicslondaResults = list()
  for(i in 1:length(data))
  {
    omicslondaResults[[i]] = omicslonda(formula = Count ~ Time, df = data[[i]], n.perm = n.perm, fit.method = "ssgaussian", points = points,
                                 text = paste("simdF_",i, sep = "" ), parall = FALSE, pvalue.threshold = pvalue_threshold,     
                                 adjust.method = "BH", col = c("blue", "green"), 
                                 prefix = paste(prefix, "_Test_F", i, sep = ""), ylabel = "NormalizedCounts", 
                                 DrawTestStatDist = FALSE)
  }
  
  
  omicslonda_adjustedPvalue = lapply(omicslondaResults, function(x) x$detailed$adjusted.pvalue)
  
  x = omicslonda_adjustedPvalue
  for(i in 1:length(x))
  {
    tmp = x[[i]]
    x[[i]][which(tmp > (pvalue_threshold/2))] = 0
    x[[i]][which(tmp <= (pvalue_threshold/2))] = 1
  }
  omicslonda_sigificance =  x
  
  
  #DA_truth = matrix(rep(c(rep(0,49), rep(1,100), rep(0,50)), length(omicslonda_sigificance)), nrow = length(omicslonda_sigificance), byrow = TRUE)
  confusion_omicslonda = data.frame()
  
  for(i in 1:length(omicslonda_sigificance)){
    len = length(omicslonda_sigificance[[i]])
    min.point = min(omicslondaResults[[i]]$detailed$points)
    max.point = max(omicslondaResults[[i]]$detailed$points)
                 
    if(pattern == "earlydiff")
    {
      DA_truth = c(rep(0, 50-min.point), rep(1,100), rep(0,max.point - 150))
    } else if(pattern == "middlediff"){
      DA_truth = c(rep(0, 100-min.point), rep(1,100), rep(0,max.point - 200))
    } else if(pattern == "latediff"){
      if(len<250){
        print("in <250")
        DA_truth = c(rep(0, 150-min.point), rep(1,max.point-150))#, rep(0,max.point - 250))
      } else {
        print("in else")
        DA_truth = c(rep(0, 150-min.point), rep(1,100), rep(0,max.point - 250))
      }
    } else if(pattern == "verylatediff"){
      if(len<300){
        print("in < 300 ")
        DA_truth = c(rep(0, 200-min.point), rep(1,max.point-200))    #, rep(0,max.point - 300))
      } else{
        print("in else 300")
        DA_truth = c(rep(0, 200-min.point), rep(1,100), rep(0,max.point - 300))
      }
    } else if(pattern == "nondiff"){
      DA_truth = c(rep(0, len))# - 150), rep(1,100), rep(0,50))
    } else{
      print("wrong choice of patters")
      stop(message)
    }
      
    
    conf = evaluation(omicslonda_sigificance[[i]], DA_truth)
    confusion_omicslonda = rbind(confusion_omicslonda, conf)
    
  }

  write.csv(confusion_omicslonda, file = paste(prefix, "_confusion_omicslonda", ".csv", sep=""), row.names=FALSE)
  return(confusion_omicslonda=confusion_omicslonda)
}





evaluation = function(predicted, truth)
{
    cat("len of predicted = ", length(predicted), "\n")
    cat("len of truth = ", length(truth), "\n")
    print(truth)
    print(predicted)
    TP = sum(predicted*truth == 1)
    TN = sum((predicted + truth) == 0)
    FP = sum((predicted * (!truth)) == 1)
    FN = sum(((!predicted) * truth) == 1)
    
    con_vec = cbind(TP = TP, FP = FP, TN = TN, FN = FN)
  
  return(con_vec)
}
