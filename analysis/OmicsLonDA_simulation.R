# install.packages("mvnfast")
# install.packages("simstudy")



library("simstudy")
library("ggplot2")

library(nlme)
library(ggplot2)
library(gss)
library(plyr)
library(pracma)
library(parallel)
library(doParallel)
library(zoo)
library(pheatmap)
library(viridis)

setwd("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")
#setwd("C:/Users/ametwall/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")


source("OmicsLonDA/R/OmicsLonDA.R")
source("OmicsLonDA/R/CurveFitting.R")
source("OmicsLonDA/R/Visualization.R")
source("OmicsLonDA/R/Permutation.R")
source("OmicsLonDA/R/Normalization.R")
source("OmicsLonDA/R/OmicsLonDA_Evaluation.R")

# https://github.com/cran/simstudy/blob/master/inst/doc/simstudy.Rmd



set.seed(72626)

sim_portion = function(nsubj = 10, corcoef = 0.4, long.formula = "10")
{
  ### Inconsistent, variable interval, correlated, longitudinal dataset
  
  # Define data for step 1
  def1 <- defData(varname = "xbase", dist = "normal", formula = 20, variance = 3)
  def1 <- defData(def1, varname = "nCount", dist = "noZeroPoisson", formula = 20)
  def1 <- defData(def1, varname = "mInterval", dist = "gamma", formula = 15, variance = 0.1)
  def1 <- defData(def1, varname = "vInterval", dist = "nonrandom", formula = 0.4)
  
  def2 <- defDataAdd(varname = "mu", dist = "nonrandom",
                     formula = long.formula)
  ## formula = "10*(time<=50) + (10+(time-50))*(time>50&time<=100) + (10 + (150-time))*(time>100&time<=150) + 10*(time>150)"
  def2 <- defDataAdd(def2, varname = "var", dist = "nonrandom", formula = 5)
  
  ## Step 1: 
  #set.seed(3827367)
  dt <- genData(nsubj, def1)
  dtPeriod <- addPeriods(dt)
  
  
  ### Step 2: Create dataset with consistency
  maxPeriod <- dtPeriod[, max(time)]
  dtx <- genData(nsubj)
  dtx <- addPeriods(dtx, nPeriod = maxPeriod) 
  setnames(dtx, "period", "time")
  dtx <- addColumns(def2, dtx)
  dtx <- addCorGen(dtOld = dtx, idvar = "id", nvars = maxPeriod, rho = corcoef, 
                   corstr = "ar1", dist = "normal", param1 = "mu", param2 = "var",
                   cnames = "ZZ")
  
  ## Step 3: merge the two datasets to get the correlated inconsisted sampling
  setkey(dtx, id, time)
  setkey(dtPeriod, id, time)
  simdf <- mergeData(dtx, dtPeriod, idvars = c("id", "time"))
  
  
  ## Check correlation
  #var(as.matrix(dcast(dtx, id ~ time, value.var = "ZZ")[, -1]))[1:10, 1:10]
  
  
  # #### Visualize
  # ggplot(data = simdf, aes(x = time, y = ZZ, group=id)) +
  #   geom_point(aes(color = factor(id))) +
  #   geom_line(aes(color = factor(id))) +
  #   xlab("Day") +
  #   theme(legend.position = "none") #+
  
  return(simdf = simdf)
}



integrate = function(finalDT, finalDT_nd)
{
  ##### Integrate Datasets
  finalDT$Group = "Disease"
  finalDT_nd$Group = "Healthy"
  
  finalDT$id = paste("SA_", finalDT$id, sep = "")
  finalDT_nd$id = paste("SB_", finalDT_nd$id, sep = "")
  
  df = rbind(finalDT, finalDT_nd)
  names(df)
  simdf = df[, c("id", "time", "ZZ", "Group")]
  names(simdf) = c("Subject", "Time", "Count", "Group")
  
  simdf$Subject = as.factor(simdf$Subject)
  simdf$Group = as.factor(simdf$Group)
  
  #simdf_subset =  simdf[sample(nrow(simdf), 100), ]
  
  ## convert to dataframe
  simdf_subset = as.data.frame(simdf)
  ## Remove duplicates
  simdf_subset = unique(simdf_subset)
  # Remove NA entries
  simdf_subset  = simdf_subset[-which(is.na(simdf_subset$Count)),]
  # remove samples at zero timepoints
  simdf_subset = simdf_subset[-which(simdf_subset$Time==0),]
  # remove samples with Count < 0
  if(length(which(simdf_subset$Count <= 0))){
    simdf_subset = simdf_subset[-which(simdf_subset$Count <= 0),]
  }

  return(simdf_subset)
}


# ##### Simulate differentially abundant/expressed dataset
# diff_1 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10*(time<=50) + (10+(time-50))*(time>50&time<=100) + (10 + (150-time))*(time>100&time<=150) + 10*(time>150)")
# diff_2 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10")
# diff.df = integrate(diff_1, diff_2)
# 
# ggplot(data = diff.df, aes(x = Time, y = Count, group = Subject)) +
#   geom_point(aes(color = factor(Subject))) +
#   geom_line(aes(color = factor(Subject))) +
#   xlab("Day") +
#   theme(legend.position = "none")
# 
# 
# 
# ##### Simulate middle differentially abundant/expressed dataset
# diff_1 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10*(time<=100) + (10+(time-100))*(time>100&time<=150) + (10 + (200-time))*(time>150&time<=200) + 10*(time>200)")
# diff_2 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10")
# diff.df = integrate(diff_1, diff_2)
# 
# ggplot(data = diff.df, aes(x = Time, y = Count, group = Subject)) +
#   geom_point(aes(color = factor(Subject))) +
#   geom_line(aes(color = factor(Subject))) +
#   xlab("Day") +
#   theme(legend.position = "none")
# 
# 
# ##### Simulate late differentially abundant/expressed dataset
# diff_1 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10*(time<=150) + (10+(time-150))*(time>150&time<=200) + (10 + (250-time))*(time>200&time<=250) + 10*(time>250)")
# diff_2 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10")
# diff.df = integrate(diff_1, diff_2)
# 
# ggplot(data = diff.df, aes(x = Time, y = Count, group = Subject)) +
#   geom_point(aes(color = factor(Subject))) +
#   geom_line(aes(color = factor(Subject))) +
#   xlab("Day") +
#   theme(legend.position = "none")
# 
# 
# 
# ##### Simulate very late differentially abundant/expressed dataset
# diff_1 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10*(time<=200) + (10+(time-200))*(time>200&time<=250) + (10 + (300-time))*(time>250&time<=300) + 10*(time>300)")
# diff_2 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10")
# diff.df = integrate(diff_1, diff_2)
# 
# ggplot(data = diff.df, aes(x = Time, y = Count, group = Subject)) +
#   geom_point(aes(color = factor(Subject))) +
#   geom_line(aes(color = factor(Subject))) +
#   xlab("Day") +
#   theme(legend.position = "none")
# 
# 
# 
# ##### Simulate NON-differentially abundant/expressed dataset
# nondiff_1 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10")
# nondiff_2 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10")
# nondiff.df = integrate(nondiff_1, nondiff_2)
# 
# ggplot(data = nondiff.df, aes(x = Time, y = Count, group = Subject)) +
#   geom_point(aes(color = factor(Subject))) +
#   geom_line(aes(color = factor(Subject))) +
#   xlab("Day") +
#   theme(legend.position = "none")






#### Simulated 10 features 
n = 1000 ### Number of features
diff_simulatedDataset = list()
nondiff_simulatedDataset = list()
middlediff_simulatedDataset = list()
latediff_simulatedDataset = list()
verylatediff_simulatedDataset = list()

for(i in 1:n)
{
  ## Simulate n differentially abundant/expressed features
  diff_1 = sim_portion(nsubj = 20, corcoef = 0.4, long.formula = "10*(time<=50) + (10+(time-50))*(time>50&time<=100) + (10 + (150-time))*(time>100&time<=150) + 10*(time>150)")
  diff_2 = sim_portion(nsubj = 20, corcoef = 0.4, long.formula = "10")
  diff.df = integrate(diff_1, diff_2)
  # diff.df = as.data.frame(diff.df)
  # 
  # # Remove NA entries
  # diff.df  = diff.df[-which(is.na(diff.df$Count)),]
  # # remove samples at zero timepoints
  # diff.df = diff.df[-which(diff.df$Time==0),]
  # # remove samples with Count < 0
  # if(length(which(diff.df$Count <= 0))){
  #   diff.df = diff.df[-which(diff.df$Count <= 0),]
  # }
  diff_simulatedDataset[[i]] = diff.df
  

  
  ## Simulate n NON differentially abundant/expressed features
  nondiff_1 = sim_portion(nsubj = 20, corcoef = 0.4, long.formula = "10")
  nondiff_2 = sim_portion(nsubj = 20, corcoef = 0.4, long.formula = "10")
  nondiff.df = integrate(nondiff_1, nondiff_2)
  # nondiff.df = as.data.frame(nondiff.df)
  # 
  # # Remove NA entries
  # nondiff.df  = nondiff.df[-which(is.na(nondiff.df$Count)),]
  # # remove samples at zero timepoints
  # nondiff.df = nondiff.df[-which(nondiff.df$Time==0),]
  # # remove samples with Count < 0
  # if(length(which(nondiff.df$Count <= 0))){
  #   nondiff.df = nondiff.df[-which(nondiff.df$Count <= 0),]
  # }
  nondiff_simulatedDataset[[i]] = nondiff.df
  
  
  
  ## Simulate n middle differentially abundant/expressed features
  middlediff_1 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10*(time<=100) + (10+(time-100))*(time>100&time<=150) + (10 + (200-time))*(time>150&time<=200) + 10*(time>200)")
  middlediff_2 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10")
  middlediff.df = integrate(middlediff_1, middlediff_2)
  # middlediff.df = as.data.frame(middlediff.df)
  # 
  # # Remove NA entries
  # middlediff.df  = middlediff.df[-which(is.na(middlediff.df$Count)),]
  # #remove samples at zero timepoints
  # middlediff.df = middlediff.df[-which(middlediff.df$Time==0),]
  # # remove samples with Count < 0
  # if(length(which(middlediff.df$Count <= 0))){
  #   middlediff.df = middlediff.df[-which(middlediff.df$Count <= 0),]
  # }
  middlediff_simulatedDataset[[i]] = middlediff.df
  
  
  # Simulate n late differentially abundant/expressed features
  latediff_1 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10*(time<=150) + (10+(time-150))*(time>150&time<=200) + (10 + (250-time))*(time>200&time<=250) + 10*(time>250)")
  latediff_2 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10")
  latediff.df = integrate(latediff_1, latediff_2)
  # latediff.df = as.data.frame(latediff.df)
  # 
  # # Remove NA entries
  # latediff.df  = latediff.df[-which(is.na(latediff.df$Count)),]
  # #remove samples at zero timepoints
  # latediff.df = latediff.df[-which(latediff.df$Time==0),]
  # # remove samples with Count < 0
  # if(length(which(latediff.df$Count <= 0))){
  #   latediff.df = latediff.df[-which(latediff.df$Count <= 0),]
  # }
  latediff_simulatedDataset[[i]] = latediff.df
  
  
  # Simulate n very late differentially abundant/expressed features
  verylatediff_1 = diff_1 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10*(time<=200) + (10+(time-200))*(time>200&time<=250) + (10 + (300-time))*(time>250&time<=300) + 10*(time>300)")
  verylatediff_2 = sim_portion(nsubj = 10, corcoef = 0.4, long.formula = "10")
  verylatediff.df = integrate(verylatediff_1, verylatediff_2)
  # verylatediff.df = as.data.frame(verylatediff.df)
  # 
  # # Remove NA entries
  # verylatediff.df  = verylatediff.df[-which(is.na(verylatediff.df$Count)),]
  # #remove samples at zero timepoints
  # verylatediff.df = verylatediff.df[-which(verylatediff.df$Time==0),]
  # # remove samples with Count < 0
  # if(length(which(verylatediff.df$Count <= 0))){
  #   verylatediff.df = verylatediff.df[-which(verylatediff.df$Count <= 0),]
  # }
  verylatediff_simulatedDataset[[i]] = verylatediff.df
}









save(diff_simulatedDataset, file = "simulatedDataset_diff_OmicsLonDA_optimized.RData")
save(nondiff_simulatedDataset, file = "simulatedDataset_nondiff_OmicsLonDA_optimized.RData")
save(middlediff_simulatedDataset, file = "simulatedDataset_middlediff_OmicsLonDA_optimized.RData")
save(latediff_simulatedDataset, file = "simulatedDataset_latediff_OmicsLonDA_optimized.RData")
save(verylatediff_simulatedDataset, file = "simulatedDataset_verylatediff_OmicsLonDA_optimized.RData")
load("simulatedDataset_diff_OmicsLonDA_optimized.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_nondiff_OmicsLonDA_optimized.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_middlediff_OmicsLonDA_optimized.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_latediff_OmicsLonDA_optimized.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_verylatediff_OmicsLonDA_optimized.RData", envir = parent.frame(), verbose = FALSE)


# 
# save(diff_simulatedDataset, file = "simulatedDataset_diff_OmicsLonDA_optimized_1000.RData")
# save(nondiff_simulatedDataset, file = "simulatedDataset_nondiff_OmicsLonDA_optimized_1000.RData")
# save(middlediff_simulatedDataset, file = "simulatedDataset_middlediff_OmicsLonDA_optimized_1000.RData")
# save(latediff_simulatedDataset, file = "simulatedDataset_latediff_OmicsLonDA_optimized_1000.RData")
# save(verylatediff_simulatedDataset, file = "simulatedDataset_verylatediff_OmicsLonDA_optimized_1000.RData")


#### Visualize raw data
## Note: Warnings here are because we remove datapoints that have Time>500
jpeg("omicslonda_simulation_diff_1.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(diff_simulatedDataset[[1]], aes(Time, Count, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = " Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) + #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()


jpeg("omicslonda_simulation_nondiff_1.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(nondiff_simulatedDataset[[1]], aes(Time, Count, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = " Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) + #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()



jpeg("omicslonda_simulation_middlediff_1.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(middlediff_simulatedDataset[[1]], aes(Time, Count, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = " Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) + #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()


jpeg("omicslonda_simulation_latediff_1.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(latediff_simulatedDataset[[1]], aes(Time, Count, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = " Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) + #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()

jpeg("omicslonda_simulation_verylatediff_1.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(verylatediff_simulatedDataset[[1]], aes(Time, Count, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = " Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) + #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()







## Normalization
normalize = function(df){
  ## Remove time point 0
  y = df #[-which(df$Time==0),]
  y$rawCount = df$Count
  subjects = unique(y$Subject)
  y$normalizedCount = 0
  for(subj in subjects)
  {
    m = min(y[y$Subject==subj, "Time"])
    y[y$Subject==subj, ]$normalizedCount = log(y[y$Subject==subj, ]$Count/y[y$Subject==subj & y$Time==m, ]$Count)
  }
  
  
  ### May add sudo code
  y$normalizedCount = y$normalizedCount + 0.000000001
  y$Count = y$normalizedCount
  return(y)
}

diff_simulatedDataset_norm = list()
for(i in 1:length(diff_simulatedDataset)){
  diff_simulatedDataset_norm[[i]] = normalize(diff_simulatedDataset[[i]])
}

nondiff_simulatedDataset_norm = list()
for( i in 1:length(nondiff_simulatedDataset)){
  nondiff_simulatedDataset_norm[[i]] = normalize(nondiff_simulatedDataset[[i]])
}


middlediff_simulatedDataset_norm = list()
for(i in 1:length(middlediff_simulatedDataset)){
  middlediff_simulatedDataset_norm[[i]] = normalize(middlediff_simulatedDataset[[i]])
}


latediff_simulatedDataset_norm = list()
for(i in 1:length(latediff_simulatedDataset)){
  latediff_simulatedDataset_norm[[i]] = normalize(latediff_simulatedDataset[[i]])
}

## Solved warning. It was because duplicates in simulation dataframe. message: In y[y$Subject == subj, ]$Count/y[y$Subject == subj & y$Time ==  : .... longer object length is not a multiple of shorter object length
verylatediff_simulatedDataset_norm = list()
for(i in 1:length(verylatediff_simulatedDataset)){
  verylatediff_simulatedDataset_norm[[i]] = normalize(verylatediff_simulatedDataset[[i]])
}







save(diff_simulatedDataset_norm, file = "simulatedDataset_diff_OmicsLonDA_optimized_normalized.RData")
save(nondiff_simulatedDataset_norm, file = "simulatedDataset_nondiff_OmicsLonDA_optimized_normalized.RData")
save(middlediff_simulatedDataset_norm, file = "simulatedDataset_middlediff_OmicsLonDA_optimized_normalized.RData")
save(latediff_simulatedDataset_norm, file = "simulatedDataset_latediff_OmicsLonDA_optimized_normalized.RData")
save(verylatediff_simulatedDataset_norm, file = "simulatedDataset_verylatediff_OmicsLonDA_optimized_normalized.RData")
load("simulatedDataset_diff_OmicsLonDA_optimized_normalized.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_nondiff_OmicsLonDA_optimized_normalized.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_middlediff_OmicsLonDA_optimized_normalized.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_latediff_OmicsLonDA_optimized_normalized.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_verylatediff_OmicsLonDA_optimized_normalized.RData", envir = parent.frame(), verbose = FALSE)


# save(diff_simulatedDataset_norm, file = "simulatedDataset_diff_OmicsLonDA_optimized_normalized_1000.RData")
# save(nondiff_simulatedDataset_norm, file = "simulatedDataset_nondiff_OmicsLonDA_optimized_normalized_1000.RData")
# save(middlediff_simulatedDataset_norm, file = "simulatedDataset_middlediff_OmicsLonDA_optimized_normalized_1000.RData")
# save(latediff_simulatedDataset_norm, file = "simulatedDataset_latediff_OmicsLonDA_optimized_normalized_1000.RData")
# save(verylatediff_simulatedDataset_norm, file = "simulatedDataset_verylatediff_OmicsLonDA_optimized_normalized_1000.RData")



load("simulatedDataset_diff_OmicsLonDA_optimized_normalized_1000.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_nondiff_OmicsLonDA_optimized_normalized_1000.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_middlediff_OmicsLonDA_optimized_normalized_1000.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_latediff_OmicsLonDA_optimized_normalized_1000.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_verylatediff_OmicsLonDA_optimized_normalized_1000.RData", envir = parent.frame(), verbose = FALSE)


#### Visualization
jpeg("omicslonda_simulation_diff_1_normalized.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(diff_simulatedDataset_norm[[1]], aes(Time, normalizedCount, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = "Adjusted Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) +  #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) +
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()



jpeg("omicslonda_simulation_nondiff_1_normalized.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(nondiff_simulatedDataset_norm[[1]], aes(Time, normalizedCount, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = "Adjusted Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) + #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) +
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()


jpeg("omicslonda_simulation_middlediff_2_normalized.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(middlediff_simulatedDataset_norm[[2]], aes(Time, normalizedCount, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = "Adjusted Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) +  #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) +
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()


jpeg("omicslonda_simulation_latediff_1_normalized.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(latediff_simulatedDataset_norm[[1]], aes(Time, normalizedCount, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = "Adjusted Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) +  #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) +
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()


jpeg("omicslonda_simulation_verylatediff_1_normalized.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(verylatediff_simulatedDataset_norm[[1]], aes(Time, normalizedCount, colour = Group, group = interaction(Group, Subject)))
p + geom_point(size = 1, alpha = 0.5) + geom_line(size = 1, alpha = 0.7) +  theme_bw() +
  #ggtitle(paste("Feature = ", text, sep = "")) + 
  labs(y = "Adjusted Molecule Level", x = "Time") +
  scale_colour_manual(values = c("blue", "green")) +  #, breaks = c("0", "1"),
  #                   labels = c(group.levels[1], group.levels[2])) +
  xlim(0,500) +
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=15, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=15, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="top") #+ scale_x_continuous(breaks = waiver())
dev.off()






############# Global Testing
library(nlme)

### Test on raw Counts
m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=diff_simulatedDataset[[1]])
anova(m1)

m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=nondiff_simulatedDataset[[1]])
anova(m1)

m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=middlediff_simulatedDataset[[1]])
anova(m1)

m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=latediff_simulatedDataset[[1]])
anova(m1)

m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=verylatediff_simulatedDataset[[1]])
anova(m1)


### Test on normalized counts
m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=diff_simulatedDataset_norm[[1]])
anova(m1)

m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=nondiff_simulatedDataset_norm[[1]])
anova(m1)

m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=middlediff_simulatedDataset_norm[[1]])
anova(m1)

m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=latediff_simulatedDataset_norm[[1]])
anova(m1)

m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=verylatediff_simulatedDataset_norm[[1]])
anova(m1)



### Test normalize and unnormalized, Time and Time:Group, for all features, for all patterns
globalTesting = function (df_normalized){
  sig_df_unnorm_TimeGroup = vector()
  sig_df_unnormTime = vector()
  
  sig_df_norm_TimeGroup = vector()
  sig_df_normTime = vector()
  
  for(i in 1:length(df_normalized))
  {
    #print(i)
    unnormalized_lme <- lme(rawCount ~ Time * Group, random=~1|Subject, data = df_normalized[[i]])
    normalized_lme <- lme(normalizedCount ~ Time * Group, random=~1|Subject, data = df_normalized[[i]])
    a1_unnormalized_lme = anova(unnormalized_lme)
    a1_normalized_lme = anova(normalized_lme)
    
    if(a1_unnormalized_lme["Time:Group", "p-value"]<0.05){
      sig_df_unnorm_TimeGroup = c(sig_df_unnorm_TimeGroup, i)
      #cat("Significant  feature = ", i, "\n")
    }
    if(a1_unnormalized_lme["Time", "p-value"]<0.05){
      sig_df_unnormTime = c(sig_df_unnormTime, i)
      #cat("Significant  feature = ", i, "\n")
    }
    
    if(a1_normalized_lme["Time:Group", "p-value"]<0.05){
      sig_df_norm_TimeGroup = c(sig_df_norm_TimeGroup, i)
      #cat("Significant  feature = ", i, "\n")
    }
    if(a1_normalized_lme["Time", "p-value"]<0.05){
      sig_df_normTime = c(sig_df_normTime, i)
      #cat("Significant  feature = ", i, "\n")
    }
  }
  
  return(list(sig_df_unnorm_TimeGroup = length(sig_df_unnorm_TimeGroup), sig_df_unnormTime = length(sig_df_unnormTime), 
              sig_df_norm_TimeGroup = length(sig_df_norm_TimeGroup), sig_df_normTime = length(sig_df_normTime)))
}

globalTesting(df_normalized = nondiff_simulatedDataset_norm)
globalTesting(df_normalized = diff_simulatedDataset_norm)
## Middle has issues 
globalTesting(df_normalized = middlediff_simulatedDataset_norm) 
globalTesting(df_normalized = latediff_simulatedDataset_norm)
globalTesting(df_normalized = verylatediff_simulatedDataset_norm)









###### Testing the interaction term only
globalTesting_interaction = function (df_normalized){
  sig_df_unnorm_TimeGroup = vector()
  #sig_df_unnormTime = vector()
  
  sig_df_norm_TimeGroup = vector()
  #sig_df_normTime = vector()
  
  for(i in 1:length(df_normalized))
  {
    #print(i)
    unnormalized_lme <- lme(rawCount ~ Time:Group, random=~1|Subject, data = df_normalized[[i]])
    normalized_lme <- lme(normalizedCount ~ Time:Group, random=~1|Subject, data = df_normalized[[i]])
    a1_unnormalized_lme = anova(unnormalized_lme)
    a1_normalized_lme = anova(normalized_lme)
    
    if(a1_unnormalized_lme["Time:Group", "p-value"]<0.05){
      sig_df_unnorm_TimeGroup = c(sig_df_unnorm_TimeGroup, i)
      #cat("Significant  feature = ", i, "\n")
    }
    # if(a1_unnormalized_lme["Time", "p-value"]<0.05){
    #   sig_df_unnormTime = c(sig_df_unnormTime, i)
    #   #cat("Significant  feature = ", i, "\n")
    # }
    # 
    if(a1_normalized_lme["Time:Group", "p-value"]<0.05){
      sig_df_norm_TimeGroup = c(sig_df_norm_TimeGroup, i)
      #cat("Significant  feature = ", i, "\n")
    }
    # if(a1_normalized_lme["Time", "p-value"]<0.05){
    #   sig_df_normTime = c(sig_df_normTime, i)
    #   #cat("Significant  feature = ", i, "\n")
    # }
  }
  
  return(list(sig_df_unnorm_TimeGroup = length(sig_df_unnorm_TimeGroup), 
              sig_df_norm_TimeGroup = length(sig_df_norm_TimeGroup)))
}

globalTesting_interaction(df_normalized = nondiff_simulatedDataset_norm)
globalTesting_interaction(df_normalized = diff_simulatedDataset_norm)
## Middle has issues 
globalTesting_interaction(df_normalized = middlediff_simulatedDataset_norm) 
globalTesting_interaction(df_normalized = latediff_simulatedDataset_norm)
globalTesting_interaction(df_normalized = verylatediff_simulatedDataset_norm)





###### OMICSLONDA
library(nlme)
library(ggplot2)
library(gss)
library(plyr)
library(pracma)
library(parallel)
library(doParallel)
library(zoo)
library(pheatmap)
library(viridis)


source("OmicsLonDA/R/OmicsLonDA.R")
source("OmicsLonDA/R/CurveFitting.R")
source("OmicsLonDA/R/Visualization.R")
source("OmicsLonDA/R/Permutation.R")
source("OmicsLonDA/R/Normalization.R")
source("OmicsLonDA/R/OmicsLonDA_Evaluation.R")


# ### Workout the evaluation dimension
points = seq(1, 500, length.out = 500)
evaluation_summary_earlydiff_norm_f1_f2 = automateEval_OmicsLonDA(data = diff_simulatedDataset_norm[1:100], n.perm = 100, points = points,
                                                           pvalue_threshold = 0.05, prefix = "Sunday_diff_normalized_1_2",
                                                           pattern = "earlydiff")

evaluation_summary_middlediff_norm_f1_f2 = automateEval_OmicsLonDA(data = middlediff_simulatedDataset_norm[1:100], n.perm = 100, points = points,
                                                             pvalue_threshold = 0.05, prefix = "Sunday_middlediff_normalized_1_2",
                                                             pattern = "middlediff")

evaluation_summary_latediff_norm_f1_f2 = automateEval_OmicsLonDA(data = latediff_simulatedDataset_norm[1:100], n.perm = 100, points = points,
                                                                   pvalue_threshold = 0.05, prefix = "Sunday_latediff_normalized_1_2",
                                                                   pattern = "latediff")


## TODO: solve warning: In (function (..., deparse.level = 1)  ... :
# number of rows of result is not a multiple of vector length (arg 1)
evaluation_summary_verylatediff_norm_f1_f2 = automateEval_OmicsLonDA(data = verylatediff_simulatedDataset_norm[1:100], n.perm = 100, points = points,
                                                                   pvalue_threshold = 0.05, prefix = "Sunday_verylatediff_normalized_1_2",
                                                                   pattern = "verylatediff")

evaluation_summary_nondiff_norm_f1_f2 = automateEval_OmicsLonDA(data = nondiff_simulatedDataset_norm[1:100], n.perm = 100, points = points,
                                                                   pvalue_threshold = 0.05, prefix = "Sunday_nondiff_normalized_1_2",
                                                                   pattern = "nondiff")







## Test omicslonda alone
points = seq(1, 500, length.out = 500)
output.omicslonda_diff_1 = omicslonda(formula = Count ~ Time, df = diff_simulatedDataset_norm[[1]], n.perm = 10, 
                                      fit.method = "ssgaussian", points = points,
                                      text = "simdf_1", parall = FALSE, pvalue.threshold = 0.05,
                                      adjust.method = "BH", col = c("blue", "green"),
                                      prefix = "TestOmicsLonDA_norm_stepbystep_f1_Sunday", ylabel = "NormalizedCounts",
                                      DrawTestStatDist = FALSE, time.unit = "days")


# points = seq(1, 200, length.out = 200)
# output.omicslonda_diff_4 = omicslonda(formula = Count ~ Time, df = diff_simulatedDataset[[4]], n.perm = 100, fit.method = "ssgaussian", points = points,
#                                       text = "simdf_1", parall = FALSE, pvalue.threshold = 0.05,
#                                       adjust.method = "BH", col = c("blue", "green"),
#                                       prefix = "TestOmicsLonDA_ssgaussian_simdatasets_NEWWWWWWWW_f4_night_new", ylabel = "NormalizedCounts",
#                                       DrawTestStatDist = FALSE)
# 
# 
# output.omicslonda_nondiff_1 = omicslonda(formula = Count ~ Time, df = nondiff_simulatedDataset[[1]], n.perm = 100, fit.method = "ssgaussian", points = points,
#                                          text = "simdf_1", parall = FALSE, pvalue.threshold = 0.05,
#                                          adjust.method = "BH", col = c("blue", "green"),
#                                          prefix = "TestOmicsLonDA_ssgaussian_simdatasets_nondiff_f1", ylabel = "NormalizedCounts",
#                                          DrawTestStatDist = FALSE)



# points = seq(1, 200, length.out = 100)
# output.omicslonda_diff_f1_norm = omicslonda(formula = normalizedCount ~ Time, df = diff_simulatedDataset_norm[[1]], 
#                                       n.perm = 100, fit.method = "ssgaussian", points = points,
#                                       text = "Simulated Omic Feature", parall = FALSE, pvalue.threshold = 0.05,
#                                       adjust.method = "BH", col = c("blue", "green"),
#                                       prefix = "TestOmicsLonDA_ssgaussian_simdatasets_f1_diff_normalized", 
#                                       ylabel = "Adjusted molecule level",
#                                       DrawTestStatDist = TRUE)




#### Test KS for studentization
normalize_natural = function(df){
  ## Remove time point 0
  y = df #[-which(df$Time==0),]
  subjects = unique(y$Subject)
  y$normalizedCount = 0
  for(subj in subjects)
  {
    m = min(y[y$Subject==subj, "Time"])
    y[y$Subject==subj, ]$normalizedCount = log(y[y$Subject==subj, ]$Count/y[y$Subject==subj & y$Time==m, ]$Count)
  }
  
  
  ### May add sudo code
  y$normalizedCount = y$normalizedCount + 0.000000001
  y$Count = y$normalizedCount
  return(y)
}



diff_simulatedDataset_norm_natural = list()
for(i in 1:length(diff_simulatedDataset)){
  diff_simulatedDataset_norm_natural[[i]] = normalize(diff_simulatedDataset[[i]])
}


## Test normality of the normalized counts
data = diff_simulatedDataset_norm_natural
ks_pvalue = numeric(length(diff_simulatedDataset_norm))
for(i in 1:length(data))
{
  hist(data[[i]]$Count)
  #res = ks.test(data[[i]]$Count, "pnorm", mean(data[[i]]$Count), sd(data[[i]]$Count))
  res = ks.test.t(data[[i]]$Count)
  ks_pvalue[i] = res$p.value
}
z = p.adjust(ks_pvalue, method = "BH")
sum(z<0.05)



set.seed(1021)
beta.true <- c(location = 0, scale = 1, df = 4)
xx <- rt(n = 1000, df = beta.true['df'])
ks.test.t(xx)
ks.test.t(xx, beta.true)



## Ks test for Student t-distribution
library(LambertW)
ks.test.t(x)




### TODO: autocorrelation figres and spectral analysis figure
x = diff_simulatedDataset_norm[[1]]

s1 = x[x$Subject == "SA_1",]

ar(s1$Count, aic = TRUE, order.max = NULL,
   method = c( "mle"))

ar(lh)
ar(lh, method = "burg")
