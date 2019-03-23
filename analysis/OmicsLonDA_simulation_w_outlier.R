library(simstudy)
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
library(lubridate)

setwd("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")
#setwd("C:/Users/ametwall/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")

set.seed(72626)

source("OmicsLonDA/R/omicslonda.R")
source("OmicsLonDA/R/omicslondaHelper.R")
source("OmicsLonDA/R/omicslondaMCPermutation.R")
source("OmicsLonDA/R/omicslondaVisualization.R")


##############################################################################################
### Simulate inconsistent, variable interval, correlated, longitudinal dataset
##############################################################################################
sim_portion = function(nsubj = 10, corcoef = 0.4, pattern.formula = "10")
{
  #nsubj = 20
  #corcoef = 0.4
  
  sim_all = data.frame()
  for (i in 1:nsubj)
  {
    baseline = toString(rnorm(1,10,10))
    print(baseline)
    # long.formula = sprintf("%s*(time<=50) + (%s+(time-50))*(time>50&time<=100) + (%s + (150-time))*(time>100&time<=150) + %s*(time>150)",
    #                   baseline, baseline, baseline, baseline)
    
    if(pattern.formula == "%s"){
      long.formula = sprintf(pattern.formula, baseline)
      print(long.formula)
    } else{
      if(i<4)
      {
        long.formula = sprintf(pattern.formula,
                               baseline, baseline, baseline, baseline)
        print(long.formula)
      } else {
        long.formula = sprintf("%s", baseline)
      }
      
    }
    
    #diff_1 = sim_portion(nsubj = 20, corcoef = 0.4, long.formula = pattern)
    
    # TODO: do for loop of 20 subjects. Make formula =1 instead of 20, then concatenate 
    # all of the measurements and rename the id to have ordered numbers
    
    # Define data for step 1
    def1 <- defData(varname = "xbase", dist = "normal", formula = 20, variance = 3)
    def1 <- defData(def1, varname = "nCount", dist = "noZeroPoisson", formula = 20)
    def1 <- defData(def1, varname = "mInterval", dist = "gamma", formula = 15, variance = 0.1)
    def1 <- defData(def1, varname = "vInterval", dist = "nonrandom", formula = 0.4)
    
    def2 <- defDataAdd(varname = "mu", dist = "nonrandom",
                       formula = long.formula)
    def2 <- defDataAdd(def2, varname = "var", dist = "nonrandom", formula = 5)
    
    ## Step 1: 
    dt <- genData(1, def1)
    dtPeriod <- addPeriods(dt)
    
    
    ### Step 2: Create dataset with consistency
    maxPeriod <- dtPeriod[, max(time)]
    dtx <- genData(1)
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
    
    simdf$id = i
    sim_all = rbind(sim_all, simdf)
  }
  
  ## Check correlation
  #var(as.matrix(dcast(ss, id ~ time, value.var = "ZZ")[, -1]))[1:10, 1:10]
  #acf(dtx$ZZ)
  return(simdf = sim_all)
}


##### Integrate Datasets
integrate = function(finalDT, finalDT_nd)
{
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
  
  ## convert to dataframe
  simdf_subset = as.data.frame(simdf)
  ## Remove duplicates
  simdf_subset = unique(simdf_subset)
  # Remove NA entries
  simdf_subset  = simdf_subset[-which(is.na(simdf_subset$Count)),]
  # remove samples at zero timepoints
  simdf_subset = simdf_subset[-which(simdf_subset$Time==0),]
  # TODO: remove samples with Count < 0 (check this simulation setting)
  if(length(which(simdf_subset$Count <= 0))){
    simdf_subset = simdf_subset[-which(simdf_subset$Count <= 0),]
  }
  return(simdf_subset)
}


## Normalization
clr_normalize = function(df){
  subjects = unique(df$Subject)
  
  ### Add psuedo count
  df$Count = df$Count + 1
  df$rawCount = df$Count
  df$normalizedCount = 0
  
  for(subj in subjects)
  {
    m = min(df[df$Subject==subj, "Time"])
    df[df$Subject==subj, ]$normalizedCount = log(df[df$Subject==subj, ]$Count/df[df$Subject==subj & df$Time==m, ]$Count)
  }
  
  return(df)
}



diff_normalize = function(df){
  subjects = unique(df$Subject)
  df$Count = df$Count
  df$rawCount = df$Count
  df$normalizedCount = 0
  
  for(subj in subjects)
  {
    m = min(df[df$Subject==subj, "Time"])
    df[df$Subject==subj, ]$normalizedCount = df[df$Subject==subj, ]$Count - df[df$Subject==subj & df$Time==m, ]$Count
  }
  
  return(df)
}




#### Simulate n features 
n = 10
earlydiff_simulatedDataset = list()
middlediff_simulatedDataset = list()
latediff_simulatedDataset = list()
verylatediff_simulatedDataset = list()
nondiff_simulatedDataset = list()

for(i in 1:n)
{
  ## Simulate n early differentially abundant/expressed features
  pattern = "%s*(time<=50) + (%s+(time-50))*(time>50&time<=100) + (%s + (150-time))*(time>100&time<=150) + %s*(time>150)"
  print(pattern)
  diff_1 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = pattern)
  diff_2 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = "%s")
  diff.df = integrate(diff_1, diff_2)
  earlydiff_simulatedDataset[[i]] = diff.df
}
{
  ## Simulate n middle differentially abundant/expressed features
  pattern = "%s*(time<=100) + (%s+(time-100))*(time>100&time<=150) + (%s + (200-time))*(time>150&time<=200) + %s*(time>200)"
  print(pattern)
  middlediff_1 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = pattern)
  middlediff_2 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = "%s")
  middlediff.df = integrate(middlediff_1, middlediff_2)
  middlediff_simulatedDataset[[i]] = middlediff.df
  
  
  # Simulate n late differentially abundant/expressed features
  pattern = "%s*(time<=150) + (%s+(time-150))*(time>150&time<=200) + (%s + (250-time))*(time>200&time<=250) + %s*(time>250)"
  print(pattern)
  latediff_1 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = pattern)
  latediff_2 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = "%s")
  latediff.df = integrate(latediff_1, latediff_2)
  latediff_simulatedDataset[[i]] = latediff.df
  
  
  # Simulate n very late differentially abundant/expressed features
  pattern = "%s*(time<=200) + (%s+(time-200))*(time>200&time<=250) + (%s + (300-time))*(time>250&time<=300) + %s*(time>300)"
  print(pattern)
  verylatediff_1 = diff_1 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = pattern)
  verylatediff_2 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = "%s")
  verylatediff.df = integrate(verylatediff_1, verylatediff_2)
  verylatediff_simulatedDataset[[i]] = verylatediff.df
  
  
  ## Simulate n NON differentially abundant/expressed features
  nondiff_1 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = "%s")
  nondiff_2 = sim_portion(nsubj = 20, corcoef = 0.4, pattern.formula = "%s")
  nondiff.df = integrate(nondiff_1, nondiff_2)
  nondiff_simulatedDataset[[i]] = nondiff.df
}




## Apply normalization
## TODO: waring is coming from duplicate records with the same minimum time. It seems that we simulate duplicate subjects
earlydiff_simulatedDataset_norm = list()
for(i in 1:length(earlydiff_simulatedDataset)){
  earlydiff_simulatedDataset_norm[[i]] = clr_normalize(earlydiff_simulatedDataset[[i]])
}


middlediff_simulatedDataset_norm = list()
for(i in 1:length(middlediff_simulatedDataset)){
  middlediff_simulatedDataset_norm[[i]] = clr_normalize(middlediff_simulatedDataset[[i]])
}


latediff_simulatedDataset_norm = list()
for(i in 1:length(latediff_simulatedDataset)){
  latediff_simulatedDataset_norm[[i]] = clr_normalize(latediff_simulatedDataset[[i]])
}


verylatediff_simulatedDataset_norm = list()
for(i in 1:length(verylatediff_simulatedDataset)){
  verylatediff_simulatedDataset_norm[[i]] = clr_normalize(verylatediff_simulatedDataset[[i]])
}


nondiff_simulatedDataset_norm = list()
for( i in 1:length(nondiff_simulatedDataset)){
  nondiff_simulatedDataset_norm[[i]] = clr_normalize(nondiff_simulatedDataset[[i]])
}





###############################################
######## Save simulated raw/normalized data
###############################################
save(earlydiff_simulatedDataset_norm, file = "simulatedDataset_earlydiff_OmicsLonDA_diffBaseline.RData")
save(middlediff_simulatedDataset_norm, file = "simulatedDataset_middlediff_OmicsLonDA_diffBaseline.RData")
save(latediff_simulatedDataset_norm, file = "simulatedDataset_latediff_OmicsLonDA_diffBaseline.RData")
save(verylatediff_simulatedDataset_norm, file = "simulatedDataset_verylatediff_OmicsLonDA_diffBaseline.RData")
save(nondiff_simulatedDataset_norm, file = "simulatedDataset_nondiff_OmicsLonDA_diffBaseline.RData")






###############################################
#### Visualize simulated raw/normalized data
###############################################

#### Visualize raw data
## Note: Warnings here are because we remove datapoints that have Time>500
jpeg("omicslonda_simulation_earlydiff_1_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(earlydiff_simulatedDataset[[1]], aes(Time, Count, colour = Group, group = interaction(Group, Subject)))
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


jpeg("omicslonda_simulation_middlediff_1_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
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


jpeg("omicslonda_simulation_latediff_1_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
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

jpeg("omicslonda_simulation_verylatediff_1_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
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


jpeg("omicslonda_simulation_nondiff_1_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
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





#### Visualize normalized data
jpeg("omicslonda_simulation_diff_1_normalized_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(earlydiff_simulatedDataset_norm[[1]], aes(Time, normalizedCount, colour = Group, group = interaction(Group, Subject)))
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


jpeg("omicslonda_simulation_middlediff_1_normalized_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
p = ggplot(middlediff_simulatedDataset_norm[[1]], aes(Time, normalizedCount, colour = Group, group = interaction(Group, Subject)))
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


jpeg("omicslonda_simulation_latediff_1_normalized_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
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


jpeg("omicslonda_simulation_verylatediff_1_normalized_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
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


jpeg("omicslonda_simulation_nondiff_1_normalized_diffBaseline.jpg", res = 300, height = 10, width = 15, units = 'cm')
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




###############################################
######## Load simulated raw/normalized data
###############################################
load("simulatedDataset_earlydiff_OmicsLonDA.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_middlediff_OmicsLonDA.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_latediff_OmicsLonDA.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_verylatediff_OmicsLonDA.RData", envir = parent.frame(), verbose = FALSE)
load("simulatedDataset_nondiff_OmicsLonDA.RData", envir = parent.frame(), verbose = FALSE)




###############################################
############# Global Testing
###############################################
### Test normalize and unnormalized, Time and Time:Group, for all features, for all patterns
globalTesting = function (df_normalized){
  sig_df_unnorm_TimeGroup = vector()
  sig_df_unnormTime = vector()
  
  sig_df_norm_TimeGroup = vector()
  sig_df_normTime = vector()
  
  for(i in 1:length(df_normalized))
  {
    unnormalized_lme <- lme(rawCount ~ Time * Group, random=~1|Subject, data = df_normalized[[i]])
    normalized_lme <- lme(normalizedCount ~ Time * Group, random=~1|Subject, data = df_normalized[[i]])
    a1_unnormalized_lme = anova(unnormalized_lme)
    a1_normalized_lme = anova(normalized_lme)
    
    if(a1_unnormalized_lme["Time:Group", "p-value"]<0.05){
      sig_df_unnorm_TimeGroup = c(sig_df_unnorm_TimeGroup, i)
    }
    if(a1_unnormalized_lme["Time", "p-value"]<0.05){
      sig_df_unnormTime = c(sig_df_unnormTime, i)
    }
    
    if(a1_normalized_lme["Time:Group", "p-value"]<0.05){
      sig_df_norm_TimeGroup = c(sig_df_norm_TimeGroup, i)
    }
    if(a1_normalized_lme["Time", "p-value"]<0.05){
      sig_df_normTime = c(sig_df_normTime, i)
    }
  }
  
  return(list(sig_df_unnorm_TimeGroup = length(sig_df_unnorm_TimeGroup), sig_df_unnormTime = length(sig_df_unnormTime), 
              sig_df_norm_TimeGroup = length(sig_df_norm_TimeGroup), sig_df_normTime = length(sig_df_normTime)))
}


globalTesting(df_normalized = earlydiff_simulatedDataset_norm)
## Sumulation of middle difference is distinct between the two groups
globalTesting(df_normalized = middlediff_simulatedDataset_norm) 
globalTesting(df_normalized = latediff_simulatedDataset_norm)
globalTesting(df_normalized = verylatediff_simulatedDataset_norm)
globalTesting(df_normalized = nondiff_simulatedDataset_norm)





###### Testing the interaction term only
globalTesting_interaction = function (df_normalized){
  sig_df_unnorm_TimeGroup = vector()
  sig_df_norm_TimeGroup = vector()
  
  for(i in 1:length(df_normalized))
  {
    unnormalized_lme <- lme(rawCount ~ Time:Group, random=~1|Subject, data = df_normalized[[i]])
    normalized_lme <- lme(normalizedCount ~ Time:Group, random=~1|Subject, data = df_normalized[[i]])
    a1_unnormalized_lme = anova(unnormalized_lme)
    a1_normalized_lme = anova(normalized_lme)
    
    if(a1_unnormalized_lme["Time:Group", "p-value"]<0.05){
      sig_df_unnorm_TimeGroup = c(sig_df_unnorm_TimeGroup, i)
    }
    if(a1_normalized_lme["Time:Group", "p-value"]<0.05){
      sig_df_norm_TimeGroup = c(sig_df_norm_TimeGroup, i)
    }
  }
  
  return(list(sig_df_unnorm_TimeGroup = length(sig_df_unnorm_TimeGroup), 
              sig_df_norm_TimeGroup = length(sig_df_norm_TimeGroup)))
}


 
globalTesting_interaction(df_normalized = earlydiff_simulatedDataset_norm)
## Sumulation of middle difference is distinct between the two groups
globalTesting_interaction(df_normalized = middlediff_simulatedDataset_norm) 
globalTesting_interaction(df_normalized = latediff_simulatedDataset_norm)
globalTesting_interaction(df_normalized = verylatediff_simulatedDataset_norm)
globalTesting_interaction(df_normalized = nondiff_simulatedDataset_norm)





###############################################
#############   Test OmicsLonDA
###############################################

## Test OmicsLonDA on one feature
points = seq(1, 500, length.out = 500)
output.omicslonda_earlydiff_1 = omicslonda(formula = normalizedCount ~ Time, df = earlydiff_simulatedDataset_norm[[1]], n.perm = 500, 
                                      fit.method = "ssgaussian", points = points,
                                      text = "sim_f1", parall = FALSE, pvalue.threshold = 0.05,
                                      adjust.method = "BH", col = c("blue", "green"),
                                      prefix = "OmicsLonDA_clr_earlydiff_f1_new", ylabel = "CLR-NormalizedCount",
                                      DrawTestStatDist = FALSE, time.unit = "days")





###############################################
#############   Evaluate OmicsLonDA
###############################################
automateEval_OmicsLonDA = function(data = data, n.perm = 100, points = points, pvalue_threshold = 0.05, prefix = "MMM",
                                   pattern = "NNN")
{
  omicslondaResults = list()
  for(i in 1:length(data))
  {
    omicslondaResults[[i]] = omicslonda(formula = normalizedCount ~ Time, df = data[[i]], n.perm = n.perm, fit.method = "ssgaussian", points = points,
                                        text = paste("Sim_F_",i, sep = "" ), parall = FALSE, pvalue.threshold = pvalue_threshold,     
                                        adjust.method = "BH", col = c("blue", "green"), 
                                        prefix = paste(prefix, "_Test_F", i, sep = ""), ylabel = "NormalizedCounts", 
                                        DrawTestStatDist = FALSE, time.unit = "days")
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
  confusion_omicslonda = data.frame()
  
  for(i in 1:length(omicslonda_sigificance)){
    len = length(omicslonda_sigificance[[i]])
    min.point = min(omicslondaResults[[i]]$detailed$points)
    max.point = max(omicslondaResults[[i]]$detailed$points)
    
    if(pattern == "earlydiff")
    {
      DA_truth = c(rep(0, 50-min.point), rep(1,100), rep(0,max.point - 150))
    } else if (pattern == "middlediff"){
      DA_truth = c(rep(0, 100-min.point), rep(1,100), rep(0,max.point - 200))
    } else if (pattern == "latediff"){
      if(len<250){
        print("in <250")
        DA_truth = c(rep(0, 150-min.point), rep(1,max.point-150))
      } else {
        print("in else")
        DA_truth = c(rep(0, 150-min.point), rep(1,100), rep(0,max.point - 250))
      }
    } else if (pattern == "verylatediff"){
      if(len<300){
        print("in < 300 ")
        DA_truth = c(rep(0, 200-min.point), rep(1,max.point-200)) 
      } else{
        print("in else 300")
        DA_truth = c(rep(0, 200-min.point), rep(1,100), rep(0,max.point - 300))
      }
    } else if(pattern == "nondiff"){
      DA_truth = c(rep(0, len))
    } else{
      print("wrong choice of simulation pattern")
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



######### Test Various patterns
points = seq(1, 500, length.out = 500)
## TODO: Specify normalizedCount, numbe of permutation, for the 1000 features


## TODO: solve warnings: 
# 1: In chol.default(wk1, pivot = TRUE) :
#   the matrix is either rank-deficient or indefinite
# 2: In (function (..., deparse.level = 1)  :
#          number of rows of result is not a multiple of vector length (arg 1)
evaluation_summary_earlydiff_norm = automateEval_OmicsLonDA(data = earlydiff_simulatedDataset_norm, n.perm = 100, points = points,
                                                           pvalue_threshold = 0.05, prefix = "evaluate_omicslonda_earlydiff_normalized_diffBaseline_outlier",
                                                           pattern = "earlydiff")

evaluation_summary_middlediff_norm = automateEval_OmicsLonDA(data = middlediff_simulatedDataset_norm[1:2], n.perm = 100, points = points,
                                                             pvalue_threshold = 0.05, prefix = "middlediff_omicslonda_normalized_1_2_diffBaseline",
                                                             pattern = "middlediff")

evaluation_summary_latediff_norm = automateEval_OmicsLonDA(data = latediff_simulatedDataset_norm[1:2], n.perm = 100, points = points,
                                                                   pvalue_threshold = 0.05, prefix = "evaluate_omicslonda_latediff_normalized_1_2_diffBaseline",
                                                                   pattern = "latediff")


evaluation_summary_verylatediff_norm = automateEval_OmicsLonDA(data = verylatediff_simulatedDataset_norm[1:2], n.perm = 100, points = points,
                                                                   pvalue_threshold = 0.05, prefix = "evaluate_omicslonda_verylatediff_normalized_1_2_diffBaseline",
                                                                   pattern = "verylatediff")

evaluation_summary_nondiff_norm = automateEval_OmicsLonDA(data = nondiff_simulatedDataset_norm[1:2], n.perm = 100, points = points,
                                                                   pvalue_threshold = 0.05, prefix = "evaluate_omicslonda_nondiff_normalized_diffBaseline",
                                                                   pattern = "nondiff")




###############################################
########## MetaLonDA
###############################################
#library(MetaLonDA)



source("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/MetaLonDA_dev/v1.1.5/MetaLonDA/R/Metalonda.R")
source("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/MetaLonDA_dev/v1.1.5/MetaLonDA/R/Visualization.R")
source("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/MetaLonDA_dev/v1.1.5/MetaLonDA/R/Permutation.R")
source("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/MetaLonDA_dev/v1.1.5/MetaLonDA/R/CurveFitting.R")
source("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/MetaLonDA_dev/v1.1.5/MetaLonDA/R/Normalization.R")


automateEval_MetaLonDA = function(data = data, n.perm = 100, points = points, pvalue_threshold = 0.05, prefix = "MMM",
                                   pattern = "NNN", method = "nbinomial")
{
  metalondaResults = list()
  for(i in 1:length(data))
  {
    
    metalondaResults[[i]] = metalonda(Count = data[[1]]$Count, 
                                    Time = data[[1]]$Time, 
                                    Group = data[[1]]$Group,
                                    ID = data[[1]]$Subject, 
                                    n.perm = n.perm, fit.method = method, points = points,
                                    text = paste(prefix, "_Test_F", i, sep = ""), parall = FALSE, 
                                    pvalue.threshold = pvalue_threshold,     
                                    adjust.method = "BH", col = c("black", "green"), prefix = prefix)
  }
  
  
  metalonda_adjustedPvalue = lapply(metalondaResults, function(x) x$detailed$adjusted.pvalue)
  
  x = metalonda_adjustedPvalue
  for(i in 1:length(x))
  {
    tmp = x[[i]]
    x[[i]][which(tmp > pvalue_threshold)] = 0
    x[[i]][which(tmp <= pvalue_threshold)] = 1
  }
  metalonda_sigificance =  x
  confusion_metalonda = data.frame()
  
  for(i in 1:length(metalonda_sigificance)){
    len = length(metalonda_sigificance[[i]])
    min.point = min(metalondaResults[[i]]$detailed$points)
    max.point = max(metalondaResults[[i]]$detailed$points)
    
    if(pattern == "earlydiff")
    {
      DA_truth = c(rep(0, 50-min.point), rep(1,100), rep(0,max.point - 150))
    } else if (pattern == "middlediff"){
      DA_truth = c(rep(0, 100-min.point), rep(1,100), rep(0,max.point - 200))
    } else if (pattern == "latediff"){
      if(len<250){
        print("in <250")
        DA_truth = c(rep(0, 150-min.point), rep(1,max.point-150))
      } else {
        print("in else")
        DA_truth = c(rep(0, 150-min.point), rep(1,100), rep(0,max.point - 250))
      }
    } else if (pattern == "verylatediff"){
      if(len<300){
        print("in < 300 ")
        DA_truth = c(rep(0, 200-min.point), rep(1,max.point-200)) 
      } else{
        print("in else 300")
        DA_truth = c(rep(0, 200-min.point), rep(1,100), rep(0,max.point - 300))
      }
    } else if(pattern == "nondiff"){
      DA_truth = c(rep(0, len))
    } else{
      print("wrong choice of simulation pattern")
      stop(message)
    }
    
    
    conf = evaluation(metalonda_sigificance[[i]], DA_truth)
    confusion_metalonda = rbind(confusion_metalonda, conf)
  }
  
  write.csv(confusion_metalonda, file = paste(prefix, "_confusion_metalonda", ".csv", sep=""), row.names=FALSE)
  return(confusion_metalonda=confusion_metalonda)
}



points = seq(1, 500, length.out = 500)
evaluation_summary_earlydiff_norm_nbinomial = automateEval_MetaLonDA(data = earlydiff_simulatedDataset_norm, n.perm = 100, points = points,
                                                           pvalue_threshold = 0.05, prefix = "evaluate_metalonda_earlydiff_normalized_diffBaseline_outlier_nbinomial",
                                                           pattern = "earlydiff", method = "nbinomial")

evaluation_summary_earlydiff_norm_lowess = automateEval_MetaLonDA(data = earlydiff_simulatedDataset_norm, n.perm = 100, points = points,
                                                           pvalue_threshold = 0.05, prefix = "evaluate_metalonda_earlydiff_normalized_diffBaseline_outlier_lowess",
                                                           pattern = "earlydiff", method = "lowess")




evaluation_summary_middlediff_norm = automateEval_MetaLonDA(data = middlediff_simulatedDataset_norm[1:2], n.perm = 100, points = points,
                                                           pvalue_threshold = 0.05, prefix = "evaluate_metalonda_middlediff_normalized_diffBaseline",
                                                           pattern = "earlydiff")

evaluation_summary_latediff_norm = automateEval_MetaLonDA(data = latediff_simulatedDataset_norm[1:2], n.perm = 100, points = points,
                                                           pvalue_threshold = 0.05, prefix = "evaluate_metalonda_latediff_normalized_diffBaseline",
                                                           pattern = "earlydiff")

evaluation_summary_verylatediff_norm = automateEval_MetaLonDA(data = verylatediff_simulatedDataset_norm[1:2], n.perm = 100, points = points,
                                                           pvalue_threshold = 0.05, prefix = "evaluate_metalonda_verylatediff_normalized_diffBaseline",
                                                           pattern = "earlydiff")

evaluation_summary_nondiff_norm = automateEval_MetaLonDA(data = nondiff_simulatedDataset_norm[1:2], n.perm = 100, points = points,
                                                           pvalue_threshold = 0.05, prefix = "evaluate_metalonda_nondiff_normalized_diffBaseline",
                                                           pattern = "earlydiff")


# View(earlydiff_simulatedDataset_norm[[1]])
# output.metalonda.f5 = metalonda(Count = earlydiff_simulatedDataset_norm[[1]]$Count, 
#                                 Time = earlydiff_simulatedDataset_norm[[1]]$Time, 
#                                 Group = earlydiff_simulatedDataset_norm[[1]]$Group,
#                                 ID = earlydiff_simulatedDataset_norm[[1]]$Subject, 
#                                 n.perm = 20, fit.method = "nbinomial", points = points,
#                                 text = "TEST_simF_5", parall = FALSE, pvalue.threshold = 0.05,     
#                                 adjust.method = "BH", col = c("black", "green"))
# 
# ### Prepare count table
# count_table = data.frame()
# for(i in 1:length(earlydiff_simulatedDataset_norm)){
#   count_table=rbind(count_table, earlydiff_simulatedDataset_norm[[i]]$Count)
# }
# colnames(count_table) = 1:ncol(count_table)
# 
# output.metalonda.all = metalondaAll(Count = count_table[1:2,], Time = earlydiff_simulatedDataset_norm[[1]]$Time, 
#                                     Group = earlydiff_simulatedDataset_norm[[1]]$Group,
#                                     ID = earlydiff_simulatedDataset_norm[[1]]$Subject, 
#                                     n.perm = 10, fit.method = "nbinomial", num.intervals = 499, 
#                                     parall = FALSE, pvalue.threshold = 0.05, adjust.method = "BH", time.unit = "hours", 
#                                     norm.method = "none", prefix = "Test_Metalonda_v1.1.3", ylabel = "Read Counts", col = c("black","green"))



###############################################
########## Test KS for studentization
###############################################
# ## Test normality of the normalized counts
# data = earlydiff_simulatedDataset_norm
# ks_pvalue = numeric(length(data))
# for(i in 1:length(data))
# {
#   hist(data[[i]]$normalizedCount)
#   #res = ks.test(data[[i]]$Count, "pnorm", mean(data[[i]]$Count), sd(data[[i]]$Count))
#   res = ks.test.t(data[[i]]$normalizedCount)
#   ks_pvalue[i] = res$p.value
# }
# z = p.adjust(ks_pvalue, method = "BH")
# sum(z<0.05)





# ## TODO: KS test for Student t-distribution
# library(LambertW)
# ks.test.t(x)
# 
# set.seed(1021)
# beta.true <- c(location = 0, scale = 1, df = 4)
# xx <- rt(n = 1000, df = beta.true['df'])
# ks.test.t(xx)
# ks.test.t(xx, beta.true)



###############################################
### Autocorrelation/spectral analysis figures 
###############################################


require(zoo)
ts <- as.zoo(c(1,2,8,1,2,2,3,2,3,2,2,1,3,2,3,1,1,2))
n <- 3
ts_n <- lag(ts, k=-n, na.pad=T)
#then compute correlation while omiting NAs

cor(ts[!is.na(ts_n)], ts_n[!is.na(ts_n)])


acf(ts)

# x = earlydiff_simulatedDataset_norm[[1]]
# s1 = x[x$Subject == "SA_1",]
# ar(s1$Count, aic = TRUE, order.max = NULL,
#    method = c( "mle"))
# ar(lh)
# ar(lh, method = "burg")
