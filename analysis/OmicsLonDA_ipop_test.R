library(zoo)
library(lubridate)
library(ggplot2)
library(gss)
library(plyr)
library(pracma)
library(parallel)
library(doParallel)
library("biomformat"); packageVersion("biomformat")
library(nlme)

rm(list=ls())
setwd("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")
load("ipop_data/Revision_MultiOmes.RData", envir = parent.frame(), verbose = FALSE)

source("R/OmicsLonDA.R")
source("R/CurveFitting.R")
source("R/Visualization.R")
source("R/Permutation.R")
source("R/Normalization.R")


## Define omicslondaAll_ipop method
omicslondaAll_ipop = function(formula = Count ~ Time, data, n.perm = 25, fit.method = "ssgaussian", 
                                   num.intervals = 100, parall = FALSE, pvalue.threshold = 0.05, 
                                   adjust.method = "BH", time.unit = "days", norm.method = "none", prefix = "TestOmicsLonDA",
                                   ylabel = "Normalized Count", col = c("blue", "firebrick"))
{
  ## Apply omicslonda for each feature
  n.features = length(data)
  detailed = list()
  summary = list()
  points = seq(1, 50, length.out = 50)
  
  
  
  ## Get Group name
  group.levels = sort(unique(data[[1]]$Group))
  if(length(group.levels) > 2){
    stop("You have more than two phenotypes.")
  }
  gr.1 = group.levels[1]
  gr.2 = group.levels[2]
  
  omicslondaResults = list()
  for (i in 1:n.features)
  {
    cat ("Feature  = ", data[[i]]$feature[1], "\n")
    omicslondaResults[[i]] = omicslonda(formula = Count ~ Time, df = data[[i]], n.perm = n.perm, fit.method = fit.method, points = points,
                                        text = unique(data[[i]]$feature), parall = parall, pvalue.threshold = pvalue.threshold,     
                                        adjust.method = adjust.method, col = col, 
                                        prefix = prefix, ylabel = ylabel, 
                                        DrawTestStatDist = FALSE)
    detailed[[i]] = omicslondaResults[[i]]$detailed
    summary[[i]] = omicslondaResults[[i]]$summary
  }
  
  summary.tmp = do.call(rbind, summary)
  summary.tmp$dominant[which(summary.tmp$dominant == 1)] = gr.1
  summary.tmp$dominant[which(summary.tmp$dominant == -1)] = gr.2
  
  ## Output table and figure that summarize the significant time intervals
  write.csv(summary.tmp, file = sprintf("%s/OmicsLonDA_TimeIntervals_%s_%s.csv", prefix, fit.method, prefix), row.names = FALSE)
  visualizeTimeIntervals(interval.details = summary.tmp, prefix, unit = time.unit, col = col, fit.method = fit.method)
  
  aggregateData = list(output.detail = detailed, output.summary = summary.tmp)
  save(aggregateData, file = sprintf("%s/OmicsLonDA_Summary_%s_%s.RData", prefix, fit.method, prefix))
  return(aggregateData)
}




#######################################
########### Metabolomics     ##########
#######################################
metabolites_seasonality_all = read.csv(file = "ipop_data/metabolites_seasonality_all.csv", row.names = 1, check.names=FALSE)

metabolites_count = metabolites_seasonality_all[-c(1:6),]
metabolites_count_table = data.frame(sapply(metabolites_count, function(x) as.numeric(as.character(x))))
rownames(metabolites_count_table) = rownames(metabolites_count)
colnames(metabolites_count_table) = colnames(metabolites_count)
metabolites_metadata_table = as.data.frame(t(metabolites_seasonality_all[c(1:6),]))
colnames(metabolites_metadata_table)[1] = "Subject"
metabolites_metadata_table$Date = as.Date(metabolites_metadata_table$Date, format = '%d-%b-%y')
metabolites_metadata_table$SampleID = rownames(metabolites_metadata_table)
metabolites_metadata_table  = merge(metabolites_metadata_table, ls, by = "SampleID")


prefix = "Infection"
xx = subset(metabolites_metadata_table, substring(state1, 1, nchar(prefix)) == prefix)
xx$Subject = factor(xx$Subject)

### TODO: Filter out subject with less than 5 infection time points. We may do the filteration based on the subject_cycle parameter later on

## Select samples from subjects with infection timepoints
subset.df = subset(metabolites_metadata_table, Subject %in% xx$Subject)
colnames(subset.df)[which(colnames(subset.df) == "CL1")] = "DaysAfterInfect"
subset.df$DaysAfterInfect = gsub("D", "", subset.df$DaysAfterInfect)

## Select Healthy and Infection time points
yy = subset(subset.df, substring(state1, 1, nchar("Infection")) == "Infection" | substring(state1, 1, nchar("Healthy")) == "Healthy")
yy = subset(yy, state1 != "Infection")


## for each datapoint with a healthy status, search for the nearest 
for(i in 1:nrow(yy)){
  if(yy[i,]$state1 == "Healthy"){
    tmp = subset(yy, Subject == yy[i,]$Subject)
    print(tmp)
    a = subset(tmp, Date>yy[i,]$Date & substring(tmp$state1, 1, nchar("Infection")) == "Infection")
    nearstTime = min(a$Date)
    yy[i,]$DaysAfterInfect = yy[i,]$Date - nearstTime
  }
}


## remove inf
remove = which(yy$DaysAfterInfect == "-Inf")
yy = yy[-remove,]


## name each infection cycle by the date of the first infection sample
yy$CycleName = ""### TODO: Name infection cycles
class(yy$CycleName) = "Date"
## Order yy
yy = yy[with(yy, order(Subject, Date)),]

write.csv(yy, file = "ipop_data/ipop_metabolites_infection_2.csv")


### TODO: rename the cycle of the block that don't have Healthy and consecutive to infection timepoints
### TODO: This is 2011-12-04
### TODO: sort df by date
tmp_date = yy[1,]$Date
tmp_subj = yy[1,]$Subject
for(i in 1:nrow(yy)){
  if(yy[i,]$state1 == "Healthy"){
    tmp = subset(yy, Subject == yy[i,]$Subject)
    #print(tmp)
    a = subset(tmp, Date>yy[i,]$Date & substring(tmp$state1, 1, nchar("Infection")) == "Infection")
    nearstTime = min(a$Date)
    yy[i,]$DaysAfterInfect = yy[i,]$Date - nearstTime
    yy[i,]$CycleName = nearstTime
    tmp_date = nearstTime
    tmp_subj = yy[i,]$Subject
  }
  else{
    tmp = subset(yy, Subject == yy[i,]$Subject)
    ## handle infection without prior healthy time points. actually we dont care. since those would be assigned NA in the CycleName by the followong line and we would remove them afterwards 
    a = subset(tmp, Date <= yy[i,]$Date & Date >= tmp_date & yy[i,]$Subject == tmp_subj & substring(tmp$state1, 1, nchar("Infection")) == "Infection")
    nearstTime = min(a$Date)
    print(as.character(nearstTime))
    yy[i,]$CycleName = nearstTime
  }
}


## TODO: remove infection samples without Healthy baseline
write.csv(yy, file = "ipop_data/ipop_metabolites_infection_3.csv")


### Read modified file
#samples_infection = read.csv("ipop_data/ipop_metabolites_infection_wCycles.csv", row.names = 1)
samples_infection = yy
samples_infection$Subject_Cycle = paste(samples_infection$Subject, as.character(samples_infection$CycleName), sep = "_")
class(samples_infection$DaysAfterInfect) = "numeric"

## Count number of unique infection cycle
length(unique(samples_infection$Subject_Cycle))


### TODO: Remove sampels without Healthy group
samples_infection_filtered = samples_infection
for(i in unique(samples_infection$Subject_Cycle)){
  block = samples_infection[which(samples_infection$Subject_Cycle == i), "state1"]
  if(!any(block == "Healthy")){
    samples_infection_filtered = samples_infection_filtered[-which(samples_infection_filtered$Subject_Cycle == i), ]
  }
}



## histiogram of BMI and age 
hist(sc$BMI)
hist(sc$Adj.age)

### Merge Gender with the df
subject_info = sc[, c("SubjectID", "Gender")]
colnames(subject_info)[1] = "Subject"
tmp = merge(samples_infection_filtered, subject_info, by = "Subject")

#### Merge
colnames(sc)[1] = "Subject"
tmp = merge(samples_infection_filtered, sc, by = "Subject")

selected_sampels = intersect(samples_infection_filtered$SampleID, tmp$SampleID)
tmp$SampleID = as.character(tmp$SampleID)
tmp2 = tmp[which(tmp$SampleID %in% selected_sampels),]


jpeg("ipop_metabolites_infection_wCycles_timepoints_22.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(tmp2, aes(x = Subject_Cycle, y = DaysAfterInfect, color = state1)) +
  theme_bw() + 
  #scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = state1), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  #scale_y_continuous(name= "Days", breaks = round(seq(0, max(timepoints.df$DaysAfterInfect) + 20, by = 1),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=12, face="bold"), legend.title = element_blank()) +
  coord_flip()
dev.off()




#### Remove all healthy timepoints except the last one
dat = tmp2
for(i in unique(dat$Subject_Cycle)){
  block = dat[which(dat$Subject_Cycle == i & dat$state1 == "Healthy"),]
  nearest = max(block$Date)
  if(block[which(block$Date == nearest),]$DaysAfterInfect < -100){
    dat = dat[-which(dat$Subject_Cycle == i),]
  } else {
    remove = block[-which(block$Date == nearest), ]
    print(dim(remove))
    if(dim(remove)[1]>0){
      dat = dat[-which(dat$SampleID %in% remove$SampleID), ]
    }
  }
}

jpeg("ipop_metabolites_infection_wCycles_timepoints_withOneHealthy.jpg", res = 300, height = 15, width = 20, units = 'cm')
ggplot(dat, aes(x = Subject_Cycle, y = DaysAfterInfect, color = state1)) +
  theme_bw() + 
  #scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject (infection episode)")+ #coord_fixed(ratio = 0.002)+
  scale_y_continuous(name = "Days after Infection") + 
  geom_point(aes(colour = state1, shape = Gender), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="right")+ scale_shape_manual(values=c(1,6)) +
  #scale_y_continuous(name= "Days", breaks = round(seq(0, max(timepoints.df$DaysAfterInfect) + 20, by = 1),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=6,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=6, face="bold"), legend.title = element_blank()) +
  coord_flip()
dev.off()




## Stats
length(dat$Subject_Cycle)
length(unique(dat$Subject_Cycle))
length(unique(dat$Subject))

table(dat$Gender)
unique(dat[, c("Subject", "Gender")])

##### Normalize and get the samples from count matrix

### Select samples 
metabolites_count_table_selected = metabolites_count_table[, tmp2$SampleID]

metabolite = list()
for(i in 1:nrow(metabolites_count_table_selected)){
  ## TODO: Check if the number matches by sample name
  print(i)
  metabolite[[i]] = cbind(tmp2, Count = t(metabolites_count_table_selected[i,]))
}



### CLR Normalization
### Need to keep the name
normalize = function(df){
  #colnames(df)[ncol(df)] = "Count"
  ## Remove time point 0
  y = df #df[-which(df$Time==0),]
  SubjectCycle = unique(y$Subject_Cycle)
  y$rawCount = df[,ncol(df)]
  y$Count = df[,ncol(df)]
  y$feature = colnames(df)[ncol(df)]
  y$normalizedCount = 0
  for(subjCyc in SubjectCycle)
  {
    max_health = max(y[y$Subject_Cycle == subjCyc & y$state1 == "Healthy", "DaysAfterInfect"])
    #m = min(y[y$Subject==subj, "Time"])
    y[y$Subject_Cycle == subjCyc, ]$normalizedCount = log(y[y$Subject_Cycle == subjCyc, ]$Count/y[y$Subject_Cycle == subjCyc & y$DaysAfterInfect == max_health, ]$Count)
  }
  
  ### May add sudo code
  #y$normalizedCount = y$normalizedCount + 0.000000001
  y$Count = y$normalizedCount
  return(y)
}


metabolite_norm = list()
for(i in 1:length(metabolite)){
  metabolite_norm[[i]] = normalize(metabolite[[i]])
}

View(metabolite_norm[[1]])



## remove healthy timepoints
metabolite_norm_filtered = list()
for(i in 1:length(metabolite_norm)){
  metabolite_norm_filtered[[i]] = metabolite_norm[[i]][-which(metabolite_norm[[i]]$state1 == "Healthy"),]
}
View(metabolite_norm_filtered[[1]])


### TODO: Remove infection occurences < 3 timepoints

jpeg("ipop_metabolites_infection_wCycles_timepoints.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(metabolite_norm_filtered[[1]], aes(x = Subject_Cycle, y = DaysAfterInfect, color = state1)) +
  theme_bw() + 
  #scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = state1, shape = Gender), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="right")+ scale_shape_manual(values=c(1,6)) +
  #scale_y_continuous(name= "Days", breaks = round(seq(0, max(timepoints.df$DaysAfterInfect) + 20, by = 1),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=8, face="bold"), legend.title = element_blank()) +
  coord_flip()
dev.off()



### Rename Variables
metabolite_norm_filtered_rename = list()
for(i in 1:length(metabolite_norm_filtered)){
  metabolite_norm_filtered_rename[[i]] = metabolite_norm_filtered[[i]][, c("Subject_Cycle", "SampleID", "Gender", "DaysAfterInfect", "Count", "feature")]
  colnames(metabolite_norm_filtered_rename[[i]]) = c("Subject", "ID", "Group", "Time", "Count", "feature")
  metabolite_norm_filtered_rename[[i]]$Subject = factor(metabolite_norm_filtered_rename[[i]]$Subject)
}







##################################
######### Run OmicsLonDA #########
##################################
points = seq(0, 50, length.out = 51)
data = metabolite_norm_filtered_rename
omicslondaResults = list()

### Test each feature
for(i in 1:length(metabolite_norm_filtered_rename))
{
  omicslondaResults[[i]] = omicslonda(formula = Count ~ Time, df = data[[i]], n.perm = 100, fit.method = "ssgaussian", points = points,
                                      text = unique(data[[i]]$feature), parall = FALSE, pvalue.threshold = 0.05,     
                                      adjust.method = "BH", col = c("pink", "black"), 
                                      prefix = paste("Feature_", unique(data[[i]]$feature), sep = ""), 
                                      ylabel = "Adjusted Metabolite Level", 
                                      DrawTestStatDist = FALSE)
}

### Test all features
omicslonda_metabolomics = omicslondaAll_ipop(formula = Count ~ Time, data = metabolite_norm_filtered_rename, n.perm = 100, fit.method = "ssgaussian", 
                                   num.intervals = 100, parall = FALSE, pvalue.threshold = 0.05, 
                                   adjust.method = "BH", time.unit = "days", norm.method = "none", 
                                   prefix = "OmicsLonDA_ipop_Metabolomics",
                                   ylabel = "Adjusted Metabolite Level", col = c("deeppink2", "black"))
### Annotate metabolites
View(metbcr.df)
View(metb.curated)

subset_omicslonda_metabolomics = omicslonda_metabolomics$output.summary
colnames(subset_omicslonda_metabolomics)[1] = "HMDB"
metabolite_subset = subset(metb.curated, HMDB %in% subset_omicslonda_metabolomics$feature)

View(merge(subset_omicslonda_metabolomics, metabolite_subset, by = "HMDB"))
annotated_metabolite_omicslonda = merge(subset_omicslonda_metabolomics, metabolite_subset, by = "HMDB")
write.csv(annotated_metabolite_omicslonda, file = "annotated_metabolite_omicslonda.csv")

## Test normality of the normalized counts
data = metabolite_norm_filtered_rename
ks_pvalue = numeric(length(data))
for(i in 1:length(data))
{
  hist(data[[i]]$Count)
  res = ks.test(data[[i]]$Count, "pnorm", mean(data[[i]]$Count), sd(data[[i]]$Count))
  ks_pvalue[i] = res$p.value
}
z = p.adjust(ks_pvalue, method = "BH")
sum(z<0.05)






# #### Test covariates
# library(nlme)
# 
# m1 <- lme(Count ~ Time * Group+BMI+Age, random=~1|Subject, data=data[[1]])
# anova(m1)
# 
# tested.df = as.matrix(gut_taxa_count_table_filtered)
# for(i in 1:nrow(tested.df))
# {
#   x = cbind(Count = as.vector(tested.df[i,]), gut_taxa_metadata_table)
#   m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=x)
#   a1 = anova(m1)
#   
#   if(a1["Time:Group", "p-value"]<0.05){
#     cat("Significant  feature = ", rownames(tested.df)[i], "\n")
#     print(a1)
#   }
# }
