library(devtools)
library(zoo)
library(lubridate)
library(ggplot2)
library(ggplot2)
library(gss)
library(plyr)
library(pracma)
library(parallel)
library(doParallel)
library(zoo)
library(pheatmap)
library(viridis)
library("biomformat"); packageVersion("biomformat")
library(nlme)

rm(list=ls())


setwd("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")
#setwd("C:/Users/ametwall/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")
#load("ipop_data/HMP_Microbiome_Count.RData", envir = parent.frame(), verbose = FALSE)
load("ipop_data/Revision_MultiOmes.RData", envir = parent.frame(), verbose = FALSE)





#######################################
########### Metabolomics     ##########
#######################################
metabolites_seasonality_all = read.csv(file = "ipop_data/metabolites_seasonality_all.csv", row.names = 1, check.names=FALSE)

####### Filter out smaples that dont have labels  #######
#metabolites_seasonality_all["IRIS"=="Unknown",]
#remov = which(metabolites_seasonality_all["IRIS",] == "Unknown" | is.na(metabolites_seasonality_all["IRIS",]))
#metabolites_seasonality_all = metabolites_seasonality_all[,-remov]
#dim(metabolites_seasonality_all)



metabolites_count = metabolites_seasonality_all[-c(1:6),]
metabolites_count_table = data.frame(sapply(metabolites_count, function(x) as.numeric(as.character(x))))
rownames(metabolites_count_table) = rownames(metabolites_count)
colnames(metabolites_count_table) = colnames(metabolites_count)

metabolites_metadata_table = as.data.frame(t(metabolites_seasonality_all[c(1:6),]))
colnames(metabolites_metadata_table)[1] = "Subject"
#colnames(metabolites_metadata_table)[4] = "Group"
#colnames(metabolites_metadata_table)[3] = "Time"


metabolites_metadata_table$Date = as.Date(metabolites_metadata_table$Date, format = '%d-%b-%y')


metabolites_metadata_table$SampleID = rownames(metabolites_metadata_table)
#metabolites_metadata_table$SampleID = gsub("X", "", metabolites_metadata_table$SampleID)
#metabolites_metadata_table$SampleID = gsub("\\.", "-", metabolites_metadata_table$SampleID)
metabolites_metadata_table  = merge(metabolites_metadata_table, ls, by = "SampleID")

#### Modify ls df
# colnames(ls)[1] = "Subject"
# colnames(ls)[2] = "SampleID"
# colnames(ls)[3] = "Time"
# metabolites_metadata_table  = merge(ls, metbcr.df, by = "SampleID")


prefix = "Infection"
xx = subset(metabolites_metadata_table, substring(state1, 1, nchar(prefix)) == prefix)
xx$Subject = factor(xx$Subject)

# ### TODO: Filter out subject with less than 5 infection time points. We may do the filteration basedon the subject_cycle parameter later on
# length(unique(xx$Subject))
# sum(table(xx$Subject)>4)
# q = table(xx$Subject)>4
# selectedSubjects = names(q[q == TRUE])

## Select samples from subjects with infection timepoints
subset.df = subset(metabolites_metadata_table, Subject %in% xx$Subject)
colnames(subset.df)[which(colnames(subset.df) == "CL1")] = "DaysAfterInfect"
subset.df$DaysAfterInfect = gsub("D", "", subset.df$DaysAfterInfect)
#subset.df = subset.df[which(subset.df$DaysAfterInfect != ""),]

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


## remove inf one
remove = which(yy$DaysAfterInfect == "-Inf")
yy = yy[-remove,]

## TODO: Commented this
# retain = c("SampleID", "Subject", "Batch", "Time", "Group", "Class", "state1", "DaysAfterInfect")
# yy = yy[,retain]



## name each infection cycle by the date of the first infection sample
yy$CycleName = ""### TODO: Name infection cycles
class(yy$CycleName) = "Date"
## Order yy
yy = yy[with(yy, order(Subject, Date)),]

write.csv(yy, file = "ipop_data/ipop_metabolites_infection_2.csv")






### TODO: rename the cycle of the block that downnot have Healthy and consecutive to infection timepoints
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





### TODO: Remove sampels without Healthy group
samples_infection_filtered = samples_infection
for(i in unique(samples_infection$Subject_Cycle)){
  block = samples_infection[which(samples_infection$Subject_Cycle == i), "state1"]
  if(!any(block == "Healthy")){
    samples_infection_filtered = samples_infection_filtered[-which(samples_infection_filtered$Subject_Cycle == i), ]
  }
}



## Add gender
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





## Statistics
length(dat$Subject_Cycle)
length(unique(dat$Subject_Cycle))
length(unique(dat$Subject))

table(dat$Gender)
unique(dat[, c("Subject", "Gender")])

##### Normalize and get the samples from count matrix
tmp2
### Select samples from 
#colnames(metabolites_count_table) = gsub("X", "", colnames(metabolites_count_table))
#colnames(metabolites_count_table) = gsub("\\.", "-", colnames(metabolites_count_table))
metabolites_count_table_selected = metabolites_count_table[, tmp2$SampleID]

metabolite = list()
for(i in 1:nrow(metabolites_count_table_selected)){
  ## TODO: Check if the number matches by sample name
  print(i)
  metabolite[[i]] = cbind(tmp2, Count = t(metabolites_count_table_selected[i,]))
}





### Normalize

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





### Rename Variable
metabolite_norm_filtered_rename = list()
for(i in 1:length(metabolite_norm_filtered)){
  metabolite_norm_filtered_rename[[i]] = metabolite_norm_filtered[[i]][, c("Subject_Cycle", "SampleID", "Gender", "DaysAfterInfect", "Count", "feature")]
  colnames(metabolite_norm_filtered_rename[[i]]) = c("Subject", "ID", "Group", "Time", "Count", "feature")
  metabolite_norm_filtered_rename[[i]]$Subject = factor(metabolite_norm_filtered_rename[[i]]$Subject)
  #metabolite_norm_filtered_rename[[i]]$Subject = factor(metabolite_norm_filtered_rename[[i]]$Subject)
}








### Rename variable and run OmicsLonda
source("OmicsLonDA/R/OmicsLonDA.R")
source("OmicsLonDA/R/CurveFitting.R")
source("OmicsLonDA/R/Visualization.R")
source("OmicsLonDA/R/Permutation.R")
source("OmicsLonDA/R/Normalization.R")
source("OmicsLonDA/R/OmicsLonDA_Evaluation.R")


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
omicslonda_metabolomics = omicslondaAll_Test_iPOP(formula = Count ~ Time, data = metabolite_norm_filtered_rename, n.perm = 100, fit.method = "ssgaussian", 
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






#### Test covariates
library(nlme)

### Add BMI and Age to this equation
m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=data[[1]])
anova(m1)

tested.df = as.matrix(gut_taxa_count_table_filtered)
for(i in 1:nrow(tested.df))
{
  x = cbind(Count = as.vector(tested.df[i,]), gut_taxa_metadata_table)
  m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=x)
  a1 = anova(m1)
  
  if(a1["Time:Group", "p-value"]<0.05){
    cat("Significant  feature = ", rownames(tested.df)[i], "\n")
    print(a1)
  }
}





## Count number of unique infection cycle
length(unique(samples_infection$Subject_Cycle))
## 27



## Do Testing Between IS & IR
length(unique(paste(samples_infection$Subject_Cycle, samples_infection$))

subset.df_IRIS = subset.df[,c("Subject", "Group")]

dim(unique(subset.df_IRIS))
View(unique(subset.df_IRIS))
## 2 IR, 7 IS
# table(unique(subset.df_IRIS))



### Time distribution

timepoints.df = subset.df[, c("Subject", "Time", "DaysAfterInfect", "state1")]
# TimePoints = data.frame(subjectID = as.factor(paste("",gut_taxa_id, sep="")), days = gut_taxa_time, 
#                         IRIS = gut_taxa_group)
# TimePoints = TimePoints[order(TimePoints$IRIS, TimePoints$subjectID),]
# TimePoints = data.frame(subjectID = factor(TimePoints$subjectID, levels = TimePoints$subjectID), Time = TimePoints$days, 
#                         Group= as.factor(TimePoints$IRIS))


#jpeg("gut_taxa_timepoints_distribution.jpg", res = 300, height = 20, width = 20, units = 'cm')
class(yy$DaysAfterInfect) = "numeric"
ggplot(yy, aes(x = Subject, y = DaysAfterInfect, color = state1)) +
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
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
# h = recordPlot(load="gridGraphics", attach=NULL)
#dev.off()





View(ls)
View(metbcr.df)
View(metb.curated)
View(metabolites_metadata_table)


metabolites_omicslondaAll = omicslondaAll(formula = Count ~ Time + state1, countTable = metabolites_count_table, 
                                          metadata = metabolites_metadata_table, 
                                          n.perm = 100, fit.method = "ssgaussian", num.intervals = 36,
                                          parall = FALSE, pvalue.threshold = 0.05,     
                                          adjust.method = "BH", time.unit = "day", norm.method = "none", 
                                          prefix = "Test_OmicsLonDA_metabolites_ipop_seasonal_all_RT", 
                                          ylabel = "Level", col = c("black", "green"))










#######################################
########### Gut Microbiome     ########
#######################################
gut_taxa_seasonality_all = read.csv(file = "ipop_data/gut_taxa_seasonality_all.csv", row.names = 1)


####### Filter out smaples that dont have labels  #######
#gut_taxa_seasonality_all["IRIS"=="Unknown",]
remov = which(gut_taxa_seasonality_all["IRIS",] == "Unknown" | is.na(gut_taxa_seasonality_all["IRIS",]))
gut_taxa_seasonality_all = gut_taxa_seasonality_all[,-remov]
dim(gut_taxa_seasonality_all)
#View(gut_taxa_seasonality_all)





gut_taxa_count = gut_taxa_seasonality_all[-c(1:8),]
gut_taxa_count_table = data.frame(sapply(gut_taxa_count, function(x) as.numeric(as.character(x))))
rownames(gut_taxa_count_table) = rownames(gut_taxa_count)
colnames(gut_taxa_count_table) = colnames(gut_taxa_count)

## filter out features that dont present in at least 25% of the samples
numLowAbundantSamples = apply(gut_taxa_count_table, 1, function(x){
  sum(x < 0.001)
})

tobeRemoved = numLowAbundantSamples > 0.25*ncol(gut_taxa_count_table)
gut_taxa_count_table_filtered = gut_taxa_count_table[!tobeRemoved,]

gut_taxa_metadata_table = as.data.frame(t(gut_taxa_seasonality_all[c(1:8),]))
colnames(gut_taxa_metadata_table)[1] = "Subject"
colnames(gut_taxa_metadata_table)[4] = "Group"
colnames(gut_taxa_metadata_table)[3] = "Time"
gut_taxa_metadata_table$Time = as.Date(gut_taxa_metadata_table$Time, format = '%d-%b-%y')
prefix = "Infection"
xx = subset(gut_taxa_metadata_table, substring(state1, 1, nchar(prefix)) == prefix)
# subset(gut_taxa_metadata_table, Class == "Diabetic")
# subset(gut_taxa_metadata_table, state1 %in% "Infection")



xx$Subject = factor(xx$Subject)
length(unique(xx$Subject))
table(xx$Subject)>5
sum(table(xx$Subject)>4)

q = table(xx$Subject)>4
selectedSubjects = names(q[q == TRUE])


subset.df = subset(gut_taxa_metadata_table, Subject %in% selectedSubjects)
write.csv(subset.df, file = "ipop_infection.csv")

## ** Select one time point to act as baseline for each subject 
## ** Select healthy timepoints

### Day of the Year
gut_taxa_metadata_table[, "Time"] = gsub("Jan", "01", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Feb", "02", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Mar", "03", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Apr", "04", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("May", "05", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Jun", "06", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Jul", "07", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Aug", "08", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Sep", "09", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Oct", "10", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Nov", "11", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = gsub("Dec", "12", gut_taxa_metadata_table[, "Time"])
gut_taxa_metadata_table[, "Time"]  = yday(as.Date(gut_taxa_metadata_table[, "Time"] , format = '%d-%m-%Y'))



#### Global Testing




library(nlme)
# m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=diff_simulatedDataset[[1]])
# anova(m1)
# 
# m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=nondiff_simulatedDataset[[1]])
# anova(m1)
# 

# countTable = gut_taxa_count_table_filtered
# metadata = gut_taxa_metadata_table

m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=nondiff_simulatedDataset[[1]])
anova(m1)

tested.df = as.matrix(gut_taxa_count_table_filtered)
for(i in 1:nrow(tested.df))
{
  x = cbind(Count = as.vector(tested.df[i,]), gut_taxa_metadata_table)
  m1 <- lme(Count ~ Time * Group, random=~1|Subject, data=x)
  a1 = anova(m1)
  
  if(a1["Time:Group", "p-value"]<0.05){
    cat("Significant  feature = ", rownames(tested.df)[i], "\n")
    print(a1)
  }
}



#View(gut_taxa_metadata_table)
#View(gut_taxa_count_table_filtered)
gut_taxa_omicslondaAll = omicslondaAll(formula = Count ~ Time + state1, countTable = gut_taxa_count_table_filtered, 
                                       metadata = gut_taxa_metadata_table, 
                                       n.perm = 100, fit.method = "ssgaussian", num.intervals = 36,
                                       parall = FALSE, pvalue.threshold = 0.05,     
                                       adjust.method = "BH", time.unit = "day", norm.method = "none", 
                                       prefix = "Test_OmicsLonDA_gut_taxa_ipop_seasonal_filtered_RT_newsampling", 
                                       ylabel = "Relative Abundance", col = c("black", "green"))






### Time distribution
TimePoints = data.frame(subjectID = as.factor(paste("",gut_taxa_id, sep="")), days = gut_taxa_time, 
                        IRIS = gut_taxa_group)
TimePoints = TimePoints[order(TimePoints$IRIS, TimePoints$subjectID),]
TimePoints = data.frame(subjectID = factor(TimePoints$subjectID, levels = TimePoints$subjectID), Time = TimePoints$days, 
                        Group= as.factor(TimePoints$IRIS))


jpeg("gut_taxa_timepoints_distribution.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(TimePoints, aes(x = subjectID, y = Time)) +
  theme_bw() + 
  scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = Group), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  scale_y_continuous(name= "Days", breaks = round(seq(0, max(TimePoints$Time) + 20, by = 100),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
# h = recordPlot(load="gridGraphics", attach=NULL)
dev.off()






#######################################
########### Clinical Data     #########
#######################################
clinical_seasonality = read.table("data/clinic.df_data_seasonality.txt", row.names = 1)
clinical_seasonality_metadata = read.table("data/clinic.df_data_seasonality_colData.txt", header = TRUE, row.names = 1)

dim(clinical_seasonality)
dim(clinical_seasonality_metadata)



t_clinical_seasonality_metadata = t(clinical_seasonality_metadata)
colnames(clinical_seasonality) = colnames(t_clinical_seasonality_metadata)
clinical_seasonality_all = rbind(t_clinical_seasonality_metadata, clinical_seasonality)
#write.csv(clinical_seasonality_all, file = "clinical_seasonality_all.csv")


####### Filter out smaples that dont have labels  #######
#clinical_seasonality_all["IRIS"=="Unknown",]
remov = which(clinical_seasonality_all["IRIS",] == "Unknown" | is.na(clinical_seasonality_all["IRIS",]))
clinical_seasonality_all = clinical_seasonality_all[,-remov]
dim(clinical_seasonality_all)


## OmicsLonDA
clinical_count = clinical_seasonality_all[-c(1:8),]
clinical_count_table = data.frame(sapply(clinical_count, function(x) as.numeric(as.character(x))))
rownames(clinical_count_table) = rownames(clinical_count)
colnames(clinical_count_table) = colnames(clinical_count)

# Filteration
numNAsamples = apply(clinical_count_table, 1, function(x){
  sum(is.na(x))
})

tobeRemoved = numNAsamples > 0.25*ncol(clinical_count_table)
clinical_count_table_filtered = clinical_count_table[!tobeRemoved,]


clinical_metadata_table = as.data.frame(t(clinical_seasonality_all[c(1:8),]))
colnames(clinical_metadata_table)[1] = "Subject"
colnames(clinical_metadata_table)[4] = "Group"
colnames(clinical_metadata_table)[3] = "Time"

clinical_metadata_table[, "Time"] = gsub("Jan", "01", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Feb", "02", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Mar", "03", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Apr", "04", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("May", "05", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Jun", "06", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Jul", "07", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Aug", "08", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Sep", "09", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Oct", "10", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Nov", "11", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = gsub("Dec", "12", clinical_metadata_table[, "Time"])
clinical_metadata_table[, "Time"]  = yday(as.Date(clinical_metadata_table[, "Time"] , format = '%d-%m-%Y'))

#View(clinical_count_table)
#View(clinical_metadata_table)

clinical_omicslondaAll = omicslondaAll(formula = Count ~ Time, countTable = clinical_count_table_filtered, 
                                       metadata = clinical_metadata_table, 
                                       n.perm = 100, fit.method = "ssgaussian", num.intervals = 36,
                                       parall = FALSE, pvalue.threshold = 0.05,     
                                       adjust.method = "BH", time.unit = "day", norm.method = "none", 
                                       prefix = "Test_OmicsLonDA_clinical_ipop_seasonal_all_RT", 
                                       ylabel = "Level", col = c("black", "green"))




### Time distribution
TimePoints = data.frame(subjectID = as.factor(paste("",clinical_id, sep="")), days = clinical_time, 
                        IRIS = clinical_group)
TimePoints = TimePoints[order(TimePoints$IRIS, TimePoints$subjectID),]
TimePoints = data.frame(subjectID = factor(TimePoints$subjectID, levels = TimePoints$subjectID), Time = TimePoints$days, 
                        Group= as.factor(TimePoints$IRIS))


jpeg("clinical_timepoints_distribution.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(TimePoints, aes(x = subjectID, y = Time)) +
  theme_bw() + 
  scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = Group), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  scale_y_continuous(name= "Days", breaks = round(seq(0, max(TimePoints$Time) + 20, by = 100),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
# h = recordPlot(load="gridGraphics", attach=NULL)
dev.off()





















#########################################
########### Nasal Microbiome     ########
#########################################
nasal_taxa_seasonality = read.table("data/nasalTaxa_data_seasonality.txt", row.names = 1)
nasal_taxa_seasonality_metadata = read.table("data/nasalTaxa_data_seasonality_colData.txt", header = TRUE, row.names = 1)

dim(nasal_taxa_seasonality)
dim(nasal_taxa_seasonality_metadata)

t_nasal_taxa_seasonality_metadata = t(nasal_taxa_seasonality_metadata)
colnames(nasal_taxa_seasonality) = colnames(t_nasal_taxa_seasonality_metadata)
nasal_taxa_seasonality_all = rbind(t_nasal_taxa_seasonality_metadata, nasal_taxa_seasonality)
write.csv(nasal_taxa_seasonality_all, file = "nasal_taxa_seasonality_all.csv")


####### Filter out smaples that dont have labels  #######
#nasal_taxa_seasonality_all["IRIS"=="Unknown",]
remov = which(nasal_taxa_seasonality_all["IRIS",] == "Unknown" | is.na(nasal_taxa_seasonality_all["IRIS",]))
nasal_taxa_seasonality_all = nasal_taxa_seasonality_all[,-remov]
dim(nasal_taxa_seasonality_all)


### OmicsLonDA
nasal_taxa_count = nasal_taxa_seasonality_all[-c(1:8),]
nasal_taxa_count_table = data.frame(sapply(nasal_taxa_count, function(x) as.numeric(as.character(x))))
rownames(nasal_taxa_count_table) = rownames(nasal_taxa_count)
colnames(nasal_taxa_count_table) = colnames(nasal_taxa_count)

## filter out features that dont present in at least 25% of the samples
numLowAbundantSamples = apply(nasal_taxa_count_table, 1, function(x){
  sum(x < 0.001)
})

tobeRemoved = numLowAbundantSamples > 0.25*ncol(nasal_taxa_count_table)
nasal_taxa_count_table_filtered = nasal_taxa_count_table[!tobeRemoved,]

dim(nasal_taxa_count_table_filtered)

nasal_taxa_metadata_table = as.data.frame(t(nasal_taxa_seasonality_all[c(1:8),]))
colnames(nasal_taxa_metadata_table)[1] = "Subject"
colnames(nasal_taxa_metadata_table)[4] = "Group"
colnames(nasal_taxa_metadata_table)[3] = "Time"

nasal_taxa_metadata_table[, "Time"] = gsub("Jan", "01", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Feb", "02", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Mar", "03", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Apr", "04", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("May", "05", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Jun", "06", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Jul", "07", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Aug", "08", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Sep", "09", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Oct", "10", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Nov", "11", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = gsub("Dec", "12", nasal_taxa_metadata_table[, "Time"])
nasal_taxa_metadata_table[, "Time"]  = yday(as.Date(nasal_taxa_metadata_table[, "Time"] , format = '%d-%m-%Y'))

View(nasal_taxa_metadata_table)
View(nasal_taxa_count_table_filtered)
nasal_taxa_omicslondaAll = omicslondaAll(formula = Count ~ Time + state1, countTable = nasal_taxa_count_table_filtered, 
                                         metadata = nasal_taxa_metadata_table, 
                                         n.perm = 100, fit.method = "ssgaussian", num.intervals = 36,
                                         parall = FALSE, pvalue.threshold = 0.05,     
                                         adjust.method = "BH", time.unit = "day", norm.method = "none", 
                                         prefix = "Test_OmicsLonDA_nasal_taxa_ipop_seasonal_all_filtered_RT", 
                                         ylabel = "Relative Abundance", col = c("black", "green"))




### Time distribution
TimePoints = data.frame(subjectID = as.factor(paste("",nasal_taxa_id, sep="")), days = nasal_taxa_time, 
                        IRIS = nasal_taxa_group)
TimePoints = TimePoints[order(TimePoints$IRIS, TimePoints$subjectID),]
TimePoints = data.frame(subjectID = factor(TimePoints$subjectID, levels = TimePoints$subjectID), Time = TimePoints$days, 
                        Group= as.factor(TimePoints$IRIS))


jpeg("nasal_taxa_timepoints_distribution.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(TimePoints, aes(x = subjectID, y = Time)) +
  theme_bw() + 
  scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = Group), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  scale_y_continuous(name= "Days", breaks = round(seq(0, max(TimePoints$Time) + 20, by = 100),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
# h = recordPlot(load="gridGraphics", attach=NULL)
dev.off()









#######################################
############# Cytokines    ############
#######################################
cytokines_seasonality = read.table("data/Cytokines_data_seasonality.txt", row.names = 1)
cytokines_seasonality_metadata = read.table("data/Cytokines_data_seasonality_colData.txt", header = TRUE, row.names = 1)

dim(cytokines_seasonality)
dim(cytokines_seasonality_metadata)


t_cytokines_seasonality_metadata = t(cytokines_seasonality_metadata)
colnames(cytokines_seasonality) = colnames(t_cytokines_seasonality_metadata)
cytokines_seasonality_all = rbind(t_cytokines_seasonality_metadata, cytokines_seasonality)
write.csv(cytokines_seasonality_all, file = "cytokines_seasonality_all.csv")

#########################################################
####### Filter out smaples that dont have labels  #######
#########################################################
#cytokines_seasonality_all["IRIS"=="Unknown",]
remov = which(cytokines_seasonality_all["IRIS",] == "Unknown" | is.na(cytokines_seasonality_all["IRIS",]))
cytokines_seasonality_all = cytokines_seasonality_all[,-remov]
dim(cytokines_seasonality_all)


### OmicsLonDA
cytokines_count = cytokines_seasonality_all[-c(1:8),]
cytokines_count_table = data.frame(sapply(cytokines_count, function(x) as.numeric(as.character(x))))
rownames(cytokines_count_table) = rownames(cytokines_count)
colnames(cytokines_count_table) = colnames(cytokines_count)

cytokines_metadata_table = as.data.frame(t(cytokines_seasonality_all[c(1:8),]))
colnames(cytokines_metadata_table)[1] = "Subject"
colnames(cytokines_metadata_table)[4] = "Group"
colnames(cytokines_metadata_table)[3] = "Time"

cytokines_metadata_table[, "Time"] = gsub("Jan", "01", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Feb", "02", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Mar", "03", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Apr", "04", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("May", "05", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Jun", "06", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Jul", "07", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Aug", "08", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Sep", "09", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Oct", "10", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Nov", "11", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = gsub("Dec", "12", cytokines_metadata_table[, "Time"])
cytokines_metadata_table[, "Time"]  = yday(as.Date(cytokines_metadata_table[, "Time"] , format = '%d-%m-%Y'))

#View(cytokines_metadata_table)
#View(cytokines_count_table)
cytokines_omicslondaAll = omicslondaAll(formula = Count ~ Time + state1, countTable = cytokines_count_table, 
                                        metadata = cytokines_metadata_table, 
                                        n.perm = 100, fit.method = "ssgaussian", num.intervals = 36,
                                        parall = FALSE, pvalue.threshold = 0.05,     
                                        adjust.method = "BH", time.unit = "day", norm.method = "none", 
                                        prefix = "Test_OmicsLonDA_cytokines_ipop_seasonal_all_RT", 
                                        ylabel = "Level", col = c("black", "green"))











### Time distribution
TimePoints = data.frame(subjectID = as.factor(paste("", cytokines_id, sep="")), days = cytokines_time, 
                        IRIS = cytokines_group)
TimePoints = TimePoints[order(TimePoints$IRIS, TimePoints$subjectID),]
TimePoints = data.frame(subjectID = factor(TimePoints$subjectID, levels = TimePoints$subjectID), Time = TimePoints$days, 
                        Group= as.factor(TimePoints$IRIS))


jpeg("cytokines_timepoints_distribution.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(TimePoints, aes(x = subjectID, y = Time)) +
  theme_bw() + 
  scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = Group), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  scale_y_continuous(name= "Days", breaks = round(seq(0, max(TimePoints$Time) + 20, by = 100),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
# h = recordPlot(load="gridGraphics", attach=NULL)
dev.off()



#######################################
############# Proteomics     ##########
#######################################
protiens_seasonality = read.table("data/Protiens_data_seasonality.txt", row.names = 1)
protiens_seasonality_metadata = read.table("data/Protiens_data_seasonality_colData.txt", header = TRUE, row.names = 1)

dim(protiens_seasonality)
dim(protiens_seasonality_metadata)

t_protiens_seasonality_metadata = t(protiens_seasonality_metadata)
colnames(protiens_seasonality) = colnames(t_protiens_seasonality_metadata)
protiens_seasonality_all = rbind(t_protiens_seasonality_metadata, protiens_seasonality)
write.csv(protiens_seasonality_all, file = "protiens_seasonality_all.csv")


####### Filter out smaples that dont have labels  #######
#protiens_seasonality_all["IRIS"=="Unknown",]
remov = which(protiens_seasonality_all["IRIS",] == "Unknown" | is.na(protiens_seasonality_all["IRIS",]))
protiens_seasonality_all = protiens_seasonality_all[,-remov]
dim(protiens_seasonality_all)


### OmicsLonDA
protiens_count = protiens_seasonality_all[-c(1:8),]
protiens_count_table = data.frame(sapply(protiens_count, function(x) as.numeric(as.character(x))))
rownames(protiens_count_table) = rownames(protiens_count)
colnames(protiens_count_table) = colnames(protiens_count)

protiens_metadata_table = as.data.frame(t(protiens_seasonality_all[c(1:8),]))
colnames(protiens_metadata_table)[1] = "Subject"
colnames(protiens_metadata_table)[4] = "Group"
colnames(protiens_metadata_table)[3] = "Time"

protiens_metadata_table[, "Time"] = gsub("Jan", "01", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Feb", "02", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Mar", "03", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Apr", "04", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("May", "05", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Jun", "06", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Jul", "07", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Aug", "08", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Sep", "09", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Oct", "10", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Nov", "11", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = gsub("Dec", "12", protiens_metadata_table[, "Time"])
protiens_metadata_table[, "Time"]  = yday(as.Date(protiens_metadata_table[, "Time"] , format = '%d-%m-%Y'))

View(protiens_metadata_table)
View(protiens_count_table)
protiens_omicslondaAll = omicslondaAll(formula = Count ~ Time + state1, countTable = protiens_count_table, 
                                       metadata = protiens_metadata_table, 
                                       n.perm = 100, fit.method = "ssgaussian", num.intervals = 36,
                                       parall = FALSE, pvalue.threshold = 0.05,     
                                       adjust.method = "BH", time.unit = "day", norm.method = "none", 
                                       prefix = "Test_OmicsLonDA_protiens_ipop_seasonal_all_RT", 
                                       ylabel = "Level", col = c("black", "green"))





### Time distribution
TimePoints = data.frame(subjectID = as.factor(protiens_metadata_table[, "Subject"]), days = protiens_metadata_table[, "Time"], 
                        IRIS = protiens_metadata_table[, "Group"])
TimePoints = TimePoints[order(TimePoints$IRIS, TimePoints$subjectID),]
TimePoints = data.frame(subjectID = factor(TimePoints$subjectID, levels = TimePoints$subjectID), Time = TimePoints$days, 
                        Group= as.factor(TimePoints$IRIS))


jpeg("proteins_timepoints_distribution.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(TimePoints, aes(x = subjectID, y = Time)) +
  theme_bw() + 
  scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = Group), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  scale_y_continuous(name= "Days", breaks = round(seq(0, max(TimePoints$Time) + 20, by = 100),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
# h = recordPlot(load="gridGraphics", attach=NULL)
dev.off()






#######################################
################# PBMC     ############
#######################################
## TODO how to read all data
## TODO: Check fill data
pbmc_seasonality = read.table("data/PBMCxCell_data_seasonality.txt", row.names = 1, check.names=FALSE, na.strings=c("","NA")) #line 2 did not have 881 elements
pbmc_seasonality_metadata = read.table("data/PBMCxCell_data_seasonality_colData.txt", header = TRUE, row.names = 1)

dim(pbmc_seasonality)
dim(pbmc_seasonality_metadata)

t_pbmc_seasonality_metadata = t(pbmc_seasonality_metadata)
colnames(pbmc_seasonality) = colnames(t_pbmc_seasonality_metadata)
pbmc_seasonality_all = rbind(t_pbmc_seasonality_metadata, pbmc_seasonality)
write.csv(pbmc_seasonality_all, file = "pbmc_seasonality_all.csv")

####### Filter out smaples that dont have labels  #######
#pbmc_seasonality_all["IRIS"=="Unknown",]
remov = which(pbmc_seasonality_all["IRIS",] == "Unknown" | is.na(pbmc_seasonality_all["IRIS",]))
pbmc_seasonality_all = pbmc_seasonality_all[,-remov]
dim(pbmc_seasonality_all)


### OmicsLonDA
metabolites_count = metabolites_seasonality_all[-c(1:8),]
metabolites_count_table = data.frame(sapply(metabolites_count, function(x) as.numeric(as.character(x))))
rownames(metabolites_count_table) = rownames(metabolites_count)
colnames(metabolites_count_table) = colnames(metabolites_count)

metabolites_metadata_table = as.data.frame(t(metabolites_seasonality_all[c(1:8),]))
colnames(metabolites_metadata_table)[1] = "Subject"
colnames(metabolites_metadata_table)[4] = "Group"
colnames(metabolites_metadata_table)[3] = "Time"

metabolites_metadata_table[, "Time"] = gsub("Jan", "01", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Feb", "02", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Mar", "03", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Apr", "04", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("May", "05", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Jun", "06", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Jul", "07", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Aug", "08", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Sep", "09", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Oct", "10", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Nov", "11", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = gsub("Dec", "12", metabolites_metadata_table[, "Time"])
metabolites_metadata_table[, "Time"]  = yday(as.Date(metabolites_metadata_table[, "Time"] , format = '%d-%m-%Y'))

View(metabolites_metadata_table)
View(metabolites_count_table)
metabolites_omicslondaAll = omicslondaAll(formula = Count ~ Time + state1, countTable = metabolites_count_table, 
                                          metadata = metabolites_metadata_table, 
                                          n.perm = 100, fit.method = "ssgaussian", num.intervals = 36,
                                          parall = FALSE, pvalue.threshold = 0.05,     
                                          adjust.method = "none", time.unit = "day", norm.method = "none", 
                                          prefix = "Test_OmicsLonDA_metabolites_ipop_seasonal_all", 
                                          ylabel = "Level", col = c("black", "green"))














#################################################
########### Pairwise Correlation    #############
#################################################
### Import OmicsLonDA RData files for data1 and data2
gut_microbiome_omicslonda = get(load("Test_OmicsLonDA_gut_taxa_ipop_seasonal_filtered_RT/OmicsLonDA_Summary_ssgaussian_Test_OmicsLonDA_gut_taxa_ipop_seasonal_filtered_RT.RData", envir = parent.frame(), verbose = FALSE))
clinical_omicslonda = get(load("Test_OmicsLonDA_clinical_ipop_seasonal_all_RT/OmicsLonDA_Summary_ssgaussian_Test_OmicsLonDA_clinical_ipop_seasonal_all_RT.RData", envir = parent.frame(), verbose = FALSE))
cytokines_omicslonda = get(load("Test_OmicsLonDA_cytokines_ipop_seasonal_all_RT/OmicsLonDA_Summary_ssgaussian_Test_OmicsLonDA_cytokines_ipop_seasonal_all_RT.RData", envir = parent.frame(), verbose = FALSE))

omicslondaCorr(data1 = gut_microbiome_omicslonda, data2 = clinical_omicslonda,
               num.intervals = 34, parall = FALSE, pvalue.threshold = 0.05,
               time.unit = "days", prefix = "OmicsLonDA_Correlation_gutTaxa_clinical",
               col = c("blue", "firebrick"))

omicslondaCorr(data1 = gut_microbiome_omicslonda, data2 = cytokines_omicslonda,
               num.intervals = 34, parall = FALSE, pvalue.threshold = 0.05,
               time.unit = "days", prefix = "OmicsLonDA_Correlation_gutTaxa_cytokines",
               col = c("blue", "firebrick"))

omicslondaCorr(data1 = clinical_omicslonda, data2 = cytokines_omicslonda,
               num.intervals = 34, parall = FALSE, pvalue.threshold = 0.05,
               time.unit = "days", prefix = "OmicsLonDA_Correlation_clinical_cytokines",
               col = c("blue", "firebrick"))





####################################################################
###########   Heatmap clustering of OmicsLonDA output  #############
####################################################################
gut_microbiome_omicslonda = get(load("Test_OmicsLonDA_gut_taxa_ipop_seasonal_filtered_RT/OmicsLonDA_Summary_ssgaussian_Test_OmicsLonDA_gut_taxa_ipop_seasonal_filtered_RT.RData", envir = parent.frame(), verbose = FALSE))
clinical_omicslonda = get(load("Test_OmicsLonDA_clinical_ipop_seasonal_all_RT/OmicsLonDA_Summary_ssgaussian_Test_OmicsLonDA_clinical_ipop_seasonal_all_RT.RData", envir = parent.frame(), verbose = FALSE))
cytokines_omicslonda = get(load("Test_OmicsLonDA_cytokines_ipop_seasonal_all_RT/OmicsLonDA_Summary_ssgaussian_Test_OmicsLonDA_cytokines_ipop_seasonal_all_RT.RData", envir = parent.frame(), verbose = FALSE))


prefix = "TEST_VIS_Seasonal"
dir.create(file.path(prefix), showWarnings = FALSE)
visualizeTimeIntervals(interval.details = gut_microbiome_omicslonda$output.summary, prefix = prefix, 
                       unit = "DAYS", col = c("blue", "green"), fit.method = "ssgaussian")

### Visulaize data as heatmap
padjust_threshold = 0.05



## z1 = FIRST OMIC
z1 = do.call(rbind, gut_microbiome_omicslonda$output.detail)
z1_padjust = do.call(rbind, z1[,"adjusted.pvalue"])
nonsig = which(z1_padjust > padjust_threshold)
sig = which(z1_padjust <= padjust_threshold)
z1_padjust[sig] = 1
z1_padjust[nonsig] = 0
z1_sign = do.call(rbind, z1[,"testStat.sign"])
z1_vis = z1_padjust*z1_sign

## add row/col names
z1_name = do.call(rbind, z1[,"feature"])
z1_features = z1_name[,1]
rownames(z1_vis) = z1_features
colnames(z1_vis) = 1:ncol(z1_vis)

## Filterout nonsignificant features
z1_vis_filtered = z1_vis[apply(z1_vis, 1, function(x) !all(x==0)),]
View(z1_vis_filtered)
View(z1_vis)



## z2 = SECOND OMIC
z2 = do.call(rbind, clinical_omicslonda$output.detail)
z2_padjust = do.call(rbind, z2[,"adjusted.pvalue"])
nonsig = which(z2_padjust > padjust_threshold)
sig = which(z2_padjust <= padjust_threshold)
z2_padjust[sig] = 1
z2_padjust[nonsig] = 0
z2_sign = do.call(rbind, z2[,"testStat.sign"])
z2_vis = z2_padjust*z2_sign

## add row/col names
z2_name = do.call(rbind, z2[,"feature"])
z2_features = z2_name[,1]
rownames(z2_vis) = z2_features
colnames(z2_vis) = 1:ncol(z2_vis)

## Filterout nonsignificant features
z2_vis_filtered = z2_vis[apply(z2_vis, 1, function(x) !all(x==0)),]
View(z2_vis_filtered)
View(z2_vis)



## z3 = THIRD OMIC
z3 = do.call(rbind, cytokines_omicslonda$output.detail)
z3_padjust = do.call(rbind, z3[,"adjusted.pvalue"])
nonsig = which(z3_padjust > padjust_threshold)
sig = which(z3_padjust <= padjust_threshold)
z3_padjust[sig] = 1
z3_padjust[nonsig] = 0
z3_sign = do.call(rbind, z3[,"testStat.sign"])
z3_vis = z3_padjust*z3_sign

## add row/col names
z3_name = do.call(rbind, z3[,"feature"])
z3_features = z3_name[,1]
rownames(z3_vis) = z3_features
colnames(z3_vis) = 1:ncol(z3_vis)

## Filterout nonsignificant features
z3_vis_filtered = z3_vis[apply(z3_vis, 1, function(x) !all(x==0)),]
View(z3_vis_filtered)
View(z3_vis)




dim(z1_vis)
dim(z2_vis)
dim(z3_vis)

z_all_vis_filtered = rbind(z1_vis_filtered, z2_vis_filtered, z3_vis_filtered)
View(z_all_vis_filtered)
dim(z_all_vis_filtered)


time = seq(1:ncol(z_all_vis_filtered)) 
group = c(rep("gut_microbiome", nrow(z1_vis_filtered)), rep("clinical", nrow(z2_vis_filtered)), rep("cytokines", nrow(z3_vis_filtered))) 
#microbiome.group = c("Bacterio", "Phy", "ss", "Bacterio", "Phy", "ss") 
#rnaseq.group = c("NA", "NAA", "brca", "Kras", "gdep", "Kar4")
annot_row =  data.frame(Group = as.factor(group))
rownames(annot_row) = rownames(z_all_vis_filtered)
annot_column = data.frame(Time_Interval = time)


## Try to have numbers as factors so that it can be drawn as discrete numbers instead of continous numbers 
# z2_vis_filteerd_df = as.data.frame(z2_vis_filteerd)
# c = apply(z2_vis_filteerd_df, 2, function(x){as.factor(x)})

cbbPalette_12 = c("#999999", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "darkorchid",
                  "mistyrose1", "sienna1")

jpeg(filename = "Visualize_Timeinterval.jpeg", res = 1200, height = 20, width = 40, units = 'cm')
pheatmap(
  mat               = z_all_vis_filtered,
  color             = c("black", "white", "green"),# inferno(15),
  border_color      = "gray",
  #breaks = c(-1,0,1),
  legend = FALSE,
  # cluster_cols      = daily_cluster_cols,
  # cluster_rows      = daily_cluster_rows,
  cluster_cols = FALSE, 
  cluster_rows = TRUE,
  #kmeans_k = 3,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_row    = annot_row,
  annotation_col    = annot_column,
  annotation_colors = list(Group = c(gut_microbiome = "#D55E00", clinical = "#F0E442", cytokines = "#0072B2"), 
                           Time_Interval = c("white", "gray")),
  #gaps_row = c(1,12,14,16,31,33), 
  drop_levels       = TRUE,
  fontsize          = 10,
  main              = ""
  
)
dev.off()





####################################
######### Quantitative    ##########
####################################


dim(gut_taxa_count_table)
dim(gut_taxa_metadata_table)
dim(gut_taxa_count_table_filtered)

length(unique(gut_taxa_metadata_table$Subject))
length(unique(gut_taxa_metadata_table$Group))
length(unique(gut_taxa_metadata_table$state1))
table(gut_taxa_metadata_table[, c("Group")])






#######################################
############## RNAseq     #############
#######################################
rna_seasonality = read.table("data/RNA_data_seasonality.txt", row.names = 1)
rna_seasonality_metadata = read.table("data/RNA_data_seasonality_colData.txt", header = TRUE, row.names = 1)

dim(rna_seasonality)
dim(rna_seasonality_metadata)


t_rna_seasonality_metadata = t(rna_seasonality_metadata)
colnames(rna_seasonality) = colnames(t_rna_seasonality_metadata)
rna_seasonality_all = rbind(t_rna_seasonality_metadata, rna_seasonality)
write.csv(rna_seasonality_all, file = "rna_seasonality_all.csv")


####### Filter out smaples that dont have labels  #######
#gut_taxa_seasonality_all["IRIS"=="Unknown",]
remov = which(rna_seasonality_all["IRIS",] == "Unknown" | is.na(rna_seasonality_all["IRIS",]) | rna_seasonality_all["Date",] == "Unknown")
rna_seasonality_all = rna_seasonality_all[,-remov]
dim(rna_seasonality_all)





rna_count = rna_seasonality_all[-c(1:8),]
rna_count_table = data.frame(sapply(rna_count, function(x) as.integer(as.character(x))))
rownames(rna_count_table) = rownames(rna_count)
colnames(rna_count_table) = colnames(rna_count)

rna_metadata_table = as.data.frame(t(rna_seasonality_all[c(1:8),]))
colnames(rna_metadata_table)[1] = "Subject"
colnames(rna_metadata_table)[4] = "Group"
colnames(rna_metadata_table)[3] = "Time"

rna_metadata_table[, "Time"] = gsub("Jan", "01", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Feb", "02", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Mar", "03", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Apr", "04", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("May", "05", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Jun", "06", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Jul", "07", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Aug", "08", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Sep", "09", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Oct", "10", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Nov", "11", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = gsub("Dec", "12", rna_metadata_table[, "Time"])
rna_metadata_table[, "Time"]  = yday(as.Date(rna_metadata_table[, "Time"] , format = '%d-%m-%Y'))


View(rna_metadata_table)
View(rna_count_table)




##### DESEQ
### TODO: Filter data based on gene count, DESEQ results, or Mixed effect model
cts = rna_count_table + 1
coldata = rna_metadata_table

dds <- DESeqDataSetFromMatrix(countData = cts[1:500,],
                              colData = coldata,
                              design= ~ Group)


dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst


dds <- DESeq(dds)

save(dds, file = "ipop_seasonal_deseq_dds.RData")

resultsNames(dds) # lists the coefficients
res <- results(dds, name="GroupIR")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")



## Mixed-effect model




#### OMICSLONDA
rna_omicslondaAll = omicslondaAll(formula = Count ~ Time + state1, countTable = rna_count_table, 
                                  metadata = rna_metadata_table, 
                                  n.perm = 100, fit.method = "ssnbinomial", num.intervals = 36,
                                  parall = FALSE, pvalue.threshold = 0.05,     
                                  adjust.method = "none", time.unit = "day", norm.method = "none", 
                                  prefix = "Test_OmicsLonDA_rna_ipop_seasonal_all", 
                                  ylabel = "Normalized Count", col = c("black", "green"))



### remove features (rows) based on any filteration criteria
###rna_count2
#################################



rna_metalonda = metalondaAll(Count = rna_count2, Time = rna_time, Group = rna_group, 
                             ID = rna_id, fit.method = "lowess", n.perm = 1000, 
                             num.intervals = 36, parall = FALSE, pvalue.threshold = 0.05, 
                             adjust.method = "none", time.unit = "days", norm.method = "none",
                             prefix = "rna_ipop_seasonal", col = c("black", "green"))



### Time distribution
TimePoints = data.frame(subjectID = as.factor(paste("", rna_id, sep="")), days = rna_time, 
                        IRIS = rna_group)
TimePoints = TimePoints[order(TimePoints$IRIS, TimePoints$subjectID),]
TimePoints = data.frame(subjectID = factor(TimePoints$subjectID, levels = TimePoints$subjectID), Time = TimePoints$days, 
                        Group= as.factor(TimePoints$IRIS))


jpeg("rna_timepoints_distribution.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(TimePoints, aes(x = subjectID, y = Time)) +
  theme_bw() + 
  scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = Group), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  scale_y_continuous(name= "Days", breaks = round(seq(0, max(TimePoints$Time) + 20, by = 100),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=15, face="bold"), legend.title = element_blank()) +
  coord_flip()
# h = recordPlot(load="gridGraphics", attach=NULL)
dev.off()





