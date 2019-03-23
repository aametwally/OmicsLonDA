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
library(compositions)
#library("biomformat"); packageVersion("biomformat")


setwd("/Users/ahmedmetwally/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")
#setwd("C:/Users/ametwall/Box Sync/Ahmed Metwally's Files/Stanford/OmicsLonDA_dev/")

set.seed(72626)

source("OmicsLonDA/R/omicslonda.R")
source("OmicsLonDA/R/omicslondaHelper.R")
source("OmicsLonDA/R/omicslondaMCPermutation.R")
source("OmicsLonDA/R/omicslondaVisualization.R")


#######################################
########### Read iPOP RData     ##########
#######################################
load("ipop_data/Revision_MultiOmes.RData", envir = parent.frame(), verbose = FALSE)
load("ipop_data/HMP_Microbiome_Count.RData")
#######################################
########### Microbiome     ##########
#######################################

### TODO: select annotated microbiome only
### TODO: Solve microbiome problem with the negative values after CLR


## Select abundant features from the Count matrix that used in iPOP analysis
ipop_microbiome = st.df
rownames(ipop_microbiome) = ipop_microbiome$HostSampleID
ipop_microbiome = ipop_microbiome[,-1]
ipop_microbiome[, c("collection.date", "CL1", "CL2", "CL3", "CL4", "SubjectID", "CollectionDate")] = list(NULL)
ipop_microbiome_selected = colnames(ipop_microbiome)


rownames(r16s.st.count.df) = r16s.st.count.df$HostSampleID
r16s.st.count.df = r16s.st.count.df[,-1]

st_selected = r16s.st.count.df[,ipop_microbiome_selected]
st_selected = st_selected + 1 ## add pesoudo count for CLR


## CLR transformation for each taxonomic
# CLR on the Phylum level
phylum_cols = which(substring(colnames(st_selected), 1, nchar("phylum")) == "phylum")
phylum_features = st_selected[,phylum_cols]
phylum_features_clr = clr(phylum_features)

# CLR on the Class level and then rmcorr
class_cols = which(substring(colnames(st_selected), 1, nchar("class")) == "class")
class_features = st_selected[,class_cols]
class_features_clr = clr(class_features)


# CLR on the Order level and then rmcorr
order_cols = which(substring(colnames(st_selected), 1, nchar("order")) == "order")
order_features = st_selected[,order_cols]
order_features_clr = clr(order_features)


# CLR on the Family level and then rmcorr
family_cols = which(substring(colnames(st_selected), 1, nchar("family")) == "family")
family_features = st_selected[,family_cols]
family_features_clr = clr(family_features)

# CLR on the Genus level and then rmcorr
genus_cols = which(substring(colnames(st_selected), 1, nchar("genus")) == "genus")
genus_features = st_selected[,genus_cols]
genus_features_clr = clr(genus_features)

### Combine each taxonomic level correlations
ipop_st_clr_abundant = cbind(phylum_features_clr, class_features_clr, order_features_clr, family_features_clr, genus_features_clr)
save(ipop_st_clr_abundant, file = "ipop_data/ipop_st_clr_abundant.RData")




#######################################
########### Metadata     ##########
#######################################
ipop_metadata = ls
ipop_metadata$status = ipop_metadata$CL3
ipop_metadata$status[which(ipop_metadata$status == "")] = ipop_metadata$CL4[which(ipop_metadata$status == "")]
colnames(ipop_metadata)[1] = "Subject"
colnames(ipop_metadata)[which(colnames(ipop_metadata) == "CL1")] = "DaysAfterInfect"
ipop_metadata = ipop_metadata[-which(ipop_metadata$CollectionDate==""),]
ipop_metadata$CollectionDate =  as.Date(ipop_metadata$CollectionDate, format = '%m/%d/%y')
subject_info = sc
colnames(subject_info)[1] = "Subject"
ipop_metadata = merge(ipop_metadata, subject_info, by = "Subject")


## Select samples that have microbiome samples
ipop_metadata = subset(ipop_metadata, SampleID %in% rownames(ipop_st_clr_abundant))







####################################################
########### Extract infection cohort     ###########
####################################################
prefix = "Infection"
xx = subset(ipop_metadata, substring(status, 1, nchar(prefix)) == prefix)
xx$Subject = factor(xx$Subject)





## Select subjects with infection timepoints
subset.df = subset(ipop_metadata, Subject %in% xx$Subject)
subset.df$DaysAfterInfect = gsub("D", "", subset.df$DaysAfterInfect)

## Select Healthy and Infection time points
subset.df = subset(subset.df, substring(status, 1, nchar("Infection")) == "Infection" | substring(status, 1, nchar("Healthy")) == "Healthy")
subset.df = subset(subset.df, status != "Infection") ## Remove infection samples without annotated infection days

# z=setdiff(subset.df$SampleID,subset.df_tmp$SampleID)
# subset(subset.df, subset.df$SampleID %in% z)


## for each datapoint with a healthy status, search for the nearest infection timepoint 
for(i in 1:nrow(subset.df)){
  if(subset.df[i,]$status == "Healthy"){
    tmp = subset(subset.df, Subject == subset.df[i,]$Subject)
    print(tmp)
    a = subset(tmp, CollectionDate>subset.df[i,]$CollectionDate & substring(tmp$status, 1, nchar("Infection")) == "Infection")
    nearstTime = min(a$CollectionDate)
    subset.df[i,]$DaysAfterInfect = subset.df[i,]$CollectionDate - nearstTime
  }
}


## remove inf
remove = which(subset.df$DaysAfterInfect == "-Inf")
subset.df = subset.df[-remove,]

## name each infection cycle by the date of the first infection sample
subset.df$CycleName = ""   ### TODO: Name infection cycles
class(subset.df$CycleName) = "Date"
## Order subset.df
subset.df = subset.df[with(subset.df, order(Subject, CollectionDate)),]
write.csv(subset.df, file = "ipop_data/ipop_infection_microbiome.csv")


### TODO: rename the cycle of the block that don't have Healthy and consecutive to infection timepoints
### TODO: This is 2011-12-04
### TODO: sort df by date
tmp_date = subset.df[1,]$CollectionDate
tmp_subj = subset.df[1,]$Subject
for(i in 1:nrow(subset.df)){
  if(subset.df[i,]$status == "Healthy"){
    tmp = subset(subset.df, Subject == subset.df[i,]$Subject)
    #print(tmp)
    a = subset(tmp, CollectionDate>subset.df[i,]$CollectionDate & substring(tmp$status, 1, nchar("Infection")) == "Infection")
    nearstTime = min(a$CollectionDate)
    subset.df[i,]$DaysAfterInfect = subset.df[i,]$CollectionDate - nearstTime
    subset.df[i,]$CycleName = nearstTime
    tmp_date = nearstTime
    tmp_subj = subset.df[i,]$Subject
  }
  else{
    tmp = subset(subset.df, Subject == subset.df[i,]$Subject)
    ## handle infection without prior healthy time points. actually we dont care. since those would be assigned NA in the CycleName by the followong line and we would remove them afterwards 
    a = subset(tmp, CollectionDate <= subset.df[i,]$CollectionDate & CollectionDate >= tmp_date & subset.df[i,]$Subject == tmp_subj & substring(tmp$status, 1, nchar("Infection")) == "Infection")
    nearstTime = min(a$CollectionDate)
    print(as.character(nearstTime))
    subset.df[i,]$CycleName = nearstTime
  }
}



## TODO: remove infection samples without Healthy baseline
write.csv(subset.df, file = "ipop_data/ipop_infection_microbiome_3.csv")


### Read modified file
#samples_infection = read.csv("ipop_data/ipop_metabolites_infection_wCycles.csv", row.names = 1)
samples_infection = subset.df
samples_infection$Subject_Cycle = paste(samples_infection$Subject, as.character(samples_infection$CycleName), sep = "_")
class(samples_infection$DaysAfterInfect) = "numeric"



### TODO: Remove sampels without Healthy group
samples_infection_filtered = samples_infection
for(i in unique(samples_infection$Subject_Cycle)){
  block = samples_infection[which(samples_infection$Subject_Cycle == i), "status"]
  if(!any(block == "Healthy")){
    samples_infection_filtered = samples_infection_filtered[-which(samples_infection_filtered$Subject_Cycle == i), ]
  }
}




#### Remove all healthy timepoints except the last one. Healthy timepoint has to exist 100 days before infection
dat = samples_infection_filtered
for(i in unique(dat$Subject_Cycle)){
  print(i)
  block = dat[which(dat$Subject_Cycle == i & dat$status == "Healthy"),]
  nearest = max(block$CollectionDate)
  if(block[which(block$CollectionDate == nearest),]$DaysAfterInfect < -100){
    dat = dat[-which(dat$Subject_Cycle == i),]
  } else {
    remove = block[-which(block$CollectionDate == nearest), ]
    print(dim(remove))
    if(dim(remove)[1]>0){
      dat = dat[-which(dat$SampleID %in% remove$SampleID), ]
    }
  }
}




##################################
### Visualization
##################################
jpeg("ipop_metabolites_infection_wCycles_timepoints_manyHealthyTimePoints.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(samples_infection_filtered, aes(x = Subject_Cycle, y = DaysAfterInfect, color = status)) +
  theme_bw() + 
  #scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = status), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="top")+ scale_shape_manual(values=c(1,6)) +
  #scale_y_continuous(name= "Days", breaks = round(seq(0, max(timepoints.df$DaysAfterInfect) + 20, by = 1),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=8,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=12, face="bold"), legend.title = element_blank()) +
  coord_flip()
dev.off()






jpeg("ipop_metabolites_infection_wCycles_timepoints_withOneHealthy.jpg", res = 300, height = 15, width = 20, units = 'cm')
ggplot(dat, aes(x = Subject_Cycle, y = DaysAfterInfect, color = status)) +
  theme_bw() + 
  #scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject (infection episode)")+ #coord_fixed(ratio = 0.002)+
  scale_y_continuous(name = "Days after Infection") + 
  geom_point(aes(colour = status, shape = Gender), stat='identity', size = 1, show.legend = TRUE) +
  theme(legend.position="right")+ scale_shape_manual(values=c(1,6)) +
  #scale_y_continuous(name= "Days", breaks = round(seq(0, max(timepoints.df$DaysAfterInfect) + 20, by = 1),1)) +
  theme(axis.text.x = element_text(colour="black",size=12,angle=45,hjust=1,vjust=1,face="plain"),
        axis.text.y = element_text(colour="black",size=6,angle=0,hjust=1,vjust=0.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0.5,face="bold"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.text=element_text(size=6, face="bold"), legend.title = element_blank()) +
  coord_flip()
dev.off()



##################################
## Stats
##################################
## Count number of unique infection cycle
length(unique(samples_infection_filtered$Subject_Cycle))
length(unique(dat$Subject_Cycle))
length(unique(dat$Subject))
table(dat$Gender)


## histogram of BMI and age 
s = unique(dat[, c("Subject", "Gender", "BMI", "Adj.age")])
hist(s$BMI)
hist(s$Adj.age)











##################################
### Select microbiome samples 
##################################
ipop_st_clr_abundant_selected = ipop_st_clr_abundant[dat$SampleID, ]

### May add psuedo count
ipop_st_clr_abundant_selected = ipop_st_clr_abundant_selected + 0.000000001

microbiome = list()
for(i in 1:ncol(ipop_st_clr_abundant_selected)){
  ## TODO: Check if the number matches by sample name
  print(i)
  microbiome[[i]] = cbind(dat, Count = ipop_st_clr_abundant_selected[,i], Feature = colnames(ipop_st_clr_abundant_selected)[i])
}


##################################
### OmicsLonDA
##################################

##############################################################################################
## Define omicslondaAll_ipop method
##############################################################################################
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
    cat ("Feature  = ", data[[i]]$Feature[1], "\n")
    omicslondaResults[[i]] = omicslonda(formula, df = data[[i]], n.perm = n.perm, fit.method = fit.method, points = points,
                                        text = unique(data[[i]]$Feature), parall = parall, pvalue.threshold = pvalue.threshold,     
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


##################################
### CLR Normalization
##################################
normalize = function(df){
  SubjectCycle = unique(df$Subject_Cycle)
  df$rawCount = df[, "Count"]
  df$normalizedCount = 0
  for(subjCyc in SubjectCycle)
  {
    max_healthy_tp = max(df[df$Subject_Cycle == subjCyc & df$status == "Healthy", "DaysAfterInfect"])
    df[df$Subject_Cycle == subjCyc, ]$normalizedCount = log(df[df$Subject_Cycle == subjCyc, ]$Count/df[df$Subject_Cycle == subjCyc & df$DaysAfterInfect == max_healthy_tp, ]$Count)
  }
  
  return(df)
}


microbiome_norm = list()
for(i in 1:length(microbiome)){
  microbiome_norm[[i]] = normalize(microbiome[[i]])
}

View(microbiome_norm[[1]])



## remove healthy timepoints
metabolite_norm_filtered = list()
for(i in 1:length(metabolite_norm)){
  metabolite_norm_filtered[[i]] = metabolite_norm[[i]][-which(metabolite_norm[[i]]$status == "Healthy"),]
}
View(metabolite_norm_filtered[[1]])



jpeg("ipop_metabolites_infection_wCycles_timepoints_normalized_fileterd_infectioOnly.jpg", res = 300, height = 20, width = 20, units = 'cm')
ggplot(metabolite_norm_filtered[[1]], aes(x = Subject_Cycle, y = DaysAfterInfect, color = status)) +
  theme_bw() + 
  #scale_colour_manual(values=c( "black", "green"),  breaks=c("IR", "IS")) +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_discrete(name ="Subject")+ #coord_fixed(ratio = 0.002)+
  geom_point(aes(colour = status, shape = Gender), stat='identity', size = 1, show.legend = TRUE) +
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
metabolite_norm_filtered_rename = metabolite_norm_filtered
for(i in 1:length(metabolite_norm_filtered_rename)){
  metabolite_norm_filtered_rename[[i]] = metabolite_norm_filtered_rename[[i]][, c("Subject_Cycle", "SampleID", "Gender", 
                                                                           "DaysAfterInfect", "Count", "rawCount", "normalizedCount", "Feature",
                                                                           "IRIS", "SSPG", "FPG", "Class", "Ethnicity", "Adj.age",
                                                                           "BMI")]
  colnames(metabolite_norm_filtered_rename[[i]]) = c("Subject", "ID", "Group", "Time", "Count", "rawCount", "normalizedCount", "Feature",
                                                     "IRIS", "SSPG", "FPG", "Class", "Ethnicity", "Adj.age",
                                                     "BMI")
  metabolite_norm_filtered_rename[[i]]$Subject = factor(metabolite_norm_filtered_rename[[i]]$Subject)
  metabolite_norm_filtered_rename[[i]]$Feature = as.character(metabolite_norm_filtered_rename[[i]]$Feature)
}

##################################
##### Global Testing & Test covariates
##################################
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




##################################
######### Run OmicsLonDA #########
##################################
points = seq(0, 50, length.out = 51)
data = metabolite_norm_filtered_rename
omicslondaResults = list()

### Test each feature
for(i in 1:length(metabolite_norm_filtered_rename))
{
  omicslondaResults[[i]] = omicslonda(formula = normalizedCount ~ Time, df = data[[i]], n.perm = 100, fit.method = "ssgaussian", points = points,
                                      text = as.character(unique(data[[i]]$Feature)), parall = FALSE, pvalue.threshold = 0.05,     
                                      adjust.method = "BH", col = c("pink", "black"), 
                                      prefix = paste("Feature_", unique(data[[i]]$Feature), sep = ""), 
                                      ylabel = "Adjusted Metabolite Level", 
                                      DrawTestStatDist = FALSE)
}


### Test all features
omicslonda_metabolomics = omicslondaAll_ipop(formula = normalizedCount ~ Time, data = metabolite_norm_filtered_rename, n.perm = 100, fit.method = "ssgaussian", 
                                   num.intervals = 100, parall = FALSE, pvalue.threshold = 0.05, 
                                   adjust.method = "BH", time.unit = "days", norm.method = "none", 
                                   prefix = "OmicsLonDA_ipop_Metabolomics",
                                   ylabel = "Adjusted Metabolite Level", col = c("deeppink2", "black"))



##################################
### Annotate metabolites
### TODO: Need to be revised
##################################
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






