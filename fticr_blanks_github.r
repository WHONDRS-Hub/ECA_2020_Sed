rm(list=ls(all=T))
library(stringr); library(devtools);library("plyr")
library("readr"); library(tidyverse); library(readxl);library(crayon)

##################################################################

study = "ECA2"

min.DI.samples = 2 # minimum number of DI process blank samples a peak is observed in to be removed from the clean data. a value of 1 will remove any peak observed even once.

input.dir = "Processed_Data/"
dat = read.csv(paste0(input.dir,"Processed_ECA_Sediment_Extractions_Data.csv"),row.names = 1)
mol = read.csv(paste0(input.dir,"Processed_ECA_Sediment_Extractions_Mol.csv"),row.names = 1)

di.dat = dat[,grep(pattern = "DI",x = colnames(dat))]
  
real.dat = dat[,-grep(pattern = "DI",x = colnames(dat))]
  
# Making data to presence absence
di.dat[di.dat>0] <- 1 
real.dat[real.dat>0] <- 1 
  
# number of times a peak is seen in the DI process blanks
row.sum.di = as.data.frame(rowSums((di.dat)));colnames(row.sum.di)= "sum"

# now drop peaks based min.DI.samples
peaks.to.drop = rownames(row.sum.di)[which(row.sum.di$sum >= min.DI.samples)]
  
# now drop DI peaks from real data
data.clean = real.dat[-which(rownames(real.dat) %in% peaks.to.drop),]
  
# a check that things went well
print(c("This should be zero:",nrow(real.dat)-nrow(data.clean)-length(peaks.to.drop)))

# generating a clean mol file
mol.clean = mol[which(rownames(mol) %in% rownames(data.clean)),]

# checking that data.clean peaks are identical to mol.clean peaks
print(c("This should be TRUE:",identical(rownames(mol.clean),rownames(data.clean))))

#  temp = merge(row.sum.di,real.dat, by = 0)
  
#for (i in 3:ncol(temp)) {
#  for (j in 1:nrow(temp)){
#    if (temp$sum[j] >= min.DI.samples){
#      temp[j,i] = 0
#    }}
#  print(j)
#}
  
#rownames(temp) = temp$Row.names 
#data.clean = temp[,3:ncol(temp)]


# some diagnostics
length(peaks.to.drop) # number of peaks dropped

real.orig.peaks.per.sample = apply(X = real.dat,MARGIN = 2,FUN = sum)
hist(real.orig.peaks.per.sample)

real.clean.peaks.per.sample = apply(X = data.clean,MARGIN = 2,FUN = sum)
hist(real.clean.peaks.per.sample)

write.csv(mol.clean, paste0(input.dir,paste0("Processed_ECA2_",min.DI.samples,"_out_of_14_Clean_Mol.csv")))
write.csv(data.clean, paste0(input.dir,paste0("Processed_ECA2_",min.DI.samples,"_out_of_14_Clean_Data.csv")))



