# this is to compile preliminary NPOC and TN data from across batches into a single file
# that integrated file will be used to help select samples for NMR, etc.

setwd("//PNL/Projects/ECA_Project/Sediment_Collection_2020/Sediment_NPOC_TN/NPOC_TN_Direct_from_Instrument")

files = list.files(pattern = "Batch")

dat.comp = NULL

for (i in files) {
  
  dat.temp = read.csv(i,stringsAsFactors = F)
  dat.comp = rbind(dat.comp,cbind(dat.temp,strsplit(x = i,split = "_")[[1]][5]))
  
}

colnames(dat.comp)[ncol(dat.comp)] = "Batch"
head(dat.comp)

# remove DI
dat.comp = dat.comp[-grep(pattern = "_DI_",x = dat.comp$Name),]

# remove samples above the standard curves
dat.comp = dat.comp[-which(dat.comp$NPOC_mg_C_per_L > 50 | dat.comp$TN_mg_N_per_L > 3),]
range(dat.comp$NPOC_mg_C_per_L)
range(dat.comp$TN_mg_N_per_L)

# if none, that's good
dat.comp$Name[which(duplicated(x = dat.comp$Name)==T)]

write.csv(dat.comp,"ECA_NPOC.TN_Compiled_Direct_from_Instrument.csv",quote = F,row.names = F)

# define unique sites, should be 48
unique.sites = unique(substr(x = dat.comp$Name,start = 1,stop = 9))
length(unique.sites)

# put together summary stats
dat.summary = NULL

for (i in unique.sites) {
  
  npoc.temp = dat.comp$NPOC_after_correcting_for_dilution[grep(pattern = i,x = dat.comp$Name)]
  tn.temp = dat.comp$TN_after_correcting_for_dilution[grep(pattern = i,x = dat.comp$Name)]
  
  dat.summary = rbind(dat.summary,c(i,mean(npoc.temp),median(npoc.temp),var(npoc.temp),mean(tn.temp),median(tn.temp),var(tn.temp)))
  
}

colnames(dat.summary) = c("Site","Mean.NPOC","Med.NPOC","Var.NPOC","Mean.TN","Med.TN","Var.TN")
dat.summary[,2:ncol(dat.summary)] = as.numeric(dat.summary[,2:ncol(dat.summary)])
dat.summary = as.data.frame(dat.summary)
head(dat.summary)

write.csv(dat.summary,"ECA_NPOC.TN_Summary_Direct_from_Instrument.csv",quote = F,row.names = F)


