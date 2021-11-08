# using this to select samples for NMR and GC-MS
# count number of samples in each site (wanting sites with 10 samples)
# quantify beta dispersion with goal of selecting sites spanning a range

rm(list=ls());graphics.off()

#### function to turn factors to character

fac.to.char.fun = function(matrix.in) {
  
  i = sapply(matrix.in,is.factor)
  matrix.in[i] = lapply(matrix.in[i],as.character)
  
  return(matrix.in)
  
}

## Switches
filter.cal = T # Remove poorly calibrated samples (i.e., low QC)
match.form = F # Only analyze data that have formulas
clean.names = F # Switch to remove excess EMSL information
month.code = "Mar" # This works in conjuction with the clean.names switch - removes month information


# ############################# #
#### Load required libraries ####
# ############################# #

require(vegan) # For broad ecology functions
library(pdftools) # to merge  pdfs
require(reshape2); require(ggplot2); require(ggthemes) # For prettier graphs

# ################## #
#### Load in data ####
# ################## #

# Set working directory
in.dir = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Processed_Data/"
out.dir = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/BetaDisp/"
poor.cal.dir = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Formularity_output/"

# Processed ICR Data
data = read.csv(paste0(in.dir,list.files(path = in.dir,pattern = "*Clean_Data.csv")), row.names = 1)
mol = read.csv(paste0(in.dir,list.files(path = in.dir,pattern = "*Clean_Mol.csv")), row.names = 1)

samples = colnames(data)[-1]
samples.site = substr(samples,1,9)
unique.sites = unique(samples.site)  

# Fixing column names if they begin with numbers
if(length(grep("^X", colnames(data))) > 0){
  colnames(data) = gsub("^X", "", colnames(data))
} # R begins integer column names with X's - this fixes that

# Removing poorly calibrated samples, if set by the flag above
if(filter.cal == T){
  poor.cal = read.csv(paste0(poor.cal.dir,list.files(path = poor.cal.dir,pattern = "*_Poorly_Calibrated_Samples.csv")))
  
  if(length(poor.cal[,1]) > 1){
    data = data[,-which(colnames(data) %in% gsub("-", ".", poor.cal$samples))]
  } else {
    stop("You've specified to remove poor calibrants, but there were none provided.")
  }
}


# #################### #
#### Error checking ####
# #################### #

# Checking row names consistency between molecular info and data
if(identical(x = row.names(data), y = row.names(mol)) == FALSE){
  stop("Something is incorrect: the mol. info and peak counts don't match")
}

# Checking to ensure "FREDA_Processing.R" was run
if(length(which(mol$C13 == 1)) > 0){
  stop("Isotopic signatures weren't removed")
}

if(length(grep("QC_SRFAII", colnames(data))) > 0){
  stop("Suwannee River standards are still in the data")
}

if(max(data) > 1){
  print("Data is not presence/absence")
  data[data > 1] = 1
}

if(clean.names){
  colnames(data) = gsub(paste0("[0-9]{2}", month.code, "[0-9]{2}.*p"), "p", colnames(data))
  colnames(data) = gsub("_1_01.*$", "", colnames(data))
}

## count up the number of samples per site
samples.per.site = numeric()
for (i in unique.sites) {
  
  samples.per.site = rbind(samples.per.site,c(i,length(grep(pattern = i,x = colnames(data)))))
  
}
samples.per.site = fac.to.char.fun(as.data.frame(samples.per.site))
colnames(samples.per.site) = c('Site',"Num.of.Samples")
samples.per.site$Num.of.Samples = as.numeric(samples.per.site$Num.of.Samples)

## calculate beta dispersion for each site
samples.per.site$betadisp = -999

for (i in unique.sites) {
  
  dat.temp = data[,grep(pattern = i,x = colnames(data))]
  zero.spp.temp = apply(dat.temp, 1 , sum); zero.spp.temp = names(zero.spp.temp)[which(zero.spp.temp == 0)]
  dat.temp = dat.temp[-which(rownames(dat.temp) %in% zero.spp.temp),]
  dist.temp = vegdist(x = t(dat.temp))
  beta.disp.temp = betadisper(d = dist.temp,group = factor(rep(1,ncol(dat.temp)),labels = "temp"))
  beta.disp.temp = mean(beta.disp.temp$distances)
  
  samples.per.site$betadisp[which(samples.per.site$Site == i)] = beta.disp.temp
  
}

beta.disp.use = samples.per.site[which(x = samples.per.site$Num.of.Samples > 0),]

write.csv(beta.disp.use,paste0(out.dir,"ECA2_FTICR_BetaDisp.csv"),quote = F,row.names = F)

# not sure what the goal of this was
#sites.to.use = c( # these are the 4 sites to use, selected to evenly span the range of betadisp
#beta.disp.use$Site[which.min(beta.disp.use$betadisp)],
#beta.disp.use$Site[which.max(beta.disp.use$betadisp)],
#beta.disp.use$Site[which.min(abs(beta.disp.use$betadisp - (min(beta.disp.use$betadisp) + (max(beta.disp.use$betadisp) - min(beta.disp.use$betadisp))/3)))],
#beta.disp.use$Site[which.min(abs(beta.disp.use$betadisp - (min(beta.disp.use$betadisp) + 2*(max(beta.disp.use$betadisp) - min(beta.disp.use$betadisp))/3)))])

#hist(beta.disp.use$betadisp)
#abline(v=c(beta.disp.use$betadisp[which(beta.disp.use$Site %in% sites.to.use)]),col=2)

  