# using this to select samples for NMR and GC-MS
# count number of samples in each site (wanting sites with 10 samples)
# quantify beta dispersion with goal of selecting sites spanning a range

rm(list=ls());graphics.off()

## Switches
filter.cal = F # Remove poorly calibrated samples (i.e., low QC)
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

# Processed ICR Data
data = read.csv(list.files(pattern = "*Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "*Mol.csv"), row.names = 1)
raw.data = data

samples = colnames(raw.data)[-1]
samples.site = substr(samples,1,9)
unique.sites = unique(samples.site)  

for (v in 1:length(unique.sites)){
  data = raw.data[grep(pattern = unique.sites[v], colnames(raw.data))]
  # Fixing column names if they begin with numbers
  if(length(grep("^X", colnames(data))) > 0){
    colnames(data) = gsub("^X", "", colnames(data))
  } # R begins integer column names with X's - this fixes that
  
  # Removing poorly calibrated samples, if set by the flag above
  if(filter.cal == T){
    poor.cal = read.csv(list.files(pattern = "*_Poorly_Calibrated_Samples.csv"))
    
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
  #if(identical(x = row.names(data), y = row.names(mol)) == FALSE){
  #stop("Something is incorrect: the mol. info and peak counts don't match")
  #}
  
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
  
  