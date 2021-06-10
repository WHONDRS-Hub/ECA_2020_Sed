### Merging null bMNTD matrices into an array and calculating bNTI
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com
# JCS 2021: James.Stegen@pnnl.gov

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013
# modified from https://github.com/danczakre/Meta-Metabolome_Ecology/blob/master/FTICR_bNTI_merge-nulls.R; Danczak et al. 2020

rm(list=ls());graphics.off()

# File Name Labels
tree_type = "MCD" # MCD or TW or TWCD

# Load in libraries/functions
library(picante)
library(reshape2)
library(abind)

acomb = function(...) abind(..., along = 3)


# ################## #
#### Data Loading ####
# ################## #

in.dir =  "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_Dendrograms/"
rand.in.dir = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_Randomizations/"
out.dir = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_bNTI_Outcomes/"

unique.sites = substr(list.files(path = in.dir,pattern = "Data.csv"),start = 1,stop = 9)

for (curr.site in unique.sites) {
  
  data = read.csv(paste0(in.dir,curr.site,"_FTICR_Data.csv"), row.names = 1) # Importing the site level data  
  tree = read.tree(paste0(in.dir,curr.site,"_MCD_UPGMA.tre")) # Importing the site specific dendrogram
  
  # ###################### #
  #### bNTI Calculation ####
  # ###################### #

  # Converting data to presence/absence
  data[data>1] = 1

  # Matching the tree to peak data
  # Any dropped tips indicate an error that needs to be fixed
  phylo = match.phylo.data(tree, data)

  # Calculating bMNTD for my samples
  coph = cophenetic(phylo$phy)

  bMNTD = as.matrix(comdistnt(t(phylo$data), coph, abundance.weighted = F, exclude.conspecifics = F))

  # Merging the separate bMNTD files
  files = list.files(path = paste(rand.in.dir,curr.site,"_",tree_type,"_Null_Results/", sep = "")
                   , pattern = "bMNTD_rep", full.names = T) # Listing files

  rand.bMNTD = NULL # Dummy object

  for(curr.file in files){
    temp = as.data.frame(read.csv(curr.file, row.names = 1))
    rand.bMNTD = c(rand.bMNTD, list(temp))
  } # Merging individual 

  rand.bMNTD = do.call(acomb, rand.bMNTD)
  rm("curr.file")

  # Calculate bNTI
  bNTI = matrix(c(NA), nrow = ncol(rand.bMNTD[,,1]), ncol = ncol(rand.bMNTD[,,1]))

  for(i in 1:(ncol(rand.bMNTD[,,1])-1)){
   for(j in (i+1):ncol(rand.bMNTD[,,1])){
      m = rand.bMNTD[j,i,] # Just setting all the randomizations for a given comparison to a matrix
     bNTI[j,i] = ((bMNTD[j,i]-mean(m))/sd(m)) # The bNTI calculation
    }
  }

  dimnames(bNTI) = dimnames(rand.bMNTD[,,1])
  rm("m", "j")

  write.csv(bNTI, paste(out.dir,curr.site, "_", tree_type, "_bNTI_", length(files), ".csv", sep = ""), quote = F)

}
