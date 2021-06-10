### Calculating null bMNTD variants
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013
# # JCS 2020: James.Stegen@pnnl.gov modified from https://github.com/danczakre/Meta-Metabolome_Ecology/blob/master/FTICR_bNTI_create-nulls.R; Danczak et al. 2020

rm(list=ls());graphics.off()

range = 1:50 # number of randomizations
Sample_Name = "Dataset_Name"
tree_type = "MCD" # MCD or TW or TWCD

#-----------------#

library(vegan)
library(picante)

###################################
#### Data Loading and cleaning ####
###################################

in.dir =  "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_Dendrograms/"
out.dir = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_Randomizations/"

unique.sites = substr(list.files(path = in.dir,pattern = "Data.csv"),start = 1,stop = 9)

for (curr.site in unique.sites) {

  data = read.csv(paste0(in.dir,curr.site,"_FTICR_Data.csv"), row.names = 1) # Importing the site level data  
  tree = read.tree(paste0(in.dir,curr.site,"_MCD_UPGMA.tre")) # Importing the site specific dendrogram

  # Creating necessary directories

  if(!dir.exists(paste0(out.dir,curr.site,"_",tree_type, "_Null_Results"))){
    dir.create(paste0(out.dir,curr.site,"_",tree_type, "_Null_Results"))
  }


  ####################################
  #### Beginning the bNTI Process ####
  ####################################

  # Converting to presence/absence
  data[data>1] = 1

  # Matching the tree to the newly rarefied OTU dataset
  # Any dropped tips indicate an error that needs to be fixed
  phylo = match.phylo.data(tree, data)

  # Running cophenetic outside of the for-loop
  coph = cophenetic(phylo$phy)

  # Calculating the bMNTD for the number of requested random distributions
  print(paste(date(), " - Start for loop"))

  for(i in range){
    bMNTD.rand = as.matrix(comdistnt(t(phylo$data), taxaShuffle(coph), abundance.weighted = F, exclude.conspecifics = F))
    write.csv(bMNTD.rand, paste(out.dir,curr.site,"_",tree_type, "_Null_Results/FTICR_",curr.site, "_", tree_type, "_bMNTD_rep", i, ".csv", sep = ""), quote = F)
    rm("bMNTD.rand")
  
    print(c(curr.site,date(),i))
  } # Performing the calculations using the data table but with randomized positions across the dendrogram

  print(paste(date(), " - End for loop"))

}
