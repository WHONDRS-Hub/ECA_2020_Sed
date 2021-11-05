### Generating the molecular characteristics dendrogram
# RED 2020; robert.danczak@pnnl.gov; danczak.6@gmail.com
# edited from RED's script to generate a dendrogram for each ECA field site

rm(list=ls());graphics.off()

options(digits = 10)

require(phangorn) # For tree based functions
#require(ggtree) # For tree visualization
require(vegan) # For vegdist
require(picante)

#### function to turn factors to character

fac.to.char.fun = function(matrix.in) {
  
  i = sapply(matrix.in,is.factor)
  matrix.in[i] = lapply(matrix.in[i],as.character)
  
  return(matrix.in)
  
}

## Switches
filter.cal = T # Remove poorly calibrated samples (i.e., low QC)
clean.names = F # Switch to remove excess EMSL information
month.code = "Mar" # This works in conjuction with the clean.names switch - removes month information


# ################## #
#### Load in data ####
# ################## #

# Load in data
in.dir =  "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Processed_Data/"
out.dir = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_Dendrograms/"
poor.cal.dir = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Formularity_output/"

data = read.csv(paste0(in.dir,list.files(path = in.dir,pattern = "*Clean_Data.csv")), row.names = 1) # load in data file with samples
mol = read.csv(paste0(in.dir,list.files(path = in.dir,pattern = "*Clean_Mol.csv")), row.names = 1) # Load in molecular data

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

# Removing peaks that have no formula assignments
mol = mol[-which(mol$MolForm %in% NA),]

#################
# subset mol to be masses observed only in a given field site
# each field site will have it's own dendrogram
#################

samples = colnames(data)[-1]
samples.site = substr(samples,1,9)
unique.sites = unique(samples.site) 

pdf(paste(out.dir,"ECA_MCD_UPGMA_All_Sites.pdf", sep = ""))

for (i in unique.sites) {
  
  dat.temp = data[,grep(pattern = i,x = colnames(data))]
  zero.spp.temp = apply(dat.temp, 1 , sum); zero.spp.temp = names(zero.spp.temp)[which(zero.spp.temp == 0)]
  dat.temp = dat.temp[-which(rownames(dat.temp) %in% zero.spp.temp),]
  
  mol.temp = mol[which(rownames(mol) %in% rownames(dat.temp)),] # this only includes peaks that are present in samples in the site of interest
  dat.temp = dat.temp[which(rownames(dat.temp) %in% rownames(mol.temp)),]
  
  # Setting objects for useful parameters
  Mol.Info = mol.temp[,c("C", "H", "O", "N", "S", "P", "DBE", "AI_Mod", "kdefect.CH2"), drop = F]
  Mol.Ratio = mol.temp[,c("OtoC_ratio", "HtoC_ratio", "NtoC_ratio", "PtoC_ratio", "NtoP_ratio")]
  
  # ##################### #
  #### Generating MCD for the specific site ####
  # ##################### #

  # Pairwise distance between peaks
  Mol.Info = as.data.frame(apply(Mol.Info, 2, scale), row.names = row.names(Mol.Info)) # Generating a distance matrix based upon the provided parameters

  # Create tree
  tree = as.phylo(hclust(vegdist(Mol.Info, "euclidean"), "average")) # Converting the distance matrix into a tree

  # check to make sure dendrogram and data align. If anything prints indicating tips being dropped, something went wrong
  phylo.temp = match.phylo.data(tree, dat.temp)
  
  # Quick visualization of the dendrogram
  plot.phylo(tree,type = "fan",show.tip.label = F,show.node.label = F); mtext(text = i,side = 3)
 
  # Writing tre
  write.tree(tree, paste(out.dir,i, "_MCD_UPGMA.tre", sep = ""))
  
  # Writing dat.temp
  write.csv(dat.temp,paste(out.dir,i, "_FTICR_Data.csv", sep = ""),quote = F,row.names = T)
  
  rm('dat.temp','mol.temp','Mol.Info','tree','phylo.temp','Mol.Ratio','zero.spp.temp')

}

dev.off()
