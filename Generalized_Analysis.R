### Generalized analysis script to detect outliers
# This is a script that doesn't require any external data (i.e., transformations, trees, etc.)
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
#setwd("~/Documents/Other 48hr Sampling Campaigns/VGC_COL20/Jan21_PB_5ppm/")

main.dir = ("/Users/gara009/OneDrive - PNNL/Desktop/ECA_ALL_VGC/FTMS_Output/")
output.dir = paste0(main.dir,"QAQC_output/")
setwd(main.dir)
#setwd(easycsv::choose_dir()) # Set work directory

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


# ###################### #
#### Misc. Processing ####
# ###################### #

# Creating ggplot themes which will reduce code repition
hori_x_theme = theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())

vert_x_theme = theme_bw()+
  theme(text = element_text(size = 10),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())

vert_vgc_theme = theme_bw()+
  theme(text = element_text(size = 8),
        axis.text.x = element_text(colour = "black", angle = 90, hjust = 1, vjust = 1),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())

# Matching data to formulas, if desired
if(match.form == T){
  data = data[!(mol$MolForm %in% NA),]
  mol = mol[!(mol$MolForm %in% NA),]
}

# ########################### #
#### Molecular Information ####
# ########################### #

### Plotting molecular characteristics
# Generating data frame of molecular data
char = data.frame(Mass.mean = rep(NA, ncol(data)), NOSC.mean = NA, AI.mean = NA,
                  DBE.mean = NA, N.mean = NA, S.mean = NA, P.mean = NA, row.names = colnames(data), 
                  stringsAsFactors = F)

for(i in 1:ncol(data)){
  temp = data[which(data[,i] > 0), i, drop = F] # Need to keep names, looking at columns
  temp = mol[row.names(temp),]
  
  char$Mass.mean[i] = mean(as.numeric(row.names(temp)), na.rm = T) # Each of these lines is averaging a given characteristic across a sample
  char$NOSC.mean[i] = mean(temp$NOSC, na.rm = T)
  char$AI.mean[i] = mean(temp$AI, na.rm = T)
  char$DBE.mean[i] = mean(temp$DBE, na.rm = T)
  char$N.mean[i] = mean(temp$N, na.rm = T)
  char$S.mean[i] = mean(temp$S, na.rm = T)
  char$P.mean[i] = mean(temp$P, na.rm = T)
  
  rm("temp")
} # I'm not sure how to do this without the for-loop, but I'm simply just finding the mean for peak stats


# Now we can melt and plot the characteristics
char = melt(data = as.matrix(char)) # Melting data to plot using ggplot
char = rbind(char, data.frame(Var1 = colnames(data), Var2 = "Peak_Count", value = colSums(data), row.names = NULL)) # Adding peak counts as a variable

#print(
  ggplot(char, aes(x = Var1, y = value))+
    geom_bar(stat="identity") + xlab(NULL)+
    facet_grid(Var2~., scales = "free_y")+
    scale_fill_stata()+
    vert_vgc_theme
  ggsave(file=paste0(output.dir, "Chracteristics_",unique.sites[v],"_",Sys.Date(),".pdf"))
  #) # Plotting characteritics

### Plotting elemental composition by sample
el.data = data[!is.na(mol$MolForm),] # Creating data and mol objects peaks assigned a molecular formula
el.mol = mol[!is.na(mol$MolForm),]
el.comp = matrix(data = 0, nrow = ncol(el.data), ncol = length(unique(el.mol$El_comp)), dimnames = list(colnames(el.data), 
                                                                                                  unique(el.mol$El_comp))) # Creating a matrix to store el comp data
for(i in 1:nrow(el.comp)){
  temp = el.mol[which(el.data[,i] > 0),] # Mol data for a given sample
  
  for(j in 1:ncol(el.comp)){
    el.comp[i,j] = length(grep(colnames(el.comp)[j], temp$El_comp))
  }
} # Counting the number of times a given elemental composition appears in a dataset

el.comp = apply(el.comp, 1, function(x) (x/sum(x))*100) # Relative abundance transformation
el.comp = melt(as.matrix(el.comp))

#print(
  ggplot(el.comp, aes(x = Var2, y = value))+
    geom_bar(stat = "identity", aes(fill = Var1))+
    xlab(NULL) + ylab("Relative Abundance (%)")+
    vert_x_theme
  ggsave(file=paste0(output.dir, "Elemental_composition_",unique.sites[v],"_",Sys.Date(),".pdf"))
  
#) # Plotting stacked bar charts for elemental composition

# ########################### #
#### Multivariate Analysis ####
# ########################### #

### Plotting PCA
pca = prcomp(x = t(data)) # Calculating PCA
pca = as.data.frame(scores(pca)) # Converting to PCA scores in order to plot using ggplot

#print(
  ggplot(data = pca, aes(x = PC1, y = PC2))+
    geom_point() + geom_label(aes(label = row.names(pca)))+
    hori_x_theme
  ggsave(file=paste0(output.dir, "PCA_",unique.sites[v],"_",Sys.Date(),".pdf"))
#) # Plotting PCA

### Beta-diversity
# Creating distance matrix
dist = vegdist(x = t(data), method = "jaccard") # Using Jaccard for historical reasons (ICR data is often analyzed using it)

# Plotting a Jaccard heatmap
dist.melt = melt(as.matrix(dist))

#print(
  ggplot(data = dist.melt, aes(x = Var1, y = Var2, fill = value))+
    geom_tile() + scale_fill_gradient2(low = "gray100", mid = "gray80", high = "darkred", midpoint = 0.4)+
    xlab(NULL) + ylab(NULL)+
    vert_x_theme
  ggsave(file=paste0(output.dir, "Heat_map_",unique.sites[v],"_",Sys.Date(),".pdf"))
#)

# Plotting Jaccard NMDS
nms = metaMDS(dist, trymax = 1000) # Determining NMDS
nms = as.data.frame(scores(nms)) # Conveting to scores

#print(
  ggplot(data = nms, aes(x = NMDS1, y = NMDS2))+
    geom_point()+geom_label(aes(label = row.names(nms)))+
    scale_color_stata()+
    hori_x_theme
  ggsave(file=paste0(output.dir, "NMDS_",unique.sites[v],"_",Sys.Date(),".pdf"))
#) # Plotting NMS graph
  
    setwd(output.dir)
    infiles = list.files(pattern=unique.sites[v])
    outfile = paste0(output.dir,"merged_plots/",unique.sites[v],"_QAQC.pdf")
    # Combine pdfs by site
    pdf_combine(infiles, output = outfile)
    

}

