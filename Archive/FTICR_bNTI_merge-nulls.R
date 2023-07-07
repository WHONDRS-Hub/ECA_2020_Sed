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
library(moments)

acomb = function(...) abind(..., along = 3)


# ################## #
#### Data Loading ####
# ################## #

in.dir =  "C:/Users/steg815/OneDrive - PNNL/Desktop/Temp_bNTI/MCD_Dendrograms/"
rand.in.dir = "C:/Users/steg815/OneDrive - PNNL/Desktop/Temp_bNTI/MCD_Randomizations/"
out.dir = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_bNTI_Outcomes/"

unique.sites = substr(list.files(path = in.dir,pattern = "Data.csv"),start = 1,stop = 9)

bnti.metrics.comp = numeric()

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
  write.csv(bMNTD, paste(out.dir,curr.site, "_", tree_type, "_bMNTD.csv", sep = ""), quote = F)
  
  jaccard = as.matrix(vegdist(x = t(phylo$data),method = "jaccard",binary = T,diag = F,upper = F))
  write.csv(jaccard, paste(out.dir,curr.site, "_", tree_type, "_Jaccard.csv", sep = ""), quote = F)
  
  if (identical(rownames(jaccard),rownames(bMNTD)) == T & identical(colnames(jaccard),colnames(bMNTD)) == T) {
    
    pdf(paste(out.dir,curr.site, "_", tree_type, "_bMNTD_v_Jaccard.pdf", sep = ""))
    par(pty="s")
    mod.to.plot = as.dist(bMNTD) ~ as.dist(jaccard)
    plot(mod.to.plot,ylab=expression(paste(beta,"MNTD Distance",sep="")),xlab="Jaccard Dissimilarity",cex.lab=2,cex.axis=1.5,main=curr.site)
    mod.for.reg = as.dist(bMNTD) ~ as.dist(jaccard)
    mod.lm = summary(lm(mod.for.reg))
    abline(mod.lm,lwd=2,col=4)
    p.val = mod.lm$coefficients[2,4]
    r.sq = round(mod.lm$r.squared,digits = 2)
    slope = round(mod.lm$coefficients[1,2],digits=4)
    if (p.val > 0.0001) {
      mtext(text = paste(" p = ",round(p.val,digits = 3),sep=""),line = -3.5,adj = 0,side = 3)
    } else { mtext(text = " p < 0.0001",line = -3.5,adj = 0,side = 3) }
    mtext(text = substitute(paste(" ", R^2," = ", r.sq),list(r.sq=r.sq)),line = -2.25,adj = 0,side = 3)
    mtext(text = paste(" Slope = ",slope,sep=""),line = -1,adj = 0,side = 3)
    dev.off()
    
  } else{
    
    print("Error: bMNTD and Jaccard have different col/row names")
    print(curr.site)
    break()
    
  }
  
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

  bNTI.dist = as.dist(bNTI)
  
  pdf(paste(out.dir,curr.site, "_", tree_type, "_bNTI_Hist.pdf", sep = ""))
  par(pty="s")
  
  hist(bNTI.dist,main=curr.site,xlab=expression(paste(beta,"NTI",sep="")),cex.lab=2,cex.axis=1.5,
       xlim=c(min(c(-2.5,min(bNTI.dist)-2)),max(c(2.5,max(bNTI.dist)+2)))
       ); abline(v=c(-2,2),lwd=2,col=2)
  
  dev.off()

  var.sel = length(which(bNTI.dist > 2)) / length(bNTI.dist)
  hom.sel = length(which(bNTI.dist < I(-2))) / length(bNTI.dist)
  stoch = length(which(abs(bNTI.dist) < 2)) / length(bNTI.dist)
  proc.tot = var.sel + hom.sel + stoch
  
  bnti.metrics.comp = rbind(bnti.metrics.comp,c(
    curr.site,
    nrow(bNTI),
    min(bNTI.dist),
    max(bNTI.dist),
    mean(bNTI.dist),
    median(bNTI.dist),
    sd(bNTI.dist),
    skewness(bNTI.dist),
    kurtosis(bNTI.dist),
    var.sel,
    hom.sel,
    stoch,
    proc.tot,
    p.val,
    r.sq,
    slope
  ))
  
  print(curr.site)
  
}

colnames(bnti.metrics.comp) = c("Site",
"Sample_Number",
"Min_bNTI",
"Max_bNTI",
"Mean_bNTI",
"Median_bNTI",
"SD_bNTI",
"Skew_bNTI",
"Kurt_bNTI",
"Var_Selection",
"Hom_Selection",
"Stochasticity",
"Process_Tot_Check",
"bMNTD_Jacc_pval",
"bMNTD_Jacc_Rsq",
"bMNTD_Jacc_Slope")

bnti.metrics.comp = as.data.frame(bnti.metrics.comp)

write.csv(bnti.metrics.comp, paste(out.dir,tree_type, "_bNTI_and_regress_metrics.csv", sep = ""), quote = F,row.names = F)


