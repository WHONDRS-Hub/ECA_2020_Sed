# combine NPOC, betadispersion, bnti and texture to help pick 5 sites to use for NMR analyses
# want high NPOC, maybe fine texture, and maximal range in bnti and/or betadispersion
# only using sites with all 10 samples

out.dir = "//PNL/Projects/ECA_Project/Sediment_Collection_2020/ECA_NMR_2020_Sediments/Picking_NMR_Samples/"

# read in beta dispersion
beta.disp = read.csv("//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/BetaDisp/ECA_BetaDisp_10_Samp_Only.csv",stringsAsFactors = F)
head(beta.disp)
str(beta.disp)

# read in npoc and TN
npoc = read.csv("//PNL/Projects/ECA_Project/Sediment_Collection_2020/Sediment_NPOC_TN/NPOC_TN_Direct_from_Instrument/ECA_NPOC.TN_Summary_Direct_from_Instrument.csv",stringsAsFactors = F)
head(npoc)
str(npoc)

# read texture
texture = read.csv("//PNL/Projects/ECA_Project/Sediment_Collection_2020/Sediment_Texture/VGC_texture.csv",stringsAsFactors = F)
colnames(texture)[which(colnames(texture) == 'Kit_ID')] = "Site"
head(texture)
str(texture)

# read bnti metrics
bnti = read.csv("//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_bNTI_Outcomes/MCD_bNTI_and_regress_metrics.csv")
bnti$Site = as.character(bnti$Site)
head(bnti)
str(bnti)

# merge everything
merged.data = merge(bnti,texture,by = 'Site')
merged.data = merge(merged.data,npoc,by = 'Site')
merged.data = merge(merged.data,beta.disp,by = 'Site')
dim(merged.data)

hist(merged.data$Median_bNTI)

# selecting NMR samples. going for sets with high NPOC, but without losing main range in bNTI
# median of 40 didn't drop the overall range in bNTI
hist(merged.data$Median_bNTI[which(merged.data$Med.NPOC > 40)])

trim.dat = merged.data[which(merged.data$Med.NPOC > 40),]
hist(trim.dat$Median_bNTI)
nmr.sites = c('ECA2_0009','ECA2_0061','ECA2_0049','ECA2_0041','ECA2_0055') # took min, max, and evenly spaced in between
abline(v = c(trim.dat$Median_bNTI[which(trim.dat$Site %in% nmr.sites)]),col=2)

# selecting GC-MS samples. using the NMR samples and adding more
# dropping the lowest NPOC sets
hist(merged.data$Median_bNTI[which(merged.data$Med.NPOC > 10)])
gc.trim.dat = merged.data[which(merged.data$Med.NPOC > 10),]
gc.sites = c('ECA2_0009','ECA2_0061','ECA2_0049','ECA2_0041','ECA2_0055','ECA2_0056','ECA2_0002','ECA2_0001','ECA2_0047','ECA2_0036','ECA2_0043','ECA2_0028','ECA2_0013')
abline(v = c(gc.trim.dat$Median_bNTI[which(gc.trim.dat$Site %in% gc.sites)]),col=2)

# metaT sites are the same as for GC-MS
# metaP sites are the min and max for median bNTI ECA2_0009 and ECA2_0055

