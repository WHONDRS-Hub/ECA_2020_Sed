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

# sed moisture
moisture = read.csv("//PNL/Projects/ECA_Project/Sediment_Collection_2020/SedimentMoisture/02_ProcessedData/merged_weights.csv")
moisture$Site = substring(moisture$ID_scanned,first = 1,last = 9)
head(moisture)

moist.sum = numeric()

for (curr.site in unique(moisture$Site)) {
  
  med.wet = median(moisture$percent_water_content_wet[which(moisture$Site == curr.site)])
  mean.wet = mean(moisture$percent_water_content_wet[which(moisture$Site == curr.site)])
  med.dry = median(moisture$percent_water_content_dry[which(moisture$Site == curr.site)])
  mean.dry = mean(moisture$percent_water_content_dry[which(moisture$Site == curr.site)])
  
  moist.sum = rbind(moist.sum,c(curr.site,med.wet,mean.wet,med.dry,mean.dry))
  
}

moist.sum = as.data.frame(moist.sum)
colnames(moist.sum) = c("Site","Med_Wet_Moist","Mean_Wet_Moist","Med_Dry_Moist","Mean_Dry_Moist")
head(moist.sum)
moist.sum$Site = as.character(moist.sum$Site)
moist.sum$Med_Wet_Moist = as.numeric(as.character(moist.sum$Med_Wet_Moist))
moist.sum$Mean_Wet_Moist = as.numeric(as.character(moist.sum$Mean_Wet_Moist))
moist.sum$Med_Dry_Moist = as.numeric(as.character(moist.sum$Med_Dry_Moist))
moist.sum$Mean_Dry_Moist = as.numeric(as.character(moist.sum$Mean_Dry_Moist))
str(moist.sum)

# merge everything
merged.data = merge(bnti,texture,by = 'Site')
merged.data = merge(merged.data,npoc,by = 'Site')
merged.data = merge(merged.data,beta.disp,by = 'Site')
merged.data = merge(merged.data,moist.sum,by = 'Site')
dim(merged.data)

## merge bnti with only moisture to maximize number of sites
bnti.moist = merge(x = bnti,y = moist.sum,by = 'Site'); dim(bnti.moist)
plot(bnti.moist$Median_bNTI ~ bnti.moist$Med_Dry_Moist)
summary(lm(bnti.moist$Median_bNTI ~ bnti.moist$Med_Dry_Moist))
plot(bnti.moist$Median_bNTI ~ bnti.moist$Med_Wet_Moist)
plot((bnti.moist$Mean_bNTI) ~ (bnti.moist$Med_Dry_Moist))
plot((bnti.moist$SD_bNTI) ~ (bnti.moist$Mean_Dry_Moist))
plot((bnti.moist$Min_bNTI) ~ (bnti.moist$Mean_Dry_Moist))
plot((bnti.moist$Max_bNTI) ~ (bnti.moist$Mean_Dry_Moist))
plot((bnti.moist$SD_bNTI) ~ (bnti.moist$Med_Dry_Moist))

# start trying constraint based regression
!! need to run through a loop of different number of break points and record the p.val and r.sq for each and make a plot of those parameters vs. number of break points
step.size = (max(bnti.moist$Med_Dry_Moist) - min(bnti.moist$Med_Dry_Moist))/10
start.val = min(bnti.moist$Med_Dry_Moist)
end.val = start.val + step.size
parsed.bnti = numeric()

for (i in 1:11) {
  
  bnti.temp = bnti.moist[which(bnti.moist$Med_Dry_Moist >= start.val & bnti.moist$Med_Dry_Moist < end.val),]
  if(nrow(bnti.temp) > 1) {
    max.val = max(bnti.temp$Median_bNTI)
    moist.val = bnti.temp$Med_Dry_Moist[which.max(bnti.temp$Median_bNTI)]
  }
  
  if(nrow(bnti.temp) == 1) {
    max.val = bnti.temp$Median_bNTI
    moist.val = bnti.temp$Med_Dry_Moist
  }
  
  if(nrow(bnti.temp) == 0) {
    max.val = NA
    moist.val = NA
  }
  
  parsed.bnti = rbind(parsed.bnti,c(max.val,moist.val))
  
  start.val = end.val
  end.val = start.val + step.size
  
}
colnames(parsed.bnti) = c("max.bnti","dry.moisture")
parsed.bnti = as.data.frame(parsed.bnti)
parsed.bnti = parsed.bnti[-which(is.na(parsed.bnti$max.bnti)),]

plot((bnti.moist$Median_bNTI) ~ (bnti.moist$Med_Dry_Moist))
mod.to.plot = parsed.bnti$max.bnti ~ parsed.bnti$dry.moisture
points(mod.to.plot,pch=19,col=4)
mod = summary(lm(mod.to.plot))
abline(mod,lwd=2,col=4)

# make bnti histogram
pdf("//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_bNTI_Outcomes/bNTI_MCD_Histogram.pdf")
par(pty="s")
hist(merged.data$Median_bNTI,xlim=c(-2.5,39),xlab=expression(paste("Within Site Median ",beta,"NTI",sep="")),cex.lab=2,cex.axis=1.5,main="")
abline(v=c(-2,2),lwd=2,col=2)
dev.off()


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

