# combine data to generate figures for manuscript on MCD null modeling of FTICR data from ECA 2020 campaign

min.num.samples = 10

out.dir = "//PNL/Projects/ECA_Project/ECA_Manuscripts/OM_Null_Modeling/"

bnti.in.path = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_bNTI_Outcomes/"

bnti.files = list.files(path = bnti.in.path,pattern = "MCD_bNTI_999")

# compile bnti values across all sites to enable a density function

bnti.comp = numeric()

for (i in bnti.files) {
  
  bnti.temp = read.csv(file = paste0(bnti.in.path,i),row.names = 1) 
  
  if (nrow(bnti.temp) >= min.num.samples) {
    
    bnti.comp = c(bnti.comp,as.numeric(as.vector(as.dist(bnti.temp))))
    
  } else {
    
    print(c(nrow(bnti.temp),i))
    
  }
  
}

bnti.comp.den = density(x = bnti.comp)

# make one histogram for each site with the among site compiliation overlaid

pdf(paste0(out.dir,"All_Sites_bNTI_Hists.pdf"),width=8)

for (i in bnti.files) {
  
  bnti.temp = read.csv(file = paste0(bnti.in.path,i),row.names = 1) 
  
  if (nrow(bnti.temp) >= min.num.samples) {
    
    par(pty="s")
    hist(as.numeric(as.vector(as.dist(bnti.temp))),xlim=c(-2.5,max(bnti.comp)),cex.lab=2,cex.axis=1.5,main=paste0("Site ",substr(x = i,start = 8,stop = 9)),xlab=expression(paste(beta,"NTI")),ylab="Within-Site Observations")
    abline(v=c(-2,2),col=2,lty=2,lwd=2)
    par(new = TRUE) 
    plot(bnti.comp.den,axes=F,lwd=2,col=4,xlab="",ylab="",main="")
    box()
    axis(side = 4, at = pretty(range(bnti.comp.den$y)),cex.axis=1.5)      # Add second axis
    mtext("Among Site Density", side = 4, line = 3,cex=2) 
    
  } 
  
}

dev.off()

# make a plot with two panels for the main ms. 
# each panel has a histogram for the site with minimum bnti and maximum bnit (9 and 55)
# label a and b
# overlay the among site density

pdf(paste0(out.dir,"Main_Min_Max_Hists.pdf"),width=15)
par(mfrow=c(1,2),pty="s")
panel.label = "a "

for (i in bnti.files) {
  
  bnti.temp = read.csv(file = paste0(bnti.in.path,i),row.names = 1) 
  
  if (nrow(bnti.temp) >= min.num.samples & substr(x = i,start = 8,stop = 9) %in% c('09','55')) {
    
    hist(as.numeric(as.vector(as.dist(bnti.temp))),xlim=c(-2.5,max(bnti.comp)),cex.lab=2,cex.axis=1.5,main="",xlab=expression(paste(beta,"NTI")),ylab="Within-Site Observations")
    abline(v=c(-2,2),col=2,lty=2,lwd=2)
    par(new = TRUE) 
    plot(bnti.comp.den,axes=F,lwd=2,col=4,xlab="",ylab="",main="")
    box()
    axis(side = 4, at = pretty(range(bnti.comp.den$y)),cex.axis=1.5)      # Add second axis
    mtext("Among Site Density", side = 4, line = 3,cex=2) 
    mtext(text = panel.label,side = 3,line = -1.5,adj = 1,cex=2)
    
    panel.label = "b "
    
  } 
  
}

dev.off()


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
merged.data = merge(bnti,texture,by = 'Site',all = T)
merged.data = merge(merged.data,npoc,by = 'Site',all = T)
merged.data = merge(merged.data,beta.disp,by = 'Site',all = T)
merged.data = merge(merged.data,moist.sum,by = 'Site',all = T)
dim(merged.data)

## merge bnti with only moisture to maximize number of sites, but need at least the min sample
bnti.moist = merge(x = bnti,y = moist.sum,by = 'Site'); 
bnti.moist = bnti.moist[which(bnti.moist$Sample_Number >= min.num.samples),]
dim(bnti.moist)

# start trying constraint based regression
#!! need to run through a loop of different number of break points and record the p.val and r.sq for each and make a plot of those parameters vs. number of break points
num.of.breaks = 10

pdf(paste(out.dir,"MCD_Med_bNTI_v_Both_Moisture.pdf",sep=""),width=15)
par(pty="s",mfrow=c(1,2))

# do moisture on dry basis first
moist.var = "Med_Dry_Moist"
step.size = (max(bnti.moist[,moist.var]) - min(bnti.moist[,moist.var]))/num.of.breaks
start.val = min(bnti.moist[,moist.var])
end.val = start.val + step.size
parsed.bnti = numeric()

for (i in 1:(num.of.breaks + 1)) {
  
  bnti.temp = bnti.moist[which(bnti.moist[,moist.var] >= start.val & bnti.moist[,moist.var] < end.val),]
  if(nrow(bnti.temp) > 1) {
    max.val = max(bnti.temp$Median_bNTI)
    moist.val = bnti.temp[which.max(bnti.temp$Median_bNTI),moist.var]
    site.temp = bnti.temp$Site[which.max(bnti.temp$Median_bNTI)]
  }
  
  if(nrow(bnti.temp) == 1) {
    max.val = bnti.temp$Median_bNTI
    moist.val = bnti.temp[,moist.var]
    site.temp = bnti.temp$Site
  }
  
  if(nrow(bnti.temp) == 0) {
    max.val = NA
    moist.val = NA
    site.temp = NA
  }
  
  parsed.bnti = rbind(parsed.bnti,c(site.temp,max.val,moist.val))
  
  start.val = end.val
  end.val = start.val + step.size
  
}
colnames(parsed.bnti) = c("Site","max.bnti","moisture")
parsed.bnti = as.data.frame(parsed.bnti)
parsed.bnti$Site = as.character(parsed.bnti$Site)
parsed.bnti$max.bnti = as.numeric(as.character(parsed.bnti$max.bnti))
parsed.bnti$moisture = as.numeric(as.character(parsed.bnti$moisture))
parsed.bnti = parsed.bnti[-which(is.na(parsed.bnti$max.bnti)),]
dim(parsed.bnti)


plot((bnti.moist$Median_bNTI) ~ (bnti.moist[,moist.var]),ylab=expression(paste("Site-Level Median ",beta,"NTI",sep="")),xlab="Site-Level Median Moisture (per dry mass)",cex.lab=2,cex.axis=1.5)
mod.to.plot = parsed.bnti$max.bnti ~ parsed.bnti$moisture
points(mod.to.plot,pch=19,col=4)
mod.lm = summary(lm(mod.to.plot))
abline(mod.lm,lwd=2,col=4)
p.val = round(mod.lm$coefficients[2,4],digits = 2)
r.sq = round(mod.lm$r.squared,digits = 2)
mtext(text = paste("p = ",p.val," ",sep=""),line = -1.5,adj = 1,side = 3,cex=1.5)
mtext(text = substitute(paste(R^2," = ", r.sq," "), list(r.sq=r.sq)),line = -3,adj = 1,side = 3,cex=1.5)
mtext(text = " a",side = 3,line = -1.5,adj = 0,cex=2)


# do moisture on wet basis
moist.var = "Med_Wet_Moist"
step.size = (max(bnti.moist[,moist.var]) - min(bnti.moist[,moist.var]))/num.of.breaks
start.val = min(bnti.moist[,moist.var])
end.val = start.val + step.size
parsed.bnti = numeric()

for (i in 1:(num.of.breaks + 1)) {
  
  bnti.temp = bnti.moist[which(bnti.moist[,moist.var] >= start.val & bnti.moist[,moist.var] < end.val),]
  if(nrow(bnti.temp) > 1) {
    max.val = max(bnti.temp$Median_bNTI)
    moist.val = bnti.temp[which.max(bnti.temp$Median_bNTI),moist.var]
    site.temp = bnti.temp$Site[which.max(bnti.temp$Median_bNTI)]
  }
  
  if(nrow(bnti.temp) == 1) {
    max.val = bnti.temp$Median_bNTI
    moist.val = bnti.temp[,moist.var]
    site.temp = bnti.temp$Site
  }
  
  if(nrow(bnti.temp) == 0) {
    max.val = NA
    moist.val = NA
    site.temp = NA
  }
  
  parsed.bnti = rbind(parsed.bnti,c(site.temp,max.val,moist.val))
  
  start.val = end.val
  end.val = start.val + step.size
  
}
colnames(parsed.bnti) = c("Site","max.bnti","moisture")
parsed.bnti = as.data.frame(parsed.bnti)
parsed.bnti$Site = as.character(parsed.bnti$Site)
parsed.bnti$max.bnti = as.numeric(as.character(parsed.bnti$max.bnti))
parsed.bnti$moisture = as.numeric(as.character(parsed.bnti$moisture))
parsed.bnti = parsed.bnti[-which(is.na(parsed.bnti$max.bnti)),]
dim(parsed.bnti)


plot((bnti.moist$Median_bNTI) ~ (bnti.moist[,moist.var]),ylab=expression(paste("Site-Level Median ",beta,"NTI",sep="")),xlab="Site-Level Median Moisture (per wet mass)",cex.lab=2,cex.axis=1.5)
mod.to.plot = parsed.bnti$max.bnti ~ parsed.bnti$moisture
points(mod.to.plot,pch=19,col=4)
mod.lm = summary(lm(mod.to.plot))
#abline(mod.lm,lwd=2,col=4)
p.val = round(mod.lm$coefficients[2,4],digits = 2)
r.sq = round(mod.lm$r.squared,digits = 2)
mtext(text = paste("p = ",p.val," ",sep=""),line = -1.5,adj = 1,side = 3,cex=1.5)
mtext(text = substitute(paste(R^2," = ", r.sq," "), list(r.sq=r.sq)),line = -3,adj = 1,side = 3,cex=1.5)
mtext(text = " b",side = 3,line = -1.5,adj = 0,cex=2)

dev.off()

!! next step is making the histogram of the median values, code just below, just need to make sure it's good


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

