rm(list=ls())

library(crayon)

background.color = "#d9d9d9" #to set background of plots
data.color = "#abd9e9" # to set color for data like bars in histogram
sig.lines.color = "#d73027" # to set color for vertical lines in plots
kern.den.color = "#fee090"  # to set color for kernal density lines
reg.dat.color = "#4575b4" # to set color for points in regressions

# combine data to generate figures for manuscript on MCD null modeling of FTICR data from ECA 2020 campaign

min.num.samples = 10

out.dir = "//PNL/Projects/ECA_Project/ECA_Manuscripts/OM_Null_Modeling/"

bnti.in.path = "//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/Null_Modeling/MCD_bNTI_Outcomes/"

bnti.files = list.files(path = bnti.in.path,pattern = "MCD_bNTI_999")

# generate list of sites that are used

sites.used = numeric()

for (i in bnti.files) {
  
  bnti.temp = read.csv(file = paste0(bnti.in.path,i),row.names = 1) 
  
  if (nrow(bnti.temp) >= min.num.samples) {
    
    sites.used = c(sites.used,substr(x = i,start = 1,stop = 9))
      
    }

}

sites.used = as.data.frame(sites.used)
colnames(sites.used) = "ECA_Site"
head(sites.used)

write.csv(sites.used,paste0(out.dir,"ECA2_Sites_Used.csv"),row.names = F,quote = F)

# compile bnti values across all sites to enable a density function

bnti.comp = numeric()
bnti.min.site = -999
bnti.max.site = -999
bnti.min.median = 999
bnti.max.median = -999

for (i in bnti.files) {
  
  bnti.temp = read.csv(file = paste0(bnti.in.path,i),row.names = 1) 
  
  if (nrow(bnti.temp) >= min.num.samples) {
    
    bnti.comp = c(bnti.comp,as.numeric(as.vector(as.dist(bnti.temp))))
    
    if (median(as.numeric(as.vector(as.dist(bnti.temp)))) < bnti.min.median ) {
      
      bnti.min.median = median(as.numeric(as.vector(as.dist(bnti.temp))))
      bnti.min.site = i
      
      print(c('min',i,bnti.min.median))
      
    }
    
    if (median(as.numeric(as.vector(as.dist(bnti.temp)))) > bnti.max.median ) {
      
      bnti.max.median = median(as.numeric(as.vector(as.dist(bnti.temp))))
      bnti.max.site = i
      
      print(c('max',i,bnti.max.median))
      
    }
    
    
  } else {
    
    print(c(nrow(bnti.temp),i))
    
  }
  
}

bnti.comp.den = density(x = bnti.comp,from = min(bnti.comp),to = max(bnti.comp))

# make one histogram for each site with the among site compiliation overlaid

pdf(paste0(out.dir,"SI_All_Sites_bNTI_Hists.pdf"),width=8)

for (i in bnti.files) {
  
  bnti.temp = read.csv(file = paste0(bnti.in.path,i),row.names = 1) 
  
  if (nrow(bnti.temp) >= min.num.samples) {
    
    par(pty="s")
    hist(as.numeric(as.vector(as.dist(bnti.temp))),xlim=c(floor(min(bnti.comp)),ceiling(max(bnti.comp))),cex.lab=2,cex.axis=1.5,main=paste0("Site ",substr(x = i,start = 8,stop = 9)),xlab=expression(paste(beta,"NTI")),ylab="Within-Site Observations")
    abline(v=floor(min(bnti.comp)),lwd=8000,col=background.color)
    par(new=T,pty="s")
    hist(as.numeric(as.vector(as.dist(bnti.temp))),xlim=c(floor(min(bnti.comp)),ceiling(max(bnti.comp))),cex.lab=2,cex.axis=1.5,main=paste0("Site ",substr(x = i,start = 8,stop = 9)),xlab=expression(paste(beta,"NTI")),ylab="Within-Site Observations",col=data.color)
    abline(v=c(-2,2),lwd=2,col=sig.lines.color,lty=2)
    axis(side = 1, at = c(floor(min(bnti.comp)),ceiling(max(bnti.comp))),cex.axis=1.5)
    par(new = TRUE) 
    plot(bnti.comp.den,axes=F,lwd=3,col=kern.den.color,xlab="",ylab="",main="")
    box()
    axis(side = 4, at = pretty(range(bnti.comp.den$y)),cex.axis=1.5)      # Add second axis
    mtext("Among Site Density", side = 4, line = 3,cex=2) 
    
  } 
  
}

dev.off()

# make a plot with two panels for the main ms. 
# each panel has a histogram for the site with minimum bnti and maximum bnit (9 and 41)
# label a and b
# overlay the among site density

pdf(paste0(out.dir,"Main_Min_Max_Hists.pdf"),width=15)
par(mfrow=c(1,2),pty="s")
panel.label = "a "

for (i in bnti.files) {
  
  bnti.temp = read.csv(file = paste0(bnti.in.path,i),row.names = 1) 
  
  if (nrow(bnti.temp) >= min.num.samples & substr(x = i,start = 8,stop = 9) %in% c('09','41')) {
    
    cat(red("Note that the min and max median sites are hard coded here, check to make sure the correct sites are called"))
    hist(as.numeric(as.vector(as.dist(bnti.temp))),xlim=c(floor(min(bnti.comp)),ceiling(max(bnti.comp))),cex.lab=2,cex.axis=1.5,main="",xlab=expression(paste(beta,"NTI")),ylab="Within-Site Observations")
    abline(v=floor(min(bnti.comp)),lwd=8000,col=background.color)
    par(new=T,pty="s")
    hist(as.numeric(as.vector(as.dist(bnti.temp))),xlim=c(floor(min(bnti.comp)),ceiling(max(bnti.comp))),cex.lab=2,cex.axis=1.5,main=paste0("Site ",substr(x = i,start = 8,stop = 9)),xlab=expression(paste(beta,"NTI")),ylab="Within-Site Observations",col=data.color)
    abline(v=c(-2,2),lwd=2,col=sig.lines.color,lty=2)
    axis(side = 1, at = c(floor(min(bnti.comp)),ceiling(max(bnti.comp))),cex.axis=1.5)
    par(new = TRUE) 
    plot(bnti.comp.den,axes=F,lwd=3,col=kern.den.color,xlab="",ylab="",main="")
    box()
    axis(side = 4, at = pretty(range(bnti.comp.den$y)),cex.axis=1.5)      # Add second axis
    mtext("Among Site Density", side = 4, line = 3,cex=2) 
    mtext(text = panel.label,side = 3,line = -1.5,adj = 1,cex=2)
    
    panel.label = "b "
    
  } 
  
}

dev.off()


# read in beta dispersion
beta.disp = read.csv("//PNL/Projects/ECA_Project/ECA_Sediment_Extraction_ICR_Data/BetaDisp/ECA2_FTICR_BetaDisp.csv",stringsAsFactors = F)
head(beta.disp)
str(beta.disp)

# read in npoc and TN
#npoc = read.csv("//PNL/Projects/ECA_Project/Sediment_Collection_2020/Sediment_NPOC_TN/NPOC_TN_Direct_from_Instrument/ECA_NPOC.TN_Summary_Direct_from_Instrument.csv",stringsAsFactors = F)
#head(npoc)
#str(npoc)

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
#merged.data = merge(merged.data,npoc,by = 'Site',all = T)
merged.data = merge(merged.data,beta.disp,by = 'Site',all = T)
merged.data = merge(merged.data,moist.sum,by = 'Site',all = T)
dim(merged.data)

## merge bnti with only moisture to maximize number of sites, but need at least the min sample
bnti.moist = merge(x = bnti,y = moist.sum,by = 'Site'); 
bnti.moist = bnti.moist[which(bnti.moist$Sample_Number >= min.num.samples),]
dim(bnti.moist)

# start trying constraint based regression
num.of.breaks = 10

pdf(paste(out.dir,"Main_Med_bNTI_v_Both_Moisture.pdf",sep=""),width=15)
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
points(mod.to.plot,pch=19,col=reg.dat.color)
mod.lm = summary(lm(mod.to.plot))
abline(mod.lm,lwd=2,col=reg.dat.color)
p.val = round(mod.lm$coefficients[2,4],digits = 2)
r.sq = round(mod.lm$r.squared,digits = 2)
mtext(text = paste("p = ",p.val," ",sep=""),line = -1.5,adj = 1,side = 3,cex=1.5)
mtext(text = substitute(paste(R^2," = ", r.sq," "), list(r.sq=r.sq)),line = -3,adj = 1,side = 3,cex=1.5)
mtext(text = "   a",side = 3,line = -1.5,adj = 0,cex=2)


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
points(mod.to.plot,pch=19,col=reg.dat.color)
mod.lm = summary(lm(mod.to.plot))
#abline(mod.lm,lwd=2,col=reg.dat.color)
p.val = round(mod.lm$coefficients[2,4],digits = 2)
r.sq = round(mod.lm$r.squared,digits = 2)
mtext(text = paste("p = ",p.val," ",sep=""),line = -1.5,adj = 1,side = 3,cex=1.5)
mtext(text = substitute(paste(R^2," = ", r.sq," "), list(r.sq=r.sq)),line = -3,adj = 1,side = 3,cex=1.5)
mtext(text = " b",side = 3,line = -1.5,adj = 0,cex=2)

dev.off()

######
# make histogram of bnti median values
pdf(paste(out.dir,"Main_Hist_of_Med_bNTI.pdf",sep=""),width = 15)
  par(pty="s",mfrow=c(1,2))
  hist(merged.data$Median_bNTI,xlim=c(floor(min(bnti.comp)),ceiling(max(bnti.comp))),xlab=expression(paste("Within Site Median ",beta,"NTI",sep="")),cex.lab=2,cex.axis=1.5,main="")
  abline(v=floor(min(bnti.comp)),lwd=8000,col=background.color)
  par(new=T,pty="s")
  hist(merged.data$Median_bNTI,xlim=c(floor(min(bnti.comp)),ceiling(max(bnti.comp))),xlab=expression(paste("Within Site Median ",beta,"NTI",sep="")),cex.lab=2,cex.axis=1.5,main="",col = data.color)
  abline(v=c(-2,2),lwd=2,col=sig.lines.color,lty=2)
  axis(side = 1, at = c(floor(min(bnti.comp)),ceiling(max(bnti.comp))),cex.axis=1.5)
  box()
dev.off()
#######


#############################
haven't edited below here

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

