# combine NPOC, betadispersion, and texture to help pick 5 sites to use for NMR analyses
# want high NPOC, fine texture, and maximal range in betadispersion
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

!! next step is to merge all the data and trim it down