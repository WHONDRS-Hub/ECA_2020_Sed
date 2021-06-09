### Processing Formularity log files
# This step identifies poorly calibrated samples

require(easycsv)

Sample_Name = "ECA_Sediment_Extractions" # Sample name for output files

setwd(easycsv::choose_dir()) # Set work directory
log = readLines(list.files(pattern = "*.log$")) # Read in the log file from Formularity

### Parsing out the relevant data
output = data.frame(samples = log[grep("Calibration of", log)], r.squared = log[grep("linear_calibration_R\\^2", log)], 
                    distinct = log[grep("distinct", log)], points = log[grep("linear_calibration=calibration", log)])

### Cleaning up the data
# Sample name
output$samples = gsub("Calibration of ", "", output$samples)
output$samples = gsub(".xml", "", output$samples)

# R-squared
output$r.squared = as.numeric(gsub("linear_calibration_R\\^2=", "", output$r.squared))

# Distinct calibrants
output$distinct = substring(output$distinct, regexpr("\\(", output$distinct))
output$distinct = as.numeric(gsub("\\(|distinct\\)", "", output$distinct)) # Two steps for distinct because the information is complicated

# Total and remove calibrants
output$points = gsub("linear.*count ", "", output$points) # Total calibrants
output$removed = as.numeric(gsub("^.*points ", "", output$points)) # Removed after interations
output$points = as.numeric(gsub(" remov.*$", "", output$points))

# Final number of calibrants (after removal)
output$final = output$points-output$removed

# Percent of calibrants remaining
output$remaining = output$final/output$points

### Determining bad calibrations
if(mean(output$r.squared) < 0.2){
  # This is the pattern that I observed for HJ Andrews data - again, this is empirical and subject to change
  bad = output[which(output$final < mean(output$final)*0.9),] # This doesn't catch PP48-12-3, so I needed to add that manually.
} else {
  # These patterns were determined empirically - I observed problems with this range of samples
  bad.1 = output[which(output$final <= 120),] # If the sample had 7 or fewer calibrants, its not great
  bad.4 = output[which(output$remaining <= 0.5),] # If more than half the calibrants needed to be removed, the calibration wasn't great
  
  bad = merge(bad.1, bad.4, all = T)
  rm("bad.1", "bad.2", "bad.3")
}

output$samples = gsub("_ReSPE", "", output$samples) 
### Writing data out
write.csv(output, paste(Sample_Name, "_All_Calibrations.csv", sep = ""), quote = F, row.names = F)
write.csv(bad, paste(Sample_Name, "_Poorly_Calibrated_Samples.csv", sep = ""), quote = F, row.names = F)

