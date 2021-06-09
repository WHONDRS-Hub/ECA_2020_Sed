### Renaming East River Data

setwd("C:/Users/gara009/OneDrive - PNNL/Documents/ECA_all_samples_4-13-2021/Final_XML")
file.list = list.files(pattern = ".xml")

new.names = gsub("Ste.*ECA", "ECA", file.list)
new.names = gsub("BAT.*p", "p", new.names)
new.names = gsub("_1_01.*$", "", new.names)
new.names = gsub("_AlderInf_144SA_", "", new.names)
new.names = gsub("_AlderInf_144SA_", "", new.names)
new.names = gsub("03Mar21", "", new.names)
new.names = gsub("15Mar21", "", new.names)
new.names = gsub("_ReSPE", "", new.names) 


if(!dir.exists("Renamed_Files")){
  dir.create("Renamed_Files")
}

for(f in 1:length(file.list)){
  
  rename = new.names[f]
  
  if(file.exists(paste0("Renamed_Files/", rename, ".xml"))){
    rename = paste0(rename, "_rep2")
  }
  
  rename = paste0(rename, ".xml")
  
  file.copy(file.list[f], paste0("Renamed_Files/", rename))
}
