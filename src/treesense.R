## remove duplicate entries; replace "," with "." and convert temperature and magnitude to numeric
files <- list.files("data/raw/TreeSense", pattern = ".csv")
for(filename in files){
  file <- read.csv(file.path("data/raw/TreeSense", filename))
  file <- file[!duplicated(file$time),]
  file$time <- as.POSIXct(file$time)
  file$temperature <- as.numeric(gsub(",", ".", file$temperature))
  file$magnitude <- as.numeric(gsub(",", ".", file$magnitude))
  colnames(file) <- c("datetime", "temp_in_degrees_c", "magnitude_in_kohm", "battery_in_percent", "solar_in_percent") #rename coloumns
  write.csv(file, paste0("data/TreeSense/", filename), row.names = F)
}




