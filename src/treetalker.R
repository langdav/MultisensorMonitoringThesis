#author: David Langenohl
#last modified: 17.08.2021
#description: calculate NDVI values from TreeTalker spectra
#NOTE: NDVI = (nir-red)/(nir+red); nir = 860 nm, red = 650 nm

rm(list = ls())
library(plyr);library(dplyr);library(ggplot2)

tt_data <- readRDS("data/TreeTalker/TT_spectraldata.RDS")
tt_data$AS7263_610 <- as.numeric(tt_data$AS7263_610)
tt_data$TT_ID <- as.integer(tt_data$TT_ID)
#head(tt_data)

## remove entries < 5
tt_data <- tt_data %>% filter(tt_data[,c(5:16)] > 5)

## convert digital numbers to corrected spectral bands
# NIR
tt_data$AS7263_610 <- -312.45+(1.6699*tt_data$AS7263_610)
tt_data$AS7263_680 <- -561.56+(1.5199*tt_data$AS7263_680)
tt_data$AS7263_730 <- -1511.2+(1.6209*tt_data$AS7263_730)
tt_data$AS7263_760 <- -1012.5+(1.4549*tt_data$AS7263_760)
tt_data$AS7263_810 <- 91.58+(0.8414*tt_data$AS7263_810)
tt_data$AS7263_860 <- 3347.88+(0.531*tt_data$AS7263_860)

# visible light
tt_data$AS7262_450 <- -212.62+(0.4562*tt_data$AS7262_450)
tt_data$AS7262_500 <- -232.13+(0.6257*tt_data$AS7262_500)
tt_data$AS7262_550 <- -842.1+(1.0546*tt_data$AS7262_550)
tt_data$AS7262_570 <- -666.72+(1.0462*tt_data$AS7262_570)
tt_data$AS7262_600 <- -328.08+(0.8654*tt_data$AS7262_600)
tt_data$AS7262_650 <- 202.77+(0.7829*tt_data$AS7262_650)

## reduce to date range between 15.03. - 01.06.
#tt_data <- tt_data %>% filter(as.Date(timestamp) > as.Date("2021-03-14") & as.Date(timestamp) < as.Date("2021-06-02"))

## merge tree ids with tree talker data
tt_tree_ids <- read.csv("data/TreeTalker/treetalker_tree_ids.csv")
tt_data <- merge(tt_data, tt_tree_ids, by.x = "TT_ID", by.y = "tt_id")

## merge mean values per date with budburst data
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

tt_data$date <- as.Date(tt_data$timestamp)
tt_data <- tt_data[with(tt_data, order(tree_id,date)),]
tt_pheno <- merge(tt_data, phenoclasses, by = c("tree_id","date"), all.x = F, all.y = F)

## remove unnecessary rows and rearrange columns
tt_long_phenoclasses <- subset(tt_pheno,select = -c(record_number, device_type, integration_time, gain, TT_ID))

tt_long_phenoclasses <- tt_long_phenoclasses %>% relocate(timestamp, .after = date)

tt_long_phenoclasses <- tt_long_phenoclasses[order(tt_long_phenoclasses$tree_id),]
tt_long_phenoclasses <- tt_long_phenoclasses[order(tt_long_phenoclasses$tree_id, tt_long_phenoclasses$date),]

tt_long_phenoclasses <- tt_long_phenoclasses[with(tt_long_phenoclasses, order(tree_id,date)),]

## save
save(tt_long_phenoclasses, file = "data/treetalker/tt_phenoclasses_all_spectra.RData")

## calculate and add NDVI
#NDVI = (nir-red)/(nir+red); nir = 860 nm, red = 650 nm
tt_pheno$ndvi <- (tt_pheno$AS7263_860-tt_pheno$AS7262_650)/(tt_pheno$AS7263_860+tt_pheno$AS7262_650)

ndvi_long_format_phenoclasses_treetalker <- subset(tt_pheno, select = c(tree_id, date, timestamp, ndvi, budburst, budburst_perc))

## save
save(ndvi_long_format_phenoclasses_treetalker, file = "out/treetalker/ndvi_all_with_phenoclasses_treetalker.RData")

## remove outliers within each trees single-day-NDVI-values, based on the IQR
for(tree in unique(ndvi_long_format_phenoclasses_treetalker$tree_id)){
  cat("Processing", tree, "\n")
  days <- unique(ndvi_long_format_phenoclasses_treetalker$date[which(ndvi_long_format_phenoclasses_treetalker$tree_id == tree)])
  for(day in 1:length(days)){
    if(length(boxplot.stats(ndvi_long_format_phenoclasses_treetalker$ndvi[which(ndvi_long_format_phenoclasses_treetalker$tree_id == tree & 
                                                  ndvi_long_format_phenoclasses_treetalker$date == days[day])])$out) != 0){
      ndvi_long_format_phenoclasses_treetalker <- ndvi_long_format_phenoclasses_treetalker[-which(ndvi_long_format_phenoclasses_treetalker$tree_id == tree & 
                                        ndvi_long_format_phenoclasses_treetalker$date == days[day] & 
                                        ndvi_long_format_phenoclasses_treetalker$ndvi %in% c(boxplot.stats(ndvi_long_format_phenoclasses_treetalker$ndvi[which(ndvi_long_format_phenoclasses_treetalker$tree_id == tree & 
                                                                                                     ndvi_long_format_phenoclasses_treetalker$date == days[day])])$out)),]
    }
  }
}

## save
save(ndvi_long_format_phenoclasses_treetalker, file = "out/treetalker/outlier_free_ndvi_mean_per_tree_with_phenoclasses_treetalker.RData")
