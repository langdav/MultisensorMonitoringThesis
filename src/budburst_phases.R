library(dplyr);library(ggplot2);library(viridis);library(tidyr)

sen1_results <- readRDS("data/satellite_data/satellite_indices/sentinel1_indices.RDS")
sen1_results$date <- as.Date(sen1_results$date)

budburst <- read.csv2("data/data_spring_phenology_mof_21.csv")
budburst$Date <- as.Date(budburst$Date, format = "%d.%m.%Y")
budburst <- budburst %>% filter(Tree_ID %in% as.character(unique(sen1_results$tree_id)))
# A - Buds, dormant, winter aspect; B - Buds swell and get round; C -  Buds have swollen and have burst, D - Buds have burst, leaves are coiled, E - Light foliage, F - Full foliage
# Budburst: C, D 
# Classification for ggplots: 0% budburst, <50% budburst , <50% budburst, 100% foliage
budburst$budburst_phase <- NA
for(i in 1:nrow(budburst)){
  if(budburst$Phase.A[i] == 100) budburst$budburst_phase[i] <- "0% budburst"
  if(budburst$Phase.A[i] < 100 & budburst$Phase.A[i] > 50) budburst$budburst_phase[i] <-  "<50% budburst"
  if(budburst$Phase.A[i] <= 50) budburst$budburst_phase[i] <-  ">50% budburst"
  if(budburst$Phase.F[i] == 100) budburst$budburst_phase[i] <-  "100% foliage"
}
budburst$budburst_phase <- factor(budburst$budburst_phase, levels = c("0% budburst","<50% budburst",">50% budburst","100% foliage"))
budburst <- budburst[,c("Date", "Tree_ID", "budburst_phase")]

## observations were made in irregular intervals; fill days between observations with the values of the last observation, so that there is a value for every day of sentinel data
for(i in as.character(unique(sen1_results$tree_id))){
  if(!exists("budburst_filled")){
    budburst_filled <- budburst %>%
      filter(Tree_ID %in% i) %>% 
      complete(Date = seq.Date(min(Date), max(Date), by="day")) %>%
      fill("budburst_phase", "Tree_ID")
  } else {
    budburst_filled <- rbind(budburst_filled,budburst %>%
                               filter(Tree_ID %in% i) %>% 
                               complete(Date = seq.Date(min(Date), max(Date), by="day")) %>%
                               fill("budburst_phase", "Tree_ID"))
  }
}
