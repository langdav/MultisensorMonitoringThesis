#author: David Langenohl
#last modified: 18.08.2021
#description: create budburst phases data frame from phenological observations
#NOTE: A - Buds, dormant, winter aspect; B - Buds swell and get round; C -  Buds have swollen and have burst, D - Buds have burst, leaves are coiled, E - Light foliage, F - Full foliage
#NOTE: no budburst = A, B, C; budburst = D, E, F

rm(list = ls())
library(dplyr);library(ggplot2);library(viridis);library(tidyr)

budburst <- read.csv2("data/budburst_data/data_spring_phenology_mof_21.csv")
budburst$date <- as.Date(budburst$Date, format = "%d.%m.%Y")
budburst <- budburst %>% filter(budburst$Tree_ID %in% unique(budburst$Tree_ID)[-c(1,2)])
#budburst <- budburst %>% filter(Tree_ID %in% as.character(unique(sen1_results$tree_id)))
#budburst <- budburst %>% filter(Tree_ID %in% c("mof_cst_00001","mof_cst_00003","mof_cst_00006","mof_cst_00013","mof_cst_00032","mof_cst_00036","mof_cst_00050", "BSF_1"))
# A - Buds, dormant, winter aspect; B - Buds swell and get round; C -  Buds have swollen and have burst, D - Buds have burst, leaves are coiled, E - Light foliage, F - Full foliage
# budburst: D, E, F
# Classification for ggplots: 0% budburst, <50% budburst , <50% budburst, 100% foliage
budburst$budburst <- NA
budburst$budburst_perc <- NA


###########################################################################
## two budburst phases; budburst or no_budburst; percentage of budburst ##
#########################################################################
for(i in 1:nrow(budburst)){
  if(budburst$Phase.D[i] > 0 || budburst$Phase.E[i] > 0 || budburst$Phase.F[i] > 0){
    budburst$budburst[i] <- T
    budburst$budburst_perc[i] <- sum(budburst[i, c(6:8)])
  } else {
    budburst$budburst[i] <- F
    budburst$budburst_perc[i] <- 0
  }
}

#budburst$budburst_phase <- factor(budburst$budburst, levels = c("no_budburst","budburst"))
budburst <- budburst[,c("date", "Tree_ID", "budburst", "budburst_perc")]

## observations were made in irregular intervals; fill days between observations with the values of the last observation, so that there is a value for every day of sentinel data
for(i in unique(budburst$Tree_ID)){
  if(!exists("budburst_filled")){
    budburst_filled <- budburst %>%
      filter(Tree_ID %in% i) %>% 
      complete(date = seq.Date(min(date), max(date), by="day")) %>%
      fill("budburst", "budburst_perc", "Tree_ID")
  } else {
    budburst_filled <- rbind(budburst_filled,budburst %>%
                               filter(Tree_ID %in% i) %>% 
                               complete(date = seq.Date(min(date), max(date), by="day")) %>%
                               fill("budburst", "budburst_perc", "Tree_ID"))
  }
}

budburst_filled <- budburst_filled[,c("Tree_ID","date","budburst","budburst_perc")]
colnames(budburst_filled)[1] <- "tree_id"

write.csv(budburst_filled, "data/budburst_data/budburst_long.csv", row.names = F)

