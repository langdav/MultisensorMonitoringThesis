## remove duplicate entries; replace "," with "." and convert temperature and magnitude to numeric
# files <- list.files("data/raw/TreeSense", pattern = ".csv")
# for(filename in files){
#   file <- read.csv(file.path("data/raw/TreeSense", filename))
#   file <- file[!duplicated(file$time),]
#   file$time <- as.POSIXct(file$time)
#   file$temperature <- as.numeric(gsub(",", ".", file$temperature))
#   file$magnitude <- as.numeric(gsub(",", ".", file$magnitude))
#   colnames(file) <- c("datetime", "temp_in_degrees_c", "magnitude_in_kohm", "battery_in_percent", "solar_in_percent") #rename coloumns
#   write.csv(file, paste0("data/TreeSense/", filename), row.names = F)
#}
#################################################
library(ggplot2);library(dplyr);library(tidyr);library(data.table)

## read in camera positions
pos <- read.csv("data/TreeSense/treesense_camera_positions.csv")

for(id in unique(pos$Tree)){
  numbers <- pos %>% filter(Tree %in% id) %>% pull(Sensor)
  
  #files <- list()
  files <- c()
  for(i in numbers){
    tmp <- read.csv(paste0("data/TreeSense/marburg_", i, ".csv"))
    tmp$datetime <- as.Date(tmp$datetime)
    assign(paste0("file_",i),tmp)
    #files <- append(files,paste0("file_",i))
    files <- c(files,as.name(paste0("file_",i)))
  }
  
  if(id == "BSF1") stuff <- list(file_1,file_2,file_3,file_4,file_5)
  else if(id == "CST0001") stuff <- list(file_6,file_7,file_8,file_9,file_10)
  else if(id == "CST0050") stuff <- list(file_11,file_12,file_13,file_14,file_15)
  
  treesense_single <- as.data.frame(data.table::rbindlist(stuff)[,lapply(.SD,mean), list(datetime)])
  treesense_single <- treesense_single[order(treesense_single$datetime),] #sort by date
  
  #add tree_id; and change the f*cking name, as those trees are named differently in every single document!!!!
  if(id == "CST0001"){
    treesense_single$tree_id <- rep("mof_cst_00001", nrow(treesense_single)) 
  } else if(id == "CST0050"){
    treesense_single$tree_id <- rep("mof_cst_00050", nrow(treesense_single)) 
  } else {
    treesense_single$tree_id <- rep("BSF_1", nrow(treesense_single)) 
  }
  
  if(!exists("treesense_all")){
    treesense_all <- treesense_single
  } else {
    treesense_all <- bind_rows(treesense_all, treesense_single)
  }
}

write.csv(treesense_all, "data/TreeSense/treesense_all.csv")

## plot
treesense_all %>%
  ggplot(aes(x=datetime, y=magnitude_in_kohm, color = tree_id)) +
  geom_line(size = 1) +
  ggtitle(id)

## load budburst data
budburst_long <- read.csv("data/budburst_data/budburst_long.csv")
budburst_long$date <- as.Date(budburst_long$date)
budburst_long <- budburst_long %>% filter(tree_id %in% treesense_all$tree_id)

## merge budburst_filled with sen1_results
treesense_budburst <- dplyr::left_join(treesense_all, budburst_long, by = c("datetime" = "date", "tree_id"))
treesense_budburst <- treesense_budburst %>% filter(datetime %in% seq.Date(min(budburst_long$date), max(budburst_long$date), by = "day"))
head(treesense_budburst)

## plots
treesense_budburst %>% 
  filter(tree_id %in% "BSF_1") %>% 
  ggplot(aes(x=datetime, y=magnitude_in_kohm)) +
  geom_point(aes(color=budburst)) +
  geom_smooth(method = "lm")

## per tree magnitudes; colored by budburst phase; linear model
treesense_budburst %>% 
  ggplot(aes(x=datetime, y=magnitude_in_kohm)) +
  geom_point(aes(color=budburst)) +
  geom_smooth(method = "lm") +
  facet_grid(tree_id~., scales = "free_y")

## boxplots; per tree distribution of magnitudes; grouped by budburst phase
treesense_budburst %>% 
  ggplot(aes(x=budburst, y=magnitude_in_kohm, color = tree_id)) +
  geom_boxplot(notch = F)

## boxplots; per budburst_perc distribution of magnitudes
treesense_budburst %>% 
  ggplot(aes(x=budburst_perc, y=magnitude_in_kohm, color = as.factor(budburst_perc))) +
  geom_boxplot(notch = F) +
  geom_jitter(aes(x=budburst_perc, y=magnitude_in_kohm))
