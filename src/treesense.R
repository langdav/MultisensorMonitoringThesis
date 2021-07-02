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
#################################################
library(ggplot2);library(dplyr);library(tidyr);library(data.table)

## read in camera positions
pos <- read.csv("data/TreeSense/treesense_camera_positions.csv")
id <- unique(pos$Tree)[2]
numbers <- pos %>% filter(Tree %in% id) %>% pull(Sensor)

for(i in numbers){
  tmp <- read.csv(paste0("data/TreeSense/marburg_", i, ".csv"))
  tmp$datetime <- as.Date(tmp$datetime)
  assign(paste0("file_",i),tmp)
}

stuff <- list(file_6,file_7,file_8,file_9,file_10)
#stuff <- list(file_1,file_2,file_3,file_4,file_5)
#stuff <- list(file_11,file_12,file_13,file_14,file_15)
test <- as.data.frame(rbindlist(stuff)[,lapply(.SD,mean), list(datetime)])
test <- test[order(test$datetime),] #sort by date

#add tree_id; and change the f*cking name, as those trees are named differently in every single document!!!!
if(id == "CST0001"){
  test$tree_id <- rep("mof_cst_00001", nrow(test)) 
} else if(id == "CST0050"){
  test$tree_id <- rep("mof_cst_00050", nrow(test)) 
} else {
  test$tree_id <- rep("BSF_1", nrow(test)) 
}



test %>%
  ggplot(aes(x=datetime, y=magnitude_in_kohm)) +
  geom_line(size = 1) +
  ggtitle(id)




budburst <- read.csv2("data/data_spring_phenology_mof_21.csv")
budburst$Date <- as.Date(budburst$Date, format = "%d.%m.%Y")
budburst <- budburst %>% filter(Tree_ID %in% c("BSF_1", "mof_cst_00001", "mof_cst_00050"))
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
for(i in as.character(unique(budburst$Tree_ID))){
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

## merge budburst_filled with sen1_results
test2 <- dplyr::left_join(test, budburst_filled, by = c("datetime" = "Date", "tree_id" = "Tree_ID"))
head(test2)

########################
test %>% ggplot(aes(x=datetime, y=magnitude_in_kohm)) +
  geom_rect(aes(xmin=datetime, xmax=budtime, ymin = min(sen1_results[,i]), ymax=max(sen1_results[,i]), fill = budburst_phase), color = NA, alpha=0.5) +
  scale_fill_brewer(palette="Greys") + #Greys = RColorBrewer Code; alternative (good) colorspaces: YlGn, Greens
  geom_line(size = 1) +
  facet_grid(tree_id~., scales = "free_y") +
  scale_color_viridis(discrete = T) + 
  theme_light() +
  xlab("Date") +
  ylab(names(sen1_results)[i]) +
  labs(fill = "Budburst Phase") +
  labs(color = "Tree ID")