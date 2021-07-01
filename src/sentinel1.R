############################################################################
## prepare data ##
##########################################################################
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

## merge budburst_filled with sen1_results
sen1_results <- dplyr::left_join(sen1_results, budburst_filled, by = c("date" = "Date", "tree_id" = "Tree_ID"))
#sen1_results$budburst_phase[is.na(sen1_results$budburst_phase)] <- 1

## add another coloumn needed for a nice ggplot 
sen1_results <- sen1_results[order(sen1_results$tree_id),]
sen1_results$budtime <- c(sen1_results$date[-1], (tail(sen1_results$date, 1)+1))
for(i in unique(sen1_results$tree_id)){
  sen1_results$budtime[which(sen1_results$date == tail(sen1_results[which(sen1_results$tree_id == i),"date"],1))] <- tail(sen1_results[which(sen1_results$tree_id == i),"date"],1)
}

############################################################################
## create ggplots; timeseries of Sentinel1 indices; per tree comparisons ##
##########################################################################
library(ggplot2);library(viridis);library(dplyr)

sen1_results <- readRDS("data/satellite_data/satellite_indices/sentinel1_indices.RDS")
sen1_results$date <- as.Date(sen1_results$date)
results <- names(sen1_results)[5:40]
per_tree_comparisons <- list()

for(i in results){
  out <- sen1_results %>%
    ggplot(aes(x=date, y=sen1_results[,i], group=tree_id, color=tree_id)) +
    geom_line(size = 1) +
    scale_color_viridis(discrete = T) + 
    theme_light()
  per_tree_comparisons[[i]] <- out
}


###########################################################################################
## create ggplots; per tree timeseries of Sentinel1 indices + time periods of budburst ##
#########################################################################################
library(dplyr);library(sf);library(raster)

sen1_results <- readRDS("data/satellite_data/satellite_indices/sentinel1_indices.RDS")
sen1_results$date <- as.Date(sen1_results$date)
products <- names(sen1_results)[-c(1:4)]

budburst <- read.csv2("data/data_spring_phenology_mof_21.csv")
budburst$Date <- as.Date(budburst$Date, format = "%d.%m.%Y")
budburst <- budburst %>% filter(Tree_ID %in% as.character(unique(sen1_results$tree_id)))
# A - Buds, dormant, winter aspect; B - Buds swell and get round; C -  Buds have swollen and have burst, D - Buds have burst, leaves are coiled, E - Light foliage, F - Full foliage
# Budburst: C, D 
# Classification for ggplots: 0% budburst, <50% budburst , <50% budburst, 100% foliage
budburst$class <- NA
for(i in 1:nrow(budburst)){
  if(budburst$Phase.A[i] == 100) budburst$class[i] <- 1 #0% budburst
  if(budburst$Phase.A[i] < 100 & budburst$Phase.A[i] > 50) budburst$class[i] <- 2 #<50% budburst
  if(budburst$Phase.A[i] <= 50) budburst$class[i] <- 3 #<50% budburst
  if(budburst$Phase.F[i] == 100) budburst$class[i] <- 4 #100% foliage
}

## create ggplots for every tree and every indices
for(i in as.character(unique(sen1_results$tree_id))){
  
  single_tree <- budburst %>% filter(Tree_ID %in% i)
  start <- single_tree$Date[!duplicated(single_tree$class)]
  end <- c(start[2:length(start)], tail(single_tree$Date, 1))
  rects <- data.frame(start=start, end=end, class=unique(single_tree$class))
  
  sen1_results_single_tree <- sen1_results %>% filter(tree_id %in% i)
  for(u in products){
    sen1_plot <- sen1_results_single_tree %>%
      ggplot(aes(x=date, y=sen1_results_single_tree[,u], group=tree_id, color=tree_id)) +
      geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(sen1_results_single_tree[,u]),
                                                   ymax=max(sen1_results_single_tree[,u]), group=class), color=NA, fill="gray35", alpha=c(0,0.4,0.6,0.8)) +
      geom_line(size = 1, color = "black") +
      theme_light()
    #alternative color scheme: c(NA, "chartreuse2", "green", "forestgreen")
    
    ggsave(paste0("out/sentinel1/per_tree_budburst_classes/",i, "_", u, ".png"), sen1_plot)
  }
}

####################################################
## create per tree timeseries of budburst phases ##
##################################################
library(dplyr);library(ggplot2);library(viridis)
sen1_results <- readRDS("data/satellite_data/satellite_indices/sentinel1_indices.RDS")
sen1_results$date <- as.Date(sen1_results$date)

budburst <- read.csv2("data/data_spring_phenology_mof_21.csv")
budburst$Date <- as.Date(budburst$Date, format = "%d.%m.%Y")
budburst <- budburst %>% filter(Tree_ID %in% as.character(unique(sen1_results$tree_id)))
# A - Buds, dormant, winter aspect; B - Buds swell and get round; C -  Buds have swollen and have burst, D - Buds have burst, leaves are coiled, E - Light foliage, F - Full foliage
# Budburst: C, D 
# Classification for ggplots: 0% budburst, <50% budburst , <50% budburst, 100% foliage
budburst$class <- NA
for(i in 1:nrow(budburst)){
  if(budburst$Phase.A[i] == 100) budburst$class[i] <- 1 #0% budburst
  if(budburst$Phase.A[i] < 100 & budburst$Phase.A[i] > 50) budburst$class[i] <- 2 #<50% budburst
  if(budburst$Phase.A[i] <= 50) budburst$class[i] <- 3 #<50% budburst
  if(budburst$Phase.F[i] == 100) budburst$class[i] <- 4 #100% foliage
}

budburst %>%
  ggplot(aes(x=Date, y=class, group=Tree_ID, color=Tree_ID)) +
  geom_line(size = 1) +
  facet_grid(Tree_ID ~ ., scales = "free_y") +
  scale_x_date(date_breaks = "1 week") +
  ylim("0% budburst","<50% budburst",">50% budburst","<100% foliage") +
  ylab("Budburst Phases") +
  labs(color = "Tree ID") +
  scale_color_viridis(discrete = T)

ggsave(paste0("out/budburst_phases.png"), last_plot(), height = 10, width = 15)

############################################################################################
## per sentinel1-indice plots; comparison of all 5 trees; budburst phases as background ##
##########################################################################################
## for-loop to create plots for all sentinel1 indices
for(i in 5:40){
  out <- sen1_results %>%
    ggplot(aes(x=date, y=sen1_results[,i], group=tree_id, color=tree_id)) +
    geom_rect( aes(xmin=date, xmax=budtime, ymin = min(sen1_results[,i]), ymax=max(sen1_results[,i]), fill = budburst_phase), color = NA, alpha=0.5) +
    scale_fill_brewer(palette="Greys") + #Greys = RColorBrewer Code; alternative (good) colorspaces: YlGn, Greens
    geom_line(size = 1) +
    facet_grid(tree_id~., scales = "free_y") +
    scale_color_viridis(discrete = T) + 
    theme_light() +
    xlab("Date") +
    ylab(names(sen1_results)[i]) +
    labs(fill = "Budburst Phase") +
    labs(color = "Tree ID")
  ggsave(paste0("out/sentinel1/sen1_indices/per_indice_timelines/",names(sen1_results)[i], ".png"), out, height = 15, width = 10)
}



############################################################################
## per sentinel1-indice plots; mean values of indices per budburst phase ##
##########################################################################
## create indice means per tree and budburst phase

for(id in as.character(unique(sen1_results$tree_id))){
  for(phase in unique(sen1_results$budburst_phase)[2:5]){
    means <- colMeans(Filter(is.numeric, sen1_results %>% 
                               filter(tree_id %in% id) %>% 
                               filter(budburst_phase %in% phase)))
    if(!exists("indice_means")){
      indice_means <- data.frame(id,phase,as.data.frame(t(means)))
    } else {
      indice_means <- rbind(indice_means, data.frame(id,phase,as.data.frame(t(means))))
    }
  }
}
rm("i", "id", "means", "phase")
indice_means$phase <- factor(indice_means$phase, levels = c("0% budburst","<50% budburst",">50% budburst","100% foliage"))

## crate colorblind-friendly palette; see http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ and https://jfly.uni-koeln.de/color/
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
greyPalette <- RColorBrewer::brewer.pal(5, "Greys")[-1] #as the first color is to light, it gets removed; this palette could be used for all ggplots

for(i in 3:38){
  out <- indice_means %>% ggplot(aes(x=phase, y=indice_means[,i])) + 
    geom_bar(stat = "identity", width=0.2, aes(fill = phase)) +
    scale_fill_manual(values= cbPalette) +
    #scale_fill_brewer(palette="Greys") + # as alternatvie to scale_fill_manual
    facet_grid(id~., scales = "free_y") +
    theme_light() +
    xlab("Budburst Phase") +
    ylab(names(indice_means)[i]) +
    labs(fill = "Budburst Phase") +
    labs(color = "Tree ID")
  ggsave(paste0("out/sentinel1/sen1_indices/per_indice_means/",names(indice_means)[i], ".png"), out, height = 15, width = 10)
}

#######################
## Misc: highlight budburst days as red dots on each line
#######################
highlight_max <- sen1_results %>% 
  filter(date %in% as.Date(c("2021-05-06","2021-04-17")))

sen1_results %>%
  ggplot(aes(x=date, y=sen1_results[,7], group=tree_id, color=tree_id)) +
  geom_line(size = 1) +
  geom_point(data=highlight_max, aes(x=date,y=highlight_max[,7]), color='red',size=3) +
  scale_color_viridis(discrete = T) + 
  theme_light()

