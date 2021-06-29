library(sf);library(raster)

sen1_results <- readRDS("data/satellite_data/satellite_indices/sentinel1_indices.RDS")
sen1_pred <- readRDS("data/satellite_data/satellite_results/prediction/sentinel1_predictors.RDS")

plot(sen1_pred[[1]][[2]]$gamma_vv)

hist(sen1_results$gamma_vv)
hist(sen1_results$gamma_vh)

boxplot(sen1_results$gamma_vv)


tree1 <- sen1_results[which(sen1_results$tree_id == unique(sen1_results$tree_id)[1]),]
tree1_ts <- ts(tree1[,c(5:39)])
plot(x = as.Date(tree1$date), y = tree1_ts[,1], type = "l")

tree2 <- sen1_results[which(sen1_results$tree_id == unique(sen1_results$tree_id)[2]),]
tree2_ts <- ts(tree2[,c(5:39)])
plot(x = as.Date(tree2$date), y = tree2_ts[,1], type = "l")

par(mfrow=c(3, 2))
for(i in 1:length(unique(sen1_results$tree_id))){
  tree <- sen1_results[which(sen1_results$tree_id == unique(sen1_results$tree_id)[i]),]
  tree_ts <- ts(tree[,c(5:39)])
  plot(x = as.Date(tree$date), y = tree_ts[,1], type = "l")
}

as.character(sen1_results$date[which(sen1_results$tree_id == unique(sen1_results$tree_id)[1])])


############################################################################
## create ggplots; timeseries of Sentinel1 results; per tree comparisons ##
##########################################################################
library(ggplot2)
library(viridis)

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


## highlight budburst days as red dots on each line
highlight_max <- sen1_results %>% 
  filter(date %in% as.Date(c("2021-05-06","2021-04-17")))

sen1_results %>%
  ggplot(aes(x=date, y=sen1_results[,i], group=tree_id, color=tree_id)) +
  geom_line(size = 1) +
  geom_point(data=highlight_max, aes(x=date,y=highlight_max[,i]), color='red',size=3) +
  scale_color_viridis(discrete = T) + 
  theme_light()




## highlight region data
# define days of budburst
start <- as.Date(c("2021-04-26", "2021-05-03", "2021-05-12")) #test dates
end <- as.Date(c("2021-04-27", "2021-05-04", "2021-05-13")) #test dates
rects <- data.frame(start=start, end=end, group=seq_along(start))

sen1_results %>%
  ggplot(aes(x=date, y=sen1_results[,i], group=tree_id, color=tree_id)) +
  geom_line(size = 1) +
  geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(sen1_results[,i]),
                                               ymax=max(sen1_results[,i]), group=group), color="transparent", fill="orange", alpha=0.3) +
  scale_color_viridis(discrete = T) + 
  theme_light()





#########################
library(dplyr)
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

single_tree <- budburst %>% filter(Tree_ID %in% as.character(unique(sen1_results$tree_id)[1]))
start <- single_tree$Date[!duplicated(single_tree$class)]
end <- c(start[2:length(start)], tail(single_tree$Date, 1))
rects <- data.frame(start=start, end=end, class=unique(single_tree$class))

sen1_results_single_tree <- sen1_results %>% filter(tree_id %in% as.character(unique(sen1_results$tree_id)[2]))

sen1_results_single_tree %>%
  ggplot(aes(x=date, y=sen1_results_single_tree[,5], group=tree_id, color=tree_id)) +
  geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(sen1_results_single_tree[,5]),
                                               ymax=max(sen1_results_single_tree[,5]), group=class), color=NA, fill="gray35", alpha=c(0,0.4,0.6,0.8)) +
  geom_line(size = 1, color = "black") +
  theme_light()


#c(NA, "chartreuse2", "green", "forestgreen")

stuff <- names(sen1_results)[-c(1:4)]

for(i in as.character(unique(sen1_results$tree_id))){
  
  single_tree <- budburst %>% filter(Tree_ID %in% i)
  start <- single_tree$Date[!duplicated(single_tree$class)]
  end <- c(start[2:length(start)], tail(single_tree$Date, 1))
  rects <- data.frame(start=start, end=end, class=unique(single_tree$class))
  
  sen1_results_single_tree <- sen1_results %>% filter(tree_id %in% i)
  for(u in stuff){
    sen1_plot <- sen1_results_single_tree %>%
      ggplot(aes(x=date, y=sen1_results_single_tree[,u], group=tree_id, color=tree_id)) +
      geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(sen1_results_single_tree[,u]),
                                                   ymax=max(sen1_results_single_tree[,u]), group=class), color=NA, fill="gray35", alpha=c(0,0.4,0.6,0.8)) +
      geom_line(size = 1, color = "black") +
      theme_light()
    
    ggsave(paste0("out/sentinel1/per_tree_budburst_classes/",i, "_", u, ".png"), sen1_plot)
  }
}
