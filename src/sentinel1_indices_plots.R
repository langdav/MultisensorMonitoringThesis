############################################################################
## prepare data ###########################################################
##########################################################################
rm(list = ls())
library(dplyr);library(ggplot2);library(viridis);library(tidyr)
sen1_results <- readRDS("data/satellite_data/satellite_indices/sentinel1_indices.RDS")
sen1_results$date <- as.Date(sen1_results$date)
sen1_results <- sen1_results[,-c(3,4)]

# phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

## merge sen1_results with phenoclasses
sen1_phenoclasses_long <- merge(sen1_results, phenoclasses, by = c("tree_id","date"), all.x = F, all.y = F) #merge with phenoclasses
write.csv(sen1_phenoclasses_long, "out/sentinel1/indicesi_long_format_phenoclasses_sentinel1.csv", row.names = F)

##########################################################################################################################################
##  create colorblind-friendly palette; can be used with scale_fill_manual(values= cbPalette) or scale_color_manual(values= cbPalette) ##
##see http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/ and https://jfly.uni-koeln.de/color/ ##########################################
#######################################################################################################################################
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
greyPalette <- RColorBrewer::brewer.pal(5, "Greys")[-1] #as the first color is to light, it gets removed; this palette could be used for all ggplots


############################################################################
## create ggplots; timeseries of Sentinel1 indices; per tree comparisons ##
##########################################################################

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

####################################################
## create per tree timeseries of budburst phases ##
##################################################

budburst %>%
  ggplot(aes(x=Date, y=budburst_phase, group=Tree_ID, color=Tree_ID)) +
  geom_line(size = 1) +
  facet_grid(Tree_ID ~ ., scales = "free_y") +
  scale_x_date(date_breaks = "1 week") +
  ylab("Budburst Phases") +
  labs(color = "Tree ID") +
  scale_color_viridis(discrete = T)

ggsave(paste0("out/sentinel1/budburst_phases.png"), last_plot(), height = 10, width = 15)

############################################################################################
## per sentinel1-indice plots; comparison of all 5 trees; budburst phases as background ##
##########################################################################################
## add another coloumn needed for geom_rects
sen1_results <- sen1_results[order(sen1_results$tree_id),]
sen1_results$budtime <- c(sen1_results$date[-1], (tail(sen1_results$date, 1)+1))
for(i in unique(sen1_results$tree_id)){
  sen1_results$budtime[which(sen1_results$date == tail(sen1_results[which(sen1_results$tree_id == i),"date"],1))] <- tail(sen1_results[which(sen1_results$tree_id == i),"date"],1)
}

## for-loop to create plots for all sentinel1 indices
for(i in 5:40){
  out <- sen1_results %>%
    ggplot(aes(x=date, y=sen1_results[,i], group=tree_id, color=tree_id)) +
    geom_rect(aes(xmin=date, xmax=budtime, ymin = min(sen1_results[,i]), ymax=max(sen1_results[,i]), fill = budburst_phase), color = NA, alpha=0.5) +
    scale_fill_brewer(palette="Greys") + #Greys = RColorBrewer Code; alternative (good) colorspaces: YlGn, Greens
    geom_line(size = 1) +
    facet_grid(tree_id~., scales = "free_y") +
    scale_color_viridis(discrete = T) + 
    theme_light() +
    xlab("Date") +
    ylab(names(sen1_results)[i]) +
    labs(fill = "Budburst Phase") +
    labs(color = "Tree ID")
  ggsave(paste0("out/sentinel1/plots/sen1_indices/per_indice_timelines/",names(sen1_results)[i], ".png"), out, height = 15, width = 10)
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
  ggsave(paste0("out/sentinel1/plots/sen1_indices/per_indice_means/",names(indice_means)[i], ".png"), out, height = 15, width = 10)
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

