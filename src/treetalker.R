rm(list = ls())
library(plyr);library(dplyr);library(ggplot2)

tt_data <- readRDS("data/TreeTalker/TT_spectraldata.RDS")
tt_data$AS7263_610 <- as.numeric(tt_data$AS7263_610)
head(tt_data)

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

## mean values per date
tt_data$date <- as.Date(tt_data$timestamp)
tt_data_mean <- aggregate(tt_data[, 5:16], list(tt_data$date), mean)
tt_data_mean <- rename(tt_data_mean, "date" = "Group.1")

## some plotting
tt_data_mean %>% 
  ggplot(aes(x=date,y=AS7263_610, group=date)) +
  geom_point() +
  geom_smooth()


######
