rm(list = ls())
library(plyr);library(dplyr);library(ggplot2)

tt_data <- readRDS("data/TreeTalker/TT_spectraldata.RDS")
tt_data$AS7263_610 <- as.numeric(tt_data$AS7263_610)
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

## mean values per date
tt_data$date <- as.Date(tt_data$timestamp)
tt_data_mean <- aggregate(tt_data[, 5:16], list(tt_data$date), mean)
tt_data_mean <- rename(tt_data_mean, "date" = "Group.1")

## some plotting
tt_data_mean %>% 
  ggplot(aes(x=date,y=AS7263_610, group=date)) +
  geom_point() +
  geom_smooth()


## merge mean values per date with budburst data
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

tt_pheno <- merge(tt_data_mean, phenoclasses, by = "date", all.x = T, all.y = F)


## find and remove outliers
# boxplot for visual detection of outliers
boxplot(tt_data$AS7263_610)

# which values are the outliers (based on Inter Quantile Range: Q25 - 1.5 * IQR  <-> Q75 + 1.5*IQR; IQR = Q75 - Q25)
quantile(tt_data$AS7263_610, 0.25) #first quartile
quantile(tt_data$AS7263_610, 0.75) #third quartile
IQR(tt_data$AS7263_610) #innerquartile range

boxplot.stats(tt_data$AS7263_610)$out # values wich lay beyond the interquartile  range
out_ind <- which(tt_data$AS7263_610 %in% c(boxplot.stats(tt_data$AS7263_610)$out)) 
out_ind #outlier rows
tt_data$AS7263_610[-out_ind]

# alternative method: find outliers according to 2.5 and 97.5 percentiles
quantile(tt_data$AS7263_610, 0.025) #lower quantile; all values < this value are outliers
quantile(tt_data$AS7263_610, 0.975) #upper quantile; all values > this value are outliers

out_ind <- which(tt_data$AS7263_610 < quantile(tt_data$AS7263_610, 0.025) | tt_data$AS7263_610 > quantile(tt_data$AS7263_610, 0.975))
out_ind #outlier rows

# log - boxplot for visual detection of outliers
boxplot(log(tt_data$AS7263_610))

# which log-values are the outliers (based on Inter Quantile Range: Q25 - 1.5 * IQR  <-> Q75 + 1.5*IQR; IQR = Q75 - Q25)
boxplot.stats(log(tt_data$AS7263_610))$out

##################################################
## histogram with normal distribution curve overlay 
data <- log(tt_data$AS7263_610)
std <- sd(data,na.rm=T) #standard deviation of data; for normal distribution
m <- mean(data,na.rm=T) #mean of data; for normal distribution
hist(data, freq=F)
curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red")

## histogram with normal distribution curve overlay for all spectra
data <- log(tt_data[,5:16])
std <- lapply(data, sd, na.rm = T)
m <- lapply(data, mean, na.rm = T)

par(mfrow=c(4,3))
for(i in 1:12) {
  hist(data[,i], freq=F, main = paste(colnames(data)[i]))
  curve(dnorm(x,mean=m[[i]],sd=std[[i]]), add=TRUE,col="red")
}
