library(readxl)
library(dplyr)
library(ggplot2)
library(broom)


#read in raw data files

data1 <- read_xls("~/Documents/GitHub/C_gigas/Polyploids/docs/20190206_FirstReadingSTRIGG.xls")
data2 <- read_xls("~/Documents/GitHub/C_gigas/Polyploids/docs/20190206_2ndReadingSTRIGG.xls")
data3 <- read_xls("~/Documents/GitHub/C_gigas/Polyploids/docs/20190206_3rdReadingSTRIGG.xls")

Feb2019 <- merge(data1[,c(3,6)],data2[,c(3,6)], by = "Well")
Feb2019 <- merge(Feb2019,data3[,c(3,6)], by = "Well")
colnames(Feb2019) <- c("Well", "Abs1", "Abs2", "Abs3")
Feb2019$mean <- apply(Feb2019[,2:4],1, mean)
Feb2019$sd <- apply(Feb2019[,2:4],1, sd)
max(Feb2019$sd)
#[1] 0.01871462


#SDs are very small between readings so go with data1 readings

#Calculate standard curve

#first read in sample names from plate map list
Feb_plateMap <- read_xlsx("~/Documents/GitHub/C_gigas/Polyploids/docs/20190131_Ronit'sSamplesforGlobalDNAMeth.xlsx", sheet = 3)

#make a list of samples to not be considered
rmv_wells <- c("no_sample", "contam", "no_ab")

Feb_data <- merge(Feb_plateMap[which(!(Feb_plateMap$Sample %in% rmv_wells)),], data1[,c(3,6)])

#simplify absorbance column name
colnames(Feb_data)[4] <- "Absorbance"

#Calculate Average absorbance and SD for each sample
Feb_data_avg <- Feb_data %>% group_by(Sample) %>% summarize(avgAbs=mean(Absorbance))
Feb_data_sd <- Feb_data  %>% group_by(Sample) %>% summarize(SD=sd(Absorbance))

Feb_data_avg <- merge(Feb_data_avg,Feb_data_sd, by = "Sample")

Feb_data_avg <- unique(merge(Feb_data[,c("Type","Sample")], Feb_data_avg, by = "Sample"))

#make new data frame with only STDs
curve <- Feb_data_avg[grep("STD", Feb_data_avg$Type),]
print(curve)

#make new column for %meth
curve$perc_meth <- c(0.1,0.2,0.5,1,5,0,10)
#order data by perc_meth
curve <- curve[order(curve$perc_meth),]
print(curve)

#plot curve
ggplot(curve, aes(perc_meth, avgAbs)) + geom_point() 

##last two points are plateauing so remove them
ggplot(curve[which(curve$perc_meth < 5.0),], aes(perc_meth, avgAbs)) + geom_point() + geom_smooth(method = "lm", se = FALSE)

#find regression line
fit <- lm(curve[which(curve$perc_meth < 5.0),"avgAbs"] ~ curve[which(curve$perc_meth < 5.0),"perc_meth"] )
#find the R-squared
rsquared <- summary(fit)$r.squared
rsquared
#[1] 0.9978712

#find the slope (https://www.cyclismo.org/tutorial/R/linearLeastSquares.html)
slope <- fit$coefficients[[2]]
slope

#plot regression line with equation and r-squared (https://www.listendata.com/2016/07/add-linear-regression-equation-and.html)
linear = function(k) {
  z <- list(xx = format(coef(k)[1], digits = 5),
            yy = format(abs(coef(k)[2]), digits = 5),
            r2 = format(summary(k)$r.squared, digits = 5));
  if (coef(k)[1] >= 0)  {
    eq <- substitute(italic(Y) == yy %*% italic(X) + xx*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(Y) == yy %*% italic(X) - xx*","~~italic(r)^2~"="~r2,z)  
  }
  as.character(as.expression(eq));              
}

x <- curve[which(curve$perc_meth < 5.0),"perc_meth"]
y <- curve[which(curve$perc_meth < 5.0),"avgAbs"]
fo = y ~ x
linplot <- ggplot(data = curve[which(curve$perc_meth < 5.0),], aes(x = perc_meth, y = avgAbs)) + geom_smooth(method = "lm", se=FALSE, color="black", formula = fo) +  geom_point() + theme_bw()
linplot1 = linplot + annotate("text", x = 0.3, y = 0.7, label = linear(lm(fo)), colour="black", size = 5, parse=TRUE)
linplot1

#calculate perc_meth for samples using equation in EpiGentek MethylFlash kit
#5-mC% = (sample OD - NC OD) / (slope x input sample DNA in ng) * 100%

NC_OD <- Feb_data_avg[grep("5 ul NC", Feb_data_avg$Sample),"avgAbs"]
slopexDNAamt <- slope * 25

Feb_data_avg$avgAbs_NC <- Feb_data_avg$avgAbs - NC_OD
Feb_data_avg$perc_meth <- Feb_data_avg$avgAbs_NC/slopexDNAamt *100

#make new data frame with only sample info
Feb_data_avg_Samples <- Feb_data_avg[which(Feb_data_avg$Type == "Sample"),]

#Add sample description to data frame
for (i in 1:nrow(Feb_data_avg_Samples)){
  if(substr(Feb_data_avg_Samples$Sample[i],1,1)=="D"){
    Feb_data_avg_Samples$ploidy[i] <- "diploid"
  }
  if(substr(Feb_data_avg_Samples$Sample[i],1,1)=="T"){
    Feb_data_avg_Samples$ploidy[i] <- "triploid"
  }
  if(as.numeric(substr(Feb_data_avg_Samples$Sample[i],2,3)) > 10){
    Feb_data_avg_Samples$treatment[i] <- "heat_stress"
  }
  if(as.numeric(substr(Feb_data_avg_Samples$Sample[i],2,3)) < 10){
    Feb_data_avg_Samples$treatment[i] <- "control"
  }
}

Feb_data_avg_Samples$condition <- paste(Feb_data_avg_Samples$ploidy, Feb_data_avg_Samples$treatment, sep = "_")

#plot absorbances x experimental group

ggplot(Feb_data_avg_Samples, aes(x = condition, y = perc_meth)) + geom_boxplot(aes(fill = condition), show.legend = F) + scale_fill_manual(values = c("dodgerblue4", "lightblue1", "darkgreen", "palegreen")) + xlab("Sample") + theme_bw() + ylab("% 5-mC") + ggtitle("Ploidy x heat stress effect on global 5-mC DNA methylation")


aov_2way <- aov(perc_meth ~ ploidy + treatment + ploidy:treatment, data = Feb_data_avg_Samples)
aov_2way_model_summary <- glance(aov_2way)

#p-value is significant at 0.05 so do a Tukeys test to see which effect is significant
tuk <- TukeyHSD(aov(perc_meth ~ ploidy + treatment + ploidy:treatment, data = Feb_data_avg_Samples))
tuk

#based on the Tukey's test results, the diploid heat stress %5-mC methylation is significantly different from 1. the diploid control 2. the triploid control and 3. the triploid heatstress. The triploid heat stress %5-mC is not different than the diploid or triploid controls. So in summary, diploid ctenidia shows a significant increase in global 5-mC methylation in response to heatstress while triploid ctenidia does not. 

