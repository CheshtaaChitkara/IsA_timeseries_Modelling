


#working directory
setwd("C:/Users/cheshtaac/Documents/ISA _sampling/Data")


# load packages
library(phyloseq)
library(vegan)
library(writexl)
library(openxlsx)
library(WriteXLS)
library(readxl)
library(tidyr)
library(lubridate)
library(dplyr)

# Import data

# Environmental data - Metadata_fix.xls
#env<-read.table("clipboard") # copied only until 'dayl'
env <- read.xlsx("metadata_fix.xlsx")
env<- env%>% mutate(Date = make_date(Year, Month, Day))

Depth<-1+(env$Depth<50)

data<-readRDS("ISAFjord_phyloseqfull_nometazoa_nosmallreads.rds")
data
otu <- otu_table(otu, taxa_are_rows = T)
env <- sample_data(env)
tax <- tax_table(data)
sample_names(otu)-> sample_names(env)

data <- phyloseq(otu, tax, env)
otu<-data.frame(otu_table(data)) # extract the OTU table of 'taxonomic' abundance

#filter
filter <- phyloseq::genefilter_sample(data, filterfun_sample(function(x) x>10), A = 0.01*nsamples(data))
#more than 10 individuals in more than 1% of the samples
ps_filtered <- prune_taxa(filter, data)

# Normalization
#with phyloseq
total = median(sample_sums(ps_filtered))
standf = function(x, t=total) round(t * (x / sum(x)))
data_norm = transform_sample_counts(ps_filtered, standf)
otu_norm <- data.frame(otu_table(data_norm))
no_seq$seq2 <- colSums(otu_norm)

env$no_of_Sequences <- no_seq$colSums.otu.
env$seq2 <- no_seq$seq2

colnames(otu_norm)<-gsub("X","",colnames(otu_norm))
# 1745(rows or otu) 279 (samples)
#how do you decide which filtering criteria to use? because in our dataset, if max stuff is removed then we will have sequences from only one season for example,
#and ignore the rest. we need to see what to remove and what not to remove.
#filtering criteria cannot be so strict


#checking if envt data is correlated or not
envn <- env[,c(1:26, 34:36)] 
envn[,27] <- as.factor(envn[,27])
str(envn)
pairs(envn)
#nitrates and phosphates. nitrates and silicates. chla10 and chla GFF. density and salinity. 



#################################################### BRAY CURTIS DISSIMILARITY - v/S TIME  beta div#######################################################
#community measurement - bray curtis dissim
#### 6.1.1) Distance from first month measurement
otu1 <- data.frame(t(otu_norm)) 
beta <- c()
s<-otu1[FALSE,]
for (i in 1:12){
  s[i,]<-colSums(otu1[env$Month==i,])/sum(as.numeric(colSums(otu1[env$Month==i,])))*100
}
for (i in 1:nrow(otu1)) {
  st <- otu1[i,]
  b <- vegdist(rbind(otu1[4,], st), method = "bray")
  # Add to the list
  beta <- rbind(beta, c(i, b))
}
colnames(beta) <- c("time", "beta")
beta <- as.data.frame(beta)

par(mfrow=c(1,1))
plot(env$Date, beta$beta, type="p", xlab="year", ylab="Bray-Curtis dissimilarity from monthly average", ylim=c(0,1))
points(env$Date, beta$beta, type="l")
abline(v=as.Date(paste0(2011:2019,"-01-01")), col="grey")
mtext(at=as.Date(paste0(2011:2019,"-06-15")), 2011:2019)


### test for a trend in dissiimalirty changes over time (is there a drift in community stucture)
require(EnvStats)
print(kendallSeasonalTrendTest(beta$beta, season=env$Month, year= env$Year))
#kendall takes into account seasonality
#no significant trend because LCL and UCL are neg and pos and not pos and pos(positive trend) or neg and neg(negative trend)
#also the p value is over 0.005 so not significant


#another way to do what we did above.
require(Kendall)
require(xts)
betadf<-data.frame(Date=env$Date, beta=beta$beta)
### convert to xts time series object
betats <- xts(betadf[,2], order.by=as.Date(betadf[,1], "%Y-%m-%d"))
plot(betats)
MannKendall(betats)
periodicity(betats)
#tau - answers the question : how strong the trend is? 0 -lowest, 1- highest. So 0.11 tau means very slight increase in trend.
#This means the community composition is changing but very slightly
#this test shows there is a slight increase but the previous test shows nothing. 
#so consider the previous test but then know that there might be a slight change in composition

### check visually for seasonality and non-linear/ non-monotonous trends
decompose.xts <-
  function (x, type = c("additive", "multiplicative"), filter = NULL) 
  {
    dts <- decompose(as.ts(x), type, filter)
    dts$x <- .xts(dts$x, .index(x))
    dts$seasonal <- .xts(dts$seasonal, .index(x))
    dts$trend <- .xts(dts$trend, .index(x))
    dts$random <- .xts(dts$random, .index(x))
    
    with(dts,
         structure(list(x = x, seasonal = seasonal, trend = trend,
                        random = if (type == "additive") x - seasonal - trend else x/seasonal/trend, 
                        figure = figure, type = type), class = "decomposed.xts"))
  }

plot.decomposed.xts <-
  function(x, ...)
  {
    xx <- x$x
    if (is.null(xx))
      xx <- with(x,
                 if (type == "additive") random + trend + seasonal
                 else random * trend * seasonal)
    p <- cbind(observed = xx, trend = x$trend, seasonal = x$seasonal, random = x$random)
    plot(p, main = paste("Decomposition of", x$type, "time series"), multi.panel = 4,
         yaxis.same = FALSE, major.ticks = "days", grid.ticks.on = "days", ...)
  }

betats<-to.monthly(betats)[,1]
plot(decompose(ts(betats, frequency = 12)))
#from this we can see that there is seasonality in the dataset, and the mankendall couldnt remove it. we dont see any trends in the dataset. maybe a slight one
#random - 

#betadec<-decompose.xts(betats)
#plot.decomposed.xts(betadec)
#par(mfrow=c(2,2))
#plot(betadec$x, main="Raw")
#plot(betadec$seasonal,main="Season")
#plot(betadec$trend, main="Trend")
#plot(betadec$random, main="Random")


#no drift in beta diversity. no drift in community composition. 
#put any figure


############################################################ ALPHA DIV ##################################################
# random subsampling
#  loading function

source("SubsampleNGS.R")

#  Subsampling
#Alpha <- SubsampleNGS(df,           # input data
  #                    n = 1000,        # repeat subsamping 100 times
   #                   sub = 10000 # subsample to the minimum library size in the dataset
#)



#df<-relSpecies_all
#df<-relGenus_all
#df<-relFamily_all
#df<-relOrder_all
#df<-relClass_all
df <- otu1
df <- data.frame(t(df))



# calculate and plot basic alpha diveristy indices
alpha1<-diversity(t(df), index = "invsimpson")
alpha1<-diversity(df, index = "shannon")
 
print(kendallSeasonalTrendTest(alpha1, season=env$Month, year= env$Year))
#no significant trend in diversity
plot(alpha1)

## loop for richness (number of genera)
alpha1<-c()

for (i in 1:ncol(t(df))){
  alpha1[i]<-length(t(df)[t(df)[,i]>0,i])
}


require(EnvStats)
print(kendallSeasonalTrendTest(alpha1, season=env$Month, year= env$Year))

plot(alpha1)
#looks like a regime shift
#################### TEST if the tipping point is real
#### is there a new stable system?

#looking for bimodality
#if we have two separate regimes then we shud have two different normal distributions for the data
#first way to look - histogram
###https://www.nature.com/articles/s41559-020-1273-8
require(diptest)
#alpha1 <- regData1$Series
plot(cut(alpha1,100))
hist(alpha1, breaks = 25)
dip.test(alpha1)
#### bi modal distribution (--> evidence for tipping point)
#there is more than one distribution - good indicator for regime shift
#a different system in the first two years and the others, but we cant be sure if the first two years are representative or theyre weird themselves necause of some reason
#remember - higher sampling frequency in first two years - remember that

####https://www.nature.com/articles/s41598-021-93843-z#Sec7
require(changepoint)
cpt.mean(alpha1)
#11/4/2013 - this is when the regime shift happened. after this date - when we see the value 91. thats the 91st data point
#after regularizing the date, then we do the same changepoint test we get april as the changing month for the regime shift
cpt.meanvar(alpha1)
#this also shows 91 - that means the data is robust
#after regularizing the regime shift shows to be around 2019, when the taxa is going back to the same curve as 2012.
cpt.var(alpha1)
# this is just the variance. all of these are different ways to look at the same thing. and we showed that we get similar results from all three.
##### statistical confirmed tipping point :)
#maybe dont call it tipping point. need more data to call it tipping point. "sudden loss of diversity" instead.


## another estimate to show that there is really a regime shift. 

#regularize the time- interpolation to see the mid of every month with the data
#we do this, because we sample every week in 2012 and 2013, so our richness could be biased that way. 
#so we get the data
require(pastecs)
alpha_new <- data.frame(Date = as.POSIXct(env$Date), alpha = alpha1)
Data2<-alpha_new[,1:2]
d<-difftime(as.Date(Data2[,1]), as.Date(Data2[1,1]), units="days")
Data2<-data.frame(as.vector(d), as.Date(Data2[,1]), as.numeric(Data2[,2])) 
str(Data2) ### Check that the values are in the right format
Data2[1,2]
tail(Data2)


 #check if the data is significantly different - the two groups before and after the regime shift
 alpha_new$log <- NA
alpha_new$log <- alpha_new[1:87,] == TRUE
write.xlsx(alpha_new, "alpha_new.xlsx")
alpha_new <- read.xlsx("alpha_new.xlsx")
alpha_new <- alpha_new[-c(276:279),]
kruskal.test(alpha~log, data=alpha_new)
#the two regimes are significantly different

#regulation of series
regData<-regul(Data2[,1], 
               Data2[,3], 
               frequency=12, 
               units="daystoyears",
               datemin="14/12/2011", 
               dateformat ="d/m/Y", tol=13,#tol=13
               n=97, method="linear")

plot(regData)

#check if there is significant difference between the two regime shifts but with regularized dataset
regDat <- regData$y
regDat$log <- NA
write.xlsx(regDat, "regDat.xlsx")
regDat <- read.xlsx("regDat.xlsx")
regDat <- regDat[-c(94:95),]
kruskal.test(Series~log, data=regDat)

## list of lost OTU's - divide the dataset into before April 2013 and after April 2013 until month 96 
#take rowsums and compare the taxa lost during that regime shift

 #dividing data into before and after regime shift
data_1213 <- subset_samples(data_norm,  Year == 2011 | Year == 2012 | Year == 2013 & Month < 4 )
data1319 <- subset_samples(data_norm, Year == 2013 & Month >= 4 | Year == 2014 | Year == 2015 | Year == 2016 | Year == 2017 | Year == 2018 | Year == 2019 & Month <= 11)


#checking if the species numbers changed from one regime to the next 
otu2 <- data.frame(otu_table(data_1213))
row <- data.frame(rowSums(otu2))
#normalise
row <- row/87
otu3 <- data.frame(otu_table(data1319))
row2 <- data.frame(rowSums(otu3))
#normalise
row2 <- row2/190
rowr <- cbind(row, row2)

 row$logical <- (row$rowSums.otu2.>0) == TRUE
sum(row$logical, na.rm = TRUE)
  
row2$logical <- (row2$rowSums.otu3.>0) == TRUE
sum(row2$logical, na.rm = TRUE)

#discuss with tobias - because otu table data not same as his, then we cant really see how many species

#nothing is lost from the row 2 to row 3, but there is still a reduction in richness anyways

################################################### ABUNDANT TAXA DINOPHYTES #####################################################

dino <- subset_taxa(ps_filtered, Division == "Ochrophyta")

otu_dino <- data.frame(otu_table(dino))
taxa_dino <- data.frame(tax_table(dino))
source("PlotAbund.R")

taxa_pooled <- taxa.pooler(data.frame(otu_dino, (taxa_dino)))
genus_all <- t(taxa_pooled$Genus)
rel_genus <- prop.table(genus_all,2)*100
fam_all <- t(taxa_pooled$Family)
rel_fam <- prop.table(fam_all, 2)*100
sp_all <- t(taxa_pooled$Species)
rel_sp <- prop.table(sp_all, 2)*100
PlotAbund(rel_sp, 2) #,colorPalette = c("put your own colors"))
PlotAbund(rel_fam, 2)
PlotAbund(rel_genus, 2) #plots top 2 genuses of every sample

######################################################  RDA : UNCONSTRAINED ANALYSIS ############################################################

#RDA 
#separate into different clusters based on cluster analysis from raul's course

#why rda?
#check for example the temp ranges - if its long enough

#months 1,2,3

clus1 <- subset_samples(data_norm, Month == 1 | Month == 2 | Month == 3 | Month == 4 & Day < 19 ) #1,2,3,4 - Temp, maybe sSal 
clus2 <- subset_samples(data_norm, Month == 10| Month == 11| Month == 12) #10,11,12 - Watermass, Temp, chla, sal
clus3 <- subset_samples(data_norm, Month == 7 | Month == 8 | Month == 9 ) #7 ,8 , 9 - Temp
clus4 <- subset_samples(data_norm, Month == 4 & Day >= 19 | Month == 5 & Day < 19)  # 4,5 - watermass and phosphate
clus5 <- subset_samples(data_norm, Month == 6 | Month == 5 & Day >= 19)  #5 , 6 - watermass



#clus1
env1 <- data.frame(sample_data(clus1))
env1 <- env1[,c(3,7,8,12,21:23,25,31)]
str(env1)
for (i in 1:8)
{
  env1[,i] <- scale(as.numeric(env1[,i]))
}



env1[,9] <- as.factor(env1[,9])
str(env1)
otu11 <- data.frame(otu_table(clus1))
otu11 <- data.frame(t(otu11))

merge1 <- cbind(env1, otu11)
merge1 <- na.omit(merge1)
str(merge1)

otu11 <- data.frame(merge1[,12:1756])
env1 <- data.frame(merge1[,1:11])

pairs(env1)
#dont use silicates and nitrates at the same time. as they are correlated. same with nitrates and phsophates. 
#but silicates and phosphates can be used togther. 
#at the end use sili and phos because they are community drivers and not nitrates - that helps in growth

#check with other envt factors as well and other periods.
a <- rda(otu11 ~   env1$Temp + env1$Sal ) #bray curtis distance with hellinger trabnsformation
anova(a, by = "terms")
plot(a, type = "t")
# strongest variable changes over time- so if Year is the strongest - that means higher the year more the community drifts

#clus2 

env2 <- data.frame(sample_data(clus2))
env2 <- env2[,c(3,7,8,12,21:23,25,31)] #didnt inlcude bacteria and virus, too many NA's
str(env2)
for (i in 1:8)
{
  env2[,i] <- scale(as.numeric(env2[,i]))
}

env2[,9] <- as.factor(env2[,9])
str(env2)
otu12 <- data.frame(otu_table(clus2))
otu12 <- data.frame(t(otu12))

merge2 <- cbind(env2, otu12)
merge2 <- na.omit(merge2)
str(merge2)

otu12 <- data.frame(merge2[,10:1754])
env2 <- data.frame(merge2[,1:9])

pairs(env2)


#check with other envt factors as well and other periods.
a2 <- rda(otu12 ~    env2$watermass_value + env2$chla.GFF + env2$watermass_value)
anova(a2, by = "terms")
plot(a2, type = "t")


#clus3 

env3 <- data.frame(sample_data(clus3))
env3 <- env3[,c(3,7,8,12,21:23,25,31)]
str(env3)
for (i in 1:8)
{
  env3[,i] <- scale(as.numeric(env3[,i]))
}

env3[,9] <- as.factor(env3[,9])
str(env3)
otu13 <- data.frame(otu_table(clus3))
otu13 <- data.frame(t(otu13))

merge3 <- cbind(env3, otu13)
merge3 <- na.omit(merge3)
str(merge3)

otu13 <- data.frame(merge3[,10:1754])
env3 <- data.frame(merge3[,1:9])

pairs(env3)
#dont use silicates and nitrates at the same time. as they are correlated. same with nitrates and phsophates. 
#but silicates and phosphates can be used togther. 
#at the end use sili and phos because they are community drivers and not nitrates - that helps in growth

#check with other envt factors as well and other periods.
a3 <- rda(otu13 ~  env3$Temp + env3$Nitrates +env3$F + env3$watermass_value )
anova(a3, by = "terms")
plot(a3, type = "t")
# strongest variable changes over time- so if Year is the strongest - that means higher the year more the community drifts


#clus4 

env4 <- data.frame(sample_data(clus4))
env4 <- env4[,c(3,7,8,12, 21:23,25,31)]
str(env4)
for (i in 1:8)
{
  env4[,i] <- scale(as.numeric(env4[,i]))
}

env4[,9] <- as.factor(env4[,9])
str(env4)
otu14 <- data.frame(otu_table(clus4))
otu14 <- data.frame(t(otu14))

merge4 <- cbind(env4, otu14)
merge4 <- na.omit(merge4)
str(merge4)

otu14 <- data.frame(merge4[,10:1754])
env4 <- data.frame(merge4[,1:9])

pairs(env4)
#dont use silicates and nitrates at the same time. as they are correlated. same with nitrates and phsophates. 
#but silicates and phosphates can be used togther. 
#at the end use sili and phos because they are community drivers and not nitrates - that helps in growth

#check with other envt factors as well and other periods.
a4 <- rda(otu14 ~  env4$watermass_value  + env4$Nitrates  + env4$chla.GFF )
anova(a4, by = "terms")
plot(a4, type = "t")
# strongest variable changes over time- so if Year is the strongest - that means higher the year more the community drifts


#clus5 

env5 <- data.frame(sample_data(clus5))
env5 <- env5[,c(3,7,8,12, 21:23,25,31)]
str(env5)
for (i in 1:8)
{
  env5[,i] <- scale(as.numeric(env5[,i]))
}

env5[,9] <- as.factor(env5[,9])
str(env5)
otu15 <- data.frame(otu_table(clus5))
otu15 <- data.frame(t(otu15))

merge5 <- cbind(env5, otu15)
merge5 <- na.omit(merge5)
str(merge5)

otu15 <- data.frame(merge5[,10:1754])
env5 <- data.frame(merge5[,1:9])

pairs(env5)
#dont use silicates and nitrates at the same time. as they are correlated. same with nitrates and phsophates. 
#but silicates and phosphates can be used togther. 
#at the end use sili and phos because they are community drivers and not nitrates - that helps in growth

#check with other envt factors as well and other periods.
a5 <- rda(otu15 ~  env5$watermass_value + env5$Temp + env5$Phosphate  + env5$chla.GFF + env5$Sal)
anova(a5, by = "terms")
plot(a5, type = "t")
# strongest variable changes over time- so if Year is the strongest - that means higher the year more the community drifts

#RDA FOR REGULARIZED DATASET
otu_top1 <- data.frame(t(abund_rel[c(1:20),]))

ps1.com.fam <- aggregate_top_taxa2(ps_filtered, "Species", top = 10)
otu_top <- data.frame(otu_table(ps1.com.fam))
otu_top <- data.frame(t(otu_top))
env_date <- data.frame(sample_data(ps_filtered))
otu_top1$date <- env_date$Date

#taking out "other"
otu_top <- otu_top[,-9]
#regularization - what species count it is supposed to be in the mid of every month. So we have one value for one month, instead of so many for example in 2012

otu_top1$d<-difftime(as.Date(otu_top1$date), as.Date(otu_top1$date[1]), units="days") ### convert dates to days since first sampling
## 2920 days --> 98 months (2920/30)
require(pastecs)
regData1<-list()
for (i in 1:ncol(otu_top1)){
  regData1[[i]]<-regul(as.vector(otu_top1$d), otu_top1[,i], frequency=12, units="daystoyears",
                      datemin="14/12/2011", 
                      dateformat ="d/m/Y", tol=8,#tol=13
                      n=98, method="spline")
  #try different methods to see if we get better red and black fits
}

regdat <- (regData1[[1]])$y
regdat$Gyrodinium_spirale

##########################################  TOP 10 SPECIES IN DIFFERENT CLUSTERS ############################################################
#Physical method to find out the top 10 species. 
otu_all <- data.frame(otu_table(data_norm))
tax_all <- data.frame(tax_table(data_norm))
taxa_pooled <- taxa.pooler(data.frame(otu_all, (tax_all)))
class_all <- t(taxa_pooled$Class)
rel_class <- prop.table(class_all,2)*100
genus_all <- t(taxa_pooled$Genus)
rel_genus <- prop.table(genus_all,2)*100
fam_all <- t(taxa_pooled$Family)
rel_fam <- prop.table(fam_all, 2)*100
sp_all <- t(taxa_pooled$Species)
rel_sp <- prop.table(sp_all, 2)*100
h <- PlotAbund(rel_sp, 2) #,colorPalette = c("put your own colors"))
PlotAbund(rel_fam, 2)
PlotAbund(rel_genus, 2)
PlotAbund(rel_class, 2)

#aggregate top taxa method
ps1.com.fam <- aggregate_top_taxa2(data_norm, "Species", top = 10)
ps1.com.fam.rel <-  microbiome::transform(ps1.com.fam, "compositional")
plot.composition.relAbun <- plot_composition(ps1.com.fam.rel,sample.sort = "Depth_new",x.label = "date")
plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
plot.composition.relAbun <- plot.composition.relAbun + scale_fill_brewer("Species", palette = "Paired") + theme_bw() 
plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Deep-Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))
print(plot.composition.relAbun)


#top 10 species in 5 clusters
#clus 1 - dino groups + gyrodinium spirale, mast, picozoa, prorocentrum
#clus 2 - dino groups but lesser + gyrodinium helveticum + gyro spirale + mast + picozoa + prorocentrum
#clus 3 - dino groups (even lesser) + gymnodinium+ gyro fusiforme+ gyro helveticum + gyro spirale + heterocapsa rotunda+ picozoa + prorocentrum
#clus 4 - one dino group + gymnodinium + gyro fusi + gyro helve+ gyro spirale + leegardiella + micromonas arctica + micromonas (other one) 
#         pentapharsodinium + strombidiidae
#clus 5 -  dino group (one) + chytriodinium_roseum + gyro dominans + gyro fusi + gyro helve + gyro spirale + heterocapsa_rotunda + pelagophyceaea
#         Stephanoecidae + woloszynskia_halophila  

##############################################################  MARSS MODEL ###########################################################################  

ps1.com.fam <- aggregate_top_taxa2(data_norm, "Species", top = 10)
ps1.com.fam.otu <- subset(otu_table(data_norm), rownames(otu_table(data_norm)) %in% c('otu0027', 'otu0033', 'otu0012', 'otu0011', 'otu0020', 'otu0009', 'otu0013',
                                                                                  'otu0018', 'otu0004', 'otu0003', 'otu0001', 'otu0008', 'otu0007', 'otu0006', 
                                                                                  'otu0014', 'otu0002', 'otu0061', 'otu0024', 'otu0021', 'otu0028'))
ps1.com.fam.tax <- subset(tax_table(data_norm), rownames(tax_table(data_norm)) %in% c('otu0027', 'otu0033', 'otu0012', 'otu0011', 'otu0020', 'otu0009', 'otu0013',
                                                                                  'otu0018', 'otu0004', 'otu0003', 'otu0001', 'otu0008', 'otu0007', 'otu0006', 
                                                                                  'otu0014', 'otu0002', 'otu0061', 'otu0024', 'otu0021', 'otu0028'))
 
ps1.com.fam <- merge_phyloseq(ps1.com.fam.otu, ps1.com.fam.tax, sample_data(data_norm))
otu_top <- data.frame(otu_table(ps1.com.fam))
 otu_top <- data.frame(t(otu_top))
env_date <- data.frame(sample_data(ps1.com.fam))
otu_top$date <- env_date$Date

#taking out "other"
otu_top <- otu_top[,-9]
 #regularization - what species count it is supposed to be in the mid of every month. So we have one value for one month, instead of so many for example in 2012
 
otu_top$d<-difftime(as.Date(otu_top$date), as.Date(otu_top$date[1]), units="days") ### convert dates to days since first sampling
## 2920 days --> 98 months (2920/30)
require(pastecs)
regData<-list()
for (i in 1:ncol(otu_top)){
  regData[[i]]<-regul(as.vector(otu_top$d), otu_top[,i], frequency=12, units="daystoyears",
                      datemin="14/12/2011", 
                      dateformat ="d/m/Y", tol=9,#tol=13
                      n=98, method="spline")
  #try different methods to see if we get better red and black fits
}

#n = number of months = 98
par(mfrow=c(3,4))
for(i in 1:12){
  plot(regData[[i]])}
par(mfrow=c(3,4))
for(i in 13:20){
  plot(regData[[i]])}
#each plot is for each species - red color predicts what it should be (interpolation), and the black is the data
#change the tol so the red and the black fit well with each other

ts0<-tseries(regData[[1]])
#converting to time series object
for (i in 2:ncol(otu_top)){
  tsi<-tseries(regData[[i]])
  tsi<-ts.union(ts0, tsi)
  ts0<-tsi
}
colnames(tsi)<-as.character(colnames(otu_top[,1:ncol(otu_top)]))
tsAll <- tsi[,1:(ncol(otu_top)-2)]

as.data.frame(tsAll)
#regularized data - per month

relGenus_reg <- prop.table(t(as.data.frame(tsAll)), 2) * 100
#proportions

require(zoo)
colnames(relGenus_reg)<-as.character(as.yearmon(time(tsAll)))


relGenus_reg<-relGenus_reg[,1:98]
#we need regularized time series objects to remove seasonality

#DECOMPOSE THIS REGULARIZED DATASET - basically taking out seasonality

#do rda on regularized data set - to see if watermass influences the community


#checking if different species show different trends.
plot(decompose(tsAll[,1])) #this just shows us the different trends...shizz...
#species 1 in tsAll[,1] - shows the regime shft that we see in richness.
#how do we check if random is actually random? - 
#cant handle 0's too well - so we might need to pool some species together - all gyrodinium to gyro spp.
#all top 10 species abundance is decreasing- so many 
#try rarefying and then finding the top 10 taxa. mayve aggregate top taxa is looking for top taxa in 2012-13 samples cos its so much

tsAll #regularized and converted to time series object
str(tsAll)
plot(tsAll[,1])
decompose(tsAll[,1])
plot(decompose(tsAll[,2]))


#### NEXT
#### data frame with x
Datts<-tsAll
Datdf<-data.frame(Date=as.Date(time(Datts)), as.data.frame(Datts))
par(mfrow=c(3,4))

#for loop i in 1:10 because i have top 20 species
for (i in 1:20){plot(tsAll[,i])}

#### data frame with deseason
require(forecast)
Datdf_deseas<-data.frame(Date=Datdf$Date)
for (i in 1:20){
  Datdf_deseas[,i+1]<-seasadj(decompose(Datts[,i]))
}
colnames(Datdf_deseas)<-colnames(Datdf)
par(mfrow=c(3,4))
for (i in 1:20){plot(seasadj(decompose(tsAll[,i])))}
#this plot is basically random + trend because seaonality is removed. So we can check for each species why we have the peaks, knowing that season isnt affecting them



#### --> export as csv
write.csv(Datdf, file="community_groups_original.csv")
write.csv(Datdf_deseas, file="community_groups_deseason.csv")

#### --> import to MARSS

#deseasonlaise envtl data
#removing non numeric stuff
env_new1 <- data.frame(env[,-c(27,29:33)])
str(env_new1)

for (i in 1:26)
{
  env_new1[,i] <- as.numeric(env_new1[,i])
}

env_new1$d<-difftime(as.Date(env_new1$Date), as.Date(env_new1$Date[1]), units="days")
env_new1 <- env_new1 %>% relocate(Date, .before = d)
require(pastecs)
regData<-list()
for (i in 1:ncol(env_new1)){
  regData[[i]]<-regul(as.vector(env_new1$d), env_new1[,i], frequency=12, units="daystoyears",
                      datemin="14/12/2011", 
                      dateformat ="d/m/Y", tol=13,#tol=13
                      n=98, method="spline")
}

#n = number of months = 98
par(mfrow=c(3,4))
for(i in 1:12){
  plot(regData[[i]])}
#each plot is for each species - red color predicts what it should be (interpolation), and the balack is the data
#change the tol so the red and the black fit well with each other

par(mfrow=c(3,4))
for(i in 13:26){
  plot(regData[[i]])}


ts0<-tseries(regData[[1]])
#converting to time series object
for (i in 2:ncol(env_new1)){
  tsi<-tseries(regData[[i]])
  tsi<-ts.union(ts0, tsi)
  ts0<-tsi
}
colnames(tsi)<-as.character(colnames(env_new1[,1:ncol(env_new1)]))
tsAll <- tsi[,1:(ncol(env_new1)-1)]


as.data.frame(tsAll)
#refularized data - per month

#DECOMPOSE THIS REGULARIZED DATASET - basically taking out seasonality
#do rda on regularized data set - to see if watermass influences the community

plot(decompose(tsAll[,10])) #this just shows us the different trends...shizz...
#species 1 in tsAll[,1] - shows the regime shft that we see in richness.
#how do we check if random is actually random? - 
#cant handle 0's too well - so we might need to pool some species together - all gyrodinium to gyro spp.
#all top 10 species abundance is decreasing- so many 
#try rarefying and then finding the top 10 taxa. mayve aggregate top taxa is looking for top taxa in 2012-13 samples cos its so much

tsAll #regularized and converted to time series object
str(tsAll)
plot(tsAll[,1])
decompose(tsAll[,1])
plot(decompose(tsAll[,2]))


#### NEXT
#### data frame with x
Datts<-tsAll
Datdf<-data.frame(Date=as.Date(time(Datts)), as.data.frame(Datts))
par(mfrow=c(3,4))

#for loop i in 1:10 because i have top 10 species
for (i in 1:12){plot(tsAll[,i])}

#### data frame with deseason
require(forecast)
Datdf_deseas<-data.frame(Date=Datdf$Date)
for (i in 1:30){
  Datdf_deseas[,i+1]<-seasadj(decompose(Datts[,i]))
}
colnames(Datdf_deseas)<-colnames(Datdf)
par(mfrow=c(3,4))
for (i in 1:12){plot(seasadj(decompose(tsAll[,i])))}
#this plot is basically random + trend because seaonality is removed. So we can check for each species why we have the peaks, knowing that season isnt affecting them



#### --> export as csv
write.csv(Datdf, file="environmental_data_original.csv")
write.csv(Datdf_deseas, file="envtl_data_deseason.csv")

