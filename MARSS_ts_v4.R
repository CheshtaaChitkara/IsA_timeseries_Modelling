###################################################################################
#### Multivariate time series analyses ###########################################
##################################################################################
#### v1.4 Tobias Vonnahme 26.10.2021 #############################################
#################################################################################
#### Input: deseasonalized & regularized covariates (except AOI/NAOI stay as raw data) 
####        Community data, main groups, regularized and deseasonalized #########
#################################################################################
### 1. Data import and merging
### 2. Trend analyses
### 3. Autocorrelation tests
### 4. MARSS model for Nutrients using phycial drivers/covariates
### 5. Use MARSS model output to fill NAs in Nutrient data
### 6. MARSS model for phytoplankton community (main groups)
####################################################################################
### required packages: readxl, Kendall, MARSS
###################################################################################

setwd("C:/Users/cheshtaac/Documents/ISA _sampling/Data") #set working directory


########################################################
#### 1) Read in the data and merge into one data frame ###
#######################################################

### covariates (simplified data to 2D, regularized, deseasonalized)
covdat<-read.csv("envtl_data_deseason.csv")
covdat<-read.csv("environmental_data_original.csv")

if(colnames(covdat)[1] == "X"){covdat<-covdat[,-1]} #remove the first column that is sometimes added via csv conversion in R

### community data (subset of main groups, regularized, deseasonalized)
comdat<-read.csv("community_groups_deseason.csv")
#comdat<-read.csv("Phytoplankton_community_groups_original.csv")

if(colnames(comdat)[1] == "X"){comdat<-comdat[,-1]} #remove the first column that is sometimes added via csv conversion in R


#### merge community and covariate data frames
## check dimensions
dim(comdat)
dim(covdat)

## check if dates are the same
identical(covdat$Date, comdat$Date)
#identical(comdat$Date, as.character(as.Date(AO$Date)))
## make dates identifcal (its a regularized datasets so the difference is negligible)
covdat$Date<-comdat$Date

### merge
Alldat<-merge(covdat, comdat, by="Date", all=TRUE)
View(Alldat)
Alldat <- Alldat[-98,] #no data so remove it




###################################
#### 2) Trend analyses
#### Mann Kendall test testing for monotonous non-linear trends
#####################################

require(Kendall)

p<-data.frame(variable=rep(NA, ncol(Alldat)-1), p=rep(NA, ncol(Alldat)-1), tau=rep(NA, ncol(Alldat)-1),
              padj=rep(NA, ncol(Alldat)-1), sign=rep(NA, ncol(Alldat)-1))

for (i in 2:ncol(Alldat)){
  print(colnames(Alldat[i]))
  print(MannKendall(Alldat[,i]))
  p[i,1]<-colnames(Alldat[i])
  p[i,2]<-MannKendall(Alldat[,i])$sl
  p[i,3]<-MannKendall(Alldat[,i])$tau
}

p<-p[-1,]

p$padj<-p.adjust(p$p, method = "fdr") #this is adjustment for multiple testing -  if we do multiple tests, then by chance if we get a good p value we need to adjust for that
#look at padj to be statistically correct


for (i in 1:nrow(p)){
  if (p$padj[i]<=0.001){p$sign[i] <- "***"}
  else if (p$padj[i]<=0.01){p$sign[i] <- "**"}
  else if (p$padj[i]<=0.05){p$sign[i] <- "*"}
  else if (p$p[i]<=0.05){p$sign[i] <- "(*)"}
  else {p$sign[i] <- "ns"}
}

p
#if tau is negative then decreasing trend. goes between -1 and +1 - no trend. And then the ns means non significant trend. so its not changing over time. 
#this is with regularized deseaonlised data.
#significantly decreasing sequecing depth - means that we need to rarefy or normalise in some way to not see this significant trend. 
#very imp table can already check for trends in this


#################################################################
### 3) Autocorrelation tests
##################################################################

# BUT looking at the plots its not a monotonous trend, but stochastic variation
#### Check for autocorrelation (lines outside the blue dotted lines show potentially significant autocorrelation at with a time lag of x)
 
par(mfrow=c(3,4))
for (i in 3:13){
  autoplot(acf(na.omit(Alldat[,i]))) 
  }
dev.off()
#No.of sequences, Year, MM, Day, Depth, Sal, Temp, Cond, Dens, O, F
par(mfrow=c(3,4))
for (i in 14:25){
  acf(na.omit(Alldat[,i]))}
dev.off()
# Dayl, picoml, nanoml, crytptoml, synml, hnfml, bacteria, virus, silicate, nitrates, phosphate, chla10µm
par(mfrow = c(3,4))
for (i in 26:31){
  acf(na.omit(Alldat[,i]))}
dev.off()
#chla.GFf, integreated chla, sequences_after_normalisation, AO_index, NAO index, date.1
par(mfrow = c(3,4))
for (i in 32:43){
  autoplot(acf(na.omit(Alldat[,i])))
}

par(mfrow = c(3,4))
for (i in 44:51){
  autoplot(acf(na.omit(Alldat[,i])))
}
#enough autocorrelation amongst all envtl and community data.
#all show autocorrel except - pelagophyceae, heterocapsa pygmea, phaeocystis, synml
#and we cannot have this for any model. so we need a model thats made for time series. 
#this shows that RDA doesnt help in us understanding our data, this is why we need these time series models. So if we have two samples that are close by in an RDA that cud mean theyre closer because of the time, and not because of the species present




#### check the orders of a potential ARIMA model - easy model that doesnt take envtl data -  baseline time series model
for (i in 2:ncol(Alldat)){
  print(auto.arima(Alldat[,i]))}


############################################################################
#### 4) Multivariate state space modelling ###################################
#### MARSS model (Temp and sal as variables and everything else as covariates)
############################################################################
Alldat <- Alldat[-97,]
Alldat$Year<-as.numeric(substr(Alldat[,1], start=1, stop=4)) #add column with the year for later subsetting
full2<-as.matrix(Alldat[,4:ncol(Alldat)]) #convert dataframe to matrix format
years<-full2[,"Year"]>=2011 & full2[,"Year"]<=2020 #option to subset the dataset to certain years

colnames(full2) # check the order of the columns
dat <- t(full2[years,colnames(full2)[c(14,15)]]) #subset data of temp and sal
dat <- dat[,-c(1:2)] #removing first two samples because in covariates we have first two samples with NA flow counts
covariates <- t(full2[years,colnames(full2)[c(5,6,19:21,23,26,27)]]) #simplest model of covariates. covariates should be the drivers of the dat. so is picoml driving temp sal?
covariates <- covariates[,-c(1:2)] # removing first two samples, because they have NA's in flow counts

dim(dat) 
row.names(covariates)

# z-score the response variables
the.mean <- apply(dat,1,mean,na.rm=TRUE)
the.sigma <- sqrt(apply(dat,1,var,na.rm=TRUE))
dat <- (dat-the.mean)*(1/the.sigma)

## ----msscov-z-score-covar-data-----------------------------------------------------------------
the.mean <- apply(covariates,1,mean,na.rm=TRUE)
the.sigma <- sqrt(apply(covariates,1,var,na.rm=TRUE))
covariates <- (covariates-the.mean)*(1/the.sigma)

## ----msscov-plank-plot, fig=TRUE, echo=FALSE, fig.cap='(ref:msscov-plank-dat)', warning=FALSE----
LWA <- ts(cbind(t(dat), t(covariates)), start=c(2011,8), freq=12) #bind the observations and covariates
LWA1 <- LWA[,1:10]
plot(LWA, main="", yax.flip=TRUE) # plot the time series
LWA2 <- LWA[,11:17]
plot(LWA2, main="", yax.flip=TRUE)

### Model testing
require(MARSS)
A <- U <- x0 <- "zero"       #### No trend/drift 
c = covariates         #### 
C <- "unequal"#"unconstrained"   #### effects of different process errors on states ,  "equal" gives lower AIC, but doesnt make sense here!
#each covariate can have a diff strength for the driver
Q <- "diagonal and unequal"#"unconstrained"  ### Test also "diagonal and unequal" "diagonal and equal" "equalvarcov" "identity" "unconstrained"
Z <- "identity"     ### only 1 time series (no different states)
B <- "diagonal and equal" ### state at t-1 affects state at t equally among both states, tried "diagonal and unequal" "identity
R <- "equalvarcov" #matrix(list(0.05,0,0,0.1),2,2) ### same effect of different observation errors (same sampling + similar measurement methods) diagonal and unequal before
tinitx <- 1             
y <- dat 

model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.1<-MARSS(y, model = model.list) #AIC: 515 - covariates - AO, NAO
                                    #AIC: 534 - covariates - temp, sal, nao, ao, silicates, phosphates, nitrates 
                                    #AIC: 530 - covariates - temp, sal, nao, ao, silicates, phosphates
                                    #AIC: 532 - covariates - temp, sal, nao, ao, silicates, phosphates, chla
                                    #AIC: 536 - covariates - temp, sal, nao, ao, silicates, nitrates, phosphates, chla

#### First model estimate

U2 <- "unconstrained"
model.list <- list(B = B, U = U2, Q = Q, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.2<-MARSS(y, model = model.list) #AIC: 519 


C2<-"unconstrained"
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, c=c, C=C2, x0=x0, tinitx=tinitx)
kem.3<-MARSS(y, model = model.list) #AIC: 515
#equal fit
#### -> same fit, unconstrained has a higher DF than unequal, also uequal effects of different covariates makes more sense
### equal is another option but makes no sense and has higher AIC  -> Discard

Q2<-"unconstrained"
model.list <- list(B = B, U = U, Q = Q2, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.4<-MARSS(y, model = model.list) #AIC: 517
#Better
#### Lower AIC and higher values for Q, meaninful with unequeal covariates and variance 
#### process errors have different variance for different taxa (var, diagonal) and year-to year changes covary/ process errors are dependent (inherent in compositionality data!!) -> Keep

B2 <- "diagonal and unequal"
model.list <- list(B = B2, U = U, Q = Q2, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.5<-MARSS(y, model = model.list) #AIC: 519


R2 <- "diagonal and unequal"
model.list <- list(B = B, U = U, Q = Q2, Z = Z, A = A, R = R2, c=c, C=C, x0=x0, tinitx=tinitx)
kem.6<-MARSS(y, model = model.list) #AIC: 525
#better but: "Warning: the  R.(NOx_int,NOx_int)  parameter value has not converged."
#### observation errors are different for the different taxa, but year-to year changes in the obs error are independent (diagonal) -> Keep


### define the best model fit
kem.t<-kem.3
kem.t
### calculate confidence intervals (if CI do not overlap with 0 thye are considered significant)
MARSSparamCIs(kem.1)

### diagnostic plots
plot(kem.t)



#######################################################################
#### 6) MARSS model for communities
##########################################################
##########################################################

### prepare dataset with observations (y) and covariates (c) in matrix format
### Matrix conversion for community
full<-as.matrix(Alldat[,2:ncol(Alldat)])
years <- full[,"Year"]>=2011 & full[,"Year"]<=2019


########################################
#### Best model fit
colnames(full)
dat <- t(full[years,colnames(full)[c(32,34,35)]]) # deseasonalized groups, subset observation data of the most abundant algae groups (Chaetoceros, Thalassiosira, Phaeocystis)
dat <- t(full[years,colnames(full)[c(31:50)]]) # deseasonalized All groups
covariates <- t(full[years,colnames(full)[c(7,8,16,25,28,29)]]) #avoid correlated stuff in covariates
covariates <- covariates[,-c(1:2)]
dat <- dat[,-c(1:2)]

dim(dat)
dim(covariates)



# z-score the response variables
the.mean <- apply(dat,1,mean,na.rm=TRUE)
the.sigma <- sqrt(apply(dat,1,var,na.rm=TRUE))
dat <- (dat-the.mean)*(1/the.sigma)
## ----msscov-z-score-covar-data-----------------------------------------------------------------
the.mean <- apply(covariates,1,mean,na.rm=TRUE)
the.sigma <- sqrt(apply(covariates,1,var,na.rm=TRUE))
covariates <- (covariates-the.mean)*(1/the.sigma)
## ----msscov-plank-plot, fig=TRUE, echo=FALSE, fig.cap='(ref:msscov-plank-dat)', warning=FALSE----
LWA <- ts(cbind(t(dat), t(covariates)), start=c(2011,8), freq=12) #bind the observations and covariates
LWA1 <- LWA[,1:10]
plot(LWA1, main="", yax.flip=TRUE) # plot the time series
#plots data after standardization. just the data



### Model testing
### own model (community)
A <- U <- x0 <- "zero"       #### No trend/drift -> justified by MK test
c = covariates         #### 
C <- "unequal"#"unconstrained"   #### effects of different process errors on states ,  "equal" gives lower AIC, but doesnt make sense here!
Q2 <- "unconstrained"#"unconstrained"  ### Test also "diagonal and unequal" "diagonal and equal" "equalvarcov" "identity" "unconstrained"
Z <- "identity"     ### only 1 time series (no different states)
B <- "diagonal and equal" ### state at t-1 affects state at t equally among both states, tried "diagonal and unequal" "identity
R2 <- "diagonal and unequal" #matrix(list(0.05,0,0,0.1),2,2) ### same effect of different observation errors (same sampling + similar measurement methods) diagonal and unequal before
tinitx <- 1              ### Needed for estimating B
y <- dat 


model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.1<-MARSS(y, model = model.list) #AIC: 5313 - covariates - temp, sal, ao, nao   
                                    #AIC: 5087 - covariates - temp, sal, nao, ao, bac, vir, nitrates, chla
                                    #AIC: 5092 - covariates - temp, sal, nao, ao, bac, vir, silicates, phosophates, chla
                                    #AIC: 5070 - covariates - temp, sal, nao, ao, picoml, cryptoml, hnfml, bac, vir, silicates, phosphates,chla -  picoml not imp, virus not imp
                                    #AIC: 5055 - covariates - temp, sal, nao, ao, cryptoml, hnfml, bac, silicates, phosphates, chla
                                    #AIC: 5052 - covariates - temp, sal, nao, ao, cryptoml, hnfml, bac, nitrates, chla - hnf not so imp, bac not so imp, nitrates not so imp 
                              #BEST #AIC: 5017 - covariates - temp, sal, nao, ao, cryptoml, chla  - maybe because dinos and ciliates feed on cryptophytes - gyrodinium selctive grazer for crypto
                                    #AIC: 5046 - covariates - temp, sal, nao, ao, cryptoml, silicates, phosphates, chla
                                    #AIC: 5061 - covariates - temp, sal, nao, ao, hnfml, silicates, phosphates, chla
                                    #AIC: 5181 - covariates - temp, sal, nao, ao, hnfml, chla


U2 <- "unconstrained" #allow for a trend or a drift
model.list <- list(B = B, U = U2, Q = Q, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.2<-MARSS(y, model = model.list) #AIC: 5057 - covariates - temp, sal, nao, ao, cryptoml, chla


C2<-"unconstrained" #any combination/strength of covariate drivers is possible. how important is each driver.
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, c=c, C=C2, x0=x0, tinitx=tinitx)
kem.3<-MARSS(y, model = model.list) #AIC: 5017 - covariates - temp, sal, nao, ao, cryptoml, chla


Q2<-"unconstrained" #
model.list <- list(B = B, U = U, Q = Q2, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.4<-MARSS(y, model = model.list) #AIC: 4611 - covariates - temp, sal, nao, ao, cryptoml, chla


B2 <- "diagonal and unequal"
model.list <- list(B = B2, U = U, Q = Q2, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.5<-MARSS(y, model = model.list) #AIC: 4542 - covariates - temp, sal, nao, ao, cryptoml, chla


R2 <- "diagonal and unequal"
model.list <- list(B = B, U = U, Q = Q2, Z = Z, A = A, R = R2, c=c, C=C, x0=x0, tinitx=tinitx)
kem.6<-MARSS(y, model = model.list) #AIC: 4615

### define the best model fit
kem.t<-kem.5
kem.t
### calculate confidence intervals (if CI do not overlap with 0 thye are considered significant)
 MARSSparamCIs(kem.t)

### diagnostic plots
par(mfrow = c(3,4))
plot(kem.t)

#########################################################################################
#### FORECAST ###########
#########################################################################################

require(forecast)
fr <- forecast(kem.5, newdata=list(y=NULL, c=kem.5$call$model$c, d=NULL), h = 15)
plot(fr)

#######################################################################################

