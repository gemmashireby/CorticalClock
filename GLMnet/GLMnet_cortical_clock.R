############ BUILDING  AN EPIGENETIC CLOCK SPECIFIC FOR THE BRAIN (NOT SPECIFIC WITHIN BRAIN REGIONS) - using 450K and EPIC data 
############ AUTHOR: GEMMA SHIREBY
############ DATE: 18-09-2019
############ USING GLMNET - PENALIZED REGRESSION METHODS 
############ based on methodology first used by S. Horvath (Horvath, S. (2013). DNA methylation age of human tissues and cell types. Genome biology, 14(10), 3156)

## LOAD LIBRARIES
#install.packages("glmnet")
library(glmnet)

### LOAD DATA 
setwd("/mnt/data1/Gemma/Brain_clock/Age_prediction/Brain_epigenetic_clock/all_sample_clock/")
load("BigBrain_LBB2_ROSMAP_Jaffe_betas_dasen_pheno.rdat") # pheno betas.dasen
dim(pheno) # 1397    6
dim(betas.dasen) # 383547   1397


## Get random sample 
# set.seed(101) # Set Seed so that same sample can be reproduced in future
# sample = sample.split(pheno$Sentrix_Full, SplitRatio = .75) # Selecting 75% of data as sample from total 'n' rows of the data
# train = subset(pheno, sample == TRUE)
# dim(train) # n = 1047    
# summary(as.factor(train$Sex))
# # F   M 
# # 362 685 
# test  = subset(pheno, sample == FALSE)
# dim(test) # n = 350   
# summary(as.factor(test$Sex))
# # F   M 
# # 144 206 
# write.csv(train, "Brain_clock_all_450K_brain_data_training_sample.csv", row.names=F)
# write.csv(test, "Brain_clock_all_450K_brain_data_testing_sample.csv", row.names=F)

################### READ IN TRAINING

train<-read.csv("Brain_clock_all_450K_brain_data_training_sample.csv")

### get some details on the sample
phenoTrain<-pheno[which(pheno$Sentrix_Full %in% train$Sentrix_Full),]
summary(as.factor(phenoTrain$BrainRegion))

summary(as.factor(phenoTrain$BrainBank))

################### READ IN TESTING

test<-read.csv("Brain_clock_all_450K_brain_data_testing_sample.csv")
## get some details on the samples 

phenoTest<-pheno[which(pheno$Sentrix_Full %in% test$Sentrix_Full ),]
summary(as.factor(phenoTest$BrainRegion))
summary(as.factor(phenoTest$BrainBank))


pdf("figures/Histogram testing ages.pdf", width=8, height=6)
hist(test$Age, xlab="Age",main="Testing Ages Distribution")
dev.off()

pdf("figures/Histogram training ages.pdf", width=8, height=6)
hist(train$Age, xlab="Age",main="Training Ages Distribution")
dev.off()

#################################################################################################################
######## TRANSFORM YOUNG AGE (to account for the logaritmic relationship from 0-20 years)  ######################
#################################################################################################################

# EQUATIONS FOR LOG TRANSFORMATION
# F(age)=log(age+1)-log(adult.age+1)age<= 20
# F(age)=(age-adult.age)/(adult.age+1) if age> 20

trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
test$AgeT=trafo(test$Age)


train$AgeT=trafo(train$Age)

## transform back at the end 

anti.trafo<-function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
test$AgeAT<-anti.trafo(test$AgeT)


train$AgeAT<-anti.trafo(train$AgeT)



#################################################################################################################
################################################ MATCH TO BETAS  ################################################
#################################################################################################################

# training
betasTrain<-betas.dasen[,match(train$Sentrix_Full, colnames(betas.dasen))]
dim(betasTrain) # 383547   1047
betasTrain[1:10,1:10]
# testing
betasTest<-betas.dasen[,match(test$Sentrix_Full, colnames(betas.dasen))]
dim(betasTest) # 383547    350
betasTest[1:10,1:10]

#################################################################################################################
################################# RUN ELASTIC NET REGRESSION TO GET CLOCK COEFS #################################
#################################################################################################################

# use 10 fold cross validation to estimate the lambda parameter 
# in the training data
alpha<-0.5
glmnet.Training.CV = cv.glmnet(t(betasTrain), train$AgeT,  nfolds=10,alpha=alpha,family="gaussian") ## CV = cross-validation
# The definition of the lambda parameter:
lambda.glmnet.Training = glmnet.Training.CV$lambda.min
# Fit the elastic net predictor to the training data
glmnet.Training = glmnet(t(betasTrain),train$AgeT, family="gaussian", alpha=0.5, nlambda=100)
# Arrive at an estimate of of DNAmAge
DNAmAgeBrainTraining=predict(glmnet.Training,t(betasTrain),type="response",s=lambda.glmnet.Training)
DNAmAgeBrainTesting=predict(glmnet.Training,t(betasTest),type="response",s=lambda.glmnet.Training)

# ## transform ages back 
DNAmAgeBrainTraining[,1]<-anti.trafo(DNAmAgeBrainTraining[,1])
DNAmAgeBrainTesting[,1]<-anti.trafo(DNAmAgeBrainTesting[,1])

# 
write.csv(DNAmAgeBrainTraining,"BB_LBB2_Jaffe_ROSMAP_anti_traf_training_age_prediction.csv", row.names=F)
write.csv(DNAmAgeBrainTesting,"BB_LBB2_Jaffe_ROSMAP_anti_traf_testing_age_prediction.csv", row.names=F)


#################################################################################################################
################################# SAVE  #########################################################################
#################################################################################################################

save(glmnet.Training,lambda.glmnet.Training,  file = "BB_LBB2_Jaffe_ROSMAP_age0plus_age0-20AgeT.rda")
lambda.glmnet.Training # 0.02138613


#################################################################################################################
################################# EXTRACT THE COEFS  ############################################################
#################################################################################################################

tmp_coeffs <- coef(glmnet.Training.CV, s = "lambda.min")
myCoeff<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
dim(myCoeff) # 348   2 (347 DNAm probes + intercept)
head(myCoeff) 
write.table(myCoeff,"Brain_clock_BB_Jaffe_LBB2_ROSMAP_coefficients.txt", row.names=F, col.names=T, quote=F)