############ BUILDING  AN EPIGENETIC CLOCK SPECIFIC FOR THE BRAIN (NOT SPECIFIC WITHIN BRAIN REGIONS) - using 450K and EPIC data 
############ AUTHOR: GEMMA SHIREBY
############ DATE: 18-09-2019
############ USING GLMNET - PENALIZED REGRESSION METHODS 
############ based on methodology first used by S. Horvath (Horvath, S. (2013). DNA methylation age of human tissues and cell types. Genome biology, 14(10), 3156)

## LOAD LIBRARIES
library(glmnet)

### LOAD DATA 
setwd("")
load("data.rdat") # pheno betas.dasen

################### READ IN TRAINING

train<-read.csv("training_sample.csv")

################### READ IN TESTING

test<-read.csv("testing_sample.csv")
## get some details on the samples 



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

#################################################################################################################
################################################ MATCH TO BETAS  ################################################
#################################################################################################################

# training
betasTrain<-betas.dasen[,match(train$Sentrix_Full, colnames(betas.dasen))]
dim(betasTrain) 
# testing
betasTest<-betas.dasen[,match(test$Sentrix_Full, colnames(betas.dasen))]
dim(betasTest) 


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
write.csv(DNAmAgeBrainTraining,"filename.csv", row.names=F)
write.csv(DNAmAgeBrainTesting,"filenamen.csv", row.names=F)


#################################################################################################################
################################# SAVE  #########################################################################
#################################################################################################################

save(glmnet.Training,lambda.glmnet.Training,  file = "file.rda")
lambda.glmnet.Training # 0.02138613


#################################################################################################################
################################# EXTRACT THE COEFS  ############################################################
#################################################################################################################

tmp_coeffs <- coef(glmnet.Training.CV, s = "lambda.min")
myCoeff<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
dim(myCoeff) # 348   2 (347 DNAm probes + intercept)
head(myCoeff) 
write.table(myCoeff,"file.txt", row.names=F, col.names=T, quote=F)
