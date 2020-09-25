#########################################################################################################################
###   This function is used for predicting DNA methylation age in cortical samples                                    ###
###   It was built using cortical DNA methylation array data (Illumina 450K)                                    ###       
###   Coefficients of the predictor are based on 1037 multi-region cortical training samples		                      ###
###	  Age predictors were identified using elastic net regression (glmnet in R) - methodology adapted from:           ###
###            Horvath, S. (2013). DNA methylation age of human tissues and cell types. Genome biology, 14(10), 3156. ###                                      
###	  Author: Gemma Shireby                                                                                           ###
###   Email: gs470@exeter.ac.uk				                                                                                ###
#########################################################################################################################



CorticalClock<-function(betas, ## betas = betas matrix (rownames=cpgs, colnames=IDs)
                        pheno, ##  pheno file = file which contains IDs that match betas IDs and  contains actual Age col 
                        dir, ## directory where the ref file and coeffecients are saved 
                        IDcol, ## ID column which matches Betas IDs for your samples
                        Agecol){  ## Age column 
                        
    ## subset pheno                    
    pheno<-pheno[c(IDcol,Agecol)]
    colnames(pheno)<-c("ID","Age")                 
 
  ## check ggplot2 is loaded - if not, load, if not installd then install
  pkgTest <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }
  
  pkgTest("ggplot2")
  
  ## ensure colnames are ID and Age 
  colnames(pheno)<-c("ID", "Age")
  
  
  #########################################################
  ############# read in cortical clock coeffs #############
  #########################################################
  
  braincoef<-read.table(paste0(dir,"CorticalClockCoefs.txt",sep=""),stringsAsFactor=F,header=T)
  braincoef$probe<-as.character(braincoef$probe)
  
  #########################################################
  ### find the overlap between the probes and your data ###
  #########################################################
  
  overlap<-braincoef[which(braincoef$probe %in% rownames(betas)),]
  if (nrow(overlap) < nrow(braincoef) ){
    print("some probes are missing, we will need to impute values here - the final predictions will be less accurate")
  } else {
    print("all the probes overlap between your data and the clock probes - no need for imputing missing values")
  }
  
  ##############################################################################################################################
  ############# add reference betas to the missing probes (average DNAm across 700 control cortical samples)                 ###
  ### imputation method adapted from:  https://github.com/qzhang314/DNAm-based-age-predictor                                 ###
  ### cite: Zhang, Q., Vallerga, C. L., Walker, R. M., Lin, T., Henders, A. K., Montgomery, G. W., ... & Pitcher, T. (2019). ###
  ### Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing.             ###
  ### Genome medicine, 11(1), 1-11.                                                                                          ###
  ############################################################################################################################## 
  
  if (length(overlap) < nrow(braincoef)) {
    ## transform betas to be cpg in col
    betas<-t(betas)
    
    ###########  read in ref data and match
    load(paste0(dir,"Ref_DNAm_brain_values.rdat",sep=""))
    ref<-ref[which(names(ref) %in% braincoef$probe) , drop=F]
    
    
    betas<-betas[,colnames(betas)%in%names(ref)]
    if(ncol(betas)<length(ref)){
      missprobe<-setdiff(names(ref),colnames(betas))
      refmiss<-ref[missprobe]
      refmiss<-matrix(refmiss,ncol=length(missprobe),nrow=nrow(betas),byrow=T)
      refmiss<-as.data.frame(refmiss)
      colnames(refmiss)<-missprobe
      rownames(refmiss)<-rownames(betas)
      betas<-cbind(betas,refmiss)
    }
    
    betas<-betas[,names(ref)]     ########### match betas
    
    #################################################################################################
    ##############  for each probe replace missing value with mean value from reference betas #######
    #################################################################################################
    
    ## impute function 
    imputeNA<-function(betas){
      betas[is.na(betas)]<-mean(betas,na.rm=T)
      return(betas)
    }
    
    ## apply function 
    betasNona<-apply(betas,2,function(x) imputeNA(x))  
    
    
    ## tranform betas - CpG in row
    betas<-t(betasNona)
    
    ##############################################################################################
    ############# Age prediciton! - weighted sum of coefficients plus the intercept ##############
    ##############################################################################################
    
    braincoef<-braincoef[match(rownames(betas), braincoef$probe),]
    brainpred<-braincoef$coef%*%betas+0.577682570446177
    
    ##############################################################################################
    ### anti transform the results (accounting for logarithmic relationship in ages 0-20)      ###
    ### Same as Horvath's: Horvath, S. (2013).                                                 ###
    ### DNA methylation age of human tissues and cell types. Genome biology, 14(10), 3156.     ###                                                              #######
    ##############################################################################################
    
    anti.trafo<-function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
    brainpred<-anti.trafo(brainpred)
    
    #################################################
    #############  Save brain predictions ###########
    #################################################
    pheno<-pheno[match(colnames(betas), pheno$ID),]
    pheno$brainpred<-as.numeric(brainpred)
    write.csv(pheno, "CorticalPred.csv")
    
    ## if not impuation for missing values:
    
  } else {
    
    ##############################################################################################
    ############# Age prediciton! - weighted sum of coefficients plus the intercept ##############
    ##############################################################################################
    
    braincoef<-braincoef[match(rownames(betas), braincoef$probe),]
    brainpred<-braincoef$coef%*%betas+0.577682570446177
    
    ##############################################################################################
    ### anti transform the results (accounting for logarithmic relationship in ages 0-20)      ###
    ### Same as Horvath's: Horvath, S. (2013).                                                 ###
    ### DNA methylation age of human tissues and cell types. Genome biology, 14(10), 3156.     ###                                                              #######
    ##############################################################################################
    
    anti.trafo<-function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
    brainpred<-anti.trafo(brainpred)
    
    #################################################
    #############  Save brain predictions ###########
    #################################################
    pheno<-pheno[match(colnames(betas), pheno$ID),]
    pheno$brainpred<-as.numeric(brainpred)
    pheno$Age<-as.numeric(pheno$Age)
    write.csv(pheno, "CorticalPred.csv")
    
  }
  
  ############################################     
  ### get some stats for the predictions #####
  ############################################
  compareStats<-function(data){
    colnames(data)<-c("ID","Age","Cortical Clock")
    data<-data[complete.cases(data$Age),]
    corr<-round(cor(data[,2],data[,3]),2) ## correlation - pearsons r
    RMSE<-function(actualAge,estimatedAge){ ## root mean squared error (years)
          sqrt(mean((actualAge-estimatedAge)^2))
            }   
    rmse<-round(RMSE(data[,2],data[,3]),2)
    mad<-round(median(abs(data[,2]-data[,3])),2) ## mean absoluate deviation (years)
    stats<-matrix(ncol=1,nrow=3)
    colnames(stats)<-c("Cortical Clock")
    rownames(stats)<-c("Correlation (r)","RMSE (years)","MAD (years)")
    stats[1,]<-corr
    stats[2,]<-rmse
    stats[3,]<-mad
    print(stats)
    write.csv(stats, "accuracy_statistics.csv")
  }
  
  compareStats(pheno)
  
  
  #####################
  ### make a plot #####
  #####################
  
  ## uses ggplot
  plotAge<-function(data){
    ggplot(data=data,aes(x=Age,y=brainpred)) +
      geom_abline(intercept=0,slope=1) +
      geom_point()+
      xlab("Chronological Age")+
      ylab("Predicted Age")+
      ggtitle("Cortical Clock") +
      theme_minimal() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title=element_text(size=16,face="bold"),
            legend.text=element_text(size=14,face="bold"),
            plot.title=element_text(size=16, face="bold.italic"),
            axis.text=element_text(size=14,face="bold"),
            strip.text.x=element_text(face="bold",size=16))
    
    ggsave("CorticalClockplot.pdf", height= 5, width=4)
  }
  
  plotAge(pheno)
  
  print("Your Cortical DNAm ages are saved in the file 'CorticalPred.csv'")
  print("Accuracy statistics comparing DNAm age and chronological age are saved in the file 'accuracy_statistics.csv'")
  print("Plot of Cortical DNAm age against age is saved in the file 'CorticalClockplot.pdf'")
  
}
