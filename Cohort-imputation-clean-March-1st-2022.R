# COVID-19 UMAP-assisted clustering analysis pipeline - Imputation

# David Greenwood - Univeristy of Birmingham 2022 - drg707@student.bham.ac.uk, wdsquared.projects@gmail.com

#----------------------------------------------------------------#

setwd("/path_to/project_folder/")

require(mice)
require(DataExplorer)

# Import formatted data
{
  # Incomplete data, wide format
  
  # Development cohort
  dev.dat<- read.csv("/path_to/project_folder/csv/development_data.csv")
  
  # Validation cohort wave 1
  w1.dat<- read.csv("/path_to/project_folder/csv/wave_1_data.csv")
  
  # Validation cohort wave 2
  w2.dat<- read.csv("/path_to/project_folder/csv/wave_2_data.csv")
  
}

# Development data 
{
  # Variables for imputation
  {
    
    imp_pred <- c("age", "charlsondementiaicd10", "cancer", "CVD", "systolicbp", "diastolicbp", "heartrate",
                  "temperature", "resps", "o2sats", "Alt", "CRP", "Hion", "HCO3", "ur2", "egfr", "Hb", "lymphocytes2",
                  "baseexc","lactate2", "cough", "fever", "delirium", "GCS", "sex", "bmi", "neutlymph_ratio",
                  "pocfio2", "chestxray_r", "frailty_cts", "COPD_SleepApnea_Asthma", "dm_WithAndWithoutComp")
    
    # Filter
    clean <- dev.dat[,imp_pred]
    
    # Quantify missingness
    profile_missing(clean)
    tmp<-data.frame(round(cbind(missing=colSums(clean == 'Missing'|clean == 'missing', na.rm=T), 
                                na=colSums(is.na(clean), na.rm=T),
                                blank=colSums('' == (clean), na.rm=T))/nrow(clean),
                          4))
    tmp
    
    # Remove >40% NA variables
    remove<- rownames(tmp[tmp$missing>0.4|tmp$na>0.4|tmp$blank>0.4,])
    imp_pred<-imp_pred[imp_pred%in%remove==F]
    
    # Specifu categorical & continuous
    cat_imp_pred <- c("charlsondementiaicd10","cancer","CVD","cough","fever","delirium","sex",
                      "chestxray_r","dm_WithAndWithoutComp","COPD_SleepApnea_Asthma")
    cts_imp_pred <- imp_pred[imp_pred%in%cat_imp_pred==F]
    
    # Factorise categorical, make continuous numeric
    clean[,colnames(clean)%in%cat_imp_pred] <- data.frame(lapply(clean[,colnames(clean)%in%cat_imp_pred],factor),stringsAsFactors = F)
    clean[,colnames(clean)%in%cts_imp_pred] <- data.frame(lapply(clean[,colnames(clean)%in%cts_imp_pred],function(x){as.numeric(as.character(x))}),stringsAsFactors = F)
    
    # Summary 
    summary(clean,exclude=F)
    str(clean)
  }
  
  # Impute 
  miDat.Ext <- mice(clean, m=5,maxit=50,meth='pmm',seed=500)
  
  # Ensure no logged events 
  is.null(miDat.Ext$loggedEvents)
  
  # Create dir
  if(dir.exists("/path_to/project_folder/imputed/")==F){
    dir.create("/path_to/project_folder/imputed/")
  }
  
  # Save 
  #saveRDS(miDat.Ext,file="/path_to/project_folder/imputed/development_miData.RDS")
}

# Validation data - Wave 1 
{
  # Variables for imputation
  {
    #cat(paste0(colnames(datToImpute),collapse = "\",\""))
    #cat(paste0(gsub("_0","",colnames(datToImpute)),collapse = "\",\""))
    
    imp_pred<- c("age","sex","ethnic_r","GCS","frailty_cts","liver_dis","charlsondementiaicd10","cancer","charlsonscopdicd10","sleep_apnoea","asthma","CVD","hypertension","dm_nocomplications","charlsonsdmcomplicationsicd10","charlsonspepticulcericd10","rheum_inflam_multisys","thyroid","breathlessness","chestpain","cough","fever","headache","malaise","newdiarrhoeaorvomiting","sputum","delirium","chestxray_r","aniongp_0","temperature_0","heartrate_0","systolicbp_0","diastolicbp_0","resps_0","glucose_0","Hb_0","o2sats_0","pCO2_0","HCO3_0","WBC_0","lymphocytes2_0","CRP_0","albumin2_0","platelets2_0","bilirubin2_0","Alt_0","MCV_0","alkP_0","ur2_0","K_0","Na_0","correctedca_0","RDW_0","monocytes2_0","eosinophil_0","baseexc_0","pocfio2_0","egfr_0","bmi_0","Hion_0","lactate2_0","HCT_0")
    
    # Filter
    clean <- w1.dat[,imp_pred]
    
    # Quantify missingness
    profile_missing(clean)
    tmp<-data.frame(round(cbind(missing=colSums(clean == 'Missing'|clean == 'missing', na.rm=T), 
                                na=colSums(is.na(clean), na.rm=T),
                                blank=colSums('' == (clean), na.rm=T))/nrow(clean),
                          4))
    tmp
    
    # Remove >40% NA variables
    #remove<- rownames(tmp[tmp$missing>0.4|tmp$na>0.4|tmp$blank>0.4,])
    #imp_pred<-imp_pred[imp_pred%in%remove==F]
    
    # Specifu categorical & continuous
    cat_imp_pred<- c("sex","ethnic_r","liver_dis","charlsondementiaicd10","cancer","charlsonscopdicd10","sleep_apnoea","asthma","CVD","hypertension","dm_nocomplications","charlsonsdmcomplicationsicd10","charlsonspepticulcericd10","rheum_inflam_multisys","thyroid","breathlessness","chestpain","cough","fever","headache","malaise","newdiarrhoeaorvomiting","sputum","delirium","chestxray_r")
    cts_imp_pred <- imp_pred[imp_pred%in%cat_imp_pred==F]
    
    # Factorise categorical, make continuous numeric
    clean[,colnames(clean)%in%cat_imp_pred] <- data.frame(lapply(clean[,colnames(clean)%in%cat_imp_pred],factor),stringsAsFactors = F)
    clean[,colnames(clean)%in%cts_imp_pred] <- data.frame(lapply(clean[,colnames(clean)%in%cts_imp_pred],function(x){as.numeric(as.character(x))}),stringsAsFactors = F)
    
    # Summary 
    summary(clean,exclude=F)
    str(clean)
  }
  
  # Impute 
  miDat.w1 <- mice(clean, m=5,maxit=50,meth='pmm',seed=500)
  
  # Ensure no logged events 
  is.null(miDat.w1$loggedEvents)
  
  # Create dir
  if(dir.exists("/path_to/project_folder/imputed/")==F){
    dir.create("/path_to/project_folder/imputed/")
  }
  
  # Save 
  #saveRDS(miDat.Ext,file="/path_to/project_folder/imputed/wave1_miData.RDS") 
}

# Validation data - Wave 2
{
  # Variables for imputation
  {
    imp_pred<-c("age","sex","ethnic_r","GCS","frailty_cts","liver_dis","charlsondementiaicd10","cancer","charlsonscopdicd10","sleep_apnoea","asthma","CVD","hypertension","dm_nocomplications","charlsonsdmcomplicationsicd10","charlsonspepticulcericd10","rheum_inflam_multisys","thyroid","breathlessness","chestpain","cough","fever","headache","malaise","newdiarrhoeaorvomiting","sputum","delirium","aniongp_0","temperature_0","heartrate_0","systolicbp_0","diastolicbp_0","resps_0","glucose_0","Hb_0","o2sats_0","pCO2_0","HCO3_0","WBC_0","lymphocytes2_0","CRP_0","albumin2_0","platelets2_0","bilirubin2_0","Alt_0","MCV_0","alkP_0","ur2_0","K_0","Na_0","correctedca_0","RDW_0","monocytes2_0","eosinophil_0","baseexc_0","pocfio2_0","egfr_0","bmi_0","Hion_0","lactate2_0","HCT_0")
    
    # Filter
    clean <- w2.dat[,imp_pred]
    
    # Quantify missingness
    tmp<-profile_missing(clean)
    tmp[order(tmp$pct_missing),]
    tmp<-data.frame(round(cbind(missing=colSums(clean == 'Missing'|clean == 'missing', na.rm=T), 
                                na=colSums(is.na(clean), na.rm=T),
                                blank=colSums('' == (clean), na.rm=T))/nrow(clean),
                          4))
    tmp
    
    # Remove >40% NA variables
    #remove<- rownames(tmp[tmp$missing>0.4|tmp$na>0.4|tmp$blank>0.4,])
    #imp_pred<-imp_pred[imp_pred%in%remove==F]
    
    # Specifu categorical & continuous
    cat_imp_pred<- c("sex","ethnic_r","liver_dis","charlsondementiaicd10","cancer","charlsonscopdicd10","sleep_apnoea","asthma","CVD","hypertension","dm_nocomplications","charlsonsdmcomplicationsicd10","charlsonspepticulcericd10","rheum_inflam_multisys","thyroid","breathlessness","chestpain","cough","fever","headache","malaise","newdiarrhoeaorvomiting","sputum","delirium","chestxray_r")
    cts_imp_pred <- imp_pred[imp_pred%in%cat_imp_pred==F]
    
    # Factorise categorical, make continuous numeric
    clean[,colnames(clean)%in%cat_imp_pred] <- data.frame(lapply(clean[,colnames(clean)%in%cat_imp_pred],factor),stringsAsFactors = F)
    clean[,colnames(clean)%in%cts_imp_pred] <- data.frame(lapply(clean[,colnames(clean)%in%cts_imp_pred],function(x){as.numeric(as.character(x))}),stringsAsFactors = F)
    
    # Summary 
    summary(clean,exclude=F)
    str(clean)
  }
  
  # Impute 
  miDat.w2 <- mice(clean, m=5,maxit=50,meth='pmm',seed=500)
  
  # Ensure no logged events 
  is.null(miDat.w2$loggedEvents)
  
  # Create dir
  if(dir.exists("/path_to/project_folder/imputed/")==F){
    dir.create("/path_to/project_folder/imputed/")
  }
  
  # Save 
  #saveRDS(miDat.Ext,file="/path_to/project_folder/imputed/wave2_miData.RDS") 
}