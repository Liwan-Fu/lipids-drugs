library(TwoSampleMR)
library(MRPRESSO)
library(LDlinkR)
library(arrow)
library(dplyr)
library(vautils)
library(xlsx)
setwd('...')
ao = read.csv('ao20230210.csv',header = T)

drug=read.xlsx('drug.xlsx',sheetName = 'drug1used')

GWASid=c('ebi-a-GCST003116',
         'ieu-a-798',
         'ebi-a-GCST006414',
         'ebi-a-GCST006906',
         'ebi-a-GCST006867',
         'ukb-b-14057',
         'ukb-b-6720',
         'ukb-d-I9_DVTANDPULM',
         'ukb-d-I9_IHD',
         'ukb-d-I9_CHD',
         'finn-b-I9_MI',
         'finn-b-I9_AF',
         'finn-b-I9_HEARTFAIL',
         'finn-b-C_STROKE',
         'finn-b-E4_DM2',
         'finn-b-I9_HYPTENS',
         'finn-b-I9_VARICVE',
         'finn-b-I9_DVTANDPULM',
         'finn-b-I9_ISCHHEART',
         'finn-b-I9_CHD',
         'ebi-a-GCST007557',
         'ebi-a-GCST90002409',
         'ieu-a-79',
         'ieu-a-67',
         'ieu-a-2',
         'ukb-b-8909',
         'ukb-b-19393',
         'ukb-b-13354',
         'ebi-a-GCST90002232',
         'ebi-a-GCST000569',
         'ebi-a-GCST90002238',
         'ieu-b-118',
         'ieu-b-38',
         'ieu-b-39')

outcomes <- filter(ao,id %in% GWASid) %>%
  arrange(desc(year), desc(nsnp))
outcomes=outcomes[match(GWASid,outcomes$id),]

outcome_dat=extract_outcome_data(drug$SNP, outcomes$id)

index=which(outcomes[,'sd']>0)
temp <- outcomes[index,] %>% select(id, sd)
outcome_dat <- merge(outcome_dat, temp, by.x="id.outcome", by.y="id", all.x=TRUE)
index <- !is.na(outcome_dat$sd)
outcome_dat$beta.outcome[index] <- outcome_dat$beta.outcome[index] / outcome_dat$sd[index]
outcome_dat$se.outcome[index] <- outcome_dat$se.outcome[index] / outcome_dat$sd[index]

dat=harmonise_data(drug,outcome_dat,action = 2)
dat=dat[!duplicated(dat),]

expos=unique(drug$exposure)
outco=unique(dat$outcome)[match(GWASid,sapply(strsplit(unique(dat$outcome), 'id:',fixed = T), '[',2))]

# IVW analysis
library(MendelianRandomization)
res_IVW=matrix(NA,length(expos)*length(outco),10)
for(i in 1:length(expos)){
  for(j in 1:length(outco)){
    temp=dat %>% filter(exposure==expos[i]&outcome==outco[j]
    ) %>% select(SNP,beta.exposure,se.exposure,
                 beta.outcome,se.outcome)
    if (nrow(temp)==0) {
      res_IVW[(i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
    } else {
      LDinfo=LDmatrix(snps = temp$SNP, 
                      pop = "EUR", r2d = "r2", 
                      token = '261fbd8cb0c9', 
                      file =FALSE)
      temp=temp[match(colnames(LDinfo)[-1],temp$SNP),]
      MRInput=mr_input(bx=temp[,'beta.exposure'],
                       bxse=temp[,'se.exposure'],
                       by=temp[,'beta.outcome'],
                       byse=temp[,'se.outcome'],
                       correlation = as.matrix(LDinfo[,-1]))
      IVW_res=mr_ivw(MRInput,model= "default",correl = T)
      res_IVW[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                        IVW_res@Estimate,IVW_res@StdError,
                                        IVW_res@CILower,IVW_res@CIUpper,
                                        IVW_res@Pvalue,IVW_res@Heter.Stat)
    }
    
  }
}
colnames(res_IVW)=c('exposure','outcome','SNPs',
                    'ES','SE','CILower','CIUpper','P value','Q','Q.P value')


# WME MR-Egger and MR-PRESSO analysis
res_WME=matrix(NA,length(expos)*length(outco),7)
res_Egger=matrix(NA,length(expos)*length(outco),8)
res_MBE=matrix(NA,length(expos)*length(outco),7)
res_PRESSO=matrix(NA,length(expos)*length(outco),9)
for(i in 1:length(expos)){
  for(j in 1:length(outco)){
    temp=dat %>% filter(exposure==expos[i]&outcome==outco[j]
    ) %>% select(SNP,beta.exposure,se.exposure,
                 beta.outcome,se.outcome)
    if (nrow(temp)<=1) {
      res_WME[(i-1)*length(outco)+j,1:3]=res_Egger[
        (i-1)*length(outco)+j,1:3]=res_PRESSO[
          (i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
    } else {
      LDinfo=LDmatrix(snps = temp$SNP, 
                      pop = "EUR", r2d = "r2", 
                      token = '261fbd8cb0c9', 
                      file =FALSE)
      temp=temp[match(colnames(LDinfo)[-1],temp$SNP),]
      MRInput=mr_input(bx=temp[,'beta.exposure'],
                       bxse=temp[,'se.exposure'],
                       by=temp[,'beta.outcome'],
                       byse=temp[,'se.outcome'],
                       correlation = as.matrix(LDinfo[,-1]))
      # WME analysis
      WME_res=mr_median(MRInput,weighting = 'penalized',distribution = 'normal',
                        alpha = 0.05,iterations = 1000)
      res_WME[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                        WME_res@Estimate,WME_res@CILower,
                                        WME_res@CIUpper,WME_res@Pvalue)
      # MR-Egger analysis
      Egger_res=mr_egger(MRInput,correl = T)
      res_Egger[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                          Egger_res@Estimate,Egger_res@CILower.Est,
                                          Egger_res@CIUpper.Est,Egger_res@Pvalue.Est,
                                          Egger_res@Pvalue.Int)
      # Mode-based estimate method
      MBE_res=mr_mbe(MRInput)
      res_MBE[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                        MBE_res@Estimate,MBE_res@CILower,
                                        MBE_res@CIUpper,MBE_res@Pvalue)
      # MR-PRESSO analysis
      datasummary=data.frame(bx=temp[,'beta.exposure'],
                             bxse=temp[,'se.exposure'],
                             by=temp[,'beta.outcome'],
                             byse=temp[,'se.outcome'])
      if (nrow(datasummary)<=3) {
        res_PRESSO[(i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
      } else {
        MR_PRESSO=mr_presso(BetaOutcome = "by", BetaExposure = "bx", 
                            SdOutcome = "byse", SdExposure = "bxse", 
                            OUTLIERtest = T, DISTORTIONtest = T,
                            data = datasummary, NbDistribution = 1000,  SignifThreshold = 0.05)
        if(is.na(MR_PRESSO[["Main MR results"]][2,'P-value'])) {
          MR_PRESSO_res=as.numeric(MR_PRESSO[["Main MR results"]][1,c('Causal Estimate','Sd','P-value')])
          res_PRESSO[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                               0,
                                               MR_PRESSO[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]],
                                               NA,
                                               MR_PRESSO_res)
        } else {
          MR_PRESSO_res=as.numeric(MR_PRESSO[["Main MR results"]][2,c('Causal Estimate','Sd','P-value')])
          res_PRESSO[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                               length(MR_PRESSO[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]),
                                               MR_PRESSO[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]],
                                               MR_PRESSO[["MR-PRESSO results"]][["Distortion Test"]][["Pvalue"]],
                                               MR_PRESSO_res)
        }
      }
    }
  }
}
colnames(res_WME)=c('exposure','outcome','SNPs',
                    'ES','CILower','CIUpper','P value')
colnames(res_Egger)=c('exposure','outcome','SNPs',
                      'ES','CILower','CIUpper','P value','P_intercept')
colnames(res_MBE)=c('exposure','outcome','SNPs',
                    'ES','CILower','CIUpper','P value')
colnames(res_PRESSO)=c('exposure','outcome','SNPs',
                       'Outliers','p_glo','p_dis','ES','SE','P value')

write.xlsx(res_IVW,'Population results.xlsx',sheetName = 'IVW',row.names = F)

write.xlsx(res_WME,'Population results.xlsx',
           sheetName = 'WME',row.names = F,append = T)

write.xlsx(res_Egger,'Population results.xlsx',
           sheetName = 'Egger',row.names = F,append = T)

write.xlsx(res_MBE,'Population results.xlsx',
           sheetName = 'MBE',row.names = F,append = T)

write.xlsx(res_PRESSO,'Population results.xlsx',
           sheetName = 'PRESSO',row.names = F,append = T)


######## lipids MR analysis ########
lipids=extract_instruments(c('ieu-a-300','ieu-a-302','ieu-a-301'))
# F statistics calculated (PMID:32001714)
for(i in 1:nrow(lipids)){
  temp=lipids[i,]
  beta=temp$beta.exposure
  MAF=temp$eaf.exposure
  N=temp$samplesize.exposure
  SE=temp$se.exposure
  R2=2*(1-MAF)*MAF*beta^2/(
    2*(1-MAF)*MAF*beta^2+2*(1-MAF)*MAF*N*SE^2)
  lipids$F[i]=R2/(1-R2)*(N-2)
}
lipids=na.omit(lipids)
lipids=lipids[lipids$F>10,]

# IVW analysis for lipids
outcome_dat=extract_outcome_data(lipids$SNP, outcomes$id)

index=which(outcomes[,'sd']>0)
temp <- outcomes[index,] %>% select(id, sd)
outcome_dat <- merge(outcome_dat, temp, by.x="id.outcome", by.y="id", all.x=TRUE)
index <- !is.na(outcome_dat$sd)
outcome_dat$beta.outcome[index] <- outcome_dat$beta.outcome[index] / outcome_dat$sd[index]
outcome_dat$se.outcome[index] <- outcome_dat$se.outcome[index] / outcome_dat$sd[index]

dat=harmonise_data(lipids,outcome_dat,action = 2)
dat=dat[!duplicated(dat),]

expos=unique(lipids$exposure)
outco=unique(dat$outcome)[match(GWASid,sapply(strsplit(unique(dat$outcome), 'id:',fixed = T), '[',2))]

# IVW analysis
library(MendelianRandomization)
lipids_IVW=matrix(NA,length(expos)*length(outco),10)
for(i in 1:length(expos)){
  for(j in 1:length(outco)){
    temp=dat %>% filter(exposure==expos[i]&outcome==outco[j]
    ) %>% select(SNP,beta.exposure,se.exposure,
                 beta.outcome,se.outcome)
    if (nrow(temp)==0) {
      lipids_IVW[(i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
    } else {
      MRInput=mr_input(bx=temp[,'beta.exposure'],
                       bxse=temp[,'se.exposure'],
                       by=temp[,'beta.outcome'],
                       byse=temp[,'se.outcome'])
      IVW_res=mr_ivw(MRInput,model= "default")
      lipids_IVW[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                           IVW_res@Estimate,IVW_res@StdError,
                                           IVW_res@CILower,IVW_res@CIUpper,
                                           IVW_res@Pvalue,IVW_res@Heter.Stat)
    }
    
  }
}
colnames(lipids_IVW)=c('exposure','outcome','SNPs',
                       'ES','SE','CILower','CIUpper','P value','Q','Q.P value')

### Sensitivity analysis for lipids
lipids_WME=matrix(NA,length(expos)*length(outco),7)
lipids_Egger=matrix(NA,length(expos)*length(outco),8)
lipids_MBE=matrix(NA,length(expos)*length(outco),7)
lipids_PRESSO=matrix(NA,length(expos)*length(outco),9)
for(i in 1:length(expos)){
  for(j in 1:length(outco)){
    temp=dat %>% filter(exposure==expos[i]&outcome==outco[j]
    ) %>% select(SNP,beta.exposure,se.exposure,
                 beta.outcome,se.outcome)
    if (nrow(temp)<=1) {
      lipids_WME[(i-1)*length(outco)+j,1:3]=lipids_Egger[
        (i-1)*length(outco)+j,1:3]=lipids_PRESSO[
          (i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
    } else {
      MRInput=mr_input(bx=temp[,'beta.exposure'],
                       bxse=temp[,'se.exposure'],
                       by=temp[,'beta.outcome'],
                       byse=temp[,'se.outcome'])
      # WME analysis
      WME_res=mr_median(MRInput,weighting = 'penalized',distribution = 'normal',
                        alpha = 0.05,iterations = 1000)
      lipids_WME[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                        WME_res@Estimate,WME_res@CILower,
                                        WME_res@CIUpper,WME_res@Pvalue)
      # MR-Egger analysis
      Egger_res=mr_egger(MRInput)
      lipids_Egger[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                          Egger_res@Estimate,Egger_res@CILower.Est,
                                          Egger_res@CIUpper.Est,Egger_res@Pvalue.Est,
                                          Egger_res@Pvalue.Int)
      # Mode-based estimate method
      MBE_res=mr_mbe(MRInput)
      lipids_MBE[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                        MBE_res@Estimate,MBE_res@CILower,
                                        MBE_res@CIUpper,MBE_res@Pvalue)
      # MR-PRESSO analysis
      datasummary=data.frame(bx=temp[,'beta.exposure'],
                             bxse=temp[,'se.exposure'],
                             by=temp[,'beta.outcome'],
                             byse=temp[,'se.outcome'])
      if (nrow(datasummary)<=3) {
        res_PRESSO[(i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
      } else {
        MR_PRESSO=mr_presso(BetaOutcome = "by", BetaExposure = "bx", 
                            SdOutcome = "byse", SdExposure = "bxse", 
                            OUTLIERtest = T, DISTORTIONtest = T,
                            data = datasummary, NbDistribution = 1000,  SignifThreshold = 0.05)
        if(is.na(MR_PRESSO[["Main MR results"]][2,'P-value'])) {
          MR_PRESSO_res=as.numeric(MR_PRESSO[["Main MR results"]][1,c('Causal Estimate','Sd','P-value')])
          lipids_PRESSO[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                               0,
                                               MR_PRESSO[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]],
                                               NA,
                                               MR_PRESSO_res)
        } else {
          MR_PRESSO_res=as.numeric(MR_PRESSO[["Main MR results"]][2,c('Causal Estimate','Sd','P-value')])
          lipids_PRESSO[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                               length(MR_PRESSO[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]),
                                               MR_PRESSO[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]],
                                               MR_PRESSO[["MR-PRESSO results"]][["Distortion Test"]][["Pvalue"]],
                                               MR_PRESSO_res)
        }
      }
    }
  }
}
colnames(lipids_WME)=c('exposure','outcome','SNPs',
                    'ES','CILower','CIUpper','P value')
colnames(lipids_Egger)=c('exposure','outcome','SNPs',
                      'ES','CILower','CIUpper','P value','P_intercept')
colnames(lipids_MBE)=c('exposure','outcome','SNPs',
                    'ES','CILower','CIUpper','P value')
colnames(lipids_PRESSO)=c('exposure','outcome','SNPs',
                       'Outliers','p_glo','p_dis','ES','SE','P value')

write.xlsx(lipids_IVW,'lipids results.xlsx',sheetName = 'IVW',row.names = F)

write.xlsx(lipids_WME,'lipids results.xlsx',
           sheetName = 'WME',row.names = F,append = T)

write.xlsx(lipids_Egger,'lipids results.xlsx',
           sheetName = 'Egger',row.names = F,append = T)

write.xlsx(lipids_MBE,'lipids results.xlsx',
           sheetName = 'MBE',row.names = F,append = T)

write.xlsx(lipids_PRESSO,'lipids results.xlsx',
           sheetName = 'PRESSO',row.names = F,append = T)


######## MR analysis for GTEx ##########
GTEx_exposure=read.xlsx('GTEx_exposure_all.xlsx',1,header = T) %>% 
  filter(p_adjust<=0.05)
exposure_dat=apply(GTEx_exposure,1,function(x) {
  if (as.numeric(x['slope'])>0) {
    c(x[1:2],x['alt'],x['ref'],x[5:6],1-as.numeric(x['maf']),-as.numeric(x['slope']),x[9:11])
  } else {
    x
  }
})
exposure_dat=data.frame(t(exposure_dat),stringsAsFactors = F)
colnames(exposure_dat)=colnames(GTEx_exposure)
exposure_dat=cbind(exposure_dat[,c(1,3:6)],
                   as.data.frame(lapply(exposure_dat[,c(2,7:11)],as.numeric)))
rm(GTEx_exposure)
exposure_dat=format_data(exposure_dat,
                         type = "exposure",
                         phenotype_col = "pheno",
                         snp_col = "rs_id_dbSNP151_GRCh38p7",
                         beta_col = "slope",
                         se_col = "slope_se",
                         effect_allele_col = "alt",
                         other_allele_col = "ref",
                         eaf_col = "maf",
                         pval_col='pval_nominal',
                         chr_col = "chr",
                         pos_col = "variant_pos")
exposure_dat$id.exposure=exposure_dat$exposure

GWASid=c('ieu-a-798',
         'ukb-d-I9_IHD',
         'ukb-d-I9_CHD')
outcomes <- filter(ao,id %in% GWASid) %>%
  arrange(desc(year), desc(nsnp))
outcomes=outcomes[match(GWASid,outcomes$id),]

outcome_dat=extract_outcome_data(exposure_dat$SNP, outcomes$id)

index=which(outcomes[,'sd']>0)
temp <- outcomes[index,] %>% select(id, sd)
outcome_dat <- merge(outcome_dat, temp, by.x="id.outcome", by.y="id", all.x=TRUE)
index <- !is.na(outcome_dat$sd)
outcome_dat$beta.outcome[index] <- outcome_dat$beta.outcome[index] / outcome_dat$sd[index]
outcome_dat$se.outcome[index] <- outcome_dat$se.outcome[index] / outcome_dat$sd[index]

dat=harmonise_data(exposure_dat,outcome_dat,action = 2)
dat=dat[!duplicated(dat),]

expos=unique(exposure_dat$exposure)
outco=unique(dat$outcome)[match(GWASid,sapply(strsplit(unique(dat$outcome), 'id:',fixed = T), '[',2))]

# IVW analysis
library(MendelianRandomization)

# r2=0.2 analysis
res_IVW=matrix(NA,length(expos)*length(outco),10)
for(i in 1:length(expos)){
  for(j in 1:length(outco)){
    temp=dat %>% filter(exposure==expos[i]&outcome==outco[j]
    ) %>% select(chr,pos,SNP,beta.exposure,se.exposure,
                 beta.outcome,se.outcome,pval.exposure,pval.outcome)
    temp=clump_data(temp,clump_r2 = 0.2)
    if (nrow(temp)==0) {
      res_IVW[(i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
    } else {
      LDinfo=LDmatrix(snps = temp$SNP, 
                      pop = "EUR", r2d = "r2", 
                      token = '261fbd8cb0c9', 
                      file =FALSE)
      temp=temp[match(colnames(LDinfo)[-1],temp$SNP),]
      MRInput=mr_input(bx=temp[,'beta.exposure'],
                       bxse=temp[,'se.exposure'],
                       by=temp[,'beta.outcome'],
                       byse=temp[,'se.outcome'],
                       correlation = as.matrix(LDinfo[,-1]))
      IVW_res=mr_ivw(MRInput,model= "default",correl = T)
      res_IVW[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                        IVW_res@Estimate,IVW_res@StdError,
                                        IVW_res@CILower,IVW_res@CIUpper,
                                        IVW_res@Pvalue,IVW_res@Heter.Stat)
    }
    
  }
}
colnames(res_IVW)=c('exposure','outcome','SNPs',
                    'ES','SE','CILower','CIUpper','P value','Q','Q.P value')
IVW0.2=res_IVW

# r2=0.01 analysis
res_IVW=matrix(NA,length(expos)*length(outco),10)
for(i in 1:length(expos)){
  for(j in 1:length(outco)){
    temp=dat %>% filter(exposure==expos[i]&outcome==outco[j]
    ) %>% select(chr,pos,SNP,beta.exposure,se.exposure,
                 beta.outcome,se.outcome,pval.exposure)
    temp=clump_data(temp,clump_r2 = 0.01)
    if (nrow(temp)==0) {
      res_IVW[(i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
    } else {
      LDinfo=LDmatrix(snps = temp$SNP, 
                      pop = "EUR", r2d = "r2", 
                      token = '261fbd8cb0c9', 
                      file =FALSE)
      temp=temp[match(colnames(LDinfo)[-1],temp$SNP),]
      MRInput=mr_input(bx=temp[,'beta.exposure'],
                       bxse=temp[,'se.exposure'],
                       by=temp[,'beta.outcome'],
                       byse=temp[,'se.outcome'],
                       correlation = as.matrix(LDinfo[,-1]))
      IVW_res=mr_ivw(MRInput,model= "default",correl = T)
      res_IVW[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                        IVW_res@Estimate,IVW_res@StdError,
                                        IVW_res@CILower,IVW_res@CIUpper,
                                        IVW_res@Pvalue,IVW_res@Heter.Stat)
    }
    
  }
}
colnames(res_IVW)=c('exposure','outcome','SNPs',
                    'ES','SE','CILower','CIUpper','P value','Q','Q.P value')
IVW0.01=res_IVW


# r2=0.001 analysis
res_IVW=matrix(NA,length(expos)*length(outco),10)
for(i in 1:length(expos)){
  for(j in 1:length(outco)){
    temp=dat %>% filter(exposure==expos[i]&outcome==outco[j]
    ) %>% select(chr,pos,SNP,beta.exposure,se.exposure,
                 beta.outcome,se.outcome,pval.exposure)
    temp=clump_data(temp,clump_r2 = 0.001)
    if (nrow(temp)==0) {
      res_IVW[(i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
    } else {
      LDinfo=LDmatrix(snps = temp$SNP, 
                      pop = "EUR", r2d = "r2", 
                      token = '261fbd8cb0c9', 
                      file =FALSE)
      temp=temp[match(colnames(LDinfo)[-1],temp$SNP),]
      MRInput=mr_input(bx=temp[,'beta.exposure'],
                       bxse=temp[,'se.exposure'],
                       by=temp[,'beta.outcome'],
                       byse=temp[,'se.outcome'],
                       correlation = as.matrix(LDinfo[,-1]))
      IVW_res=mr_ivw(MRInput,model= "default",correl = T)
      res_IVW[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                        IVW_res@Estimate,IVW_res@StdError,
                                        IVW_res@CILower,IVW_res@CIUpper,
                                        IVW_res@Pvalue,IVW_res@Heter.Stat)
    }
    
  }
}
colnames(res_IVW)=c('exposure','outcome','SNPs',
                    'ES','SE','CILower','CIUpper','P value','Q','Q.P value')
IVW0.001=res_IVW

write.xlsx(IVW0.2,'GTEx IVW results.xlsx',sheetName = 'IVW0.2',row.names = F)

write.xlsx(IVW0.01,'GTEx IVW results.xlsx',
           sheetName = 'IVW0.01',row.names = F,append = T)

write.xlsx(IVW0.001,'GTEx IVW results.xlsx',
           sheetName = 'IVW0.001',row.names = F,append = T)


