###Download###
if(!require("remotes"))
  install.packages("remotes") # if necessary
library(remotes)
install_github("chr1swallace/coloc@main",build_vignettes=T)


library(TwoSampleMR)
library(coloc)
library(locuscomparer)
library(dplyr)
library(xlsx)
library(data.table)
library(ldscr)
library(ggplot2)
setwd('...')
ao = read.csv('ao20230210.csv',header = T)

GTEx_exposure=fread('GTEx_exposure_all.txt',header = T) %>% filter(pheno=='LPL_blood')
exposure_dat=format_data(GTEx_exposure,
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

FI=VariantAnnotation::readVcf("ieu-a-798.vcf.gz")
FI<-gwasvcf::vcf_to_tibble(FI,id='ieu-a-798')


FI<-FI[match(exposure_dat$SNP,FI$rsid),]
FI$pval.outcome=10^(-FI$LP)
outcome_dat=format_data(FI,
                        type = "outcome",
                        snp_col = "rsid",
                        beta_col = "ES",
                        se_col = "SE",
                        effect_allele_col = "ALT",
                        other_allele_col = "REF",
                        eaf_col = "AF",
                        pval_col='pval.outcome',
                        chr_col = "seqnames",
                        pos_col = "start",
                        samplesize_col = "SS")

rm(GTEx_exposure,FI)
dat=harmonise_data(exposure_dat,outcome_dat,action = 2)
dat=dat[!duplicated(dat),]

result=coloc.abf(dataset1 = list(beta=dat$beta.exposure,
                                 varbeta=(dat$se.exposure)^2,
                                 pvalues=dat$pval.exposure,
                                 MAF=dat$eaf.exposure,
                                 snp=dat$SNP,
                                 N=670,
                                 type="quant",
                                 position=dat$pos.outcome),
                 dataset2 = list(beta=dat$beta.outcome,
                                 varbeta=(dat$se.outcome)^2,
                                 pvalues=dat$pval.outcome,
                                 MAF=dat$eaf.outcome,
                                 snp=dat$SNP,
                                 N=171875,
                                 type="quant",
                                 position=dat$pos.outcome))

snp_result=result$results %>% filter(SNP.PP.H4>0.75)
o <- order(result$results$SNP.PP.H4,decreasing=TRUE)
cs <- cumsum(result$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]  
result$results[o,][1:w,]$snp

PPH=result[["summary"]]

### Draw colocalisation plot

eqtl=dat[,c('SNP','pval.exposure')] %>% rename(rsid=SNP,pval=pval.exposure)
gwas=dat[,c('SNP','pval.outcome')] %>% rename(rsid=SNP,pval=pval.outcome)
write.table(eqtl,'.../eqtl.tsv',row.names = F,col.names=T)
write.table(gwas,'.../gwas.tsv',row.names = F,col.names=T)
in_fn1='.../eqtl.tsv'
in_fn2='.../gwas.tsv'
jpeg('coloca(LPL(blood)_MI).jpeg',height = 1600,width = 3000,res= 300)
locuscompare(in_fn1 = in_fn1, title1 = 'eQTL(LPL-blood)',
             in_fn2 = in_fn2, title2 = 'GWAS(Myocardial infarction)')
dev.off()

# Colocalazation,Mediation Analysis,ldsc.
exposure=c("APOB_adipo","LPL_adipo","LPL_blood")
expo_label=c('APOB-adiposity','LPL-adiposity','LPL-blood')
expo_N=c(663,663,755)
outcome=c('ieu-a-798.vcf.gz','ukb-d-I9_IHD.vcf.gz','ukb-d-I9_CHD.vcf.gz')
outco_label=c('Myocardial infarction','Ischemic heart diseases','Major CHD event')
outco_N=c(171875,361194,361194)
PPH=list()
for(i in 1:3) {
  for(j in 1:3) {
    GTEx_exposure=fread('GTEx_exposure_all.txt',header = T) %>% filter(pheno==exposure[i])
    exposure_dat=format_data(GTEx_exposure,
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
    
    FI=VariantAnnotation::readVcf(outcome[j])
    FI<-gwasvcf::vcf_to_tibble(FI)
    
    
    FI<-FI[match(exposure_dat$SNP,FI$rsid),]
    FI$pval.outcome=10^(-FI$LP)
    outcome_dat=format_data(FI,
                            type = "outcome",
                            snp_col = "rsid",
                            beta_col = "ES",
                            se_col = "SE",
                            effect_allele_col = "ALT",
                            other_allele_col = "REF",
                            eaf_col = "AF",
                            pval_col='pval.outcome',
                            chr_col = "seqnames",
                            pos_col = "start",
                            samplesize_col = "SS")
    
    rm(GTEx_exposure,FI)
    dat=harmonise_data(exposure_dat,outcome_dat,action = 2)
    dat=dat[!duplicated(dat),]
    
    result=coloc.abf(dataset1 = list(beta=dat$beta.exposure,
                                     varbeta=(dat$se.exposure)^2,
                                     pvalues=dat$pval.exposure,
                                     MAF=dat$eaf.exposure,
                                     snp=dat$SNP,
                                     N=expo_N[i],
                                     type="quant",
                                     position=dat$pos.outcome),
                     dataset2 = list(beta=dat$beta.outcome,
                                     varbeta=(dat$se.outcome)^2,
                                     pvalues=dat$pval.outcome,
                                     MAF=dat$eaf.outcome,
                                     snp=dat$SNP,
                                     N=outco_N[j],
                                     type="quant",
                                     position=dat$pos.outcome))
    PPH[[j+3*(i-1)]]=result[["summary"]]
    
    eqtl=dat[,c('SNP','pval.exposure')] %>% rename(rsid=SNP,pval=pval.exposure)
    gwas=dat[,c('SNP','pval.outcome')] %>% rename(rsid=SNP,pval=pval.outcome)
    write.table(eqtl,'.../eqtl.tsv',row.names = F,col.names=T)
    write.table(gwas,'.../gwas.tsv',row.names = F,col.names=T)
    in_fn1='.../eqtl.tsv'
    in_fn2='.../gwas.tsv'
    
    jpeg(paste0('coloca(',expo_label[i],'_',outco_label[j],').jpeg'),
         height = 1600,width = 3000,res= 300)
    print(locuscompare(in_fn1 = in_fn1, title1 = paste0('eQTL(',expo_label[i],')'),
                 in_fn2 = in_fn2, title2 = paste0('GWAS(',outco_label[j],')')))
    dev.off()
  }
}
PPH=do.call('rbind',PPH)


### Two step MR for mediation 

drug=read.xlsx('drug.xlsx',sheetName = 'drug1used') %>% filter(exposure==c('APOC3','LPL'))
mediator=c('ebi-a-GCST90002232','ieu-b-118','ieu-b-38','ieu-b-39')
outcome=c('ieu-a-798','ukb-d-I9_CHD','ukb-d-I9_IHD')

# X --> M
outcomes <- filter(ao,id %in% mediator) %>%
  arrange(desc(year), desc(nsnp))
outcomes=outcomes[match(mediator,outcomes$id),]
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
outco=unique(dat$outcome)[match(mediator,sapply(strsplit(unique(dat$outcome), 'id:',fixed = T), '[',2))]

library(MendelianRandomization)
library(LDlinkR)
res_IVW=matrix(NA,length(expos)*length(outco),8)
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
                                        IVW_res@Pvalue)
    }
    
  }
}
colnames(res_IVW)=c('exposure','outcome','SNPs',
                    'ES','SE','CILower','CIUpper','P value')
XtoM=res_IVW


# X --> Y
outcomes <- filter(ao,id %in% outcome) %>%
  arrange(desc(year), desc(nsnp))
outcomes=outcomes[match(outcome,outcomes$id),]
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
outco=unique(dat$outcome)[match(outcome,sapply(strsplit(unique(dat$outcome), 'id:',fixed = T), '[',2))]

library(MendelianRandomization)
library(LDlinkR)
res_IVW=matrix(NA,length(expos)*length(outco),8)
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
                                        IVW_res@Pvalue)
    }
    
  }
}
colnames(res_IVW)=c('exposure','outcome','SNPs',
                    'ES','SE','CILower','CIUpper','P value')
XtoY=res_IVW


# M --> Y
mediators=extract_instruments(mediator)
for(i in 1:nrow(mediators)){
  temp=mediators[i,]
  beta=temp$beta.exposure
  MAF=temp$eaf.exposure
  N=temp$samplesize.exposure
  SE=temp$se.exposure
  R2=2*(1-MAF)*MAF*beta^2/(
    2*(1-MAF)*MAF*beta^2+2*(1-MAF)*MAF*N*SE^2)
  mediators$F[i]=R2/(1-R2)*(N-2)
}
mediators=na.omit(mediators)
mediators=mediators[mediators$F>10,]

outcomes <- filter(ao,id %in% outcome) %>%
  arrange(desc(year), desc(nsnp))
outcomes=outcomes[match(outcome,outcomes$id),]
outcome_dat=extract_outcome_data(mediators$SNP, outcomes$id)
index=which(outcomes[,'sd']>0)
temp <- outcomes[index,] %>% select(id, sd)
outcome_dat <- merge(outcome_dat, temp, by.x="id.outcome", by.y="id", all.x=TRUE)
index <- !is.na(outcome_dat$sd)
outcome_dat$beta.outcome[index] <- outcome_dat$beta.outcome[index] / outcome_dat$sd[index]
outcome_dat$se.outcome[index] <- outcome_dat$se.outcome[index] / outcome_dat$sd[index]

dat=harmonise_data(mediators,outcome_dat,action = 2)
dat=dat[!duplicated(dat),]

expos=unique(mediators$exposure)
outco=unique(dat$outcome)[match(outcome,sapply(strsplit(unique(dat$outcome), 'id:',fixed = T), '[',2))]

library(MendelianRandomization)
library(LDlinkR)
res_IVW=matrix(NA,length(expos)*length(outco),8)
for(i in 1:length(expos)){
  for(j in 1:length(outco)){
    temp=dat %>% filter(exposure==expos[i]&outcome==outco[j]
    ) %>% select(SNP,beta.exposure,se.exposure,
                 beta.outcome,se.outcome)
    if (nrow(temp)==0) {
      res_IVW[(i-1)*length(outco)+j,1:3]=c(expos[i],outco[j],length(temp$SNP))
    } else {
      MRInput=mr_input(bx=temp[,'beta.exposure'],
                       bxse=temp[,'se.exposure'],
                       by=temp[,'beta.outcome'],
                       byse=temp[,'se.outcome'])
      IVW_res=mr_ivw(MRInput,model= "default")
      res_IVW[(i-1)*length(outco)+j,]=c(expos[i],outco[j],length(temp$SNP),
                                        IVW_res@Estimate,IVW_res@StdError,
                                        IVW_res@CILower,IVW_res@CIUpper,
                                        IVW_res@Pvalue)
    }
    
  }
}
colnames(res_IVW)=c('exposure','outcome','SNPs',
                    'ES','SE','CILower','CIUpper','P value')
MtoY=res_IVW



expos=unique(XtoM[,'exposure'])
mediator=unique(MtoY[,'exposure'])
outcome=unique(XtoY[,'outcome'])
XtoY=data.frame(XtoY)
XtoM=data.frame(XtoM)
MtoY=data.frame(MtoY)

# Inferring the XtoY direct effect
library(TwoStepCisMR)
direct_effect=matrix(NA,length(expos)*length(mediator)*length(outcome),7)
colnames(direct_effect)=c('exposure','mediator','outcome',
                          'direct_XtoY','direct_lower','direct_upper','direct_P')

for(i in 1:length(expos)){
  for(j in 1:length(mediator)){
    for(k in 1:length(outcome)){
      go=subset(XtoY,exposure==expos[i]&outcome==outcome[k])
      gc=subset(XtoM,exposure==expos[i]&outcome==mediator[j])
      co=subset(MtoY,exposure==mediator[j]&outcome==outcome[k])
      temp=TSCMR(as.numeric(go$ES),as.numeric(go$SE),
                 as.numeric(gc$ES),as.numeric(gc$SE),
                 as.numeric(co$ES),as.numeric(co$SE))
      direct_effect[k+(j-1)*length(outcome)+(i-1)*length(mediator)*length(outcome),
                    ]=c(expos[i],mediator[j],outcome[k],
                        temp$Bgo,temp$Bgo-1.96*temp$BSSE,temp$Bgo+1.96*temp$BSSE,
                        2*pnorm(abs(temp$Bgo/temp$BSSE),lower.tail = F))
    }
  }
}


# Inferring indirect effect
indirect=matrix(NA,length(expos)*length(mediator)*length(outcome),18)
colnames(indirect)=c('exposure','mediator','outcome',
                     'XtoM_ES','XtoM_lower','XtoM_upper','XtoM_P',
                     'MtoY_ES','MtoY_lower','MtoY_upper','MtoY_P',
                     'indirect_ES','indirect_lower','indirect_upper','indirect_P',
                     'Mediated proportion','lower_%','upper_%')

for(i in 1:length(expos)){
  for(j in 1:length(mediator)){
    for(k in 1:length(outcome)){
      go=subset(XtoY,exposure==expos[i]&outcome==outcome[k])
      gc=subset(XtoM,exposure==expos[i]&outcome==mediator[j])
      co=subset(MtoY,exposure==mediator[j]&outcome==outcome[k])
      indirect_ES=as.numeric(gc$ES)*as.numeric(co$ES)
      SE_ind=sqrt(as.numeric(gc$ES)^2*as.numeric(co$SE)^2+
                    as.numeric(co$ES)^2*as.numeric(gc$SE)^2)
      pro=indirect_ES/as.numeric(go$ES)
      temp1=rnorm(1e06,indirect_ES,SE_ind)
      temp2=rnorm(1e06,as.numeric(go$ES),as.numeric(go$SE))
      pro_CI=quantile(temp1/temp2,c(0.025,0.975))
      rm(temp1,temp2)
      indirect[k+(j-1)*length(outcome)+(i-1)*length(mediator)*length(outcome),
      ]=c(expos[i],mediator[j],outcome[k],
          gc$ES,gc$CILower,gc$CIUpper,gc$P.value,
          co$ES,co$CILower,co$CIUpper,co$P.value,
          indirect_ES,indirect_ES-1.96*SE_ind,indirect_ES+1.96*SE_ind,
          2*pnorm(abs(indirect_ES/SE_ind),lower.tail = F),
          pro,pro_CI)
    }
  }
}

res_mediation=cbind(direct_effect,indirect[,-c(1:3)])
write.xlsx(res_mediation,'Mediation results.xlsx',row.names = F)

### LDSC analysis
library(ldscr)
df <- sumstats_munged_example(example = "BMI")
h2_res <- ldsc_h2(munged_sumstats = df, ancestry = "EUR")

rg_res <- ldsc_rg(
  munged_sumstats = list(
    "APOB" = sumstats_munged_example(example = "APOB"),
    "LDL" = sumstats_munged_example(example = "LDL")
  ),
  ancestry = "EUR"
)

library(ggplot2)
autoplot(rg_res)


######## Analysis @@@@@@@@@@@@@@@@@@@

LDL<-VariantAnnotation::readVcf("ieu-a-300.vcf.gz") %>% 
  gwasvcf::vcf_to_tibble(id='ieu-a-300') %>% 
  transform(Z=ES/SE) %>% 
  select('rsid','SS','Z','ALT','REF') %>% 
  rename(SNP='rsid',N='SS',A1='ALT',A2='REF')

TC<-VariantAnnotation::readVcf("ieu-a-301.vcf.gz") %>% 
  gwasvcf::vcf_to_tibble(id='ieu-a-301') %>% 
  transform(Z=ES/SE) %>% 
  select('rsid','SS','Z','ALT','REF') %>% 
  rename(SNP='rsid',N='SS',A1='ALT',A2='REF')

TG<-VariantAnnotation::readVcf("ieu-a-302.vcf.gz") %>% 
  gwasvcf::vcf_to_tibble(id='ieu-a-302') %>% 
  transform(Z=ES/SE) %>% 
  select('rsid','SS','Z','ALT','REF') %>% 
  rename(SNP='rsid',N='SS',A1='ALT',A2='REF')

MI<-VariantAnnotation::readVcf("ieu-a-798.vcf.gz") %>% 
  gwasvcf::vcf_to_tibble(id='ieu-a-798') %>% 
  transform(Z=ES/SE) %>% 
  select('rsid','SS','Z','ALT','REF') %>% 
  rename(SNP='rsid',N='SS',A1='ALT',A2='REF')

IHD<-VariantAnnotation::readVcf("ukb-d-I9_IHD.vcf.gz") %>% 
  gwasvcf::vcf_to_tibble(id='ukb-d-I9_IHD') %>% 
  transform(Z=ES/SE) %>% 
  select('rsid','SS','Z','ALT','REF') %>% 
  rename(SNP='rsid',N='SS',A1='ALT',A2='REF')

CHD<-VariantAnnotation::readVcf("ukb-d-I9_CHD.vcf.gz") %>% 
  gwasvcf::vcf_to_tibble(id='ukb-d-I9_CHD') %>% 
  transform(Z=ES/SE) %>% 
  select('rsid','SS','Z','ALT','REF') %>% 
  rename(SNP='rsid',N='SS',A1='ALT',A2='REF')

rg_res <- ldsc_rg(
  munged_sumstats = list(
    'LDL'=LDL,'TC'=TC,'TG'=TG,
    'MI'=MI,'IHD'=IHD,'CHD'=CHD
  ),
  ancestry = "EUR"
)

jpeg('Figure 6A.jpeg',height = 1400,width = 1800,res= 300)
print(autoplot(rg_res))
dev.off()


####### GO and PPI @@@@@@@@@@@@@@@@@@@
library(tidyverse)
KEGG=read.xlsx('KEGG and  GO.xlsx',sheetName = 'KEGG')
KEGG=separate(KEGG,Term,sep = ':',into = c('ID','Description'))

jpeg('Figure 6B.jpeg',height = 1400,width = 1800,res= 300)
print(ggplot(KEGG,aes(y=Description,x=Gene_ratio))+
        geom_point(aes(size=Count,color=PValue))+
        scale_color_gradient(low='red',high='blue')+
        labs(color=expression(PValue,size='Count'),
             x='Gene ratio',y='Pathways',title='KEGG Pathway Enrichment')+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5)))
dev.off()


GO=read.xlsx('KEGG and  GO.xlsx',sheetName = 'GO')
GO=separate(GO,Term,sep = '~',into = c('ID','Description'))

for(i in 1:nrow(GO)){
  description_splite=strsplit(GO$Description[i],split = ' ')
  description_collapse=paste(description_splite[[1]][1:5],collapse=' ')#选择前五个单词GO term名字
  GO$Description[i]=description_collapse
  GO$Description=gsub(pattern = 'NA','',GO$Description) 
}
GO_term_order=factor(as.integer(rownames(GO)),labels = GO$Description)
m=GO$Category
m=gsub("TERM","",m)
m=gsub("_DIRECT","",m)
GO$Category=m

jpeg('Figure 6C.jpeg',height = 1800,width = 3500,res= 300)
print(ggplot(GO,aes(x=GO_term_order,y=Count,fill=Category))+
        geom_bar(stat='identity',width=0.8)+
        scale_fill_manual(values = c('#6666FF','#33CC33','#FF6666'))+
        xlab('GO term')+
        ylab('Gene Number')+
        labs(title='GO Terms Enrich')+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(family = 'sans',face='bold',
                                         color = 'gray50',angle = 70,
                                         vjust=1,hjust=1)))
dev.off()






