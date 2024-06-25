library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(xlsx)
setwd('...')

######### plot forest plot ########
library(ggforestplot)
library(tidyverse)

Fig2A <- read.xlsx('Fig 2.xlsx',sheetName = '2A',header = T)

jpeg('Figure 2A.jpeg',height = 1500,width = 1400,res= 300)
print(forestplot(
  df = Fig2A,
  estimate = beta,
  logodds = TRUE,
  colour = Lipids,
  pvalue = pvalue,
  psignif = 0.005,
  xlab = "OR (95% CI) in discovery dataset"
)+theme(legend.position = "top"))
dev.off()


Fig2B <- read.xlsx('Fig 2.xlsx',sheetName = '2B',header = T)

jpeg('Figure 2B.jpeg',height = 1500,width = 1400,res= 300)
print(forestplot(
  df = Fig2B,
  estimate = beta,
  logodds = TRUE,
  colour = Lipids,
  pvalue = pvalue,
  psignif = 0.05,
  xlab = "OR (95% CI) in replication dataset"
)+theme(legend.position = "top"))
dev.off()

Fig2C <- read.xlsx('Fig 2.xlsx',sheetName = '2C',header = T)

jpeg('Figure 2C.jpeg',height = 1400,width = 1400,res= 300)
print(forestplot(
  df = Fig2C,
  estimate = beta,
  logodds = TRUE,
  colour = Targets,
  pvalue = pvalue,
  psignif = 0.005,
  xlab = "OR (95% CI) in discovery dataset"
)+theme(legend.position = "none"))
dev.off()


Fig2D <- read.xlsx('Fig 2.xlsx',sheetName = '2D',header = T)

jpeg('Figure 2D.jpeg',height = 1400,width = 1400,res= 300)
print(forestplot(
  df = Fig2D,
  estimate = beta,
  logodds = TRUE,
  colour = Targets,
  pvalue = pvalue,
  psignif = 0.05,
  xlab = "OR (95% CI) in replication dataset"
)+theme(legend.position = "none"))
dev.off()


Fig3A <- read.xlsx('Fig 3.xlsx',sheetName = '3A',header = T)

jpeg('Figure 3A.jpeg',height = 1680,width = 1400,res= 300)
print(forestplot(
  df = Fig3A,
  estimate = beta,
  logodds = F,
  colour = Lipids,
  pvalue = pvalue,
  psignif = 0.0036,
  xlab = "Beta (95% CI) for mediators"
)+theme(legend.position = "top"))
dev.off()


Fig3B <- read.xlsx('Fig 3.xlsx',sheetName = '3B',header = T)

jpeg('Figure 3B.jpeg',height = 1500,width = 1400,res= 300)
print(forestplot(
  df = Fig3B,
  estimate = beta,
  logodds = F,
  colour = Targets,
  pvalue = pvalue,
  psignif = 0.0036,
  xlab = "Beta (95% CI) for mediators"
)+theme(legend.position = "none"))
dev.off()

