
options(stringsAsFactors = F)
library(fitdistrplus)
library(extraDistr)
library(umap)
library(ggplot2)
library(Biobase)
library(pvca)
library(BiocParallel)
library(tidyr)
library(reshape2)
library(openxlsx)
source("../src/MyPVCA.R")
source("../src/MyCombat.R")
source("../src/MyPriorDraw.R")
source("../test/functionsFromsva2.R")
drawUMAP<-function(myd,label){
  #myd[is.na(myd)]<-0
  myumap<-umap(myd,n_neighbors=10)
  mydf<-data.frame(myumap$layout)
  mydf$label<-label
  p<-ggplot(mydf,aes(x=X1, y=X2, colour=label)) + geom_point(size=3)+
    theme(  #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      #axis.text.x = element_text(vjust = 1,angle = 45),
      axis.line.x = element_line(color="black", size = 0.5),
      axis.line.y = element_line(color="black", size = 0.5),
      panel.background = element_blank())
  print(p)
}

#library(BatchQC)

analysisBatch<-function(dat,batch,type,sample_info,out_pdf){
  pdf(out_pdf)
  #pvcaobj<-pvcaBF(t(dat),sample_info,c("batch","type"),0.1)
  #pieDraw(pvcaobj)
  drawUMAP(t(dat),as.factor(type))
  drawUMAP(t(dat),as.factor(batch))
  
  modcombat = model.matrix(~as.factor(type), data=sample_info)
  
  for (param in c("auto","parameter","noparameter")) {
    t<-system.time(dat.combat<-combat(as.matrix(dat),batch,mod = modcombat,par.prior=param))
    message(param)
    print(t)
    dat2<-dat.combat$bayesdata
    #dat2[is.na(dat2)]<-0
    combat_pvcaobj<-pvcaBF(t(dat2),sample_info,c("batch","type"),0.1)
    pieDraw(combat_pvcaobj)
    drawUMAP(t(dat2),as.factor(type))
    drawUMAP(t(dat2),as.factor(batch))
    if(param == "auto"){
      print(dat.combat$additiondata$passTest)
      for(b in unique(batch))
        drawPrior(dat.combat$additiondata,b)
    }
    
  }
  dev.off()
}
#########################################################


#proteome p

prot<-read.xlsx("1-s2.0-S0092867420306279-mmc2.xlsx",sheet='Proteomics_proteins_training',rowNames =T )

samples<-read.xlsx("1-s2.0-S0092867420306279-mmc1.xlsx",sheet='Clinical_information')
#samples<-samples[,1:5]
#samples<-samples[nchar(samples$MS.ID.b)>3,]
#rownames(samples)<-samples$MS.ID.b
#protsSample<-samples[,c("Group.d","MS.ID.b")]
protsSample<-samples
rownames(protsSample)<-samples$sample
protsSample<-protsSample[,-1]
colnames(protsSample)<-c("batch","type")
protsSample<-protsSample[colnames(prot),]
type<-protsSample$type
batch<-as.character(protsSample$batch)
rna<-apply(prot,1,function(d){
  sum(is.na(d))
})
#prot<-prot[rna<(ncol(prot)*0.5),]
prot[is.na(prot)]<-0#.1*min(prot,na.rm = T)
prot[prot>5]<-5
analysisBatch(prot,batch,type,protsSample,"cellTMT20200915b.pdf")
