
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
library(readxl)
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
  pvcaobj<-pvcaBF(t(dat),sample_info,c("batch","type"),0.1)
  pieDraw(pvcaobj)
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


#proteome plasma
# install.packages("rJava")
# library(rJava)
# prot<-gdata::read.xls("1508plasma/322-Diogenes-Precursor_Report.xls",1)
# library(openxlsx)
# prot<-read.xls(xls="1508plasma/322-Diogenes-Precursor_Report.xls",sheet=1)
prot<-read.table("1508plasma/322-Diogenes-Precursor_Report.xls",header = T,sep = "\t",stringsAsFactors=FALSE,quote = "",fill = TRUE)
prot.samples<-unique(prot$R.FileName)
prot.prots<-unique(prot$PG.ProteinAccessions)
prot3<-data.frame(matrix(NA,ncol = length(prot.prots),nrow=length(prot.samples)))
colnames(prot3)<-prot.prots
rownames(prot3)<-prot.samples
prot2<-prot[,c('R.FileName','PG.ProteinAccessions','PG.Quantity')]

for (i in 1:nrow(prot2)) {
  prot3[prot2[i,"R.FileName"],prot2[i,"PG.ProteinAccessions"]]<-prot2[i,'PG.Quantity']
}

# prot.w<-dcast(prot2,R.FileName~prot2$PG.ProteinNames,value.var='PEP.Quantity')
# prot.w2<-spread(prot2,key=PG.ProteinNames,value=PEP.Quantity)
prot3.t<-prot3

prot3<-t(prot3)
#sample annotation did not have "D_D170317_S322-1014402757-P1-E9_MHRM_R01_T0"
prot3<-prot3[,-c(which(colnames(prot3)=="D_D170317_S322-1014402757-P1-E9_MHRM_R01_T0"))]

protSample<-read.table("1508plasma/DiOGenes Raw File Annotation.txt",header=T,sep="\t")

protSample$RawFile<-gsub(".raw","",protSample$RawFile,fixed = T)
rownames(protSample)<-protSample$RawFile
protSample<-protSample[colnames(prot3),]
prot_batch<-as.character(protSample$Plate)
prot_type<-as.character(protSample$Condition)
protSample2<-data.frame(type=prot_type,batch=prot_batch)
rownames(protSample2)<-rownames(protSample)
prot4<-prot3
prot4[is.na(prot4)]<-1#.1*min(prot3,na.rm = T)
prot4<-log2(prot4)
##
prot4<-prot4[,-which(nchar(protSample2$type)<1)]
protSample3<-protSample2[-which(nchar(protSample2$type)<1),]
prot_batch<-prot_batch[-which(nchar(protSample2$type)<1)]
prot_type<-prot_type[-which(nchar(protSample2$type)<1)]
###keep type 1,2,3
type123Index<-which(protSample3$type==1|protSample3$type==2|protSample3$type==3|protSample3$type==4)
protSample3<-protSample3[type123Index,]
prot_batch<-prot_batch[type123Index]
prot_type<-prot_type[type123Index]
prot4<-prot4[,rownames(protSample3)]
prot4<-prot4[which(!grepl(";",rownames(prot4))),]
analysisBatch(prot4,prot_batch,prot_type,protSample3,"plasmaProteome20200916.pdf")
