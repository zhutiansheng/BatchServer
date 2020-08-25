
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
source("src/MyPVCA.R")
source("src/MyCombat.R")
source("src/MyPriorDraw.R")
source("test/functionsFromsva2.R")
drawUMAP<-function(myd,label){
  myd[is.na(myd)]<-0
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
    dat2[is.na(dat2)]<-0
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
#prot<-read.table("F:/复旦大学/毕业相关/BatchServer/testData/1508plasma/DiOGenes Protein-Level-DIA_Report.txt",header = T,sep = "\t")

prot<-read.table("F:/复旦大学/毕业相关/BatchServer/testData/1508plasma/DiOGenes Glycation DIA_Report.txt",header = T,sep = "\t")
prot.samples<-unique(prot$R.FileName)
prot.prots<-unique(prot$PG.ProteinNames)
prot3<-data.frame(matrix(NA,ncol = length(prot.prots),nrow=length(prot.samples)))
colnames(prot3)<-prot.prots
rownames(prot3)<-prot.samples
prot2<-prot[,c('R.FileName','PG.ProteinNames','PEP.Quantity')]

for (i in 1:nrow(prot2)) {
  prot3[prot2[i,"R.FileName"],prot2[i,"PG.ProteinNames"]]<-prot2[i,'PEP.Quantity']
}

# prot.w<-dcast(prot2,R.FileName~prot2$PG.ProteinNames,value.var='PEP.Quantity')
# prot.w2<-spread(prot2,key=PG.ProteinNames,value=PEP.Quantity)

prot3<-t(prot3)
#sample annotation did not have "D_D170317_S322-1014402757-P1-E9_MHRM_R01_T0"
prot3<-prot3[,-c(which(colnames(prot3)=="D_D170317_S322-1014402757-P1-E9_MHRM_R01_T0"))]

protSample<-read.table("F:/复旦大学/毕业相关/BatchServer/testData/1508plasma/DiOGenes Raw File Annotation.txt",header=T,sep="\t")

protSample$RawFile<-gsub(".raw","",protSample$RawFile,fixed = T)
rownames(protSample)<-protSample$RawFile
protSample<-protSample[colnames(prot3),]
prot_batch<-as.character(protSample$Plate)
prot_type<-as.character(protSample$Condition)
protSample2<-data.frame(type=prot_type,batch=prot_batch)
rownames(protSample2)<-rownames(protSample)
prot3[is.na(prot3)]<-1
prot3<-log2(prot3)
analysisBatch(prot3,prot_batch,prot_type,protSample2,"plasmaProteome.pdf")
