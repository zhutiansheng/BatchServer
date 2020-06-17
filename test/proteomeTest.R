<<<<<<< HEAD
=======
options(stringsAsFactors = F)
>>>>>>> 0f50c8dcffa8a8dfad94cbc4300ea6f134dba506
library(fitdistrplus)
library(extraDistr)
library(umap)
library(ggplot2)
library(Biobase)
library(pvca)
library(BiocParallel)
library(openxlsx)
source("src/MyPVCA.R")
<<<<<<< HEAD
source("test/MyCombat2.R")
=======
#source("test/MyCombat2.R")
source("src/MyCombat.R")
>>>>>>> 0f50c8dcffa8a8dfad94cbc4300ea6f134dba506
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
#dataset 1
data(example_batchqc_data)
batch <- batch_indicator$V1
condition <- batch_indicator$V2
sample_info<-data.frame(batch,type=condition)
rownames(sample_info)<-colnames(signature_data)
analysisBatch(signature_data,batch,condition,sample_info,"signature.pdf")
####
<<<<<<< HEAD
#proteome

prot<-read.xlsx("G:\\proteomeTest\\1-s2.0-S2589004219304407-mmc2.xlsx",sheet='E')
prot2<-read.xlsx("G:\\proteomeTest\\1-s2.0-S2589004219304407-mmc2.xlsx",sheet='D')
rownames(prot2)<-prot2$protein.accession.number
prot2<-prot2[,-c(1:8)]
protSample<-read.xlsx("G:\\proteomeTest\\1-s2.0-S2589004219304407-mmc2.xlsx",sheet='A')
rownames(protSample)<-protSample$sample_name
prot_batch<-protSample$PCT_batch_id
prot_type<-protSample$tumor_type
colnames(protSample)[7]<-"type"
colnames(protSample)[8]<-"batch"
analysisBatch(prot2,prot_batch,prot_type,protSample,"G:\\proteomeTest\\nci60Proteome.pdf")
=======
#proteome nci60
#prot<-read.xlsx("test/1-s2.0-S2589004219304407-mmc2.xlsx",sheet='E')
prot2<-read.xlsx("test/1-s2.0-S2589004219304407-mmc2.xlsx",sheet='D')
rownames(prot2)<-prot2$protein.accession.number
prot2<-prot2[,-c(1:8)]
protSample<-read.xlsx("test/1-s2.0-S2589004219304407-mmc2.xlsx",sheet='A')
rownames(protSample)<-protSample$sample_name
protSample<-protSample[colnames(prot2),]
prot_batch<-as.character(protSample$PCT_batch_id)
prot_type<-as.character(protSample$tumor_type)
#colnames(protSample)[7]<-"type"
#colnames(protSample)[8]<-"batch"
protSample2<-data.frame(type=prot_type,batch=prot_batch)
rownames(protSample2)<-rownames(protSample)
analysisBatch(prot2,prot_batch,prot_type,protSample2,"nci60Proteome2.pdf")
#proteome tpd
prot2<-read.xlsx("test/sTable3 Protein Matrix_20200106_v3.xlsx",sheet=1)
rownames(prot2)<-prot2[,1]
prot2<-t(prot2)
prot2<-prot2[-1,]
prot2[is.na(prot2)]<-0
prot2<-apply(prot2,2,as.numeric)
prot2<-data.frame(prot2)
rna<-apply(prot2, 1,function(x){sum(x==0)})
prot2<-prot2[rna<ncol(prot2)*0.8,]
protSample<-read.xlsx("test/sTable3 Protein Matrix_20200106_v3.xlsx",sheet=4)
rownames(protSample)<-protSample$MS_file_name
protSample<-protSample[colnames(prot2),]
prot_batch<-as.character(substring(protSample$MS_file_name,1,1))
prot_type<-as.character(protSample$Histopathology_type)
protSample2<-data.frame(type=prot_type,batch=prot_batch)
rownames(protSample2)<-rownames(protSample)
analysisBatch(prot2,prot_batch,prot_type,protSample2,"tpdProteome.pdf")
>>>>>>> 0f50c8dcffa8a8dfad94cbc4300ea6f134dba506
#####
#dataset 2
#####
library(bladderbatch)
data(bladderdata)
sample_info <- pData(bladderEset)
edata <- exprs(bladderEset)
batch <- sample_info$batch  
type <- as.vector(sample_info$cancer)
colnames(sample_info)[4]<-"type"
analysisBatch(edata,batch,type,sample_info,"bladder.pdf")

#data set 3
data(protein_example_data)
dat<-protein_data
sample_info<-protein_sample_info[,c("samplename", "Batch","category")]
sample_info$samplename<-paste0("X",sample_info$samplename)
rownames(sample_info)<-sample_info[,1]
colnames(sample_info)<-c("samplename", "batch","type")
sample_info<-sample_info[colnames(dat),]
type<-as.vector(sample_info$type)
batch<-sample_info$batch
analysisBatch(dat,batch,type,sample_info,"protein.pdf")
