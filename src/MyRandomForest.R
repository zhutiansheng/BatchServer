myRF<-function(d,ntree,nodesize){
  label_index<-grep("^label$",colnames(d))
  header<-colnames(d[-label_index])
  header_name<-paste0("index_",1:length(header))
  names(header)<-header_name
  colnames(d)[-label_index]<-header_name
  RandomForest <- randomForest(label ~ . ,data=d,importance=T,ntree=ntree,nodesize=nodesize,na.action=na.exclude)
  imps  <- data.frame(importance(RandomForest));
  impScore <- imps$MeanDecreaseAccuracy
  imps <- imps[order(impScore,decreasing=T),]
  orderedFeatures <- rownames(imps)
  orderedFeatures<-header[orderedFeatures]
  return(list(features=orderedFeatures,mod=RandomForest))
}
