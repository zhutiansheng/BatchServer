myRF<-function(d,ntree,nodesize){
  colnames(d)<-paste0("a",colnames(d))
  RandomForest <- randomForest(alabel ~ . ,data=d,importance=T,ntree=ntree,nodesize=nodesize,na.action=na.exclude)
  
  imps  <- data.frame(importance(RandomForest));
  impScore <- imps$MeanDecreaseAccuracy
  #return(impScore)
  imps <- imps[order(impScore,decreasing=T),]
  orderedFeatures <- rownames(imps)
  orderedFeatures<-gsub("^a","",orderedFeatures)
  return(list(features=orderedFeatures,mod=RandomForest))
}
