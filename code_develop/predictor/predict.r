library(FSelector)

## generate prediction model (adapting from old pathipred) 2015-10
get.prediction.model <- function(prettyways.matrix, exp.design, type, k, filter_paths=F) {
  ## pre-process data
  cat("Removing paths without variability...\n")
  m <- t(remove.novar.paths(prettyways.matrix))
  ed <- sapply(rownames(m), function(x) add.exp.info(x, exp.design, 1, 2, type))
  m <- data.frame(m, ed)
  
  ## filter features with Correlation-based Feature Selection (CFS)
  if(filter_paths) {    
    cat("Correlation-based Feature Selection...\n")
    m.filtered <- cfs(ed ~ ., data=m)
    m <- m[, c(m.filtered, "ed")]
  }
  
  ## get prediction model
  cat("Generating best model...\n")
  if (type == "categorical") {
    svm.type = "C-classification"
  } else if (type == "continuous") {
    svm.type = "eps-regression"
  }  
  bestmod <- get.best.model(m, svm.type, k)
  return(list("model"=bestmod, "input"=m))
}

## remove paths without variability -> equal value across samples
remove.novar.paths <- function(m) {
  equalPaths <- apply(m, 1, function(x) { all(x == x[1]) })
  varPaths <- m[!equalPaths, ]
  return(varPaths)
}

## add experimental design info to activation values matrix
add.exp.info <- function(id, df, sample.indx, class.indx, type) {
  c <- df[, sample.indx] == id
  res <- df[c, ]
  if (type == "categorical") {
    return(as.character(res[class.indx]))
  } else if (type == "continuous") {
    return(as.numeric(res[class.indx]))
  }  
}

## generate best.model
get.best.model <- function(m, svm.type, k){
  tune.svm.model <- tune.svm(ed ~ ., data=m, scale=T, type=svm.type, kernel="radial", cost=10^(1:2), gamma=10^(-6:-3), cross=k)
  bestGamma <- tune.svm.model$best.parameters[[1]]
  bestC <- tune.svm.model$best.parameters[[2]]
  best.model <- svm(ed ~ ., data=m, scale=T, type=svm.type, kernel="radial", cost=bestC, gamma=bestGamma, cross=k)
  return(best.model)
}

## Compute confusion matrix
get.confusion.mat <- function(svm.pred, best.model) {
  acc <- table(pred = svm.pred, true = best.model$input[,"ed"])
  acc.export <- as.data.frame.matrix(acc)
  acc.export <- cbind(c("", rownames(acc.export)), c(colnames(acc.export)[1],acc.export[,1]), c(colnames(acc.export)[2],acc.export[,2]))
  diag <- classAgreement(acc)$diag
  kappa <- classAgreement(acc)$kappa
  rand <- classAgreement(acc)$rand
  crand <- classAgreement(acc)$crand
  return(list("acc"=acc, "diag"=diag, "kappa"=kappa, "rand"=rand, "crand"=crand))  
}

get.stats <- function(svm.pred, best.model, edtype){
  conf.mat <- get.confusion.mat(svm.pred, best.model) 
  
  if (edtype == "categorical") {
    sensitivity <- conf.mat[["acc"]][1,1]/sum(conf.mat[["acc"]][,1]) ## True Positives / Positive
    specificity <- conf.mat[["acc"]][2,2]/sum(conf.mat[["acc"]][,2]) ## True Negatives / Negatives
    PPV <- conf.mat[["acc"]][1,1]/sum(conf.mat[["acc"]][1,]) ## True Positives / Test Outcome Positive
    NPV <- conf.mat[["acc"]][2,2]/sum(conf.mat[["acc"]][2,]) ## True Negatives / Test Outcome Negatives
    FPR <- 1 - specificity
    FNR <- 1 - sensitivity
    LRP <- sensitivity/(1-specificity) 
    LRN <- (1 - sensitivity)/specificity
    
    all.stats <- as.data.frame(list(statistic=c("Sensitivity", "Specificity", "Positive Predictive Value", 
                                                  "Negative Predictive Value", "False Positive Rate", "False Negative Rate", "Likelihood Ratio Positive", 
                                                  "Likelihood Ratio Negative", "Percentage of data points in the main diagonal", 
                                                  "Percentage of data points in the main diagonal corrected for agreement by chance", "Rand index", 
                                                  "Rand index corrected for agreement by chance", "Total Accuracy"),
                                      value=c(sensitivity, specificity, PPV, NPV, FPR, FNR, LRP, LRN, conf.mat[["diag"]], conf.mat[["kappa"]], conf.mat[["rand"]], conf.mat[["crand"]], 
                                              best.model$model$tot.accuracy)))      
  } else if (edtype == "continuous") {
      all.stats <- as.data.frame(list(statistic=c("Percentage of data points in the main diagonal", 
                                                  "Percentage of data points in the main diagonal corrected for agreement by chance", 
                                                  "Rand index", "Rand index corrected for agreement by chance", "Total Mean Squared Error"),
                                      value=c(conf.mat[["diag"]], conf.mat[["kappa"]], conf.mat[["rand"]], conf.mat[["crand"]], best.model$model$tot.MSE)))  
  }
  return(all.stats)
}

## best.model statistics
get.predmod.stats <- function(best.model, edtype) {
  svm.pred <- predict(best.model[["model"]], best.model[["input"]], decision.values=F, probability=F, na.action=na.omit)
  conf.mat <- get.confusion.mat(svm.pred, best.model) 
  stats <- get.stats(svm.pred, best.model, edtype)
  return(stats)
}

predict.newdataset <- function(model, m, exp.design, type) {
  ed <- sapply(rownames(m), function(x) add.exp.info(x, exp.design, 1, 2, type))
  m <- data.frame(m, ed)
  svm.pred <- predict(model, m, decision.values=FALSE, probability=FALSE, na.action=na.omit)
  pred.res <- as.data.frame(cbind(names(svm.pred), as.character(svm.pred)))
  names(pred.res) <- c("sample", "prediction")
  return(list("svm_pred"=svm.pred, "pred_res"=pred.res)) 
}

save.pred.res <- function(model,model_stats,pretty_results,path_vals,pathigraphs,output.folder,filter.paths, effector, conf=0.05){
    
  if(!file.exists(output.folder)){
    dir.create(output.folder)
  }
  # pathipred
  write.table(path_vals,file=paste0(output.folder,"/paths_vals.txt"),col.names=T,row.names=T,quote=F,sep="\t")
  save(model, file=paste0(output.folder,"/model.RData"))
  if(!is.null(model_stats)) write.table(model_stats,file=paste0(output.folder,"/model_stats.txt"),col.names=T,row.names=F,quote=F,sep="\t")
  # CFS + cellmaps
  if (filter.paths) {
    features <- colnames(model[["input"]])[-length(colnames(model[["input"]]))]
    write.table(features,file=paste0(output.folder,"/filtered_features_cfs.txt"),col.names=F,row.names=F,quote=F,sep="\t")  
    # CELLMAPS
    if(!exists(paste0(output.folder, "/sifs4CellMaps/")))
      dir.create(paste0(output.folder, "/sifs4CellMaps/"))
    
    colors <- c(rep("UP", times=length(features)), rep("DOWN", times=length(rownames(path_vals))-length(features))) 
    pvals <- as.numeric(c(rep(0.0001, times=length(features)), rep(0.6, times=length(rownames(path_vals))-length(features))))
    mywt <- data.frame(as.numeric(rep(0.3, times=length(rownames(path_vals)))), as.character(colors), as.numeric(rep(0.3, times=length(rownames(path_vals)))), pvals,stringsAsFactors=F)
    colnames(mywt) <- c("p.value", "UP/DOWN", "statistic", "FDRp.value")
    touse <- rownames(path_vals) %in% features
    rownames(mywt) <- c(features, rownames(path_vals)[!touse])
    
    comp_names <- sapply(strsplit(rownames(mywt),"__"),"[[",1)
    
    for(pathway in names(results$by.path)){
      this_comp <- mywt[which(comp_names==pathway),]
      print(pathway)
      write.attributes(this_comp, pathway, pathigraphs, paste0(output.folder, "/sifs4CellMaps/", pathway), moreatts=NULL, effector=effector,conf=conf)
    }
    
    dir.create(paste0(output.folder,"/report/"))
    create.html.report2(fpathigraphs,mywt,pretty_home,output.folder,effector=(decompose==F),template_name = "report_template.html",output_name = "network.html")
  }      
}

save.prednewdataset.res <- function(prediction,prediction_stats,pretty_results,path_vals,pathigraphs,output.folder){
  
  if(!file.exists(output.folder)){
    dir.create(output.folder)
  }
  write.table(path.vals,file=paste0(output.folder,"/path_vals.txt"),col.names=T,row.names=T,quote=F,sep="\t")
  write.table(prediction,file=paste0(output.folder,"/prediction_results.txt"),col.names=T,row.names=F,quote=F,sep="\t")
  write.table(prediction_stats,file=paste0(output.folder,"/prediction_stats.txt"),col.names=T,row.names=F,quote=F,sep="\t")
  
}