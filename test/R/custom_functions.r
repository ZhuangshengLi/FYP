
# Libraries ----
library(glmnet)
library(randomForest)
library(caret)
library(varSelRF)
library(VennDiagram)
library(limma)
library(psych)
library(reshape)
library(RColorBrewer)
library(readr)
library(ggplot2)
library(ggpubr)
library(Boruta)
library(ggrepel)
library(gridExtra)


# PCA_plot <- function(data, group = NULL, label = "", title = "", top_genes = 10, n_components = 10) {
#   
#   # 标准化数据并执行PCA
#   scaled_data <- scale(data, center = TRUE, scale = TRUE)
#   Data_pca <- prcomp(scaled_data)
#   
#   # 获取贡献度
#   contributions <- summary(Data_pca)$importance[2, ] * 100
#   dataScores <- data.frame(Data_pca$x[, 1:2])
#   
#   # 创建PCA散点图
#   p <- ggplot(dataScores, aes(x = PC1, y = PC2)) +
#     {
#       if (!is.null(group)) aes(col = group) else aes()
#     } +
#     geom_point(size = 4) +
#     xlab(paste0("PC1 (", round(contributions[1], 2), "%)")) +
#     ylab(paste0("PC2 (", round(contributions[2], 2), "%)")) +
#     theme_bw(base_size = 16) +
#     ggtitle(title) +
#     theme(
#       plot.title = element_text(hjust = 0.5),
#       legend.position = "right"
#     )
#   
#   # 生成横向子表格数据
#   valid_pcs <- min(n_components, length(contributions))
#   gene_names <- colnames(data)
#   
#   # 筛选有效主成分（贡献度 > 0）
#   valid_indices <- which(contributions[1:valid_pcs] > 0)
#   
#   # 为每个有效主成分创建子表格
#   table_grobs <- lapply(valid_indices, function(i) {
#     pc_loadings <- Data_pca$rotation[, i]
#     sorted_indices <- order(abs(pc_loadings), decreasing = TRUE)
#     top_indices <- sorted_indices[1:top_genes]
#     
#     # 构建子表格
#     sub_table <- data.frame(
#       Gene = gene_names[top_indices],
#       Loading = sprintf("%.4f", round(pc_loadings[top_indices], 4)),
#       stringsAsFactors = FALSE
#     )
#     colnames(sub_table)[2] <- paste0("PC", i)
#     
#     # 创建带标题的表格组合
#     title_grob <- textGrob(
#       paste0("PC", i, " (", round(contributions[i], 2), "%)"),
#       gp = gpar(fontsize = 10, fontface = "bold")
#     )
#     
#     table_grob <- tableGrob(
#       sub_table,
#       rows = NULL,
#       theme = ttheme_minimal(
#         base_size = 8,
#         padding = unit(c(3, 4), "mm"),
#         core = list(
#           bg_params = list(fill = c("grey95", "white"))
#         ))
#       )
#       
#       # 垂直组合标题和表格
#       arrangeGrob(
#         title_grob,
#         table_grob,
#         heights = c(0.3, 1),
#         ncol = 1
#       )
#   })
#     
#     # 横向排列所有子表格（仅当存在有效表格时）
#     if (length(table_grobs) > 0) {
#       combined_tables <- do.call(arrangeGrob, c(table_grobs, nrow = 1, padding = unit(0.5, "line")))
#       final_plot <- grid.arrange(
#         p,
#         combined_tables,
#         nrow = 2,
#         heights = c(2, 0.8),
#         top = textGrob(label, gp = gpar(fontsize = 12, fontface = "bold"))
#       )
#     } else {
#       final_plot <- grid.arrange(
#         p,
#         nrow = 1,
#         top = textGrob(label, gp = gpar(fontsize = 12, fontface = "bold"))
#       )
#     }
#     
#     return(invisible(final_plot))
# }
PCA_plot <- function(data, group = NULL, title = "", n_components = 10, show_legend = TRUE) {
  
  # 1. 标准化数据并执行 PCA
  scaled_data <- scale(data, center = TRUE, scale = TRUE)
  Data_pca <- prcomp(scaled_data)
  
  # 2. 获取主成分的贡献度（方差占比）
  contributions <- summary(Data_pca)$importance[2, ] * 100
  
  # 3. 只绘制前 2 个主成分在散点图上
  dataScores <- data.frame(Data_pca$x[, 1:2])
  
  # 4. 根据是否有 group，动态生成映射
  #    如果 group = NULL，则不使用颜色映射
  mapping <- if (!is.null(group)) {
    aes(x = PC1, y = PC2, color = group)
  } else {
    aes(x = PC1, y = PC2)
  }
  
  # 5. 准备散点图
  p <- ggplot(data = dataScores, mapping) +
    geom_point(size = 4) +
    xlab(paste0("PC1 (", round(contributions[1], 2), "%)")) +
    ylab(paste0("PC2 (", round(contributions[2], 2), "%)")) +
    theme_bw(base_size = 16) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "left"  # 默认放在左侧
    )
  
  # 如果不需要显示图例，则隐藏
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  # 6. 构建要显示的文本（前 n_components 个主成分）
  label_text <- paste0(
    "PC", 1:n_components, ": ", 
    round(contributions[1:n_components], 2), "%", 
    collapse = "\n"
  )
  
  # 7. 在图右上角放置注释文本
  p <- p +
    coord_cartesian(clip = "off") +
    theme(
      # top, right, bottom, left 的顺序
      plot.margin = unit(c(1, 4, 1, 1), "lines")
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = label_text,
      hjust = -0.1,  # 负值表示往左移动文本块
      vjust = 1,     # 顶部对齐
      size = 3
    )
  
  return(p)
}


OLEdetect <- function(data, group) {
  
  # Perform PCA on the input data
  pca <- prcomp(data, scale. = TRUE, rank. = 10)
  U <- pca$x
  
  # Identify outliers using robust method
  ind.out <- apply(U, 2, function(x) which((abs(x - median(x)) / mad(x)) > 6 )) %>%
    Reduce(union, .)
  
  # Visualize outliers in PCA plot
  col <- rep("black", nrow(U))
  col[ind.out] <- "red"
  qplot(U[, 1], U[, 3], color = I(col), size = I(2)) + coord_equal()
  
  # Remove outliers from the data
  newdata <- data[-ind.out, ]
  newgroup <- group[-ind.out, ]
  
  # Visualize PCA plots before and after removing outliers
  pca_plot_before <- PCA_plot(data = data, group = as.factor(group$CLASS))
  pca_plot_after <- PCA_plot(data = newdata, group = as.factor(newgroup$CLASS))
  
  # Arrange PCA plots side by side
  plot_arrangement <- ggarrange(pca_plot_before, pca_plot_after,
                                labels = c("Before", "After"),
                                ncol = 2, nrow = 1)
  
  # Print the number of outliers and their corresponding sample names
  print(paste("Number of outliers:", length(ind.out)))
  print("Outlier Sample Names:")
  print(row.names(data[ind.out, ]))
  
  # Return a list containing the arranged plots, new data, new group, and outlier indices
  return(list(plot_arrangement, newdata, newgroup, ind.out))
}

modelMetrExtract<-function(res){
  out<-data.frame()
  vars<-data.frame()
  
  for(i in 1:length(res)){
    
    if(is.null(res[[i]])){
      next
    }else{
      
      out<-rbind(out,res[[i]]$rfmodelRes$PredStats_byClass)
      vars<-rbind(vars,paste(res[[i]]$varSel,collapse = " "))
    }
  }
  colnames(out)<-names(res[[1]]$rfmodelRes$PredStats_byClass)
  finalRes<-data.frame("vars"=vars,out)
  colnames(finalRes)[1]<-"vars"
  
  return(finalRes)
}


#function to do boxplots of selected variables


boxplotSubset<-function(soma=soma, selVars,metadata=metina,tableIDs){
  
  somaNew<-data.frame(metadata,soma[,colnames(soma)%in%selVars])
  somaNew_melt<-melt(data=somaNew,id.vars=colnames(metadata))
  somaNew_melt$GeneID<-paste(tableIDs[match(somaNew_melt$variable,tableIDs[,1]),2],somaNew_melt$variable,sep="_")
  
  somaNew_melt$CLASS<-factor(somaNew_melt$CLASS,ordered=T,levels=c("non_metastasis","metastasis"))
  
  
  pp<-ggplot(somaNew_melt)+
    geom_boxplot(aes(x=CLASS,y=value,fill=CLASS))+
    facet_wrap(~GeneID, scales="free_y",)+
    theme_bw()+
    xlab("")+ylab("log2 normCounts")+
    scale_fill_manual(values=c("white","black")) + 
    theme(axis.text.x=element_text(angle=45,vjust = 0.6,hjust=0.6),legend.position = "bottom",legend.title = element_blank())
  
  print(pp)
  
  return(pp)
}

#function to do different splits of data and store them on a list ----
#the idea is to apply this function to the train data to get slightly different slices of data

applyPrepData <- function(data_train, prepData=prepData,factor_var,split=0.9, n=10) {
  result <- list()  # Create an empty list to store the iterations
  
  for (i in 1:n) {
    iteration <- prepData(data_train,factor_var,split,setSeed = "notReproduce")  # Call the prepData function
    result[[i]] <- iteration  # Store the iteration in the result list
  }
  
  return(result)
}




#function to apply one of the variable selection methods to different splits of data and store them as well as some statistics ----

applyMultiVarSelectionForCrossVal_withModel <- function(inputListwSplits, functionSelection ,thresMod=80,ntree=2000) {
  
  #inputListwSplits assumes an input that will have n elements with each element containing
  # train (train data), test (test data), factor_train, factor_test and idx
  
  #functionSelection is a function to undertake variable selection a number of times. The outputs of these functions
  #must be a list that will have at least one matrix called finalSelVars with 2 columns (probes and frequencies)
  
  #ntree is for the RF model
  
  result <- list()  # Create an empty list to store the iterations
  
  for (i in 1:length(inputListwSplits)) {
    
    train<-as.matrix(inputListwSplits[[i]]$train)
    test<-as.matrix(inputListwSplits[[i]]$test)
    factor_train<-inputListwSplits[[i]]$factor_train
    factor_test<-inputListwSplits[[i]]$factor_test
    
    VarSelec<-functionSelection(train,factor_train,100,25)
    #100 is the number of runs
    # 25 is the min of variables to keep
    #print(VarSelec$)
    
    #fit a model with those Vars (at a given selection thresh) and report the stats
    varsForModel<-VarSelec$finalSelVars[VarSelec$finalSelVars[,2]>thresMod,1]
    print(varsForModel)
    if(length(varsForModel)<2){
      next
    }
    
    rfmodelRes<-RandomForestWithStats(train[,colnames(train)%in%varsForModel],test[,colnames(test)%in%varsForModel],factor_train,factor_test,ntree=ntree)
    
    
    result[[i]] <- list("data"=inputListwSplits[[i]], "varSel"=varsForModel,"VarsModelling"=varsForModel,"rfmodelRes"=rfmodelRes,"finalSelVars"=VarSelec$finalSelVars)
    #in the line above there is a problem and i should have reported as well all the variables seleced by adding as output VarSelec$allVarsFreq or at least VarSelec$finalSelVars, now added as an after thought
  }
  
  
  
  
  return(result)
}



#extract results

concatenateVarSelec<-function(res){
  out<-c()
  for(i in 1:length(res)){
    out<-c(out,res[[i]]$varSel)
  }
  return(out)
}


#function to do 2-class Lasso selection iteratively a number of times and count the number of times each variable is selected ----

TwoClassLasso<-function(train_data=train_data,factor_var_train=factor_var_train,n_runs=100,min_appearances=25){
  
  #train_data: proportion of the data to do the variable selection with
  # factor_var_train - the multimodal factor
  # n_runs - number of times we run Lasso
  # min_appearances =25 - minimum frequency to keep a var
  
  # Initialize variable appearance counter
  var_counts <- data.frame("names"=colnames(train_data), "counts"=rep(0, ncol(train_data)) )
  
  #start a list to store results
  out<-list()
  
  #find lambda
  lambda_seq <- 10^seq(2, -2, by = -.1)
  v_output <- cv.glmnet(train_data, factor_var_train,
                        alpha = 1, lambda = lambda_seq, 
                        nfolds = 5, family = "binomial")
  best_lam <- v_output$lambda.min
  
  # Run Lasso variable selection on training data n_runs times
  for (i in 1:n_runs) {
    # Sample data for Lasso run, slightly different samples in each run
    set.seed(i)
    train_xid <- sample(seq_len(nrow(train_data)), size = floor(0.9*nrow(train_data)))
    train_x <- train_data[train_xid, ]
    factor_var_train_x<-factor_var_train[train_xid]
    
    # Fit Lasso model
    lasso_best <- glmnet(train_x, factor_var_train_x, alpha = 1, lambda = best_lam,family="binomial")
    #pred <- predict(lasso_best, s = best_lam, newx = x_test)
    
    #we select the variables and put them in a matrix
    varsLass<-coef(object = lasso_best,s = best_lam);
    
    varsLass<-as.matrix(varsLass)
    #select the ones that are not zero penalised
    varsLass_no0<-varsLass[varsLass[,1]!=0,];
    #remove intercept
    varsLass_no0<-varsLass_no0[-1]
    
    #update our counter
    var_counts[var_counts[,1]%in%names(varsLass_no0),2]<-var_counts[var_counts[,1]%in%names(varsLass_no0),2]+1
    
    
    #store vars in a list just in case
    out[[i]]<-varsLass
  }
  
  # Select variables that appear at least min_appearances times
  final_selected_vars <- var_counts[var_counts[,2]>min_appearances,]
  
  xout<-list("finalSelVars"=final_selected_vars,"allVarsFreq"=var_counts ,"allRes"=out)
  
  return(xout)
  
}



#function to do backwards elimination selection iteratively a number of times and count the number of times each variable is selected
BackSelec_wTimes<-function(train_data=train_data,factor_var_train=factor_var_train,n_runs=100,min_appearances=25){
  
  #train_data: proportion of the data to do the variable selection with
  # factor_var_train - the multimodal factor
  # n_runs - number of times we run Lasso
  # min_appearances =25 - minimum frequency to keep a var
  
  # Initialize variable appearance counter
  var_counts <- data.frame("names"=colnames(train_data), "counts"=rep(0, ncol(train_data)) )
  
  #start a list to store results
  out<-list()
  
  # Run Lasso variable selection on training data n_runs times
  for (i in 1:n_runs) {
    # Sample data for Lasso run, slightly different samples in each run
    set.seed(i)
    train_xid <- sample(seq_len(nrow(train_data)), size = floor(0.9*nrow(train_data)))
    train_x <- train_data[train_xid, ]
    factor_var_train_x<-factor_var_train[train_xid]
    
    # Fit Lasso model
    rfModel_backselVars<-varSelRF(xdata=train_x,Class=factor_var_train_x, ntree=2000, ntreeIterat=500,vars.drop.num = NULL,vars.drop.frac = 0.1)
    
    
    
    #we select the variables and put them in a matrix
    varsBack<-rfModel_backselVars$selected.vars
    
    
    # Update variable appearance counter
    
    var_counts[var_counts[,1]%in%varsBack,2] <- var_counts[var_counts[,1]%in%varsBack,2]+1
    
    #store vars in a list just in case
    out[[i]]<-varsBack
  }
  
  # Select variables that appear at least min_appearances times
  final_selected_vars <- var_counts[var_counts[,2]>min_appearances,]
  
  xout<-list("finalSelVars"=final_selected_vars,"allVarsFreq"=var_counts ,"allRes"=out)
  
  return(xout)
  
}



Boruta_wTimes <- function(train_data = train_data, factor_var_train = factor_var_train, n_runs = 100, min_appearances = 25) {
  # Initialize variable appearance counter
  var_counts <- data.frame("names" = colnames(train_data), "counts" = rep(0, ncol(train_data)))
  
  #start a list to store results
  out <- list()
  
  # Run Boruta variable selection on training data n_runs times
  for (i in 1:n_runs) {
    # Sample data for each run, slightly different samples in each run
    set.seed(i)
    train_xid <- sample(seq_len(nrow(train_data)), size = floor(0.9 * nrow(train_data)))
    train_x <- train_data[train_xid, ]
    factor_var_train_x <- factor_var_train[train_xid]
    
    # Create Boruta object
    boruta_obj <- Boruta(train_x, factor_var_train_x)
    
    # Perform Boruta feature selection
    boruta_res <- attStats(boruta_obj)
    
    # Get selected variables
    varsBack <- rownames(boruta_res[boruta_res$decision == "Confirmed", ])
    
    # Update variable appearance counter
    var_counts[var_counts[, 1] %in% varsBack, 2] <- var_counts[var_counts[, 1] %in% varsBack, 2] + 1
    
    #store vars in a list just in case
    out[[i]] <- varsBack
  }
  
  # Select variables that appear at least min_appearances times
  final_selected_vars <- var_counts[var_counts[, 2] > min_appearances, ]
  
  xout <- list("finalSelVars" = final_selected_vars, "allVarsFreq" = var_counts, "allRes" = out)
  
  return(xout)
}

Model_RF <- function(pdac,geneset,nround) { 
  #Run 50 times and find average 
  PC0_RF <- 0 ; PC1_RF <- 0;RC0_RF <- 0;RC1_RF <- 0;F10_RF <- 0;F11_RF <- 0;MF1_RF<-0 ; PR_AUC_RF <- 0;AUC_RF <- 0;MPC_RF <-0 ; MRC_RF <-0 ; MCC_RF <-0
  for (i in 1:nround) {
    print(i)
    set.seed(i)
    index <- createDataPartition(pdac$CLASS, p=0.9 , list = F ,times = 1)
    #train and test set 
    train_df <- pdac[index, ]
    test_df <-  pdac[-index, ]
    
    #convert CLASS to factor 
    train_df$CLASS <- as.factor(train_df$CLASS )
    test_df$CLASS <- as.factor(test_df$CLASS )
    
    train.data  <- train_df[ ,which( colnames(train_df) %in% geneset ) ] 
    dim(train.data)
    
    test.data  <- test_df[ ,which( colnames(test_df) %in% geneset) ] 
    dim(test.data)
    
    train.data  <- cbind(train_df$CLASS , train.data)
    names(train.data)[1] <- "CLASS"
    
    test.data  <- cbind(test_df$CLASS , test.data)
    names(test.data)[1] <- "CLASS"
    
    #cvIndex <- createFolds(factor(imbal_train$CLASS), 5, returnTrain = T)
    ctrlspecs <- trainControl(method = "cv" , number = 5,
                              savePredictions = "all",#"all" (for true) or "none" (for false).
                              classProbs = TRUE) #,index = cvIndex)
    
    
    #We use ADASYN 
    genData_ADAS = ADAS(X= train.data[,-1], target = train.data$CLASS ,K=5  )
    train_df_adas <- genData_ADAS[["data"]]  
    # keep CLASS as a df , last column
    CLASS_df <- train_df_adas[ncol(train_df_adas)]
    # move CLASS to first column 
    train_df_adas <- cbind(train_df_adas[ncol(train_df_adas)], train_df_adas[,-ncol(train_df_adas)] )
    names(train_df_adas)[1] <- "CLASS"
    train.data <- train_df_adas
    
    orig_fit <- caret::train(CLASS ~ .,
                             data = train.data,
                             method = "ranger",
                             verbose = FALSE,
                             trControl = ctrlspecs,
                             metric = 'Accuracy',#scale_pos_weight = 3,
                             tuneLength = 10)
    # Prediction 
    prob_RF <- predict(orig_fit, test.data, type="prob")
    pred_RF <- predict(orig_fit, test.data)
    res_RF <- Res_mod(prob = prob_RF , pred = pred_RF ,test_data = test.data)
    
    
    #Original
    PC0_RF[i] <- res_RF[[1]] ; PC1_RF[i] <- res_RF[[2]]
    RC0_RF[i] <- res_RF[[3]] ;RC1_RF[i] <- res_RF[[4]]
    F10_RF[i] <- res_RF[[5]] ; F11_RF[i] <- res_RF[[6]]
    MPC_RF[i] <-  res_RF[[7]] ; MRC_RF[i] <-  res_RF[[8]]
    MF1_RF[i] <- res_RF[[9]] ; PR_AUC_RF[i]<- res_RF[[10]]
    AUC_RF[i] <- res_RF[[11]];MCC_RF[i] <- res_RF[[12]]
    
    # Original
    df_res_RF <- data.frame(PC0 = PC0_RF, 
                            PC1= PC1_RF ,
                            RC0 = RC0_RF, 
                            RC1= RC1_RF ,
                            F10= F10_RF,
                            F11=F11_RF, 
                            MPC= MPC_RF,
                            MRC= MRC_RF ,
                            MF1= MF1_RF,
                            ROC_AUC = AUC_RF ,
                            PR_AUC= PR_AUC_RF,
                            MCC= MCC_RF )
    
    df_res_RF <- as.data.frame(t(df_res_RF))
    df_res_RF[is.na(df_res_RF)] = 0
    df_res_RF$mean <-  apply(df_res_RF, 1, mean)
    df_res_RF$sd <-  apply(df_res_RF, 1, sd)
    df_res_RF$CI95 <- 2*sqrt( (df_res_RF$mean * (1 - df_res_RF$mean)) / nround)
    df_res_RF <- round(df_res_RF,3)
    #View(df_res_RF)
  }
  return(df_res_RF)
}

Res_mod <- function(prob,pred,test_data) {
  #AUC 
  roc_c <- pROC::auc(test_data$CLASS, prob[,2])
  ROC_AUC <- roc_c[1]
  PR_AUC <- prauc(truth = test_data$CLASS , prob=prob[, 2],positive= "non_metastasis")
  
  #confusion matrix
  actual = as.factor(test_data$CLASS)
  predicted = as.factor(pred)
  cm = as.matrix(table(Actual = actual, Predicted = predicted)) # create the confusion matrix
  
  n = sum(cm) # number of instances
  nc = nrow(cm) # number of CLASSes
  diag = diag(cm) # number of correctly CLASSified instances per CLASS 
  rowsums = apply(cm, 1, sum) # number of instances per CLASS
  colsums = apply(cm, 2, sum) # number of predictions per CLASS
  p = rowsums / n # distribution of instances over the actual CLASSes
  q = colsums / n # distribution of instances over the predicted CLASSes
  
  precision = diag / colsums 
  recall = diag / rowsums 
  f1 = 2 * precision * recall / (precision + recall) 
  #data.frame(precision, recall, f1) 
  macroPrecision = mean(precision)
  macroRecall = mean(recall)
  macroF1 = mean(f1)
  mcc = mltools::mcc(preds = predicted,actuals = actual)
  
  #All results 
  Res <- matrix(c(precision,recall,f1,macroPrecision,macroRecall, macroF1, PR_AUC, ROC_AUC,mcc))
  return(Res) #we tell the function what to output
}


#function to split data in train/test----
prepData<-function(data,factor_var,split=0.8,setSeed="reproduce"){
  # for reproducibility
  if(setSeed=="reproduce"){
    set.seed(123)
    
  }
  
  
  train_idx <- sample(seq_len(nrow(data)), size = floor(split*nrow(data)))
  train_data <- data[train_idx, ]
  test_data <- data[-train_idx, ]
  factor_var_train<-factor(factor_var[train_idx])
  factor_var_test<-factor(factor_var[-train_idx])
  print(factor_var_train)
  #modified the function a posteriori to also return the indeces as that was a mistake on my part!
  out<-list("train"=as.matrix(train_data),"test"=as.matrix(test_data),"factor_train"=factor_var_train,"factor_test"=factor_var_test,"idx"=train_idx)
  
}
