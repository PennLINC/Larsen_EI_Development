## Load Libraries
require(ggplot2)
require(dplyr)
library(e1071)
# library(ROCR)
library(matrixTests)
library(doParallel)
theme_set(theme_classic())
setwd("/cbica/projects/alpraz_EI/scripts/")

extract_matrix2 <- function(subid,sesid,atlasname,data_type,subdivide=FALSE, subdivision=NULL){
  # This function loads the connectivity matrix for a subj, filters features (if needed; subdivide=T and specify subdivision), and converts to a vector
  if (data_type == "Alpraz") {
    fname <- sprintf("/cbica/projects/alpraz_EI/data/TASK_GSR/%s/%s/%s_CC_GSR_000.netcc",subid,sesid,atlasname)
  } else if (data_type == "PNC") {
    fname <- sprintf("/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/%s/%s/%s_CC_000.netcc",subid,sesid,atlasname)
  }
  print(fname)
  thisCC <- as.matrix(read.table(fname,header = T,skip = 4))
  
  # Select specific features if desired.
  if (subdivide==T) {
    # We want to subidivide the features and select a subset of them.
    templateCC <- thisCC
    
    # Get labels
    labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_labels.csv',atlasname),header=F)
    labels$nums = as.numeric(labels$V1)
    if (atlasname %in% c("schaefer400x7_aal","schaefer200x7_aal")) {
      community_labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_node_labels.csv',
                                           gsub(x = atlasname,pattern = "_aal",replacement = ""))
                                   ,header=T)
      if (atlasname == "schaefer200x7_aal") {
        trans_labels <- read.csv('/cbica/projects/alpraz_EI/input/atlases/Schaef200_transmodality7.csv',header = F)
      } else{
        trans_labels <- read.csv('/cbica/projects/alpraz_EI/input/atlases/Schaef400_transmodality7.csv',header = F)
      }
      
      # Just doing some re-naming and joining
      colnames(trans_labels) <- c('name','transmodality')
      trans_labels$name <- gsub(trans_labels$name,pattern = "7Networks_",replacement = "")
      trans_labels$nums = as.numeric(rownames(trans_labels))
      
      
      mat_labels <- data.frame(nums = colnames(templateCC))
      mat_labels$nums <- as.numeric(gsub(mat_labels$nums,pattern="X",replace=""))
      mat_labels <- mat_labels %>% left_join(labels,by = "nums")
      mat_labels <- mat_labels %>% left_join(community_labels,by = c("nums"="index"))
      mat_labels <- mat_labels %>%
        mutate(community_name = as.character(community_name),V2 = as.character(V2),node_names=as.character(node_names))
      mat_labels$community_name[is.na(mat_labels$community_name)]="AAL_subcortex"
      mat_labels$node_names[is.na(mat_labels$node_names)]= mat_labels$V2[is.na(mat_labels$node_names)]
      
      # Check if we are using transmodality scores
      if (subdivision %in% c("transmodal25","unimodal25")) {
        #add transmodality scores
        mat_labels <- mat_labels %>%
          left_join(trans_labels,by = "nums")
        
        #cut transmodality into quantiles
        cuts <- quantile(mat_labels$transmodality,c(0,1/4,3/4,1),na.rm = T)
        mat_labels <- mat_labels %>% 
          mutate(trans_category=cut(transmodality, breaks=cuts, labels=c("unimodal","middle","transmodal"),include.lowest = T))
        # Set names to transmodality names
        colnames(thisCC)=mat_labels$trans_category
        rownames(thisCC)=mat_labels$trans_category
        
        if (grepl(x=subdivision,pattern = "transmodal")) {
          communities="transmodal"
          commCC <- thisCC[rownames(thisCC) %in% communities,]
          vals <- as.vector(commCC)
        } else if (grepl(x=subdivision,pattern = "unimodal")) {
          communities="unimodal"
          commCC <- thisCC[rownames(thisCC) %in% communities,]
          vals <- as.vector(commCC)
          length(vals)
        } 
      } 
    }
  } 
  else{
    thisCC_vector = thisCC[upper.tri(thisCC,diag = F)]
    vals <-thisCC_vector
  }
  return(vals)
}

CV_function <- function(sublist,df,num_folds){
  ## Function for cross-validation
  foldIdxs <- data.frame(subid=sublist)
  foldIdxs$foldID <- row_number(foldIdxs$subid)
  foldIdxs$foldID <- ntile(foldIdxs$foldID,num_folds)
  fold_output<-vector(mode = "list",length = num_folds)
  for (fold in 1:num_folds) {
    trainingIDs <- as.matrix(foldIdxs %>% filter(foldID != fold) %>% select(subid))
    trainingIndex <- df$subid %in% trainingIDs # indices for training subs
    trainingData <- df[trainingIndex, 3:dim(df)[2] ]
    testData <- df[!trainingIndex, 4:dim(df)[2]] # test data. Take columns 4:end (Removes subid sesid drug).
    testLabels <- data.frame(df[!trainingIndex,c(1:3) ])
    # svm
    x <- as.matrix(trainingData[, 2:dim(trainingData)[2]])
    y <- as.factor(as.matrix(trainingData[,1]))
    svm.model <- svm(x =x, y = y, 
                     cost = 1, kernel = "linear",type = "C-classification",scale = F)
    svm.pred <- predict(svm.model, as.matrix(testData))
    
    decisionValues <- w %*% t(testData)-svm.model$rho
    distance <- decisionValues/norm(w)
    testLabels$decisionValues <- t(decisionValues)
    testLabels$distance <- t(distance)
    testLabels$model.pred = svm.pred
    fold_output[[fold]]<-testLabels
  }
  pred_out <- data.table::rbindlist(fold_output)
  return(pred_out)
}

SVM_2class <- function(df,folds,feature_selection = F,feature_proportion = .1,num_repetitions = 100){
  #SVM 2-class classifier
  # folds: number of folds for cv.
  # feature_selection: optional flag for data-driven feature selection. Selects the proportion of features indicated by feature_proportion. Not used by default.
  # num_repetitions: How many times to you want to repeat the cv process using different random splits of the data. This is just extra caution against a randomly good or bad split.
  
  cat('\nRunning SVM models.....')
  # Set up folds for CV
  if (folds == "LOO") {
    # This is for leave-one-out cross-validation
    num_folds = length(unique(df$subid))
    num_repetitions <- 1
  } else {
    num_folds = folds
    num_repetitions <- num_repetitions
  }

  svm_output <- vector(mode = "list",length = num_repetitions) #set up output object.
  
  # Create the folds
  unique_IDs <- unique(df$subid)
  subid_folds <- replicate(num_repetitions,sample(unique_IDs,size = length(unique_IDs))) # create sets of random draws.
  foldIdxs <- data.frame(subid=unique_IDs)
  foldIdxs$foldID <- row_number(foldIdxs$subid)
  foldIdxs$foldID <- ntile(foldIdxs$foldID,num_folds)
  # cat('Sending data to CV')
  for (r in 1:num_repetitions) {
    foldIdxs$subid <- subid_folds[,r] # Grab a random split for folds.
    fold_output<-vector(mode = "list",length = num_folds)
    cat(sprintf('\nrepetition  %d.... ',r))
    for (fold in 1:num_folds) {
      # cat(sprintf('\nfold  %d.... ',fold))
      trainingIDs <- as.matrix(foldIdxs %>% filter(foldID != fold) %>% select(subid))
      trainingIndex <- df$subid %in% trainingIDs # indices for training subs
      trainingData <- df[trainingIndex, 3:dim(df)[2] ] %>% arrange(drug) #this is important because libsvm automatically makes the first observation class 1, so drug must be first every time. Placebo will be class -1.
      testData <- df[!trainingIndex, 4:dim(df)[2]] # test data. Take columns 4:end (Removes subid sesid drug).
      testLabels <- data.frame(df[!trainingIndex,c(1:3) ])
      # svm
      x <- as.matrix(trainingData[, 2:dim(trainingData)[2]])
      y <- as.factor(as.matrix(trainingData[,1]))
      svm.model <- svm(x =x, y = y, 
                       cost = 1, kernel = "linear",type = "C-classification",scale = F)
      svm.pred <- predict(svm.model, as.matrix(testData))
      
      w <- t(svm.model$coefs) %*% svm.model$SV #calculate feature weights.
      # num_features <- dim(x)[2]
      decisionValues <- w %*% t(testData)-svm.model$rho # Get decision valus.
      distance <- decisionValues/norm(w) #calculate distance from the classification hyperplane
      # just adding the results to the dataframe.
      testLabels$decisionValues <- t(decisionValues)
      testLabels$distance <- t(distance)
      testLabels$model.pred = svm.pred
      fold_output[[fold]]<-testLabels
    }
    svm_output[[r]] <- data.table::rbindlist(fold_output) # saving output for the repetition.
    cat('complete\n')
  }
  # output_list <- apply(subid_folds,CV_function,df=df,num_folds = num_folds,MARGIN = 2)
  
  # Now train a model using all the data. This is for the estimation of the feature weights using the most possible data.
  final_data<-df[, 3:dim(df)[2] ]
  x <- as.matrix(final_data[, 2:dim(final_data)[2]])
  y <- as.factor(as.matrix(final_data[,1]))
  svm.model <- svm(x = x, y = y, 
                   cost = 1, kernel = "linear",type = "C-classification",scale = F)
  w <- t(svm.model$coefs) %*% svm.model$SV
  svm_results <- list(svm_output,w,svm.model)
  cat('\nFinished SVM models\n')
  
  return(svm_results)
}

AUC <- function(DecisionValues, labels){
  # This function creates an ROC curve and then calculates AUC.
  # Decision values is nx1
  # Labels is nx1
  
  # N.B.
  # Drug (class 0) is assigned as +1 by libsvm and placebo (class 1) is assigned as 0 by libsvm default. 
  # Adjust the true labels and predicted labels to match that here so that the decision values make sense.
  labels <- labels*-1+1
  
  P <- sum(labels == 1)
  N <- sum(labels == 0)
  Sorted_DecisionValues <- sort(unique(DecisionValues), decreasing = FALSE)
  numDecisionValues <- length(Sorted_DecisionValues)
  
  TP_Array <- vector(mode = "numeric",length = numDecisionValues)
  FP_Array <- vector(mode = "numeric",length = numDecisionValues)
  Accuracy_Array = vector(mode = "numeric",length = numDecisionValues)
  for (i in 1:numDecisionValues){
    thisCutoff <- Sorted_DecisionValues[i]
    thisPredictedLabels <- as.numeric(DecisionValues>thisCutoff)
    detections <- thisPredictedLabels==1
    
    TP <- sum(labels[detections] == thisPredictedLabels[detections])
    TPR <- TP/P
    FP <- sum(labels[detections]!=thisPredictedLabels[detections])
    FPR <- FP/N
    
    TP_Array[i] <- TPR
    FP_Array[i] <- FPR
    
    Accuracy_Array[i] = (TP + N - FP) / (P + N)
  }
  
  ROC_output <- data.frame(TPR =TP_Array,FPR=FP_Array,Accuracy = Accuracy_Array)
  ROC_output <- ROC_output%>%arrange(TPR,FPR)
  
  #AUC
  dFPR <- c(0,diff(ROC_output$FPR))
  dTPR <- c(0,diff(ROC_output$TPR))
  AUC <- sum(ROC_output$TPR * dFPR) + sum(dTPR * dFPR)/2
  return(AUC)
}

run_model <- function(df,folds,feature_selection = F,feature_proportion = .1,permutation_test = F, num_permutations = 10000,type = "svm"){
  ## This sends the data to the desired classifier and other functions
  if (feature_selection == T) {
    # Not done for this study #
    cat(sprintf("\nPerforming feature extraction: extracting the top %1.3f features\n",feature_proportion))
  }
  if (type == "svm") {
    print("Performing classification using SVM")
    model_results <- SVM_2class(df,folds,feature_selection = feature_selection,feature_proportion = feature_proportion,num_repetitions = 100)
  } # No other types for now
  
  prediction_output <- model_results[[1]]
  # Now we just organize the output and return the important info.
  W<-model_results[[2]]
  svm.model <- model_results[[3]]
  num_features <- sum(!is.na(W[1,]))
  if (feature_selection==T) {
    feature_counts <- colSums(!is.na(W))
    W <- rbind(W,feature_counts)
  } else{
    W <- W
  }
  accuracy_fun <- function(x) sum(x$model.pred==x$drug)/dim(x)[1] # function to calculate the accuracy. 
  accuracies <- sapply(prediction_output,accuracy_fun) #Apply the function to each repetition of the cross-validation.
  
  accuracy <- mean(accuracies)
  num_obs <- length(prediction_output[[1]]$model.pred)
  num_correct <- round(accuracy*num_obs)
  b <- binom.test(num_correct,num_obs,.5) #This is a binomial test. The p-value is not used (permutation test is used instead), but this function summarizes the data nicely.
  b$pred_data <- prediction_output
  b$accuracy <- accuracy
  if (permutation_test == T) {
    print(sprintf("Permuting %d times...",num_permutations))
    ptm = proc.time()
    # This function is set up to run in parallel, but the cluster SGE is not handling this well or consistently, so it is set to be serial now.
    # cl<-makeCluster(20)
    # registerDoParallel(cl)
    # nw <- getDoParWorkers()
    nw=1
    # # perms <- shuffleSet(n = dim(df)[1],nset = num_permutations)
    perm_acc <- matrix(nrow = num_permutations)
    perm_W <- matrix(nrow = num_permutations,ncol = num_features)
    perm_list <- list()
    perm_list = foreach(perm_chunk = idiv(num_permutations,chunks = nw),
                        .combine=c,
                        .export = c("featureExtraction","SVM_2class"),
                        .packages = c("dplyr","e1071","matrixTests")) %do% {  # This must be `dopar` to be parallel.
                          # pacc = numeric(perm_chunk)
                          # pW = matrix(nrow = perm_chunk,ncol = num_features)
                          perm_result_list=list()
                          for (p in 1:perm_chunk){
                            # thisPerm <- perms[p,]
                            cat(sprintf("\nPermutation %d",p))
                            thisDF <- df
                            # thisDF$drug <- df$drug[thisPerm]
                            permuted <- df %>% select(subid,drug) %>% group_by(subid) %>% mutate(perm_drug=sample(drug))
                            thisDF$drug <- permuted$perm_drug
                            perm_pred_result <- SVM_2class(thisDF,folds,
                                                           feature_selection = feature_selection,
                                                           feature_proportion = feature_proportion,
                                                           num_repetitions = 1)
                            # perm_pred_data <- perm_pred_result[[1]][[1]] #$pred_data[[1]]
                            # perm_W <- perm_pred_result[[2]]
                            # perm_num_features <- sum(!is.na(perm_W))
                            # pacc<-sum(perm_pred_data$model.pred==perm_pred_data$drug)/dim(perm_pred_data)[1]
                            # pdf <- data.frame(acc=pacc,pred_data=perm_pred_data)  #make a data.frame to hold the output
                            # names(pdf) = c("acc",names(perm_pred_data)) # fix names
                            # pdf <- pdf %>% group_by(acc) %>%nest(pred_data = c(names(perm_pred_data))) # nest the pred_data df so we have one line per permutation
                            # perm_result_list[p]=list(pred_output = pdf,W=perm_W) #add the data.frame to the list
                            perm_result_list[[p]] = list(perm_pred_result)
                          }
                          perm_result_list
                        }
    print("done")
    # stopCluster(cl)
    print(proc.time()-ptm)
    
    #organizing the output
    pred_data_list <- lapply(perm_list, function(x) x[[1]][[1]])
    perm_Ws <- lapply(perm_list,function(x) x[[1]][[2]])
    
    # Calculating permutation accuracies and AUCs
    perm_acc_distribution <- sapply(pred_data_list, function(x) sum(x[[1]]$model.pred==x[[1]]$drug)/length(x[[1]]$drug))
    perm_p <- sum(perm_acc_distribution>accuracy)/length(perm_acc_distribution)
    perm_auc_distribution <- sapply(pred_data_list, function(x) AUC(DecisionValues=x[[1]]$decisionValues,labels=x[[1]]$drug))

    cat(sprintf("\nPermutation p-value =  %1.3f\n",perm_p))
    W_test <- sapply(perm_Ws, FUN = function(x) {
      abs(x)>abs(W)
    })
    W_sig <- rowMeans(W_test)
    
    b$perm_p <- perm_p
    b$perm_W_sig <- W_sig
    # b$perm_Ws <- perm_Ws #this takes up way too much space
    b$perm_accs = perm_acc_distribution
    b$perm_aucs = perm_auc_distribution
    print(b$perm_p)
  }
  
  print(sprintf("Overall Accuracy: %1.3f; p = %.5f\n\n",accuracy,b$p.value))
  
  return(list(b,W,num_features,svm.model))
}

featureExtraction <- function(trainingData, feature_proportion = .1){
  # This is not used, but it will perform data driven feature selection if desired.
  features <- as.matrix(trainingData %>% select(-drug))
  drug <- features[trainingData$drug==0,]
  placebo <- features[trainingData$drug==1,]
  
  tstats<-col_t_paired(placebo,drug)
  newdata=t(tstats[order(-abs(tstats$statistic)),])
  keep_features = colnames(newdata)[1:round(feature_proportion*dim(newdata)[2])]
  
  newTrainingData <- trainingData %>% select(drug,keep_features)
  return(list(newTrainingData,keep_features))
}

#########################
####### Load Data #######
########################

FD_thresh = .5 #matching what was used in Wolf et al.
subData <- read.csv('/cbica/projects/alpraz_EI/input/alpraz_sublist_FD2.csv')
subInfo <- subData %>% filter(exists==1 & motion_pass==1)

PNCData <- read.csv('/cbica/projects/alpraz_EI/input/PNC_sublist_FD.csv')
PNCInfo <- PNCData %>% filter(exists==1 & FD < FD_thresh) # Filtered here, but this isn't actually necessary since we convervatively threshold during the developmental analyses at 0.3 mm

#####################
#### SET OPTIONS ####
#####################
classifier = "svm"
data_type ="Alpraz"

# Permutation test?
# Do we want to use a permutation test for significance testing?
# "permute_on" = yes
# "permute_off" = no
perm_test="permute_on"
num_permutations = 1000

# Do we want to pull out only certain brain areas?
# Set subdivide = TRUE if this is desired
# Set subdivision desired:
## "transmodal25" = top 25% most transmodal regions
## "unimodal25" = top 25% most unimodal regions
## "all" = all regions
## "regional" = perform classification separately for each region in the atlas.
subdivide = TRUE
subdivision = "unimodal25"
cat(sprintf("\nsubdivision = %s\n",subdivision))

# Atlas and FE  
## Create a list of atlases that we want to use for classification. We are using only schaefer 400 for primary analyses, but we also replicate with "gordon333_aal","glasser360_aal","aal116", "schaefer200x7_aal". 
atlas_list = c("schaefer400x7_aal") 
fe_list=c(1)
###############

atlas_acc <- matrix(ncol = 6)
colnames(atlas_acc)=c("accuracy", "p.value", "fe","perm.p","atlas","num_features")

## Loop over atlases and run classification
for (atlasname in atlas_list){
  cat(atlasname)
  ## Load the data
  if (file.exists(sprintf("/cbica/projects/alpraz_EI/input/CorMats/%s_%s_%s.rds",atlasname,data_type,subdivision))) {
    # The data has already been compiled into a data.frame, just load it.
    df <- readRDS(sprintf("/cbica/projects/alpraz_EI/input/CorMats/%s_%s_%s.rds",atlasname,data_type,subdivision))
  } else {
    # Compile the data if we haven't already.
    gm<- subInfo %>%
      group_by(subid, sesid) %>%
      mutate(mat = list(extract_matrix2(subid,sesid,atlasname,data_type,subdivide=subdivide,subdivision = subdivision)))
    
    if (subdivision=="regional") {
      mat <- do.call(rbind, lapply(gm$mat, function(x) data.frame(x)))
      df <- gm %>% full_join(mat,by=c("subid","sesid"))%>%select(-mat)
    } else {
      mat <- data.frame(matrix(unlist(gm$mat), nrow=dim(gm[1]), byrow=T))
      df <- cbind(gm %>%
                    ungroup() %>%
                    select(subid,sesid,drug)
                  ,mat)
    }
    
    print("saving...")
    saveRDS(df,file = sprintf("/cbica/projects/alpraz_EI/input/CorMats/%s_%s_%s.rds",atlasname,data_type,subdivision))
    print("saved")
    rm(gm)
  }
  
  ## Now run the SVM for a set of feature extraction levels or regions.
  # df<- df %>% arrange(drug)
  acc = matrix(nrow = length(fe_list),ncol = 5)
  for (n in 1:length(fe_list)) {
    num <- fe_list[n]
    if (file.exists(
      sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_%s_%s_%s_%s_results.rds",
              atlasname,data_type,subdivision,classifier,as.character(fe_list[n]),perm_test
      ))
    ) {
      results <- readRDS(
        file = sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_%s_%s_%s_%s_results.rds",
                       atlasname,data_type,subdivision,classifier,as.character(fe_list[n]),perm_test))
    } else{
      if (num == 1) {
        feature_selection = FALSE
      } else {
        feature_selection =TRUE
      }
      
      if (perm_test == "permute_on") {
        permutation_test = TRUE
      } else {
        permutation_test = FALSE
      }
      
      results <- run_model(df = df,folds = 10,
                           feature_selection = feature_selection,feature_proportion = num,
                           permutation_test = permutation_test,num_permutations = num_permutations,
                           type = classifier)
      results[[5]]=atlasname
      cat('\nFinished SVM\n')
      
      cat(
        sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_%s_%s_%s_%s_results.rds\n",
                atlasname,data_type,subdivision,classifier,as.character(fe_list[n]),perm_test)
        )
      saveRDS(results,
              file = sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_%s_%s_%s_%s_results.rds",
                             atlasname,data_type,subdivision,classifier,as.character(fe_list[n]),perm_test))
      cat('saved\n')
    }
  }
  
  ## Now we use the trained model on the PNC data.
  ## Load the PNC data
  cat("\nLoading PNC data\n")
  if (file.exists(sprintf("/cbica/projects/alpraz_EI/input/CorMats/PNC/%s_%s.rds",atlasname,subdivision))) {
    df <- readRDS(sprintf("/cbica/projects/alpraz_EI/input/CorMats/PNC/%s_%s.rds",atlasname,subdivision))
  } else {
    gm<- PNCInfo %>% 
      group_by(subid,sesid) %>%
      mutate(mat = list(extract_matrix2(subid,sesid,atlasname,"PNC",subdivide=subdivide,subdivision = subdivision)))
    if (atlasname=="gordon333_aal") {
      feat_count=gm%>%summarise(n=length(unlist(mat)))
      good_bbl=feat_count%>%filter(n==19701)%>%select(subid)
      gm <- gm %>% filter(subid%in%good_bbl$subid)
    }
    mat <- data.frame(matrix(unlist(gm$mat), nrow=dim(gm[1]), byrow=T))
    df <- cbind(gm %>% 
                  ungroup() %>% 
                  select(subid,sesid)
                ,mat)
    print("saving...")
    saveRDS(df,file = sprintf("/cbica/projects/alpraz_EI/input/CorMats/PNC/%s_%s.rds",atlasname,subdivision))
    print("saved")
    rm(gm)
  }
  cat("Done\n")
  
  #Apply the trained model
  svm.model <- results[[4]]
  testData <- as.matrix(df[, 3:dim(df)[2]])
  svm.pred <- predict(svm.model, as.matrix(testData))
  
  # Organize the outputs
  w <- t(svm.model$coefs) %*% svm.model$SV
  decisionValues <- w %*% t(as.matrix(testData))-svm.model$rho
  distance <- decisionValues/norm(w)
  df$decisionValues <- t(decisionValues)
  df$distance <- t(distance)
  df$pred <-svm.pred
  #save the output
  cat("saving PNC output....")
  saveRDS(df%>%select(subid,sesid,pred,distance,decisionValues),file = sprintf("/cbica/projects/alpraz_EI/output/PNC_predictions/%s_%s.rds",atlasname,subdivision))
  cat("Done!\n")
}

