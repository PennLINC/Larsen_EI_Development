# Varius functions that organize and plot data from the study.
# library(tidyverse)
# library(oro.nifti)
library(neurobase)
library(cowplot)
library(ggridges)
library(pracma)
# library(purrr)
# library(mgcv)
library(gratia)
# library(broom)
library(scales)
# library(kableExtra)
# library(RColorBrewer)
setwd("/cbica/projects/alpraz_EI/SubmissionGithub/")
source('tools/RainCloudPlots/tutorial_R/R_rainclouds.R')
font_size <- 16
unimodal_color = "#1c405eff"
transmodal_color = "#285c27ff"
theme_set(theme_classic(base_family = "sans",base_size = font_size))
theme_replace(axis.text=element_text(colour = "black",size = font_size))
line_size <- 1.5
point_size <- 2

feat2mat <- function(W,atlasname,community_summary = F,plots = T,templateCC = NULL,returnRankingsOnly = F,GABA=T,transmodality=T,surface=F){
  if (is.null(templateCC)) {
    fname <- sprintf("/cbica/projects/alpraz_EI/data/TASK_GSR/13783/1986/%s_CC_GSR_000.netcc",atlasname)
    templateCC <- as.matrix(read.table(fname,header = T,skip = 4))
  }
  
  #load template matrix (picked 13783 1986 at random)
  
  templateCC_vector = templateCC[upper.tri(templateCC,diag = F)]
  if (length(templateCC_vector == length(W))) {
    templateCC_vector = W
  } else {
    stop("Length of W does not match the supplied atlas")
  }
  newCC <- templateCC
  newCC[upper.tri(newCC,diag = T)] = 0
  newCC[upper.tri(newCC,diag = F)] = templateCC_vector
  newCC[lower.tri(newCC)] <- t(newCC)[lower.tri(newCC)]
  
  # Get labels
  
  labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_labels.csv',atlasname),header=F)
  mat_labels <- data.frame(nums = colnames(templateCC))
  mat_labels$nums <- as.numeric(gsub(mat_labels$nums,pattern="X",replace=""))
  mat_labels <- mat_labels %>% left_join(labels,by = c("nums"="V1"))
  mat_labels <- mat_labels %>%
    mutate(node_names = as.character(V2))
  
  # Set names to node names
  colnames(newCC)=mat_labels$node_names
  rownames(newCC)=mat_labels$node_names
  
  rvals_abs <- rowSums(abs(newCC)) #get the sum of absolute values
  rvals <- rowSums(newCC) #get the sum of signed values
  rval_df <- data.frame(rvals = rvals,rvals_abs = rvals_abs, names = as.character(names(rvals)))
  rval_df <- rval_df %>% 
    left_join(mat_labels,by = c("names" = "node_names"))
  
  if (GABA == T) {
    cat("GABRA gene analysis....\n")
    corr_vec <- vector(mode="double",length = 6)
    for (n in c(1,2,3,4,5,6)) {
      GABA_vals <- read.table(sprintf(
        '/cbica/projects/alpraz_EI/input/atlases/GABRA%d/GABRA%d_%s_threshold_0.95_vals.csv',n,n,atlasname),
        header = T)
      
      GABA_mat_labels <- rval_df %>% left_join(GABA_vals,by = c("nums"="label"))
      GABA_mat_labels$node_names <- GABA_mat_labels$names
      GABA_mat_labels$GABRA <- GABA_mat_labels$value
      colnames(GABA_mat_labels)<- gsub(colnames(GABA_mat_labels),pattern="GABRA",replacement = sprintf("GABRA%d",n))
      GABA_mat_labels$SVM_weights <- GABA_mat_labels$rvals
      GABA_mat_labels$absolute_SVM_weights <- GABA_mat_labels$rvals_abs
      brainsmash(GABA_mat_labels,c(sprintf("GABRA%d",n),"SVM_weights"),atlasname=atlasname)
      brainsmash(GABA_mat_labels,c(sprintf("GABRA%d",n),"absolute_SVM_weights"),atlasname=atlasname)

      if (atlasname=="schaefer400x7_aal") {
        cortex = GABA_mat_labels %>% filter(nums<4000)
      } else {
        cortex = GABA_mat_labels %>% filter(nums<4000)
      }
      corr_vec[n] <- cor(cortex$SVM_weights,cortex$GABRA)
    }
    print(corr_vec)
    if (file.exists(sprintf('surrogate_brainmap_corrs_%s_GABRA3_map.txt',atlasname))) {
      x=data.table::rbindlist(lapply(list.files(pattern = glob2rx("surr*brainmap*schaefer400*GAB*.txt")), data.table::fread),idcol="GABRA")
      names(x) <- c("GABRA", "r")
      x<-x %>%
        mutate(GABRA = factor(GABRA,levels = c(1,2,3,4,5,6),labels = c("GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6")),
               R2=r^2,
               source = "BrainSMASH null")
      corr_df <- data.frame(GABRA = factor(c("GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRA6")),
                            r = corr_vec,
                            R2 = corr_vec^2,
                            source = "Observed")
      print(corr_df)
      # x <- rbind(x,corr_df)
      x <- x %>%
        mutate(benzo_sensitivity = factor(GABRA %in% c("GABRA1","GABRA2","GABRA3","GABRA5"),
                                          levels=c("TRUE","FALSE"),
                                          labels=c("BZD sensitive","BZD insensitive"))
        )
      rc_plot <- ggplot(data=x,aes(x = GABRA,y = R2)) + 
        geom_violin(show.legend = TRUE,alpha = .9,aes(fill=benzo_sensitivity)) +
        scale_fill_manual(values = c("#984ea3", "#ff7f00"),guide = guide_legend(override.aes = list(color=NULL)))+
        ylab(bquote("Spatial relationship" ~ (italic(R)^2)))+
        theme(legend.title = element_blank(),axis.title.x = element_blank(),axis.text.x = element_text(angle = 45,hjust=1))
      
      rc_plot <- rc_plot + geom_point(data=corr_df,aes(x=as.factor(GABRA),y=R2,color="black"),size=4,show.legend = TRUE) +
        scale_color_manual(values = c("black"),labels = c("Observed")) 
      print(rc_plot)
      ggsave(filename = "figs/GABRA.svg",plot = rc_plot,device = "svg",width = 7,height = 4.5,units = "in")
      
    }
    cat("GABRA gene analysis complete \n")
  }

  if (transmodality == T){
    trans_labels <- read.csv('/cbica/projects/alpraz_EI/input/atlases/Schaef400_transmodality7.csv',header = F)
    
    colnames(trans_labels) <- c('name','transmodality')
    trans_labels$name <- gsub(trans_labels$name,pattern = "7Networks_",replacement = "")
    trans_labels$nums = as.numeric(rownames(trans_labels))

    mat_labels <- mat_labels %>%
      left_join(trans_labels,by = "nums")
    
    rval_df_trans <- rval_df %>% select(-nums)%>%
      left_join(mat_labels,by = c("names" = "node_names"),keep=TRUE)

    #Brainsmash
    rval_df_trans$SVM_weights <- rval_df_trans$rvals
    rval_df_trans$absolute_SVM_weights <- rval_df_trans$rvals_abs
    brainsmash(rval_df_trans,c("transmodality","SVM_weights"))
    brainsmash(rval_df_trans,c("transmodality","absolute_SVM_weights"))
    
    colnames(newCC)=mat_labels$node_names
    rownames(newCC)=mat_labels$node_names
    
    r = cor.test(rval_df_trans$transmodality,rval_df_trans$rvals)
    trans_R2 <- r$estimate^2
    dotplot <- ggplot(data=rval_df_trans,
                      aes(x=transmodality,
                          y=rvals,
                          color = community_name)) +
      geom_point(color="black")+
      geom_smooth(method = "lm",aes(x=transmodality,y=rvals),color="black") +
      annotate(x=mean(rval_df_trans$transmodality,na.rm=T),y=mean(rval_df_trans$rvals),geom = "text",label=sprintf("r = %1.2f,p = %1.3f",r$estimate,r$p.value))+
      ylab("Nodal SVM weight")
    tlegend <- get_legend(dotplot +theme(legend.position = "bottom",legend.title = element_blank(),legend.background = element_blank()))
    dotplot <- dotplot + theme(legend.position = "none")
    print(dotplot)
    ggsave(file='figs/transmodality_svmW_plot.svg',plot=dotplot,device='svg',height=4.5,width=5,units = "in")
    
    brainsmash_distribution = read.table('surrogate_brainmap_corrs_schaefer400x7_aal_Transmodality_map.txt',col.names = "brainSMASH")%>%
      mutate(R2 = brainSMASH^2, varname = "Transmodality")
    brainsmash_plot <- ggplot(data=brainsmash_distribution,aes(x =varname,y = R2)) + 
      geom_violin(show.legend = TRUE,alpha = .9,fill="grey50") +
      ylab(bquote("Spatial relationship" ~ (italic(R)^2)))+
      theme(legend.title = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank()) +
      geom_point(aes(x=varname,y=trans_R2),color="black",size=4,show.legend = FALSE) 
    print(brainsmash_plot)
    ggsave(filename = "figs/brainsmash_plot.svg",plot = brainsmash_plot,device = "svg",width = 2,height = 3.5,units = "in")

  }
      
    return_obj <- list(feat_mat = newCC)

  return(return_obj)
}

mat2brain <- function(CC,atlasname,filename){
  
  #load atlas
  fname <- sprintf("/cbica/projects/alpraz_EI/input/atlases/%s_threshold_0.95.nii.gz",atlasname)
  atlas <- readnii(fname)
  atlas_vec <- c(atlas)
  
  # Get labels
  fname <- sprintf("/cbica/projects/alpraz_EI/data/TASK_GSR/13783/1986/%s_CC_GSR_000.netcc",atlasname)
  template_CC <- as.matrix(read.table(fname,header = T,skip = 4))
  labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_labels.csv',atlasname),header=F)
  if (atlasname %in% c("schaefer200x7_aal","schaefer400x7_aal","gordon333_aal")) {
    community_labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_node_labels.csv',
                                         gsub(x = atlasname,pattern = "_aal",replacement = ""))
                                 ,header=T)
    

    mat_labels <- data.frame(nums = colnames(template_CC))
    mat_labels$nums <- as.numeric(gsub(mat_labels$nums,pattern="X",replace=""))
    mat_labels <- mat_labels %>% left_join(labels,by = c("nums"="V1"))
    mat_labels <- mat_labels %>% left_join(community_labels,by = c("nums"="index"))
    mat_labels <- mat_labels %>%
      mutate(community_name = as.character(community_name),V2 = as.character(V2))
    # mat_labels$community_name[is.na(mat_labels$community_name)]=mat_labels$V2[is.na(mat_labels$community_name)]
    mat_labels$community_name[is.na(mat_labels$community_name)]="AAL_subcortex"
    
    # Calculate values to map on brain
    vals <-rowSums(CC)
    vals <- data.frame(vals = vals,names = names(vals))
    ## Get key column
    idx <- sapply(mat_labels, function(x) sum(x%in%vals$names)/length(vals$names))
    matching_col<-colnames(mat_labels)[idx==T][1] # get first match if multiple
    mat_labels$names <- mat_labels[,matching_col]
    label_df <- mat_labels %>% 
      left_join(vals, by = "names")
    
    # Map values on brain
    voxel_df <- as.data.frame(atlas_vec) %>% left_join(label_df,by = c("atlas_vec"="nums"))
    label_vec <- as.matrix(voxel_df$vals)
    label_matrix <- array(label_vec, dim = dim(atlas))
    nifti_image <- niftiarr(atlas,label_matrix)
    cat(sprintf("\n Writing /cbica/projects/alpraz_EI/output/drug_classification/%s.nii.gz ......\n",filename))
    writenii(nifti_image,filename = sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s.nii.gz",filename))
    cat("done\n")
    
  } else {
    mat_labels <- data.frame(nums = colnames(template_CC))
    mat_labels$nums <- as.numeric(gsub(mat_labels$nums,pattern="X",replace=""))
    mat_labels <- mat_labels %>% left_join(labels,by = c("nums"="V1"))
    colnames(newCC)=mat_labels$V2
    rownames(newCC)=mat_labels$V2
  }
  
}

write_brainsmash_files <- function(df,distance_matrix,idx,fname,var_names,atlasname){
  this_distance_matrix <- distance_matrix[idx,idx]
  write.table(this_distance_matrix,file = sprintf("/cbica/projects/alpraz_EI/output/brainsmash/input_files/%s_threshold_0.95_%s_distMat.csv",atlasname,fname),
              sep = ",",row.names = F,col.names = F)
  brain_map1 <- df[idx,var_names[1]]
  write.table(brain_map1,file = sprintf("/cbica/projects/alpraz_EI/output/brainsmash/input_files/%s_threshold_0.95_%s_%s.csv",atlasname,fname,var_names[1]),
              sep = ",",row.names = F,col.names = F)
  brain_map2 <- df[idx,var_names[2]]
  write.table(brain_map2,file = sprintf("/cbica/projects/alpraz_EI/output/brainsmash/input_files/%s_threshold_0.95_%s_%s.csv",atlasname,fname,var_names[2]),
              sep = ",",row.names = F,col.names = F)
  
}

brainsmash <- function(df,var_names,atlasname="schaefer400x7_aal"){
  distance_matrix <- read.table(sprintf("/cbica/projects/alpraz_EI/input/atlases/%s_threshold_0.95_distMat.txt",atlasname))
  distance_matrix_expanded <- squareform(distance_matrix$V1)
  
  distance_matrix_labels <- read.table(sprintf("/cbica/projects/alpraz_EI/input/atlases/%s_threshold_0.95_roi_order.csv",atlasname),
                                       header = T,sep = ",",row.names = 1)
  colnames(distance_matrix_labels) <- "labels"
  labeled_brain_map <- left_join(distance_matrix_labels,df,by = c("labels"="nums"),keep=TRUE) 
  
  if (atlasname=="schaefer400x7_aal"){
    left_parcel_numbers <- labeled_brain_map$nums[grepl(x=labeled_brain_map$node_names,pattern = "LH")]
    right_parcel_numbers <- labeled_brain_map$nums[grepl(x=labeled_brain_map$node_names,pattern = "RH")]
  } else {
    left_parcel_numbers <- labeled_brain_map$nums[grepl(x=labeled_brain_map$node_names,pattern = "Left")]
    right_parcel_numbers <- labeled_brain_map$nums[grepl(x=labeled_brain_map$node_names,pattern = "Right")]
  }

  # Left cortex
  
  left_parcel_idx <- match(left_parcel_numbers,labeled_brain_map$labels)#this step is probably not necessary but is for abundance of caution for indexing.
  write_brainsmash_files(df = labeled_brain_map,distance_matrix = distance_matrix_expanded,idx = left_parcel_idx,fname = "left",var_names = var_names, atlasname = atlasname)
  
  # Right cortex
  
  right_parcel_idx <- match(right_parcel_numbers,labeled_brain_map$labels)#this step is probably not necessary but is for abundance of caution for indexing.
  write_brainsmash_files(df = labeled_brain_map,distance_matrix = distance_matrix_expanded,idx = right_parcel_idx,fname = "right",var_names = var_names, atlasname = atlasname)
  
  # All cortex
  cortex_parcel_numbers <- labeled_brain_map$nums[labeled_brain_map$nums<1000]
  cortex_parcel_idx <- match(cortex_parcel_numbers,labeled_brain_map$labels)#this step is probably not necessary but is for abundance of caution for indexing.
  
  write_brainsmash_files(df = labeled_brain_map,distance_matrix = distance_matrix_expanded,idx = cortex_parcel_idx,fname = "cortex",var_names = var_names, atlasname = atlasname)
}

ROC_curve <- function(DecisionValues, labels,perm_auc_distribution = NULL){
  # Decision values is nx1
  # Labels is nx1
  # perm_auc_distribution = permutation auc distribution. This can be added if we want to compare our AUC to a null AUC distribution.
  
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
  
  # LargestAccuracy = max(Accuracy_Array)
  ROC_output <- data.frame(TPR =TP_Array,FPR=FP_Array,Accuracy = Accuracy_Array)
  ROC_output <- ROC_output%>%arrange(TPR,FPR)
  
  #AUC
  dFPR <- c(0,diff(ROC_output$FPR))
  dTPR <- c(0,diff(ROC_output$TPR))
  AUC <- sum(ROC_output$TPR * dFPR) + sum(dTPR * dFPR)/2
  
  # Plot 
  # MaxAccuracyIdx <- which.max(ROC_output$Accuracy)
  roc_plot <- ggplot(ROC_output,aes(x=FPR,y=TPR)) +
    geom_point() + geom_line() +
    geom_abline(slope = 1, intercept = 0,linetype = "dashed") +
    # geom_point(x=ROC_output$FPR[MaxAccuracyIdx],y=ROC_output$TPR[MaxAccuracyIdx],size=3,aes(color="red"))+
    # scale_color_manual(values="red",labels=sprintf("Max Accuracy = %1.3f",LargestAccuracy))+
    # theme(legend.title = element_blank(),legend.position = c(.5,.25),legend.background = element_blank(),legend.justification = c("left","top"))+
    annotate("text",x=.58,y=.2,label = sprintf("AUC = %1.3f",AUC),hjust=0,size= theme_get()$text[["size"]]/4) +
    xlab('False positive rate (1-specificity)')+ylab("True positive rate (sensitivity)") +
    theme(axis.line.y.right = element_line(),axis.line.x.top = element_line())
  # ggsave(file='figs/ROC_plot.svg',plot=roc_plot,device='svg',height=4,width=4.5,units = "in")
  print(roc_plot)
  
  return(list(auc=AUC,roc_plot = roc_plot))
}

display_results <- function(atlasname,GSR="GSR",classifier="svm",perm_results=F,result_fname=NULL,results=NULL,subdivision="all"){
  cat(sprintf("\n## Displaying classification results for %s atlas (classifier = %s)\n",atlasname,classifier))

  if (is.null(results)) {
    if (!is.null(result_fname)) {
      results <- readRDS(result_fname)
    } else {
      if (perm_results == T){
        results <- readRDS(sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_all_%s_1_permute_on_FD5results.rds",
                                   atlasname,GSR,classifier))
      }else {
        results <- readRDS(sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_all_%s_1_permute_off_FD5results.rds",
                                   atlasname,GSR,classifier))
      }
    }
  }
  
  b <- results[[1]]
  cat("Classifier Results:")
  print(b)
  cat(sprintf('\naccuracy = %1.3f\n',b$accuracy))

  # Get ROC curve and AUC
  pred_data <- b$pred_data[[1]] #Grabbing first one just to get the drug labels, subid, sesid.
  dec_folds <- data.table::rbindlist(b$pred_data,idcol = "fold")%>%
    pivot_wider(names_from = "fold","values_from"="decisionValues",id_cols = c("subid","sesid")) %>% 
    group_by(subid,sesid)
  pred_data$mean_dec_vals <- rowMeans(dec_folds[,3:dim(dec_folds)[2]])
  AUC_obj <- ROC_curve(pred_data$mean_dec_vals,pred_data$drug)
  AUC <- AUC_obj$auc
  roc_plot <- AUC_obj$roc_plot
  ggsave(file=sprintf('figs/ROC_plot_%s.svg',subdivision),plot=roc_plot,device='svg',height=4,width=4.5,units = "in")
  cat(sprintf("\nAUC = %1.4f\n",AUC))
  
  # Look at W coefs
  W <- results[[2]]
  if (atlasname=="schaefer400x7_aal") {
    feat_mat_obj <- feat2mat(W,atlasname,community_summary = T,GABA=T,transmodality = T,surface = T)
    Wmap <- feat_mat_obj$feat_mat
    mat2brain(CC = abs(Wmap),atlasname = atlasname, filename = sprintf("%s_%s_%s_abs_weights",atlasname,GSR,classifier))
  } else{
    # If not, we don't need all the plotting of the feature weights
    feat_mat_obj <- feat2mat(W,atlasname,community_summary = F,GABA=F,transmodality = F,plots=F)
  }

  if ("perm_p" %in% names(b)) {
    cat(sprintf("\n Number of features: %d\n",results[[3]]))
    cat(sprintf("Accuracy = %1.4f",b$accuracy))
    cat(sprintf("Permutation p = %1.3f\n",b$perm_p))
    perm_acc_plot <- ggplot(data = data.frame(perm_acc_distribution=b$perm_accs),aes(x = perm_acc_distribution)) +
      geom_histogram()+geom_vline(xintercept = b$accuracy) + 
      geom_label(x = b$accuracy,y = 100,label=paste0("Observed\n p = ",as.character(round(b$perm_p,digits = 4))))+
      xlab("Classification Accuracy")+ylab("Number of Draws")+
      ggtitle("Permutation Test")+coord_cartesian(clip = "off")
    
    perm_auc_p <-sum(b$perm_aucs>AUC)/length(b$perm_aucs)

    perm_auc_plot <- ggplot(data = data.frame(perm_auc_distribution=b$perm_aucs),aes(x = perm_auc_distribution)) +
      geom_histogram()+geom_vline(xintercept = AUC) + 
      geom_label(x = AUC,y = 10,label=paste0("Observed\n p = ",perm_auc_p))+
      xlab("AUC")+ylab("Number of Draws")+
      ggtitle("Permutation Test")+coord_cartesian(clip = "off")
    print(perm_auc_plot)

    ggsave(file=sprintf('figs/perm_acc_plot_%s.svg',subdivision),plot=perm_acc_plot,device='svg',height=2.5,width=3,units = "in")
    ggsave(file=sprintf('figs/perm_auc_plot_%s.svg',subdivision),plot=perm_auc_plot,device='svg',height=2.5,width=3,units = "in")

  }
}


### function to extract derivative, confidence interval, significance, and plot for GAMs ###

get_derivs_and_plot <- function(modobj,smooth_var,low_color=NULL,hi_color=NULL){
  this_font_size = font_size
  if (is.null(low_color)){low_color = "white"}
  if (is.null(hi_color)){hi_color = "grey20"}
  derv<-derivatives(modobj,term=smooth_var)
  derv<- derv %>%
    mutate(sig = !(0 >lower & 0 < upper))
  derv$sig_deriv = derv$derivative*derv$sig
  cat(sprintf("\nSig change: %1.2f - %1.2f\n",min(derv$data[derv$sig==T]),max(derv$data[derv$sig==T])))
  d1<- ggplot(data=derv) + geom_tile(aes(x = data, y = .5, fill = sig_deriv))
  
  # Set the gradient colors
  if (min(derv$derivative)>0) {
    d1 <- d1 + scale_fill_gradient(low = low_color, high = hi_color,limits = c(0,max(derv$derivative)))
    # If you want all plots to have the same scaling, this code can be used instead-- This is desirable if you have a factor-smooth model.
    ## max_val = .5
    ## scale_fill_gradient(low = low_color,high = hi_color,limits = c(0,max_val),oob = squish)
  } else if (min(derv$derivative)<0 & max(derv$derivative)<0) {
    d1 <- d1 + scale_fill_gradient(low = hi_color, high = low_color,limits = c(min(derv$derivative),0))
  }else {
    d1 <- d1 +
      scale_fill_gradient2(low = "steelblue", midpoint = 0, mid = "white",
                           high = "firebrick",limits = c(min(derv$derivative),max(derv$derivative)))
  }

  d1 <- d1 + 
    labs(x = smooth_var,fill = sprintf("\u0394%s",smooth_var)) + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = this_font_size),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=this_font_size),
          legend.text = element_text(size = this_font_size),
          axis.title = element_text(size = this_font_size),
          legend.key.width = unit(1,"cm"),
          legend.position = "bottom",
          legend.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    guides(fill = guide_colorbar(reverse = F,direction = "horizontal",title.position = "left")) +
    geom_rect(aes(ymin=0,ymax=1,xmin=min(data),xmax=max(data)),color="black",fill="white",alpha = 0)
  return(d1)
}

# Func to visualize GAM model outputs
visualize_model <- function(modobj,smooth_var, int_var = NULL ,group_var = NULL, plabels = NULL,check_diagnostics = F,derivative_plot = F,show.data=TRUE, line_color = "black"){
  this_font_size = font_size
  if (any(class(modobj)=="gam")) {
    model <- modobj
  } else if (class(modobj$gam)=="gam") {
    model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  s<-summary(model)
  df <- model$model
  ## Generate custom line plot
  np <- 10000 #number of predicted values

  theseVars <- attr(model$terms,"term.labels")
  varClasses <- attr(model$terms,"dataClasses")
  thisResp <- as.character(model$terms[[2]])
  
  if (!is.null(int_var)) {
    # We will produce and interaction plot
    if (!any(grepl(x=as.character(model$formula),pattern = int_var))) {
      warning("int_var not recognized in model formula!")
      return()
    }
    switch (varClasses[int_var],
            "numeric" = {
              q <- quantile(df[,int_var],probs = c(.05,.95)) #pick 10% and 90% to plot
              bigq <- q[[2]]
              smallq <- q[[1]]
              values <- c(bigq,smallq)
              labs <- c(sprintf("high (%1.2f)",bigq),sprintf("low (%1.2f)",smallq))
              
              q <-quantile(rescale(df[,int_var],c(0,1)),probs = c(0,.5,1))
              limit_values <- c(q[[1]],q[[length(q)]])
              midpoint_val <- unname(q[[2]])
              cbar_vals <- unname(q)
              
              theseLabs = rep(values,each = np)
              grad_fill = T
            },
            "factor" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = rep(values,each = np)
              grad_fill = F
            },
            "ordered" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = ordered(rep(values,each = np),levels = values)
              grad_fill = F
            }
    )
    labPred <- data.frame(init = rep(0,np*length(labs)))
    labPred[,int_var] = theseLabs
    labPred$lab = rep(labs,each = np)
    labPred <- labPred[,names(labPred) !="init"]
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else if (v == int_var) {
        next
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    
    thisPred <- thisPred %>% select(-init)
    thisPred <- do.call("rbind", replicate(length(labs), thisPred, simplify = FALSE))
    
    pred <- cbind(labPred,thisPred)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    if (!is.null(group_var)) {
      pred[,group_var] = NA #these columns have to exist in the dataframe for plotting
    }
    pred[,thisResp] = 1 #these columns have to exist in the dataframe for plotting
    
    low_color = "#91bfdb"
    high_color = "#fc8d59"
    high_line = "#f46d43"
    low_line = "#4575b4"
    
    if (grad_fill == T) {
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var)) +
        geom_point(alpha = 0.65,stroke = 0, size = point_size) 
      if (!is.null(group_var)) {
        cat("adding lines")
        p1<- p1 + geom_line(aes_string(group = group_var),alpha = .5)
      }
      p1 <- p1 +
        scale_color_gradientn(colors = c(low_line,"grey90",high_line), values = cbar_vals,name = "") +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = "lab"),alpha = .3, linetype = 0) +
        scale_fill_manual(values = c(high_color,low_color)) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",group = "lab"),size = line_size) +
        labs(title = plabels)
    } else {
      black_color = "#1c405eff"
      green_color = "#285c27ff"
      
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var))
      if (show.data==TRUE) {
        print('show data is on')
        p1<- p1 +  
          geom_point(alpha = .35,stroke = 0, size = point_size,show.legend = FALSE)
          # geom_hex(color=NULL) 
      } 
      
      if (!is.null(group_var)) {
        p1<- p1 + geom_line(aes_string(group = group_var),alpha = .3,color="black")
      } 
      p1 <- p1 +
        # scale_color_brewer(type = "qual",palette = "Set1",direction = 1) +
        scale_color_manual(values = c(black_color,green_color)) +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = int_var),alpha = .5, linetype = 0,show.legend=FALSE) +
        # scale_fill_brewer(type = "qual",palette = "Set1",direction = 1) +
        scale_fill_manual(values = c(black_color,green_color)) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",color = int_var),size = line_size,show.legend = FALSE) +
        labs(title = plabels)
    }
  } else {
    
    # No interaction variable, just produce a single line plot
    int_var = "" # This is a lazy solution to making the existing code workable with no int_var.
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    pred <- thisPred %>% select(-init)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    pred[,group_var] = NA
    pred[,thisResp] = 1
    
    df <- df %>%
      gratia::add_partial_residuals(model)
    df$partial_resids <- unlist(df[,grep(x=names(df),pattern = "s(",fixed = T)])

    p1 <- ggplot(data = df, aes_string(x = smooth_var,y = "partial_resids"))
    if (show.data==TRUE) {
      p1<- p1 +  
        geom_point(alpha = .3,stroke = 0, size = point_size,color = line_color)
        # geom_hex(show.legend = TRUE) + scale_fill_gradient(low="white",high=line_color,limits = c(1, 9), oob = scales::squish)
    } 
    if (!is.null(group_var)) {
      cat("adding lines")
      p1<- p1 + geom_line(aes_string(group = group_var),alpha = .5)
    }
    p1 <- p1 + geom_ribbon(data = pred,aes_string(x = smooth_var ,y=thisResp, ymin = "selo",ymax = "sehi"),fill = line_color, alpha = .5, linetype = 0) +
      geom_line(data = pred,aes_string(x = smooth_var, y = "fit"),size = line_size,color=line_color) +
      labs(title = plabels) #+ ylim(NA,30)
  }
  
  if (derivative_plot == T) {
    # We will add a bar that shows where the derivative is significant.
    # First make some adjustments to the line plot.
    p1<- p1+theme(text = element_text(size=this_font_size),
                  axis.text = element_text(size = this_font_size),
                  axis.title.y = element_text(size = this_font_size),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.text = element_text(size = this_font_size),
                  legend.title = element_text(size = this_font_size),
                  axis.title = element_text(size = this_font_size),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "transparent",colour = NA),
                  plot.background = element_rect(fill = "transparent",colour = NA),
                  plot.margin = unit(c(.2, .2, 0, .2), "cm")) #Top, left,Bottom, right
    scatter <- list(p1)
    
    # Now add the plots using the derivative plotting function
    if (any(grepl(x = row.names(s$s.table),pattern =  ":") & grepl(x=row.names(s$s.table),pattern = int_var))) {
      # Factor levels separately if there is an interaction in the model.
      f<-formula(model) # current formula
      fterms <- terms(f)
      fac <- attr(fterms, "factors")
      idx <- which(as.logical(colSums(fac[grep(x=row.names(fac),pattern = int_var),])))
      new_terms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
      new_formula <- formula(new_terms) # Formula without any interaction terms in the model.
      
      #add derivative gradients for each level of the factor
      num_levels <- length(levels(df[,int_var]))
      level_colors <- suppressWarnings(RColorBrewer::brewer.pal(num_levels,"Set1")) #use the same palette as the line plot
      plotlist = vector(mode = "list",length = num_levels+1) # we will be building a list of plots
      plotlist[1] = scatter # first the scatter plot
      level_colors=c(black_color,green_color)
      for (fcount in 1:num_levels) {
        this_level <- levels(df[,int_var])[fcount]
        df$subset <- df[,int_var] == this_level
        this_mod <- gam(formula = new_formula,data = df,subset = subset)
        # this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        
        if (fcount != num_levels & fcount != 1){
          # get rid of redundant junk
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
          this_d$theme$legend.text = element_blank()
          this_d$theme$legend.position="none"
        }
        if (fcount == 1) {
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.text = element_blank()
          this_d$theme$legend.position="none"
        }
        if (fcount == num_levels) {
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
          legend <- get_legend(this_d)
        }
        this_d$labels$fill=NULL
        plotlist[fcount+1] = list(this_d+theme(legend.position = "none"))
      }

      pg<-cowplot::plot_grid(rel_heights = c(16*num_levels,rep(num_levels,num_levels-1),2*num_levels),plotlist = plotlist,align = "v",axis = "lr",ncol = 1)
      final_plot <- cowplot::plot_grid(pg,legend,rel_heights = c(1,.15),ncol = 1)
      # final_plot <- pg
      print(final_plot)
    } else {
      # No need to split
      d1 <- get_derivs_and_plot(modobj = modobj,smooth_var = smooth_var)
      scatter <- list(p1)
      bar <- list(d1+theme(legend.position = "none"))
      legend <- get_legend(d1+theme(legend.position = "bottom",legend.direction = "horizontal"))
      allplots <- c(scatter,bar)
      pg<-cowplot::plot_grid(rel_heights = c(1,.35),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
      final_plot <- cowplot::plot_grid(pg,legend,rel_heights = c(1,.15),ncol = 1)
      print(final_plot)
      # final_plot=pg
      
      # print(final_plot)
    }
    
  }    else {
    # No derivative plot
    p1<- p1+theme(text = element_text(size=font_size),
                  axis.text = element_text(size = font_size),
                  legend.text = element_text(size = font_size),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank())
    final_plot <- p1
    print(final_plot)
  }
  
  if (check_diagnostics == T) {
    cp <- check(b,
                a.qq = list(method = "tnorm",
                            a.cipoly = list(fill = "light blue")),
                a.respoi = list(size = 0.5),
                a.hist = list(bins = 10))
    print(cp)
  }
  return(final_plot)
}
