library(mlr3)
library(mlr3learners)
library(dplyr)   
library(mlr3tuning)
library(paradox)

all.combined.df <- all.combined.log.list$four.gp

## teminator setting
evals10 <- trm("evals", n_evals = 10)
## internal sampling setting
resamp5 = rsmp("cv", folds = 5)
tuner <- tnr("grid_search", resolution = 10)

generate_design_grid(tune_ps, 10)

multiclass.loso.tuning <- function(a = data, 
                            metadata = metadata, 
                            project = pro,
                            learner.method = "classif.ranger",
                            measure.method = "classif.acc",
                            evals = evals10,
                            resamp = resamp5, 
                            tuner = tuner
                            #tune_ps = tune_ps,
                            ){
    library(mlr3)
    library(mlr3learners)
    library(dplyr)   
    library(mlr3tuning)
    library(paradox)
                         
    ## construct learners
    train_sample <- metadata %>% filter(Project %in% project) %>% pull(Sample_ID)
    ## feat.data except the specified project
    a.2 <- a[!rownames(a) %in% train_sample, ]

    mtry.def <- sqrt(ncol(a.2) - 1)

    tune_ps <- ParamSet$new(list(
    ParamInt$new("num.trees", lower = 500, upper = 1000, default = 500), 
    ParamInt$new("mtry", lower = round(mtry.def/2), upper = round(mtry.def*2), default = round(mtry.def))))
    tune_ps

    ## normalization
    train.data.mean <- colMeans(a.2[,! names(a.2) %in% "Group"])
    train.data.sd <- apply(a.2[,! names(a.2) %in% "Group"], 2, sd)
    train.data.quantile <- quantile(train.data.sd, 0.1, names = F)

    train.data.norm <- a.2
    for(i in names(train.data.mean)){
        train.data.norm[,i] <- (a.2[,i]-train.data.mean[i])/(train.data.sd[i] + train.data.quantile)
    }

    ## frozen normalization
    test.data <- a[rownames(a) %in% train_sample, ]
    test.data.norm <- test.data
    for(i in names(train.data.mean)){
       test.data.norm[,i] <- (test.data[,i]-train.data.mean[i])/(train.data.sd[i] + train.data.quantile)
    }

    task_fourg_test <- TaskClassif$new(id = "test", 
                                       backend = test.data.norm, 
                                       target = "Group")

    ## custom sampling
    train_ind <- caret::createMultiFolds(train.data.norm$Group, k = 10, times = 10)
    test_ind <- lapply(train_ind, function(x) setdiff(1:nrow(train.data.norm), x))             

    ## tuning list
    opt.param <- list()

    ## internal validation 
    models.list <- list()
    ## split data
    prediction.internal.list <- list()
    confusion.internal.list <- list()
    acc.internal.list <- list()

    ## 
    for(i in 1:length(train_ind)){
    #train_sample <- subset(metadata, Project %in% pro)$Sample_ID
    train_set <- train_ind[[i]]
    test_set <- test_ind[[i]]

     ## construct task for training data and testing data from nesting resampling
    a.3 <- train.data.norm[train_set, ]

    task_inner_train <- TaskClassif$new(id = "fourgp_train", ## arbitrary identifier
                                        backend = a.3, ## data
                                        target="Group")
    task_inner_test <- TaskClassif$new(id = "fourgp_test", ## arbitrary identifier
                                        backend = train.data.norm[test_set,], ## data
                                        target="Group")

    # tuning the parameters
    models.list[[i]] <- lrn(learner.method)

    models.list[[i]]$predict_type = "prob"
    instance <- TuningInstanceSingleCrit$new(task = task_inner_train,
                                             learner = models.list[[i]],
                                             resampling = resamp5,
                                             measure = msr(measure.method),
                                             search_space = tune_ps,
                                             terminator = evals10)
		tuner$optimize(instance)
		
    opt.param[[i]] <- instance$result_learner_param_vals
    models.list[[i]]$param_set$values <- instance$result_learner_param_vals
		
    models.list[[i]]$train(task_inner_train)

    # predict data
    prediction <- models.list[[i]]$predict(task_inner_test)
    prediction.internal.list[[i]] <- prediction
    # calculate performance

    confusion.internal.list[[i]] <- prediction$confusion
    acc.internal.list[[i]] <- prediction$score(msr(measure.method))
    }

    ## external validation 
    prediction.external.list <- list()
    confusion.external.list <- list()
    acc.external.list <- list()
    
    for(i in 1:length(train_ind)){
    #train_sample <- subset(metadata, Project %in% pro)$Sample_ID
    #train_set <- train_ind[[i]]

    # predict data
    prediction <- models.list[[i]]$predict(task_fourg_test)
    prediction.external.list[[i]] <- prediction
    # calculate performance
    confusion.external.list[[i]] <- prediction$confusion
    acc.external.list[[i]] <- prediction$score(msr(measure.method))
    }

    return(list(train_ind = train_ind, 
                models.list = models.list, 
                opt.param = opt.param, 

                prediction.internal.list = prediction.internal.list, 
                acc.internal.list = acc.internal.list, 
                confusion.internal.list = confusion.internal.list,
                prediction.external.list = prediction.external.list, 
                confusion.external.list = confusion.external.list, 
                acc.external.list = acc.external.list))
}

## multiclass for Dif- features
multiclass.loso.result <- list()

data <- all.combined.df[, names(all.combined.df) %in% names(kw.combined.list$four.gp)]
data$Group <- factor(data$Group, levels = c("CTR", "UC", "CD", "CRC"))

for(i in unique(metadata$Project)){
  multiclass.loso.result[[i]] <- multiclass.loso.tuning(a = data, 
                                                        metadata = metadata, 
                                                        project = i,
                                                        learner.method = "classif.ranger",
                                                        measure.method = "classif.acc",
                                                        evals = evals10,
                                                        resamp = resamp5, 
                                                        tuner = tuner)
}

multiclass.loso.result.kw <- multiclass.loso.result
save(multiclass.loso.result.kw, file = "mlr3_combined_LOSO_kw.RData") 

## Binary models
data <- all.combined.df[, names(all.combined.df) %in% names(kw.combined.list$four.gp)]
data$Group <- as.character(data$Group)
data[-which(data$Group %in% "CTR"), "Group"] <- "Case"
data$Group <- factor(data$Group, levels = c("CTR","Case"))

binary.loso.result <- list()

for(i in unique(metadata$Project)){
  binary.loso.result[[i]] <- multiclass.loso.tuning(a = data, 
                                                        metadata = metadata, 
                                                        project = i,
                                                        learner.method = "classif.ranger",
                                                        measure.method = "classif.auc",
                                                        evals = evals10,
                                                        resamp = resamp5, 
                                                        tuner = tuner)
}
binary.loso.result.kw <- binary.loso.result
save(binary.loso.result.kw, file = "mlr3_combined_LOSO_binary_kw.RData") 


## Three-classes models
data <- all.combined.df[, names(all.combined.df) %in% names(kw.combined.list$four.gp)]

data$Group <- as.character(data$Group)
data <- data[-which(data$Group %in% "CTR"), ]
data$Group <- factor(data$Group, levels = c("UC", "CD", "CRC"))

three.loso.result <- list()

for(i in unique(metadata$Project)){
  three.loso.result[[i]] <- multiclass.loso.tuning(a = data, 
                                                        metadata = metadata, 
                                                        project = i,
                                                        learner.method = "classif.ranger",
                                                        measure.method = "classif.acc",
                                                        evals = evals10,
                                                        resamp = resamp5, 
                                                        tuner = tuner)
}

three.loso.result.kw <- three.loso.result
save(three.loso.result.kw, file =  "mlr3_combined_LOSO_kw_three.RData")

## CTR IBD and CRC groups
## Three-classes models

data <- all.combined.df
data <- all.combined.df[, names(all.combined.df) %in% names(kw.combined.list$four.gp)]

data$Group <- as.character(data$Group)
data[-which(data$Group %in% c("CTR", "CRC")), "Group"] <- "IBD"
data$Group <- factor(data$Group, levels = c("CTR", "IBD", "CRC"))


three.all.loso.result <- list()
for(i in unique(metadata$Project)){
  three.all.loso.result[[i]] <- multiclass.loso.tuning(a = data, 
                                                        metadata = metadata, 
                                                        project = i,
                                                        learner.method = "classif.ranger",
                                                        measure.method = "classif.acc",
                                                        evals = evals10,
                                                        resamp = resamp5, 
                                                        tuner = tuner)
}

three.all.loso.result.kw <- three.all.loso.result
save(three.all.loso.result.kw, file = "mlr3_combined_LOSO_three_all_kw.RData")


## testing results for LODO
## for examples, the lodo results of four-class model based on combined profile
multiclass.loso.prob.list <- lapply(multiclass.loso.result, 
                               function(data){
                                 y.list <- list()
                                 for (i in 1:100) {
                                   y.list[[i]] <- data$prediction.external.list[[i]];
                                   y.list[[i]] <- as.data.table(y.list[[i]]);
                                   y.list[[i]] <- data.frame(y.list[[i]])
                                 }
                                 y.df <- y.list %>% reduce(rbind) %>% 
                                   group_by(row_id) %>%
                                   mutate(mean_CTR = mean(prob.CTR),
                                          mean_UC = mean(prob.UC),
                                          mean_CD = mean(prob.CD),
                                          mean_CRC = mean(prob.CRC));
                                 y.df <- y.df %>%
                                   distinct(row_id, truth, .keep_all = T)
                                 return(y.df)
                               })
for (study in names(multiclass.loso.prob.list)) {
  multiclass.loso.prob.list[[study]]$Project <- study
  multiclass.loso.prob.list[[study]]$mean_response <- c("mean_CTR", "mean_UC", "mean_CD", "mean_CRC")[apply(multiclass.loso.prob.list[[study]][, c("mean_CTR", "mean_UC", "mean_CD", "mean_CRC")], 1, which.max)]
}
multiclass.loso.prob.df <- multiclass.loso.prob.list %>%
  reduce(rbind)

multiclass.loso.prob.df$response <- lapply(strsplit(multiclass.loso.prob.df$mean_response, split = "_", fixed = T), 
                                      function(data){ y <- data[2]}) %>% unlist
multiclass.loso.prob.df$response <- factor(multiclass.loso.prob.df$response, 
                                      levels = c("CTR", "UC", "CD", "CRC"))
multiclass.loso.prob.df$truth <- factor(multiclass.loso.prob.df$truth, 
                                   levels = c("CTR", "UC", "CD", "CRC"))

## trainng results for four-class model 
multiclass.loso.prob.internal.list <- list()
a <- all.combined.df

for (study in names(multiclass.loso.result)) {
  train_sample <- metadata %>% 
    filter(Project %in% study) %>% 
    pull(Sample_ID)
  ## feat.data except the specified project
  a.2 <- a[!rownames(a) %in% train_sample, ]
  
  data <- multiclass.loso.result[[study]]
  
  y.list <- list()
  for (i in 1:100) {
    y.list[[i]] <- data$prediction.internal.list[[i]];
    y.list[[i]] <- as.data.table(y.list[[i]]);
    y.list[[i]] <- data.frame(y.list[[i]])
    y.list[[i]]$Sample_ID <- setdiff(seq(1, nrow(a.2)), data$train_ind[[i]])
  }
  y.df <- y.list %>% reduce(rbind) %>% 
    group_by(Sample_ID) %>%
    mutate(mean_CTR = mean(prob.CTR),
           mean_UC = mean(prob.UC),
           mean_CD = mean(prob.CD),
           mean_CRC = mean(prob.CRC));
  y.df <- y.df %>%
    distinct(Sample_ID, truth, .keep_all = T)
  y.df$Project <- study
  y.df$mean_response <- c("mean_CTR", "mean_UC", "mean_CD", "mean_CRC")[apply(y.df[, c("mean_CTR", "mean_UC", "mean_CD", "mean_CRC")], 1, which.max)]
  
  multiclass.loso.prob.internal.list[[study]] <- y.df
}

multiclass.loso.prob.inter.df <- multiclass.loso.prob.internal.list %>%
  reduce(rbind)

multiclass.loso.prob.inter.df$response <- lapply(strsplit(multiclass.loso.prob.inter.df$mean_response, split = "_", fixed = T), 
                                           function(data){ y <- data[2]}) %>% unlist
multiclass.loso.prob.inter.df$response <- factor(multiclass.loso.prob.inter.df$response, 
                                           levels = c("CTR", "UC", "CD", "CRC"))
multiclass.loso.prob.inter.df$truth <- factor(multiclass.loso.prob.inter.df$truth, 
                                        levels = c("CTR", "UC", "CD", "CRC"))

## three-class models
three.all.loso.prob.list <- lapply(three.all.loso.result, 
                                    function(data){
                                      y.list <- list()
                                      for (i in 1:100) {
                                        y.list[[i]] <- data$prediction.external.list[[i]];
                                        y.list[[i]] <- as.data.table(y.list[[i]]);
                                        y.list[[i]] <- data.frame(y.list[[i]])
                                      }
                                      y.df <- y.list %>% reduce(rbind) %>% 
                                        group_by(row_id) %>%
                                        mutate(mean_CTR = mean(prob.CTR),
                                               mean_IBD = mean(prob.IBD),
                                               mean_CRC = mean(prob.CRC));
                                      y.df <- y.df %>%
                                        distinct(row_id, truth, .keep_all = T)
                                      return(y.df)
                                    })
for (study in names(three.all.loso.prob.list)) {
  three.all.loso.prob.list[[study]]$Project <- study
  three.all.loso.prob.list[[study]]$mean_response <- c("mean_CTR", "mean_IBD", "mean_CRC")[apply(three.all.loso.prob.list[[study]][, c("mean_CTR", "mean_IBD", "mean_CRC")], 1, which.max)]
}
three.all.loso.prob.df <- three.all.loso.prob.list %>%
  reduce(rbind)

three.all.loso.prob.df$response <- lapply(strsplit(three.all.loso.prob.df$mean_response, split = "_", fixed = T), 
                                           function(data){ y <- data[2]}) %>% unlist
three.all.loso.prob.df$response <- factor(three.all.loso.prob.df$response, 
                                           levels = c("CTR", "IBD", "CRC"))
three.all.loso.prob.df$truth <- factor(three.all.loso.prob.df$truth, 
                                        levels = c("CTR", "IBD", "CRC"))

## training results
a <- all.combined.df

a$Group <- as.character(a$Group)
a[-which(a$Group %in% c("CTR", "CRC")), "Group"] <- "IBD"
a$Group <- factor(a$Group, levels = c("CTR", "IBD", "CRC"))

three.all.loso.prob.internal.list <- list()

for (study in names(three.all.loso.result)) {
  train_sample <- metadata %>% 
    filter(Project %in% study) %>% 
    pull(Sample_ID)
  ## feat.data except the specified project
  a.2 <- a[!rownames(a) %in% train_sample, ]
  
  data <- three.all.loso.result[[study]]
  
  y.list <- list()
  for (i in 1:100) {
    y.list[[i]] <- data$prediction.internal.list[[i]];
    y.list[[i]] <- as.data.table(y.list[[i]]);
    y.list[[i]] <- data.frame(y.list[[i]])
    y.list[[i]]$Sample_ID <- setdiff(seq(1, nrow(a.2)), data$train_ind[[i]])
  }
  y.df <- y.list %>% reduce(rbind) %>% 
    group_by(Sample_ID) %>%
    mutate(mean_CTR = mean(prob.CTR),
           mean_IBD = mean(prob.IBD),
           mean_CRC = mean(prob.CRC));
  y.df <- y.df %>%
    distinct(Sample_ID, truth, .keep_all = T)
  y.df$Project <- study
  y.df$mean_response <- c("mean_CTR", "mean_IBD", "mean_CRC")[apply(y.df[, c("mean_CTR", "mean_IBD", "mean_CRC")], 1, which.max)]
  
  three.all.loso.prob.internal.list[[study]] <- y.df
}

three.all.loso.prob.inter.df<- three.all.loso.prob.internal.list %>%
  reduce(rbind)

three.all.loso.prob.inter.df$response <- lapply(strsplit(three.all.loso.prob.inter.df$mean_response, split = "_", fixed = T), 
                                                 function(data){ y <- data[2]}) %>% unlist
three.all.loso.prob.inter.df$response <- factor(three.all.loso.prob.inter.df$response, 
                                                 levels = c("CTR", "IBD", "CRC"))
three.all.loso.prob.inter.df$truth <- factor(three.all.loso.prob.inter.df$truth, 
                                              levels = c("CTR", "IBD", "CRC"))
                                              
# cases model
three.loso.prob.list <- lapply(three.loso.result, 
                            function(data){
                              y.list <- list()
                              for (i in 1:100) {
                                y.list[[i]] <- data$prediction.external.list[[i]];
                                y.list[[i]] <- as.data.table(y.list[[i]]);
                                y.list[[i]] <- data.frame(y.list[[i]])
                              }
                              y.df <- y.list %>% reduce(rbind) %>% 
                                group_by(row_id) %>%
                                mutate(mean_UC = mean(prob.UC),
                                       mean_CD = mean(prob.CD),
                                       mean_CRC = mean(prob.CRC));
                              y.df <- y.df %>%
                                distinct(row_id, truth, .keep_all = T)
                              return(y.df)
                            })
for (study in names(three.loso.prob.list)) {
  three.loso.prob.list[[study]]$Project <- study
  three.loso.prob.list[[study]]$mean_response <- c("mean_UC", "mean_CD", "mean_CRC")[apply(three.loso.prob.list[[study]][, c("mean_UC", "mean_CD", "mean_CRC")], 1, which.max)]
}
three.loso.prob.df <- three.loso.prob.list %>%
  reduce(rbind)
three.loso.prob.df$response <- lapply(strsplit(three.loso.prob.df$mean_response, split = "_", fixed = T), 
                                      function(data){ y <- data[2]}) %>% unlist
three.loso.prob.df$response <- factor(three.loso.prob.df$response, 
                                      levels = c("UC", "CD", "CRC"))
three.loso.prob.df$truth <- factor(three.loso.prob.df$truth, 
                                      levels = c("UC", "CD", "CRC"))

##training 
a <- all.combined.df
a$Group <- as.character(a$Group)
a <- a[-which(a$Group %in% "CTR"), ]
a$Group <- factor(a$Group, levels = c("UC", "CD", "CRC"))

three.loso.prob.internal.list <- list()

for (study in names(three.loso.result)) {
  train_sample <- metadata %>% 
    filter(Project %in% study) %>% 
    pull(Sample_ID)
  ## feat.data except the specified project
  a.2 <- a[!rownames(a) %in% train_sample, ]
  
  data <- three.loso.result[[study]]
  
  y.list <- list()
  for (i in 1:100) {
    y.list[[i]] <- data$prediction.internal.list[[i]];
    y.list[[i]] <- as.data.table(y.list[[i]]);
    y.list[[i]] <- data.frame(y.list[[i]])
    y.list[[i]]$Sample_ID <- setdiff(seq(1, nrow(a.2)), data$train_ind[[i]])
  }
  y.df <- y.list %>% reduce(rbind) %>% 
    group_by(Sample_ID) %>%
    mutate(mean_UC = mean(prob.UC),
           mean_CD = mean(prob.CD),
           mean_CRC = mean(prob.CRC));
  y.df <- y.df %>%
    distinct(Sample_ID, truth, .keep_all = T)
  y.df$Project <- study
  y.df$mean_response <- c("mean_UC", "mean_CD", "mean_CRC")[apply(y.df[, c("mean_UC", "mean_CD", "mean_CRC")], 1, which.max)]
  
  three.loso.prob.internal.list[[study]] <- y.df
}

three.loso.prob.inter.df <- three.loso.prob.internal.list %>%
  reduce(rbind)

three.loso.prob.inter.df$response <- lapply(strsplit(three.loso.prob.inter.df$mean_response, split = "_", fixed = T), 
                                                function(data){ y <- data[2]}) %>% unlist
three.loso.prob.inter.df$response <- factor(three.loso.prob.inter.df$response, 
                                                levels = c("UC", "CD", "CRC"))
three.loso.prob.inter.df$truth <- factor(three.loso.prob.inter.df$truth, 
                                             levels = c("UC", "CD", "CRC"))

## Binary models
binary.loso.prob.list <- lapply(binary.loso.result, 
                                   function(data){
                                     y.list <- list()
                                     for (i in 1:100) {
                                       y.list[[i]] <- data$prediction.external.list[[i]];
                                       y.list[[i]] <- as.data.table(y.list[[i]]);
                                       y.list[[i]] <- data.frame(y.list[[i]])
                                     }
                                     y.df <- y.list %>% reduce(rbind) %>% 
                                       group_by(row_id) %>%
                                       mutate(mean_CTR = mean(prob.CTR),
                                              mean_Case = mean(prob.Case));
                                     y.df <- y.df %>%
                                       distinct(row_id, truth, .keep_all = T)
                                     return(y.df)
                                   })
for (study in names(binary.loso.prob.list)) {
  binary.loso.prob.list[[study]]$Project <- study
  binary.loso.prob.list[[study]]$mean_response <- c("mean_CTR",  "mean_Case")[apply(binary.loso.prob.list[[study]][, c("mean_CTR",  "mean_Case")], 1, which.max)]
}
binary.loso.prob.df <- binary.loso.prob.list %>%
  reduce(rbind)

binary.loso.prob.df$response <- lapply(strsplit(binary.loso.prob.df$mean_response, split = "_", fixed = T), 
                                          function(data){ y <- data[2]}) %>% unlist
binary.loso.prob.df$response <- factor(binary.loso.prob.df$response, 
                                          levels = c("CTR", "Case"))
binary.loso.prob.df$truth <- factor(binary.loso.prob.df$truth, 
                                       levels = c("CTR", "Case"))

## AUC for training 
a <- all.combined.df 
## Binary models 
a$Group <- as.character(a$Group)
a[-which(a$Group %in% "CTR"), "Group"] <- "Case"
a$Group <- factor(a$Group, levels = c("CTR","Case"))

binary.loso.prob.internal.list <- list()

for (study in names(binary.loso.result)) {
  train_sample <- metadata %>% 
    filter(Project %in% study) %>% 
    pull(Sample_ID)
  ## feat.data except the specified project
  a.2 <- a[!rownames(a) %in% train_sample, ]
  
  data <- binary.loso.result[[study]]
  
  y.list <- list()
  for (i in 1:100) {
    y.list[[i]] <- data$prediction.internal.list[[i]];
    y.list[[i]] <- as.data.table(y.list[[i]]);
    y.list[[i]] <- data.frame(y.list[[i]])
    y.list[[i]]$Sample_ID <- setdiff(seq(1, nrow(a.2)), data$train_ind[[i]])
  }
  y.df <- y.list %>% reduce(rbind) %>% 
    group_by(Sample_ID) %>%
    mutate(mean_CTR = mean(prob.CTR), 
           mean_Case = mean(prob.Case));
  y.df <- y.df %>%
    distinct(Sample_ID, truth, .keep_all = T)
  y.df$Project <- study
  y.df$mean_response <- c("mean_CTR", "mean_Case")[apply(y.df[, c("mean_CTR", "mean_Case")], 1, which.max)]
  
  binary.loso.prob.internal.list[[study]] <- y.df
}

binary.loso.prob.inter.df <- binary.loso.prob.internal.list %>%
  reduce(rbind)

binary.loso.prob.inter.df$response <- lapply(strsplit(binary.loso.prob.inter.df$mean_response, split = "_", fixed = T), 
                                            function(data){ y <- data[2]}) %>% unlist
binary.loso.prob.inter.df$response <- factor(binary.loso.prob.inter.df$response, 
                                            levels = c("CTR", "Case"))
binary.loso.prob.inter.df$truth <- factor(binary.loso.prob.inter.df$truth, 
                                         levels = c("CTR", "Case"))

## AUC

binary.loso.prob.inter.auc <- list()
for (study in unique(binary.loso.prob.inter.df$Project)) {
  a <- binary.loso.prob.inter.df %>% 
    filter(Project %in% study) 
  forestpred <- prediction(a$mean_Case, a$truth, 
                           label.ordering = c("CTR", "Case"));
  forestperf <- performance(forestpred, "tpr", "fpr");
  forest.auc <- performance(forestpred, "auc")@y.values[[1]]
  data.roc.df <- data.frame(fpr = forestperf@x.values[[1]], 
                            tpr = forestperf@y.values[[1]],
                            Project = study)
  binary.loso.prob.inter.auc[[study]]$prediction <- forestpred
  binary.loso.prob.inter.auc[[study]]$performance <- forestperf
  binary.loso.prob.inter.auc[[study]]$AUC <- forest.auc
  binary.loso.prob.inter.auc[[study]]$roc.df <- data.roc.df
}
lapply(binary.loso.prob.inter.auc, 
       function(data){ y<- data$AUC}) %>% unlist

## summarize the LODO training results to a dataframe
## for example, summarize the results of binary models

binary.LODO.train.result <- list(MultBinary.all = (lapply(binary.loso.prob.inter.auc, function(data){data$AUC}) %>% unlist %>% data.frame),
                                 MultBinary.kw = (lapply(kw.binary.loso.prob.inter.auc, function(data){data$AUC}) %>% unlist %>% data.frame))
for (i in names(binary.LODO.train.result)) {
  binary.LODO.train.result[[i]]$Data <- i
  binary.LODO.train.result[[i]]$Project <- rownames(binary.LODO.train.result[[i]])
  
}
binary.LODO.train.df <- binary.LODO.train.result %>%
  reduce(rbind)
names(binary.LODO.train.df)[1] <- "AUC"

col.lodo.heatmap <- c("#ffe600","#e25d6e")

binary.LODO.train.result.plot <- list()                               
binary.LODO.train.result.plot[[1]] <- binary.LODO.train.df %>%  
  mutate(Data = factor(Data, levels = c("MultBinary.all", "MultBinary.kw"))) %>%
  mutate(Project = factor(Project, levels = unique(metadata$Project))) %>%
  ggplot(aes(x=Data, y=Project, fill=AUC)) + 
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUC, digits=2)"), col='white', size=6) +
  scale_fill_gradientn(colours = col.lodo.heatmap, limits=c(0.5, 1), 
                       guide=FALSE)  + 
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(vjust = 0, hjust = 0.2, 
                                   angle = 30, size = 13)) + 
  ylab('Leave-out Dataset') + xlab('')

binary.LODO.train.result.plot[[2]] <- binary.LODO.train.df %>%  
  mutate(Data = factor(Data, levels = c("MultBinary.all", "MultBinary.kw"))) %>%
  mutate(Project = factor(Project, levels = unique(metadata$Project))) %>% 
  group_by(Data) %>% 
  summarise(AUC.mean = mean(AUC)) %>%
  ggplot(aes(x=Data, y=1, fill = AUC.mean)) + 
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUC.mean, digits=2)"), col='white', size=6)+
  scale_fill_gradientn(colours = col.lodo.heatmap, limits=c(0.5, 1), 
                       guide=FALSE)  + 
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text = element_blank()) + 
  ylab('Model Average') + xlab('')

binary.LODO.test.df <- list(Binary.all = data.frame(AUC = roc(binary.loso.prob.df$truth,binary.loso.prob.df$mean_Case)$auc %>% as.numeric,
                                                        Data = "MultBinary.all"),
                            Binary.kw = data.frame(AUC = roc(kw.binary.loso.prob.df$truth, kw.binary.loso.prob.df$mean_Case)$auc %>% as.numeric, 
                                                       Data = "MultBinary.kw")) %>% reduce(rbind)

binary.LODO.test.result.plot <- binary.LODO.test.df %>%  
  mutate(Data = factor(Data, levels = c("MultBinary.all", "MultBinary.kw"))) %>% 
  group_by(Data) %>%  
  ggplot(aes(x = Data, y = 1, fill = AUC)) + 
  geom_tile() + theme_bw() +
  geom_text(aes_string(label="format(AUC, digits=2)"), col='white', size=6) +
  scale_fill_gradientn(colours = col.lodo.heatmap, limits=c(0.5, 1), 
                       guide=FALSE)  + 
  scale_x_discrete(position='top') + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text = element_blank()) + 
  ylab('LODO validation') + xlab('')


lodo.plot <- plot_grid(LODO.train.result.plot[[1]],
                       LODO.train.result.plot[[2]],
                       LODO.test.result.plot, 
                       rel_heights = c(12/16, 2/16, 2/16), 
                       align = 'v', 
                       ncol = 1)

binary.lodo.plot <- plot_grid(binary.LODO.train.result.plot[[1]],
                              binary.LODO.train.result.plot[[2]],
                              binary.LODO.test.result.plot, 
                              rel_heights = c(12/16, 2/16, 2/16), 
                              align = 'v', 
                              ncol = 1)

lodo.plot.grid <- plot_grid(lodo.plot, binary.lodo.plot, 
                            rel_widths = c(12/16, 4/16),
                            labels = "AUTO", ncol = 2, 
                            align = 'h', axis = '1')
pdf("../pdf/Figure4_SF_LODO.pdf",
    width = 12, 
    height = 8)
lodo.plot.grid
dev.off()

save(LODO.train.result.df, 
     LODO.test.result.df, 
     LODO.train.result.plot, 
     LODO.test.result.plot, 
     lodo.plot, 
     binary.LODO.train.df, 
     binary.LODO.test.df, 
     binary.LODO.train.result.plot, 
     binary.LODO.test.result.plot, 
     binary.lodo.plot, 
     lodo.plot.grid, 
     file = "SuppleFig_LODO.RData")
