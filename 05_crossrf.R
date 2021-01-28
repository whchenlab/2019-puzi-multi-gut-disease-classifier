library(purrr)
library(dplyr)
library(pROC)
library(ROCR)
library(randomForest)
library(caret)

## two types of kw.spe
kw.spe <- unique(kw.spe.df$feature)


##
crc.kw.spe <- kw.spe.df %>%
  filter(exposure %in% "CRC") %>%
  pull(feature) %>%
  unique

cd.kw.spe <- kw.spe.df %>%
  filter(exposure %in% "CD") %>%
  pull(feature) %>%
  unique

uc.kw.spe <- kw.spe.df %>%
  filter(exposure %in% "UC") %>%
  pull(feature) %>%
  unique


abund.s.list.filter.df <- merge(abund.s.list.filter[["CRC"]], 
                                abund.s.list.filter[["CD"]], 
                                   by = "features", all = T)
abund.s.list.filter.df <- merge(abund.s.list.filter.df, 
                                abund.s.list.filter[["UC"]][, !names(abund.s.list.filter[["UC"]]) %in% cd.meta$Sample_ID],
                                   by = "features", all = T)
abund.s.list.filter.df[is.na(abund.s.list.filter.df)] <- 0
rownames(abund.s.list.filter.df) <- abund.s.list.filter.df$features

abund.s.list.filter.df.t <- t(abund.s.list.filter.df[,-1])

abund.s.list.rf <- data.frame(abund.s.list.filter.df.t)
abund.s.list.rf$Group <- metadata[rownames(abund.s.list.rf), "Group"]

## For cross rf
all.spe.list <- list(four.gp = abund.s.list.rf, 
                     three.gp = abund.s.list.rf[-which(abund.s.list.rf$Group %in% "CTR"), ],
                     binary = abund.s.list.rf)

all.spe.list$binary[-which(all.spe.list$binary$Group %in% "CTR"), "Group"] <- "Case"

kw.spe.list <- lapply(all.spe.list, 
                      function(data){
                        y <- data[, names(data) %in% c(unique(kw.spe.df$feature), "Group")];
                        return(y)
                      })

save(all.spe.list, 
     kw.spe.list, 
     metadata,
     file = "../RData/species_list_rf.Rdata")

## Pathways data ------------
kw.path.df <- read.delim("../SuppleData/metafor_pathways_coef_signif.txt",
                         header = T, sep = "\t", as.is = T)

abund.path.list.filter.df <- merge(path.abund.fil.list[["CRC"]], 
                                   path.abund.fil.list[["CD"]], 
                                by = "features", all = T)
abund.path.list.filter.df <- merge(abund.path.list.filter.df, 
                                   path.abund.fil.list[["UC"]][, !names(path.abund.fil.list[["UC"]]) %in% cd.meta$Sample_ID],
                                by = "features", all = T)

abund.path.list.filter.df[is.na(abund.path.list.filter.df)] <- 0
rownames(abund.path.list.filter.df) <- abund.path.list.filter.df$features

length(intersect(names(abund.path.list.filter.df),
                 metadata$Sample_ID))

abund.path.list.filter.df.t <- t(abund.path.list.filter.df[,-1])

pathways.map <- data.frame(Pathways = dimnames(abund.path.list.filter.df.t)[[2]],
                           stringsAsFactors = F)
pathways.map$Abb <- lapply(strsplit(pathways.map$Pathways, 
                             split = ":", fixed = T), function(data){ y <- data[1]}) %>% unlist
pathways.map$Abb <- gsub(pathways.map$Abb, pattern = "+", replacement = "_", fixed = T)
pathways.map$Abb <- gsub(pathways.map$Abb, pattern = "-", replacement = "_", fixed = T)
pathways.map$Abb <- paste0( "Pathway_", pathways.map$Abb )

rownames(pathways.map) <- pathways.map$Pathways

abund.path.list.rf <- data.frame(abund.path.list.filter.df.t)
names(abund.path.list.rf) <- pathways.map$Abb
abund.path.list.rf$Group <- metadata[rownames(abund.path.list.rf), "Group"]

## merge kw.path.df
kw.path.df.merge <- merge(kw.path.df, 
                          pathways.map, 
                          by.x = "feature", 
                          by.y = "Pathways",
                          all.x = T)

kw.path  <- kw.path.df.merge %>% filter(qval.fdr < 0.05) %>%
  pull(Abb) %>%
  unique


## For cross rf
all.path.list <- list(four.gp = abund.path.list.rf, 
                     three.gp = abund.path.list.rf[-which(abund.path.list.rf$Group %in% "CTR"), ],
                     binary = abund.path.list.rf)

all.path.list$binary[-which(all.path.list$binary$Group %in% "CTR"), "Group"] <- "Case"

kw.path.list <- lapply(all.path.list, 
                      function(data){
                        y <- data[, names(data) %in% c(unique(kw.path.df.merge$Abb), "Group")];
                        return(y)
                      })

save(all.path.list, 
     kw.path.list, 
     metadata,
     file = "../RData/pathways_list_rf.Rdata")

## Create combined data
abund.combined.list.rf <- cbind(abund.s.list.rf,
                                abund.path.list.rf[rownames(abund.s.list.rf), !names(abund.path.list.rf) %in% "Group"])
all.combined.list <- list(four.gp = abund.combined.list.rf, 
                          three.gp = abund.combined.list.rf[-which(abund.combined.list.rf$Group %in% "CTR"), ],
                          binary = abund.combined.list.rf)
all.combined.list$binary[-which(all.combined.list$binary$Group %in% "CTR"), "Group"] <- "Case"

kw.combined.list <- lapply(all.combined.list, 
                       function(data){
                         y <- data[, names(data) %in% c(unique(kw.path.df.merge$Abb), unique(kw.spe.df$feature), "Group")];
                         return(y)
                       })

save(all.combined.list, 
     kw.combined.list, 
     metadata,
     file = "../RData/combined_list_rf.Rdata")

## Rscript allfeat_list01.R ../RData/combined_list_rf.Rdata combined_crossrf_result
## Rscript allfeat_list.R ../RData/pathways_list_rf.Rdata pathways_crossrf_result
## -------------------------- ten times ten-fold cross-validation
repeated.rf.confusion <- function(a, times = 10, k = 10){
  library(ROCR)
  library(dplyr)
  library(randomForest)
  library(caret)
  library(purrr)
  library(pROC) 
  library(ggplot2)
  library(reshape2)

  ##require group name is "Group"
  sample.split <- createMultiFolds(a$Group, times = times, k=k)
  a$Group <- as.factor(a$Group)
  ## normalization
    train.data.mean <- colMeans(a[,! names(a) %in% "Group"])
    train.data.sd <- apply(a[,! names(a) %in% "Group"], 2, sd)
    train.data.quantile <- quantile(train.data.sd, 0.1, names = F)

    train.data.norm <- a
    for(i in names(train.data.mean)){
        train.data.norm[,i] <- (a[,i]-train.data.mean[i])/(train.data.sd[i] + train.data.quantile)
    }
 
  a <- train.data.norm

  rf.models <- list()
  pred.matrix <- list()
  confusion.result <- list()
  for (iter in names(sample.split)) {
    train.data  <- a[sample.split[[iter]], ]
    
    test.data <- a[-sample.split[[iter]], ]
    rf.models[[iter]] <- randomForest(Group ~ .,train.data, proximity = T, importance = T, ntree = 500)
    predictions <- predict(rf.models[[iter]], test.data, type = "prob")
    pred.matrix[[iter]] <- predictions

    predictions.class <- predict(rf.models[[iter]], test.data, type= "response")
    confusion.result[[iter]] <- predictions.class
  }
  pred.matrix.list <- list()
  pred.df <- data.frame()
  pred.df.2 <- data.frame()

  fold.str <- paste("Fold", gsub(" ", "0", format(1:k)), sep = "")
  rep.str <- paste("Rep", gsub(" ", "0", format(1:times)), sep = "")
  for (rep in rep.str) {
    for (fold in fold.str) {
      iter <- paste(fold, rep, sep = ".")
      pred.df.2 <- data.frame(pred.matrix[[iter]], stringsAsFactors = F)
      pred.df <- rbind(pred.df, pred.df.2)
    }
    for (grp in names(pred.df)) {
      y <- data.frame(pred.df[,grp], stringsAsFactors = F);
      names(y) <- grp;
      y$Sample <- rownames(pred.df)
      pred.matrix.list[[grp]][[rep]] <- y
    }
    pred.df <- data.frame();
  }
  ## --- importance scores --- ##
  importance <- lapply(rf.models, function(data){ imp <- importance(data, type=1)})
  importance.df <- reduce(importance, cbind)
  importance.df <- data.frame(importance.df, stringsAsFactors = F)
  names(importance.df) <- names(rf.models)
  return(list(sample.split = sample.split,
              rf.models = rf.models,
              pred.matrix = pred.matrix,
              confusion.result = confusion.result,
              importance.df = importance.df,
              pred.matrix.list = pred.matrix.list
              #roc.df = roc.df,
              #roc.plot= roc.plot,
              #forest.auc = forest.auc,
  ))
}
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

## for example
kw.combined.log.rf.result  <- parLapply(cl, kw.combined.log.list, repeated.rf.confusion)
kw.spe.log.rf.result  <- parLapply(cl, kw.spe.log.list, repeated.rf.confusion)
kw.path.log.rf.result <- parLapply(cl, kw.path.log.list, repeated.rf.confusion)

save(kw.combined.log.rf.result, kw.spe.log.rf.result, kw.path.log.rf.result,
     file = "species_crossrf_result.RData")
stopCluster(cl)

## --------------------------
load("all_normData_rf_results.RData")
load("kw_normData_rf_results.RData")
## dif- data
data.df.group.list <- list(four.gp = data.frame(Group = kw.combined.list$four.gp$Group,
                                                stringsAsFactors = FALSE),
                           three.gp = data.frame(Group = kw.combined.list$three.gp$Group,
                                                stringsAsFactors = FALSE),
                           binary = data.frame(Group = kw.combined.list$binary$Group,
                                                 stringsAsFactors = FALSE))
rownames(data.df.group.list$four.gp) <- rownames(kw.combined.list$four.gp)
rownames(data.df.group.list$three.gp) <- rownames(kw.combined.list$three.gp)
rownames(data.df.group.list$binary) <- rownames(kw.combined.list$binary)

confusion.plot <- function(data, data.df.group){
  data.sig.result.pred <- lapply(data$pred.matrix.list,function(data){ y <- reduce(data, rbind);name <- names(y);names(y) <- c("Group","Sample"); y.2 <- y %>% group_by(Sample) %>% summarise(mean = mean(Group)); names(y.2) <- rev(name); return(y.2)})
  data.sig.result.pred.df <- data.sig.result.pred %>% reduce(merge, by = "Sample", all = T)
  data.sig.result.pred.df$max <- apply(data.sig.result.pred.df, 1, function(data){y <- names(which.max(data))})
  rownames(data.sig.result.pred.df) <- data.sig.result.pred.df$Sample
  data.sig.result.pred.df <- data.sig.result.pred.df[rownames(data.df.group), ]
  data.sig.result.pred.df$max <- factor(data.sig.result.pred.df$max)
  data.sig.result.pred.df$true <- data.df.group$Group
  data.sig.result.pred.df$true <- factor(data.sig.result.pred.df$true)
  
  data.sig.confusion <- confusionMatrix(data.sig.result.pred.df$max, data.sig.result.pred.df$true)
  
  data.sig.confusion.table <- data.frame(data.sig.confusion$table)
  data.sig.confusion.table <- data.sig.confusion.table %>% group_by(Reference) %>% mutate(sum = sum(Freq), percentage = Freq/sum) 
  data.sig.confusion.table <- data.frame(data.sig.confusion.table, stringsAsFactors = F)
  data.sig.confusion.plot <- data.sig.confusion.table %>% 
  mutate(Reference = factor(Reference, levels = c( "CRC", "CD", "UC","CTR"))) %>%
  mutate(Prediction = factor(Prediction, levels = c( "CTR", "UC", "CD","CRC"))) %>%
  ggplot(aes(x = Prediction, y = Reference)) + 
  geom_tile(aes(fill = percentage)) + theme_bw() +
  geom_text(aes_string(label="Freq"), col='black', size=7)+
  scale_fill_gradient(low = "white", high = "#ff7f7f") + 
  theme(axis.line=element_blank(), 
  axis.ticks = element_blank(),  
  strip.background = element_blank(), 
  strip.text = element_blank(),
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 14),
  panel.border=element_blank(),
  plot.title = element_text(vjust = .5, hjust = .5, size = 14)) +
    labs(title = paste0("Accuracy: ", format(data.sig.confusion$overall[1],digits = 2),"\r\n",
  " 95%CI: ", format(data.sig.confusion$overall[3],digits = 2),
  "-",format(data.sig.confusion$overall[4],digits = 2)),")")
  return(list(pred.df =data.sig.result.pred.df,
              data.sig.confusion = data.sig.confusion,
              data.sig.confusion.plot = data.sig.confusion.plot,
              data.sig.confusion.table = data.sig.confusion.table))
}

kw.combined.rf.confusion <- list()
kw.combined.rf.confusion$four.gp <- confusion.plot(kw.combined.log.rf.result$four.grp, 
                                                   data.df.group.list$four.gp)
kw.combined.rf.confusion$three.gp <- confusion.plot(kw.combined.log.rf.result$three.grp, 
                                                    data.df.group.list$three.gp)
kw.combined.rf.confusion$binary <- roc.plot(kw.combined.log.rf.result$binary,
                                            data.df.group.list$binary)
kw.combined.rf.fourgp.pred <- kw.combined.rf.confusion$four.gp$pred.df
kw.combined.rf.fourgp.pred <- merge(kw.combined.rf.fourgp.pred, 
                                    metadata[, c("Sample_ID", "Project")],
                                    by.x = "Sample", 
                                    by.y  = "Sample_ID",
                                    all = T)
kw.combined.rf.fourgp.pred %>% group_by(Project) %>%
  summarise(count = n(),
         right = sum(max ==true),
         accuracy = right/count)


kw.spe.rf.confusion <- list()
kw.spe.rf.confusion$four.gp <- confusion.plot(kw.spe.log.rf.result$four.grp, 
                                              data.df.group.list$four.gp)
kw.spe.rf.confusion$three.gp <- confusion.plot(kw.spe.log.rf.result$three.grp, 
                                               data.df.group.list$three.gp)
kw.spe.rf.confusion$binary <- roc.plot(kw.spe.log.rf.result$binary,
                                      data.df.group.list$binary)

kw.path.rf.confusion <- list()
kw.path.rf.confusion$four.gp <- confusion.plot(kw.path.log.rf.result$four.grp, 
                                              data.df.group.list$four.gp)
kw.path.rf.confusion$three.gp <- confusion.plot(kw.path.log.rf.result$three.grp, 
                                               data.df.group.list$three.gp)
kw.path.rf.confusion$binary <- roc.plot(kw.path.log.rf.result$binary,
                                       data.df.group.list$binary)

## all features
load("all_normData_rf_results.RData") 

all.combined.rf.confusion <- list()
all.combined.rf.confusion$four.gp <- confusion.plot(all.combined.log.rf.result$four.grp, 
                                                   data.df.group.list$four.gp)
all.combined.rf.confusion$three.gp <- confusion.plot(all.combined.log.rf.result$three.grp, 
                                                    data.df.group.list$three.gp)
all.combined.rf.confusion$binary <- roc.plot(all.combined.log.rf.result$binary,
                                            data.df.group.list$binary)

all.spe.rf.confusion <- list()
all.spe.rf.confusion$four.gp <- confusion.plot(all.spe.log.rf.result$four.grp, 
                                              data.df.group.list$four.gp)
all.spe.rf.confusion$three.gp <- confusion.plot(all.spe.log.rf.result$three.grp, 
                                               data.df.group.list$three.gp)
all.spe.rf.confusion$binary <- roc.plot(all.spe.log.rf.result$binary,
                                       data.df.group.list$binary)

all.path.rf.confusion <- list()
all.path.rf.confusion$four.gp <- confusion.plot(all.path.log.rf.result$four.grp, 
                                               data.df.group.list$four.gp)
all.path.rf.confusion$three.gp <- confusion.plot(all.path.log.rf.result$three.grp, 
                                                data.df.group.list$three.gp)
all.path.rf.confusion$binary <- roc.plot(all.path.log.rf.result$binary,
                                        data.df.group.list$binary)

## CTR/IBD/CRC models
load("normData_threegrp_rf_results.RData") 
## three.all.log.rf.result
data.df.group.list$three.all.gp <- data.df.group.list$four.gp
data.df.group.list$three.all.gp[which(data.df.group.list$three.all.gp$Group %in% c("UC", "CD")), "Group"] <- "IBD"

confusion.plot.2 <- function(data, data.df.group){
  data.sig.result.pred <- lapply(data$pred.matrix.list,function(data){ y <- reduce(data, rbind);name <- names(y);names(y) <- c("Group","Sample"); y.2 <- y %>% group_by(Sample) %>% summarise(mean = mean(Group)); names(y.2) <- rev(name); return(y.2)})
  data.sig.result.pred.df <- data.sig.result.pred %>% reduce(merge, by = "Sample", all = T)
  data.sig.result.pred.df$max <- apply(data.sig.result.pred.df, 1, function(data){y <- names(which.max(data))})
  rownames(data.sig.result.pred.df) <- data.sig.result.pred.df$Sample
  data.sig.result.pred.df <- data.sig.result.pred.df[rownames(data.df.group), ]
  data.sig.result.pred.df$max <- factor(data.sig.result.pred.df$max)
  data.sig.result.pred.df$true <- data.df.group$Group
  data.sig.result.pred.df$true <- factor(data.sig.result.pred.df$true)
  
  data.sig.confusion <- confusionMatrix(data.sig.result.pred.df$max, data.sig.result.pred.df$true)
  
  data.sig.confusion.table <- data.frame(data.sig.confusion$table)
  data.sig.confusion.table <- data.sig.confusion.table %>% group_by(Reference) %>% mutate(sum = sum(Freq), percentage = Freq/sum) 
  data.sig.confusion.table <- data.frame(data.sig.confusion.table, stringsAsFactors = F)
  data.sig.confusion.plot <- data.sig.confusion.table %>% 
  mutate(Reference = factor(Reference, levels = c( "CRC", "IBD","CTR"))) %>%
  mutate(Prediction = factor(Prediction, levels = c( "CTR", "IBD","CRC"))) %>%
  ggplot(aes(x = Prediction, y = Reference)) + 
  geom_tile(aes(fill = percentage)) + theme_bw() +
  geom_text(aes_string(label="Freq"), col='black', size=7)+
  scale_fill_gradient(low = "white", high = "#ff7f7f") + 
  theme(axis.line=element_blank(), 
  axis.ticks = element_blank(),  
  strip.background = element_blank(), 
  strip.text = element_blank(),
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 14),
  panel.border=element_blank(),
  plot.title = element_text(vjust = .5, hjust = .5, size = 14)) +
    labs(title = paste0("Accuracy: ", format(data.sig.confusion$overall[1],digits = 2),"\r\n",
  " 95%CI: ", format(data.sig.confusion$overall[3],digits = 2),
  "-",format(data.sig.confusion$overall[4],digits = 2)),")")
  return(list(pred.df =data.sig.result.pred.df,
              data.sig.confusion = data.sig.confusion,
              data.sig.confusion.plot = data.sig.confusion.plot,
              data.sig.confusion.table = data.sig.confusion.table))
}

three.all.cv.rf.confusion <- lapply(three.all.log.rf.result, 
                                    function(data){
                                      y <- confusion.plot.2(data, data.df.group.list$three.all.gp)
                                    })

##
save(data.df.group.list, 
     kw.combined.rf.confusion, 
     kw.spe.rf.confusion, 
     kw.path.rf.confusion, 
     all.combined.rf.confusion, 
     all.spe.rf.confusion, 
     all.path.rf.confusion, 
     three.all.cv.rf.confusion, 
     file = "rep10_10fold/models_confusion_plot_20201201.RData" )

load("rep10_10fold/models_confusion_plot_20201201.RData")

## three all groups
for (pro in names(three.all.cv.rf.confusion)) {
  three.all.cv.rf.confusion[[pro]]$data.sig.confustion.freq <- three.all.cv.rf.confusion[[pro]]$data.sig.confusion.table %>% 
    mutate(Reference = factor(Reference, levels = c( "CRC", "IBD","CTR"))) %>%
    mutate(Prediction = factor(Prediction, levels = c( "CTR", "IBD","CRC"))) %>%
    filter(Prediction == Reference) %>% 
    ggplot(aes(x = 1, y = Reference)) + 
    geom_tile(aes(fill = percentage)) + theme_bw()+
    geom_text(aes_string(label="format(percentage, digits=2)"), col='black', size=7)+
    scale_fill_gradient(low = "white", high = "#ff7f7f",limits=c(0,1)) + 
    scale_x_discrete(position='top') + 
    theme(axis.line=element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          panel.grid=element_blank(), 
          panel.border=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank(), 
          axis.title = element_text(size = 14)) + 
    xlab('TPR') + ylab('')  + 
    theme(legend.position = "none")
  
  three.all.cv.rf.confusion[[pro]]$data.sig.confusion.plot <- three.all.cv.rf.confusion[[pro]]$data.sig.confusion.plot +
    scale_fill_gradient(low = "white", high = "#ff7f7f",limits=c(0,1)) + 
    theme(legend.position = "none")
  
}


three.all.cv.rf.plot <- list()
for (pro in names(three.all.cv.rf.confusion)) {
  three.all.cv.rf.plot[[pro]] <- plot_grid(three.all.cv.rf.confusion[[pro]]$data.sig.confusion.plot,
                                           three.all.cv.rf.confusion[[pro]]$data.sig.confustion.freq,
                                           align = "h", ncol = 2,
                                           rel_widths = c(2/3, 1/3))
}


pdf("../pdf/Figure4_threeall_cv.pdf",
    width = 7, 
    height = 5)
three.all.cv.rf.plot
dev.off()

## four groups cv plot
multiclass.cv.rf.confusion <- list(all.combined = all.combined.rf.confusion$four.gp,
                                   all.spe = all.spe.rf.confusion$four.gp,
                                   all.path = all.path.rf.confusion$four.gp,
                                   kw.combined = kw.combined.rf.confusion$four.gp,
                                   kw.spe = kw.spe.rf.confusion$four.gp,
                                   kw.path = kw.path.rf.confusion$four.gp)

for (pro in names(multiclass.cv.rf.confusion)) {
  multiclass.cv.rf.confusion[[pro]]$data.sig.confustion.freq <- multiclass.cv.rf.confusion[[pro]]$data.sig.confusion.table %>% 
    mutate(Reference = factor(Reference, levels = c( "CRC", "CD", "UC","CTR"))) %>%
    mutate(Prediction = factor(Prediction, levels = c( "CTR", "UC", "CD", "CRC"))) %>%
    mutate(percentage = round(percentage, digits = 2)) %>%
    filter(Prediction == Reference) %>% 
    ggplot(aes(x = 1, y = Reference)) + 
    geom_tile(aes(fill = percentage)) + theme_bw()+
    geom_text(aes_string(label="format(percentage, digits=2)"), col='black', size=7)+
    scale_fill_gradient(low = "white", high = "#ff7f7f",limits=c(0,1)) + 
    scale_x_discrete(position='top') + 
    theme(axis.line=element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          panel.grid=element_blank(), 
          panel.border=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank(), 
          axis.title = element_text(size = 14)) + 
    xlab('TPR') + ylab('')  + 
    theme(legend.position = "none")
  
  multiclass.cv.rf.confusion[[pro]]$data.sig.confusion.plot <- multiclass.cv.rf.confusion[[pro]]$data.sig.confusion.plot +
    scale_fill_gradient(low = "white", high = "#ff7f7f",limits=c(0,1)) + 
    theme(legend.position = "none")
  
}


multiclass.cv.rf.plot <- list()
for (pro in names(multiclass.cv.rf.confusion)) {
  multiclass.cv.rf.plot[[pro]] <- plot_grid(multiclass.cv.rf.confusion[[pro]]$data.sig.confusion.plot,
                                           multiclass.cv.rf.confusion[[pro]]$data.sig.confustion.freq,
                                           align = "h", ncol = 2,
                                           rel_widths = c(2/3, 1/3))
}


pdf("../pdf/Figure4_multiclass_cv.pdf",
    width = 7, 
    height = 5)
multiclass.cv.rf.plot
dev.off()


## three groups cv plot
three.cv.rf.confusion <- list(all.combined = all.combined.rf.confusion$three.gp,
                                   all.spe = all.spe.rf.confusion$three.gp,
                                   all.path = all.path.rf.confusion$three.gp,
                                   kw.combined = kw.combined.rf.confusion$three.gp,
                                   kw.spe = kw.spe.rf.confusion$three.gp,
                                   kw.path = kw.path.rf.confusion$three.gp)

for (pro in names(three.cv.rf.confusion)) {
  three.cv.rf.confusion[[pro]]$data.sig.confustion.freq <- three.cv.rf.confusion[[pro]]$data.sig.confusion.table %>% 
    mutate(Reference = factor(Reference, levels = c( "CRC", "CD", "UC"))) %>%
    mutate(Prediction = factor(Prediction, levels = c( "UC", "CD", "CRC"))) %>%
    mutate(percentage = round(percentage, digits = 2)) %>%
    filter(Prediction == Reference) %>% 
    ggplot(aes(x = 1, y = Reference)) + 
    geom_tile(aes(fill = percentage)) + theme_bw()+
    geom_text(aes_string(label="format(percentage, digits=2)"), col='black', size=7)+
    scale_fill_gradient(low = "white", high = "#ff7f7f",limits=c(0,1)) + 
    scale_x_discrete(position='top') + 
    theme(axis.line=element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          panel.grid=element_blank(), 
          panel.border=element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_blank(), 
          axis.title = element_text(size = 14)) + 
    xlab('TPR') + ylab('')  + 
    theme(legend.position = "none")
  
  three.cv.rf.confusion[[pro]]$data.sig.confusion.plot <- three.cv.rf.confusion[[pro]]$data.sig.confusion.plot +
    scale_fill_gradient(low = "white", high = "#ff7f7f",limits=c(0,1)) + 
    theme(legend.position = "none")
  
}


three.cv.rf.plot <- list()
for (pro in names(three.cv.rf.confusion)) {
  three.cv.rf.plot[[pro]] <- plot_grid(three.cv.rf.confusion[[pro]]$data.sig.confusion.plot,
                                       three.cv.rf.confusion[[pro]]$data.sig.confustion.freq,
                                            align = "h", ncol = 2,
                                            rel_widths = c(2/3, 1/3))
}


pdf("../pdf/Figure4_three_cv.pdf",
    width = 7, 
    height = 5)
three.cv.rf.plot
dev.off()

## Binary cv plot

binary.cv.rf.confusion <- list(all.combined = all.combined.rf.confusion$binary,
                              all.spe = all.spe.rf.confusion$binary,
                              all.path = all.path.rf.confusion$binary,
                              kw.combined = kw.combined.rf.confusion$binary,
                              kw.spe = kw.spe.rf.confusion$binary,
                              kw.path = kw.path.rf.confusion$binary)

for (pro in names(binary.cv.rf.confusion)) {
  binary.cv.rf.confusion[[pro]]$data.roc.df$Type <- paste(pro, 
                                                          round(binary.cv.rf.confusion[[pro]]$forest.auc,
                                                                digits = 2), 
                                                          sep = ":")
  binary.cv.rf.confusion[[pro]]$data.roc.df$Group <- lapply(strsplit(pro, split = ".", fixed = T), 
                                                           function(data){ y <- data[2]}) %>% unlist
}
binary.cv.rf.confusion.df <-  lapply(binary.cv.rf.confusion,
                                     function(data){ y <- data$data.roc.df}) %>%
  reduce(rbind)

binary.cv.rf.confusion.df$Type2 <- "Profiles"
binary.cv.rf.confusion.df[which(binary.cv.rf.confusion.df$Type %in% c("kw.combined:0.87", 
                                                                      "kw.spe:0.84",
                                                                      "kw.path:0.8")), "Type2"] <- "Markers"

binary.cv.rf.plot <- binary.cv.rf.confusion.df %>%
  mutate(Group = factor(Group, 
                        levels = c("combined", "spe", "path"))) %>%
  mutate(Type2 = factor(Type2, 
                        levels = c("Profiles", "Markers"))) %>%
  ggplot(aes(x = fpr, y = tpr, color = Group)) + 
  geom_line(size = 0.9, aes(linetype = Type2))  +
  theme_light() + 
  labs(x = "False Positve Rate", 
       y = "True Positive Rate", 
       title = "ROC plot")  +
  theme(plot.title = element_text(hjust = .5, 
                                  face = "bold", size = 11), 
        axis.text = element_text(size = 11)) +
  scale_color_manual(values = c("#FB8072", "#B3DE69", 
                                "#80B1D3"))

pdf("../pdf/Figure4_binary_cv.pdf",
    width = 5, 
    height = 4)
binary.cv.rf.plot
dev.off()
