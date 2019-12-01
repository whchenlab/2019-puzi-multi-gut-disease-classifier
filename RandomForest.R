repeated.rf.confusion <- function(a, times = 10, k = 10){
  ##require group name is "Group"
  sample.split <- createMultiFolds(a$Group, times = times, k=k)
  a$Group <- as.factor(a$Group)
  rf.models <- list()
  pred.matrix <- list()
  confusion.result <- list()
  for (iter in names(sample.split)) {
    train.data <- a[sample.split[[iter]], ]
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

## ---------for multiple-class models ------##
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
    mutate(Reference = factor(Reference, levels = c( "CRC", "CD","UC","Control"))) %>%
    mutate(Prediction = factor(Prediction, levels = c( "Control", "UC","CD","CRC"))) %>%
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
              data.sig.confusion.table = data.sig.confusion.table,
              data.sig.confusion.plot = data.sig.confusion.plot))
}

## ---------for binary models --------------##
roc.plot <- function(data, data.df.group){
  data.sig.result.pred <- lapply(data$pred.matrix.list,function(data){ y <- reduce(data, rbind);name <- names(y);names(y) <- c("Group","Sample"); y.2 <- y %>% group_by(Sample) %>% summarise(mean = mean(Group)); names(y.2) <- rev(name); return(y.2)})
  data.sig.result.pred.df <- data.sig.result.pred %>% reduce(merge, by = "Sample", all = T)
  data.sig.result.pred.df$max <- apply(data.sig.result.pred.df, 1, function(data){y <- names(which.max(data))})
  rownames(data.sig.result.pred.df) <- data.sig.result.pred.df$Sample
  data.sig.result.pred.df <- data.sig.result.pred.df[rownames(data.df.group), ]
  data.sig.result.pred.df$max <- factor(data.sig.result.pred.df$max)
  data.sig.result.pred.df$true <- data.df.group$Group
  data.sig.result.pred.df$true <- factor(data.sig.result.pred.df$true)
  
  forestpred <- prediction(data.sig.result.pred.df[,3], data.sig.result.pred.df$true);
  forestperf <- performance(forestpred, "tpr", "fpr");
  forest.auc <- performance(forestpred, "auc")@y.values[[1]]
  data.roc.df <- data.frame(fpr = forestperf@x.values[[1]], 
                            tpr = forestperf@y.values[[1]])
  
  data.roc.plot <- ggplot(data.roc.df, aes(x = fpr, y = tpr)) + geom_line() + theme_bw() + 
    labs(x = "False Positve Rate", y = "True Positive Rate", title = "ROC plot")  + 
    annotate( "text", x = .75, y = .75, label = paste0("AUC:",round(forest.auc, 2)), size = 4) +
    theme(plot.title = element_text(hjust = .5, face = "bold", size = 11), 
          axis.text = element_text(size = 9))
  return(list(pred.df =data.sig.result.pred.df,
              data.roc.df = data.roc.df,
              forest.auc = forest.auc,
              data.roc.plot = data.roc.plot))
}

## ---------give an example ----------------##

load("data_all_rfmodel.RData")

data.df.group <- data.frame(data.list$data.df.t.nohdc$Group)
rownames(data.df.group) <- rownames(data.list$data.df.t.nohdc)
names(data.df.group)[1] <- "Group"

data.df.case <- subset(data.df.group, !Group %in% "Control")

## -----------data.noctr.rf.result -----------##
## -----------repeated randomforest for three-class model with species ---------##

confusion.result <- list()
confusion.result$spe.all <- confusion.plot(data.noctr.rf.result$data.df.t.nohdc, 
                                   data.df.case)
confusion.result$spe.kw <- confusion.plot(data.noctr.rf.result$data.df.t.kw.nohdc, 
                                           data.df.case)

for (pro in names(confusion.result)) {
  confusion.result[[pro]]$data.sig.confustion.freq <- confusion.result[[pro]]$data.sig.confusion.table %>% 
    mutate(Reference = factor(Reference, levels = c( "CRC", "CD","UC","Control"))) %>%
    mutate(Prediction = factor(Prediction, levels = c( "Control", "UC","CD","CRC"))) %>%
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
          strip.text = element_blank()) + 
    xlab('True positive rate') + ylab('')
}

for (pro in names(confusion.result)) {
  confusion.result[[pro]]$data.sig.confusion.plot <- confusion.result[[pro]]$data.sig.confusion.plot +
    scale_fill_gradient(low = "white", high = "#ff7f7f",limits=c(0,1)) + 
    theme(legend.position = "none")
}

confusion.all.plot.1 <- plot_grid(confusion.result$spe.all$data.sig.confusion.plot,
          confusion.result$spe.all$data.sig.confustion.freq,
          align = "h", ncol = 2,
          rel_widths = c(3/5, 2/5))

confusion.all.plot.2 <- plot_grid(confusion.result$spe.kw$data.sig.confusion.plot,
              confusion.result$spe.kw$data.sig.confustion.freq,
              align = "h", ncol = 2,
              rel_widths = c(3/5, 2/5))


