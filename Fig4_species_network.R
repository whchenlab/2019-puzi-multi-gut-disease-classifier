setwd("/mnt/raid5/puzi/CDdata/result/")
load("network_rawlist.Rdata")
test.kw.spearman.cor <- list()
test.kw.spearman.pvalue <- list()
for(pro in names(test.kw.list)){
  a <- rcorr(as.matrix(test.kw.list[[pro]]), type = "spearman")
  a.triu <- a$r
  a.triu[upper.tri(a.triu)] <- NA
  a.triu.df <- data.frame(a.triu, stringsAsFactors = F)
  a.triu.df$features <- rownames(a.triu.df) 
  a.triu.df <- melt(a.triu.df, id.vars = "features")
  test.kw.spearman.cor[[pro]] <- a.triu.df
  
  a.triu <- a$P
  a.triu[upper.tri(a.triu)] <- NA
  a.triu.df <- data.frame(a.triu, stringsAsFactors = F)
  a.triu.df$features <- rownames(a.triu.df) 
  a.triu.df <- melt(a.triu.df, id.vars = "features")
  test.kw.spearman.pvalue[[pro]] <- a.triu.df
}
test.kw.spearman.cor.filter <- lapply(test.kw.spearman.cor, 
                                      function(data){
                                        y <- data %>% filter(!features == variable) %>% filter(!is.na(value))
                                        return(y)
                                      })
test.kw.spearman.pvalue.filter <- lapply(test.kw.spearman.pvalue, 
                                         function(data){
                                           y <- data %>% filter(!features == variable) %>% 
                                             filter(!is.na(value))
                                           return(y)
                                         })

for (f in names(test.kw.spearman.cor.filter)) {
  test.kw.spearman.cor.filter[[f]]$disease <- strsplit(f, split = "_",fixed = T)[[1]][2]
  test.kw.spearman.cor.filter[[f]]$project <- strsplit(f, split = "_",fixed = T)[[1]][1]
  test.kw.spearman.cor.filter[[f]]$target <- paste(test.kw.spearman.cor.filter[[f]]$features, test.kw.spearman.cor.filter[[f]]$variable, sep = "_")
}
test.kw.spearman.cor.filter$PRJNA447983_cohort1_CTR$disease <- "CTR"
test.kw.spearman.cor.filter$PRJNA447983_cohort1_CRC$disease <- "CRC"
test.kw.spearman.cor.filter.df <- test.kw.spearman.cor.filter %>% reduce(rbind) 

for (f in names(test.kw.spearman.pvalue.filter)) {
  test.kw.spearman.pvalue.filter[[f]]$disease <- strsplit(f, split = "_",fixed = T)[[1]][2]
  test.kw.spearman.pvalue.filter[[f]]$project <- strsplit(f, split = "_",fixed = T)[[1]][1]
  test.kw.spearman.pvalue.filter[[f]]$target <- paste(test.kw.spearman.pvalue.filter[[f]]$features, test.kw.spearman.pvalue.filter[[f]]$variable, sep = "_")
}
test.kw.spearman.pvalue.filter$PRJNA447983_cohort1_CTR$disease <- "CTR"
test.kw.spearman.pvalue.filter$PRJNA447983_cohort1_CRC$disease <- "CRC"
test.kw.spearman.pvalue.filter.df <- test.kw.spearman.pvalue.filter %>% reduce(rbind) 

test.kw.spearman.pvalue.filter.df <- test.kw.spearman.pvalue.filter.df %>% 
  group_by(project, disease) %>% arrange(value)  %>% 
  mutate(fdr = p.adjust(p = value, method = "fdr"))

## ------------仅取关联后 Fdr值<0.05的关联 --------------##
test.kw.spearman.pvalue.filter.df.2 <- subset(test.kw.spearman.pvalue.filter.df, fdr <0.05)
names(test.kw.spearman.pvalue.filter.df.2)[3] <- "pvalue"
test.kw.spearman.pvalue.filter.df.2 <- data.frame(test.kw.spearman.pvalue.filter.df.2, 
                                                  stringsAsFactors = F)
test.kw.spearman.merge <- merge(test.kw.spearman.pvalue.filter.df.2,
                                test.kw.spearman.cor.filter.df, 
                                by = c("features","variable","disease","project","target"),
                                all.x = T)
test.kw.spearman.merge.2 <- test.kw.spearman.merge %>% filter(!is.na(value))

## --------------------prepare metadata---------------------------------##
metadata <- data.frame(nSamples = unlist(lapply(test.kw.list, nrow)),
                       stringsAsFactors = F)
metadata$project <- unlist(lapply(strsplit(names(test.kw.list), 
                                    split = "_", fixed = T), function(data){
                                      y <- data[[1]];
                                      return(y)
                                    }))
metadata$disease <- unlist(lapply(strsplit(names(test.kw.list), 
                                           split = "_", fixed = T), function(data){
                                             y <- data[[2]];
                                             return(y)
                                           }))
metadata["PRJNA447983_cohort1_CTR", "disease"] <- "CTR"
metadata["PRJNA447983_cohort1_CRC", "disease"] <- "CRC"

test.kw.spearman.merge.3 <- merge(test.kw.spearman.merge.2, metadata,
                                  by = c("disease","project"),
                                  all.x = T)
test.kw.spearman.merge.summary <- test.kw.spearman.merge.3 %>% 
	group_by(disease, target) %>% 
  summarise(count = n())
  
## ----
## -------------在疾病中关联起码在两个项目中存在---------------##

test.kw.spearman.cor.correct <- list()

for (grp in c("CTR","CRC","CD","UC")) {
  targetspe <- test.kw.spearman.merge.3 %>% filter(disease %in% grp) %>%
    pull(target)
  for (comp in targetspe) {
    sparcc.sub <- subset(test.kw.spearman.merge.3, disease %in% grp & target %in% comp)
    if (nrow(sparcc.sub)>1) {
      sub.zscore <- metacor(value, nSamples, project, 
                                   data = sparcc.sub, sm = "ZCOR")
      sub.correct.r <- (exp(2* sub.zscore$TE.random)-1)/(exp(2* sub.zscore$TE.random)+1)
      
      test.kw.spearman.cor.correct[[paste(comp, grp, sep = "_")]]$meta_z <- sub.zscore
      test.kw.spearman.cor.correct[[paste(comp, grp, sep = "_")]]$correct_r <- sub.correct.r
    }
  }
}

## 疾病数据之间分开看 CTR的关联 ----------- ## 
test.kw.spearman.cor.correct.2 <- list()
targetspe <- test.kw.spearman.merge.3 %>% 
  filter(disease %in% "CTR") %>% 
  filter(project %in% c("PRJDB4176","PRJEB10878","PRJEB12449",
                        "PRJEB27928","PRJEB6070",  
                        "PRJEB7774","PRJNA447983")) %>%
  pull(target)

for (comp in targetspe) {
  sparcc.sub <- test.kw.spearman.merge.3 %>% filter(disease %in% "CTR") %>%
    filter(target %in% comp) %>% 
    filter(project %in% c("PRJDB4176","PRJEB10878","PRJEB12449",
                          "PRJEB27928","PRJEB6070",  
                          "PRJEB7774","PRJNA447983"))
  
  if (nrow(sparcc.sub)>1) {
      sub.zscore <- metacor(value, nSamples, project, 
                                   data = sparcc.sub, sm = "ZCOR")
      sub.correct.r <- (exp(2* sub.zscore$TE.random)-1)/(exp(2* sub.zscore$TE.random)+1)
    
    test.kw.spearman.cor.correct.2[[paste(comp, "CRC" ,sep = "_")]]$meta_z <- sub.zscore
    test.kw.spearman.cor.correct.2[[paste(comp,"CRC", sep = "_")]]$correct_r <- sub.correct.r
  }
}

targetspe <- test.kw.spearman.merge.3 %>% 
  filter(disease %in% "CTR") %>% 
  filter(project %in% c("PRJNA389280","PRJNA400072", "SRP057027")) %>%
  pull(target)

for (comp in targetspe) {
  sparcc.sub <- test.kw.spearman.merge.3 %>% 
    filter(disease %in% "CTR") %>%
    filter(target %in% comp) %>% 
    filter(project %in% c("PRJNA389280","PRJNA400072", "SRP057027"))
  if (nrow(sparcc.sub)>1) {
    sub.zscore <- metacor(value, nSamples, project, 
                                 data = sparcc.sub, sm = "ZCOR")
    sub.correct.r <- (exp(2* sub.zscore$TE.random)-1)/(exp(2* sub.zscore$TE.random)+1)
    
    test.kw.spearman.cor.correct.2[[paste(comp, "CD" ,sep = "_")]]$meta_z <- sub.zscore
    test.kw.spearman.cor.correct.2[[paste(comp,"CD", sep = "_")]]$correct_r <- sub.correct.r
  }
}

targetspe <- test.kw.spearman.merge.3 %>% 
  filter(disease %in% "CTR") %>% 
  filter(project %in% c("PRJNA389280","PRJNA400072", "PRJEB1220")) %>%
  pull(target)

for (comp in targetspe) {
  sparcc.sub <- test.kw.spearman.merge.3 %>% filter(disease %in% "CTR") %>%
    filter(target %in% comp) %>% 
    filter(project %in% c("PRJNA389280","PRJNA400072", "PRJEB1220"))
  if (nrow(sparcc.sub)>1) {
    sub.zscore <- metacor(value, nSamples, project, 
                                 data = sparcc.sub, sm = "ZCOR")
    sub.correct.r <- (exp(2* sub.zscore$TE.random)-1)/(exp(2* sub.zscore$TE.random)+1)
    
    
    test.kw.spearman.cor.correct.2[[paste(comp, "UC" ,sep = "_")]]$meta_z <- sub.zscore
    test.kw.spearman.cor.correct.2[[paste(comp,"UC", sep = "_")]]$correct_r <- sub.correct.r
  }
}

## -----------提取corrected.r, fdr ------##
meta.correct.func <- function(test.kw.spearman.cor.correct = data){
  test.kw.spearman.cor.correct.zpvalue <- lapply(test.kw.spearman.cor.correct, function(data){
    y <- data$meta_z$pval.random
    return(y)
  })
  test.kw.spearman.cor.correct.zpvalue.df <- test.kw.spearman.cor.correct.zpvalue %>% reduce(rbind)
  test.kw.spearman.cor.correct.zpvalue.df <- data.frame(test.kw.spearman.cor.correct.zpvalue.df,
                                                        stringsAsFactors = F)
  names(test.kw.spearman.cor.correct.zpvalue.df)[1] <- "random.pvalue"
  test.kw.spearman.cor.correct.zpvalue.df$target <- names(test.kw.spearman.cor.correct.zpvalue)
  
  test.kw.spearman.cor.correct.zpvalue.df$disease <- "CTR"
  test.kw.spearman.cor.correct.zpvalue.df[grep("CD$", test.kw.spearman.cor.correct.zpvalue.df$target), "disease"] <- "CD"
  test.kw.spearman.cor.correct.zpvalue.df[grep("UC$", test.kw.spearman.cor.correct.zpvalue.df$target), "disease"] <- "UC"
  test.kw.spearman.cor.correct.zpvalue.df[grep("CRC$", test.kw.spearman.cor.correct.zpvalue.df$target), "disease"] <- "CRC"
  
  test.kw.spearman.cor.correct.zpvalue.df <- test.kw.spearman.cor.correct.zpvalue.df %>% filter(!is.na(random.pvalue))
  test.kw.spearman.cor.correct.zfdr.df <- test.kw.spearman.cor.correct.zpvalue.df %>% group_by(disease) %>%
    arrange(random.pvalue) %>% 
    mutate(fdr = p.adjust(random.pvalue, method = "fdr"))
  test.kw.spearman.cor.correct.zfdr.df <- data.frame(test.kw.spearman.cor.correct.zfdr.df,
                                                     stringsAsFactors = F)
  
  ########## I2 assessment -----------##
  test.kw.spearman.cor.correct.I2.df <- lapply(test.kw.spearman.cor.correct, function(data){
    y <- data$meta_z$I2
    return(y)
  }) %>% reduce(rbind)
  test.kw.spearman.cor.correct.I2.df <- data.frame(test.kw.spearman.cor.correct.I2.df,
                                                   stringsAsFactors = F)
  test.kw.spearman.cor.correct.I2.df$target <- names(test.kw.spearman.cor.correct)
  
  
  ## ------------corrected coefficient ----------##
  test.kw.spearman.cor.correct.r.df <- lapply(test.kw.spearman.cor.correct, function(data){
    y <- data$correct_r
    return(y)
  })  %>% reduce(rbind)
  
  test.kw.spearman.cor.correct.r.df <- data.frame(test.kw.spearman.cor.correct.r.df,
                                                  stringsAsFactors = F)
  names(test.kw.spearman.cor.correct.r.df)[1] <- "random.r"
  test.kw.spearman.cor.correct.r.df$target <- names(test.kw.spearman.cor.correct)
  
  test.kw.spearman.cor.correct.r.df$disease <- "CTR"
  test.kw.spearman.cor.correct.r.df[grep("CD$", test.kw.spearman.cor.correct.r.df$target), "disease"] <- "CD"
  test.kw.spearman.cor.correct.r.df[grep("UC$", test.kw.spearman.cor.correct.r.df$target), "disease"] <- "UC"
  test.kw.spearman.cor.correct.r.df[grep("CRC$", test.kw.spearman.cor.correct.r.df$target), "disease"] <- "CRC"
  
  test.kw.spearman.cor.correct.r.df <- test.kw.spearman.cor.correct.r.df %>% filter(!is.na(random.r))
  
  test.kw.spearman.cor.merge <- merge(test.kw.spearman.cor.correct.r.df,
                                      test.kw.spearman.cor.correct.zfdr.df,
                                      by = c("target","disease"), all = T)
  
  test.kw.spearman.cor.merge$target_real <- test.kw.spearman.cor.merge$target
  test.kw.spearman.cor.merge$target_real <- gsub("_CTR", "",test.kw.spearman.cor.merge$target_real)
  test.kw.spearman.cor.merge$target_real <- gsub("_CRC", "",test.kw.spearman.cor.merge$target_real)
  test.kw.spearman.cor.merge$target_real <- gsub("_CD", "",test.kw.spearman.cor.merge$target_real)
  test.kw.spearman.cor.merge$target_real <- gsub("_UC", "",test.kw.spearman.cor.merge$target_real)
  return(list(correct.zpvalue.df =  test.kw.spearman.cor.correct.zpvalue.df,
              correct.fdr.df = test.kw.spearman.cor.correct.zfdr.df,
              correct.I2.df =test.kw.spearman.cor.correct.I2.df,
              correct.r2.df = test.kw.spearman.cor.correct.r.df,
              cor.merge = test.kw.spearman.cor.merge
  ))
}

## -----------疾病数据中疾病的网络结果-----##
test.kw.spearman.cor.merge.result <- meta.correct.func(test.kw.spearman.cor.correct)

test.kw.spearman.cor.merge.filter <- subset(test.kw.spearman.cor.merge.result$cor.merge,
                                            fdr <0.05)
test.kw.spearman.cor.merge.filter.2 <- merge(test.kw.spearman.cor.merge.filter,
                                             test.kw.spearman.cor.filter.df[,c("features","variable","target")] %>% distinct(target, .keep_all = T),
                                             by.x = "target_real", by.y = "target",
                                             all.x = T)
                                             
## correlation strength quantile
cut_edge <- function(r) {
  out <- cut(abs(r), breaks = c(0, 0.4,0.6,1), include.lowest = F, 
             labels = c("1", "2", "3"))
  return(out)
}

test.kw.spearman.cor.merge.filter.2$edgestrength <- cut_edge(test.kw.spearman.cor.merge.filter.2$random.r)
test.kw.spearman.cor.merge.filter.2$edgetrend <- 1
test.kw.spearman.cor.merge.filter.2[which(test.kw.spearman.cor.merge.filter.2$random.r <0), "edgetrend"] <- -1

##------------疾病数据中正常人结果--------##
test.kw.spear.correct.2.result <- meta.correct.func(test.kw.spearman.cor.correct.2)
test.kw.spear.correct.2.filter <- subset(test.kw.spear.correct.2.result$cor.merge, 
                                         fdr < 0.05)
test.kw.spear.correct.2.filter.2 <- merge(test.kw.spear.correct.2.filter,
                                          test.kw.spearman.cor.filter.df[,c("features","variable","target")] %>% 
                                            distinct(target, .keep_all = T),
                                          by.x = "target_real", by.y = "target",
                                          all.x = T)

## 0<r<=0.4, 0.4<r<=0.6, r>0.6
test.kw.spear.correct.2.filter.2$edgestrength <- cut_edge(test.kw.spear.correct.2.filter.2$random.r)
test.kw.spear.correct.2.filter.2$edgetrend <- 1
test.kw.spear.correct.2.filter.2[which(test.kw.spear.correct.2.filter.2$random.r <0), "edgetrend"] <- -1

test.kw.spear.correct.2.filter.2$variable <- as.character(test.kw.spear.correct.2.filter.2$variable)
## ------
test.kw.spear.correct.2.filter.3 <- list()
test.kw.spear.correct.2.filter.3[["CRC"]] <- test.kw.spear.correct.2.filter.2 %>%
  filter(disease %in% "CRC")
test.kw.spear.correct.2.filter.3[["CD"]] <- test.kw.spear.correct.2.filter.2 %>%
  filter(disease %in% "CD")
test.kw.spear.correct.2.filter.3[["UC"]] <- test.kw.spear.correct.2.filter.2 %>%
  filter(disease %in% "UC")

test.kw.spear.correct.2.filter.3 <- lapply(test.kw.spear.correct.2.filter.3,
                                           function(data){
                                             data$disease <- "CTR";
                                             return(data)
                                           })
                                           
##-----只提取各个疾病的相关菌
## ----------------
kw.sum.df <- read.delim("kw_taxa.txt",
                        header = T, sep = "\t", as.is = T)
kw.spe.crc <- kw.sum.df %>% filter(Disease %in% "CRC") %>%
  pull(features)
kw.spe.cd <- kw.sum.df %>% filter(Disease %in% "CD") %>%
  pull(features)
kw.spe.uc <- kw.sum.df %>% filter(Disease %in% "UC") %>%
  pull(features)

## -----
## --------不在正常人网络中的点的prevalence ------------##
spear.net.ctr.add <- list()
spear.net.ctr.add[["CRC"]] <- test.kw.spear.correct.2.filter.3[["CRC"]] %>%
  filter(features %in% kw.spe.crc) %>% 
  filter(variable %in% kw.spe.crc)

add.spe <- kw.spe.crc[!kw.spe.crc %in% unique(c(spear.net.ctr.add[["CRC"]]$features, 
                                                spear.net.ctr.add[["CRC"]]$variable))]
for (add in add.spe) {
  spear.net.ctr.add[["CRC"]][nrow(spear.net.ctr.add[["CRC"]])+1, 
                                     ] <- c(NA, NA, "CRC",
                                            0, 0, 0, 
                                            add, add, 0, 0)
}

spear.net.ctr.add[["CD"]] <- test.kw.spear.correct.2.filter.3[["CD"]] %>%
  filter(features %in% kw.spe.cd) %>% 
  filter(variable %in% kw.spe.cd)

add.spe <- kw.spe.cd[!kw.spe.cd %in% unique(c(spear.net.ctr.add[["CD"]]$features, 
                                                spear.net.ctr.add[["CD"]]$variable))]
for (add in add.spe) {
  spear.net.ctr.add[["CD"]][nrow(spear.net.ctr.add[["CD"]])+1, 
                             ] <- c(NA, NA, "CD",
                                    0, 0, 0, 
                                    add, add, 0, 0)
}

spear.net.ctr.add[["UC"]] <- test.kw.spear.correct.2.filter.3[["UC"]] %>%
  filter(features %in% kw.spe.uc) %>% 
  filter(variable %in% kw.spe.uc)
add.spe <- kw.spe.uc[!kw.spe.uc %in% unique(c(spear.net.ctr.add[["UC"]]$features, 
                                              spear.net.ctr.add[["UC"]]$variable))]
for (add in add.spe) {
  spear.net.ctr.add[["UC"]][nrow(spear.net.ctr.add[["UC"]])+1, 
                            ] <- c(NA, NA, "UC",
                                   0, 0, 0, 
                                   add, add, 0, 0)
}

write.table(spear.net.ctr.add.abb[["CRC"]], 
           file = "spear_net_meta_CRC_CTRresult.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.table(spear.net.ctr.add.abb[["CD"]], 
            file = "spear_net_meta_CD_CTRresult.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.table(spear.net.ctr.add.abb[["UC"]], 
            file = "spear_net_meta_UC_CTRresult.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)

## --------不在疾病网络中的点的prevalence ------------##
spear.net.disease.add <- list()
test.kw.spearman.cor.merge.filter.2$variable <- as.character(test.kw.spearman.cor.merge.filter.2$variable)

spear.net.disease.add[["CRC_CRC"]] <- test.kw.spearman.cor.merge.filter.2 %>% 
  filter(disease %in% "CRC")  %>%
  filter(features %in% kw.spe.crc) %>% 
  filter(variable %in% kw.spe.crc)
add.spe <- kw.spe.crc[!kw.spe.crc %in% unique(c(spear.net.disease.add[["CRC_CRC"]]$features, 
                                                spear.net.disease.add[["CRC_CRC"]]$variable))]
for (add in add.spe) {
  spear.net.disease.add[["CRC_CRC"]][nrow(spear.net.disease.add[["CRC_CRC"]])+1, 
                                     ] <- c(NA, NA, "CRC",
                                            0, 0, 0, 
                                            add, add, 0, 0)
}


spear.net.disease.add[["CD_CD"]] <- test.kw.spearman.cor.merge.filter.2 %>% 
  filter(disease %in% "CD") %>%
  filter(features %in% kw.spe.cd) %>% 
  filter(variable %in% kw.spe.cd)

add.spe <- kw.spe.cd[!kw.spe.cd %in% unique(c(spear.net.disease.add[["CD_CD"]]$features, 
                                                spear.net.disease.add[["CD_CD"]]$variable))]
for (add in add.spe) {
  spear.net.disease.add[["CD_CD"]][nrow(spear.net.disease.add[["CD_CD"]])+1, 
                                     ] <- c(NA, NA, "CD",
                                            0, 0, 0, 
                                            add, add, 0, 0)
}

spear.net.disease.add[["UC_UC"]] <- test.kw.spearman.cor.merge.filter.2 %>% 
  filter(disease %in% "UC") %>%
  filter(features %in% kw.spe.uc) %>% 
  filter(variable %in% kw.spe.uc)

add.spe <- kw.spe.uc[!kw.spe.uc %in% unique(c(spear.net.disease.add[["UC_UC"]]$features, 
                                              spear.net.disease.add[["UC_UC"]]$variable))]
for (add in add.spe) {
  spear.net.disease.add[["UC_UC"]][nrow(spear.net.disease.add[["UC_UC"]])+1, 
                                   ] <- c(NA, NA, "UC",
                                          0, 0, 0, 
                                          add, add, 0, 0)
}

write.table(spear.net.disease.add.abb[["CD_CD"]], 
            file = "spear_net_meta_CD_result_3.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.table(spear.net.disease.add.abb[["UC_UC"]], 
            file = "spear_net_meta_UC_result_3.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)
write.table(spear.net.disease.add.abb[["CRC_CRC"]], 
            file = "spear_net_meta_CRC_result_3.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)

## --------------network centrality ----------------##
network.func <- function(data.frame){
  ## -----data as network form 
  a <- data.frame %>%
    filter(!features == variable) %>%
    select(features, variable, random.r)
  a$random.r <- as.numeric(a$random.r)
  a$random.r <- abs(a$random.r)
  names(a) <- c("from", "to", "weight")
  data <- graph.data.frame(a, directed = FALSE)
  E(data)$weight <- a$weight
  
  ## -----degree weighted
  y <- graph.strength(data,
                      vids = V(data),
                      weights = abs(E(data)$weight))
  y <- data.frame(y, stringsAsFactors = F)
  names(y) <- "nodesize"
  y$node <- rownames(y)
  ## -----eigen centrality --------------##
  eigen <- eigen_centrality(graph = data,
                        directed = FALSE,
                        weights = E(data)$weight)
  ## ------betweenness -------------------
  betweenness <- betweenness(graph = data,
                   directed = FALSE,
                   weights = E(data)$weight)
  
  return(list(network.data = data,
              weighted.degree = y,
              eigen.data = eigen,
              betweenness = betweenness))
}

network.eigen.plot <- function(spear.net.eigen.df){
  spear.net.eigen.order <- spear.net.eigen.df %>%
    group_by(Species) %>% 
    summarise(sum = sum(Eigen)) %>%
    arrange(desc(sum)) %>%
    pull(Species)
  spear.net.eigen.plot <- spear.net.eigen.df %>% 
    mutate(Species = factor(Species,
                            levels = spear.net.eigen.order)) %>%
    ggplot(aes(x = Species, y = Eigen, fill = type2)) +
    geom_bar(stat = "identity") +
    theme_bw()+
    coord_flip() +
    labs(y = "Eigenvector Centrality Scores") + 
    scale_fill_manual(values = c("#56B4E9","firebrick3"))
  return(list(spear.net.eigen.order = spear.net.eigen.order,
              spear.net.eigen.plot = spear.net.eigen.plot))
}

kw.spearman.cor.ctr.ig <- lapply(spear.net.ctr.add, network.func)
kw.spearman.cor.ig  <- lapply(spear.net.disease.add, network.func)

## --------------整合betweenness与eigen的结果
spear.net.betweenness <- lapply(c(kw.spearman.cor.ctr.ig,
                                  kw.spearman.cor.ig), function(data){
                                    y <- data.frame(data$betweenness)
                                    names(y) <- "Betweenness"
                                    y$Species <- rownames(y)
                                    return(y)
                                  })
spear.net.eigen <- lapply(c(kw.spearman.cor.ctr.ig,
                            kw.spearman.cor.ig), function(data){
                              y <- data.frame(data$eigen.data$vector)
                              names(y) <- "Eigen"
                              y$Species <- rownames(y)
                              return(y)
                            })
spear.net.degree <- lapply()
for (grp in names(spear.net.eigen)) {
  spear.net.eigen[[grp]]$type <- grp
  spear.net.betweenness[[grp]]$type <- grp
}
spear.net.eigen.df <- spear.net.eigen %>% reduce(rbind)
spear.net.eigen.df$type2 <- spear.net.eigen.df$type
spear.net.eigen.df[-grep(spear.net.eigen.df$type, 
                         pattern = "_", fixed = T), 
                   "type2"] <- paste0(spear.net.eigen.df[-grep(spear.net.eigen.df$type,
                                                               pattern = "_", fixed = T), 
                                                         "type"], "_CTR")
spear.net.eigen.df$type2 <- factor(spear.net.eigen.df$type2,
                                   levels = c( "CD_CTR", "CD_CD",
                                               "CRC_CTR", "CRC_CRC",
                                               "UC_CTR","UC_UC"))
spear.net.eigen.plot <- list()
spear.net.eigen.plot[["CRC"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.eigen.df, 
                                                                                type2 %in% c("CRC_CTR","CRC_CRC")))
spear.net.eigen.plot[["CD"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.eigen.df, 
                                                                               type2 %in% c("CD_CTR","CD_CD")))
spear.net.eigen.plot[["UC"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.eigen.df, 
                                                                               type2 %in% c("UC_CTR","UC_UC")))

spear.net.betweenness.df <- spear.net.betweenness %>% reduce(rbind)
spear.net.betweenness.df$type2 <- spear.net.betweenness.df$type
spear.net.betweenness.df[-grep(spear.net.betweenness.df$type, 
                               pattern = "_", fixed = T), 
                         "type2"] <- paste0(spear.net.betweenness.df[-grep(spear.net.betweenness.df$type,
                                                                           pattern = "_", fixed = T), 
                                                                     "type"], "_CTR")
spear.net.betweenness.df$type2 <- factor(spear.net.betweenness.df$type2,
                                         levels = c( "CD_CTR", "CD_CD",
                                                     "CRC_CTR", "CRC_CRC",
                                                     "UC_CTR","UC_UC"))
names(spear.net.betweenness.df)[1] <- "Eigen"
spear.net.bet.plot <- list() 
spear.net.bet.plot[["CRC"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.betweenness.df, 
                                                                              type2 %in% c("CRC_CTR","CRC_CRC")))
spear.net.bet.plot[["CRC"]]$spear.net.eigen.plot <- spear.net.bet.plot[["CRC"]]$spear.net.eigen.plot +
  labs(y = "Betweenness centrality")
spear.net.bet.plot[["CD"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.betweenness.df, 
                                                                             type2 %in% c("CD_CTR","CD_CD")))
spear.net.bet.plot[["CD"]]$spear.net.eigen.plot <- spear.net.bet.plot[["CD"]]$spear.net.eigen.plot +
  labs(y = "Betweenness centrality")

spear.net.bet.plot[["UC"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.betweenness.df, 
                                                                             type2 %in% c("UC_CTR","UC_UC")))
spear.net.bet.plot[["UC"]]$spear.net.eigen.plot <- spear.net.bet.plot[["UC"]]$spear.net.eigen.plot +
  labs(y = "Betweenness centrality")

## -----网络中连接数目与连接强度的变化 ------------##
## -----numbers of correlations ------------##
spear.net.cornum <- lapply(c(spear.net.ctr.add, 
                             spear.net.disease.add), 
                           function(data){
                             y <- data %>% filter(!features == variable); 
                             posi.num <- nrow(subset(y, edgetrend %in% "1"));
                             nega.num <- nrow(subset(y, !edgetrend %in% "1"));
                             return(c(posi.num, nega.num))
                           })
spear.net.cornum.df <- spear.net.cornum %>%
  reduce(rbind)
spear.net.cornum.df <- data.frame(spear.net.cornum.df, 
                                  stringsAsFactors = F)
spear.net.cornum.df$type2 <- names(spear.net.cornum)
spear.net.cornum.df[1:3,"type2"] <- paste0(spear.net.cornum.df[1:3,"type2"], 
                                            "_CTR")
spear.net.cornum.df$type <- c("CTR", "CTR", "CTR",
                               "Case","Case","Case")
spear.net.cornum.df$data <- c("CRC" ,"CD", "UC",
                               "CRC" ,"CD", "UC")
names(spear.net.cornum.df)[1:2] <- c("posi.num", "nega.num")

spear.net.cornum.df$posi.basic <- c(1,1,1, 
                                    18/12,225/15,9/5)

spear.net.cornum.plot.1 <- spear.net.cornum.df %>%  
  mutate(data = factor(data, levels = c("UC","CD","CRC"))) %>%
  mutate(type = factor(type, levels = c("CTR", "Case"))) %>%
  ggplot(aes(x = data, y = posi.basic, fill = type)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#56B4E9","firebrick3")) +
  labs(x = "Disease", 
       y = "Case-to-control ratio of positive correlation numbers") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) + 
  coord_cartesian(ylim = c(0,3)) +
  scale_y_continuous(breaks = c(0,1,2,3)) +
  geom_text(aes(label = posi.num, y = posi.basic -0.1), 
            position = position_dodge(0.9), size = 5,
            color = "white")

spear.net.cornum.plot.2 <- spear.net.cornum.df %>%  
  mutate(data = factor(data, levels = c("UC","CD","CRC"))) %>%
  mutate(type = factor(type, levels = c("CTR", "Case"))) %>%
  ggplot(aes(x = data, y = posi.basic, fill = type)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#56B4E9","firebrick3")) +
  labs(x = "Disease", 
       y = "Network density") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_cartesian(ylim = c(14,15)) +
  scale_y_continuous(breaks = c(14,15))  +
  geom_text(aes(label = posi.num, y = posi.basic -0.1), 
            position = position_dodge(0.9), size = 5,
            color = "white")

## ----------正相关个数的分布-----------------##
spear.net.cornum.plot.all <- ggarrange(spear.net.cornum.plot.2,
          spear.net.cornum.plot.1,
          heights=c(1/5, 4/5),ncol = 1, 
          nrow = 2,common.legend = TRUE,
          legend="right",align = "v") 
pdf("Figure5e.pdf",
    width = 5,
    height = 5)
spear.net.cornum.plot.all
dev.off()