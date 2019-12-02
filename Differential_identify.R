crc.kw <- list()

for (project in names(crc.abund.filter.list)) {
  y <- crc.abund.filter.list[[project]]
  x <- data.frame(t(y), stringsAsFactors = F)
  print(rownames(x))
  x.meta <- subset(metadata, Sample_ID %in% rownames(x))
  rownames(x.meta) <- x.meta$Sample_ID
  x <- x[rownames(x.meta),]
  x$Group <- x.meta$Group
  x$Group <- factor(x$Group, levels = c("CTR","CRC"))
  x.result <- data.frame(features = names(x),
                         p.value = NA,
                         greater= NA,
                         less = NA)
  features <- names(x)
  rownames(x.result) <- features
  features <- features[!features %in% "Group"]
  for (i in features) {
    a <- wilcox.test(x[,i] ~ x$Group, x)
    x.result[i, "p.value"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "greater")
    x.result[i, "greater"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "less")
    x.result[i, "less"] <- a$p.value
  }
  x.result$Project <- project
  crc.kw[[project]] <- x.result
}

crc.kw.fdr <- crc.kw %>% reduce(rbind) %>% filter(!features %in% "Group")
crc.kw.fdr <- crc.kw.fdr %>% group_by(Project) %>% arrange(p.value) %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr"))

crc.kw.fdr$trend <- "Control"
crc.kw.fdr[which(crc.kw.fdr$less<0.05), "trend"] <- "Case"

crc.kw.sum <- crc.kw.fdr %>% filter(fdr <0.05) %>% group_by(features, trend) %>% summarise(count = n())
crc.kw.CTR <- subset(crc.kw.sum, trend %in% "Control")
crc.kw.CRC <- subset(crc.kw.sum, trend %in% "Case")
delete.spe <- intersect(crc.kw.CTR$features, crc.kw.CRC$features)
##rm(crc.kw.sum.2)
crc.kw.sum.3 <- subset(crc.kw.sum, count >1)
crc.kw.sum.3 <- data.frame(crc.kw.sum.3, stringsAsFactors = F)
crc.kw.sum.3$features <- as.character(crc.kw.sum.3$features)
crc.kw.sum.3$Disease <- "CRC"
crc.kw.spe <- crc.kw.sum.3$features 

CD.kw.list <- list()
for (project in names(cd.abund.filter.list)) {
  y <- cd.abund.filter.list[[project]]
  x <- data.frame(t(y), stringsAsFactors = F)
  print(rownames(x))
  x.meta <- subset(metadata, Sample_ID %in% rownames(x))
  rownames(x.meta) <- x.meta$Sample_ID
  x <- x[rownames(x.meta),]
  x$Group <- x.meta$Group
  x$Group <- factor(x$Group, levels = c("CTR","CD"))
  x.result <- data.frame(features = names(x),
                         p.value = NA,
                         greater= NA,
                         less = NA)
  features <- names(x)
  rownames(x.result) <- features
  features <- features[!features %in% "Group"]
  for (i in features) {
    a <- wilcox.test(x[,i] ~ x$Group, x)
    x.result[i, "p.value"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "greater")
    x.result[i, "greater"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "less")
    x.result[i, "less"] <- a$p.value
  }
  x.result$Project <- project
  CD.kw.list[[project]] <- x.result
}

## ------CD-altered species identification -----------##
CD.kw.list.fdr <- CD.kw.list %>% reduce(rbind) %>% filter(!features %in% "Group")
CD.kw.list.fdr <- CD.kw.list.fdr %>% group_by(Project) %>% arrange(p.value) %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr"))

CD.kw.list.fdr$trend <- "Control"
CD.kw.list.fdr[which(CD.kw.list.fdr$less<0.05), "trend"] <- "Case"

CD.kw.sum <- CD.kw.list.fdr %>% filter(fdr <0.05) %>% group_by(features, trend) %>% summarise(count = n())
CD.kw.CTR <- subset(CD.kw.sum, trend %in% "Control")
CD.kw.CD <- subset(CD.kw.sum, trend %in% "Case")
delete.spe <- intersect(CD.kw.CTR$features, CD.kw.CD$features)
CD.kw.sum.2 <- subset(CD.kw.sum, !features %in% delete.spe)

CD.kw.sum.3 <- subset(CD.kw.sum.2, count >1)
CD.kw.sum.3 <- data.frame(CD.kw.sum.3, stringsAsFactors = F)
CD.kw.sum.3$features <- as.character(CD.kw.sum.3$features)
CD.kw.sum.3$Disease <- "CD"
CD.kw.spe <- CD.kw.sum.3$features 

UC.kw.list <- list()
for (project in names(uc.abund.filter.list)) {
  y <- uc.abund.filter.list[[project]]
  x <- data.frame(t(y), stringsAsFactors = F)
  print(rownames(x))
  x.meta <- subset(metadata, Sample_ID %in% rownames(x))
  rownames(x.meta) <- x.meta$Sample_ID
  x <- x[rownames(x.meta),]
  x$Group <- x.meta$Group
  x$Group <- factor(x$Group, levels = c("CTR","UC"))
  x.result <- data.frame(features = names(x),
                         p.value = NA,
                         greater= NA,
                         less = NA)
  features <- names(x)
  rownames(x.result) <- features
  features <- features[!features %in% "Group"]
  for (i in features) {
    a <- wilcox.test(x[,i] ~ x$Group, x)
    x.result[i, "p.value"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "greater")
    x.result[i, "greater"] <- a$p.value
    a <- wilcox.test(x[,i] ~ x$Group, alternative = "less")
    x.result[i, "less"] <- a$p.value
  }
  x.result$Project <- project
  UC.kw.list[[project]] <- x.result
}


## ------UC-altered species identification -----------##
UC.kw.list.fdr <- UC.kw.list  %>% reduce(rbind) %>% filter(!features %in% "Group")
UC.kw.list.fdr <- UC.kw.list.fdr %>% group_by(Project) %>% arrange(p.value) %>% 
  mutate(fdr = p.adjust(p.value, method = "fdr"))

UC.kw.list.fdr$trend <- "Control"
UC.kw.list.fdr[which(UC.kw.list.fdr$less<0.05), "trend"] <- "Case"

UC.kw.sum <- UC.kw.list.fdr %>% filter(fdr <0.05) %>% group_by(features, trend) %>% summarise(count = n())
UC.kw.CTR <- subset(UC.kw.sum, trend %in% "Control")
UC.kw.UC <- subset(UC.kw.sum, trend %in% "Case")
delete.spe <- intersect(UC.kw.CTR$features, UC.kw.UC $features)
UC.kw.sum.2 <- subset(UC.kw.sum, !features %in% delete.spe)

UC.kw.sum.3 <- subset(UC.kw.sum.2, count >1)
UC.kw.sum.3 <- data.frame(UC.kw.sum.3, stringsAsFactors = F)
UC.kw.sum.3$features <- as.character(UC.kw.sum.3$features)
UC.kw.sum.3$Disease <- "UC"
UC.kw.spe <- UC.kw.sum.3$features 

save(UC.kw.sum.3, CD.kw.sum.3, crc.kw.sum.3, 
     UC.kw.list, CD.kw.list, crc.kw,
     file = "Differential_Species.RData")
