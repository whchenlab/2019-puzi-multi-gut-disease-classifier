## filtering out pathways with low abundances -----------##

load("pathways.Rdata")
pathabund.filter.adenoma2 <- lapply(pathabund.filter.adenoma, function(x){
  y <- noise.removal.path(x, percent = 1e-06)
  y <- y[!rownames(y) %in% c("UNMAPPED","UNINTEGRATED"),]
})
lapply(pathabund.filter.adenoma2, dim)

for (names in names(pathabund.filter.adenoma2)) {
  pathabund.filter.adenoma2[[names]]$pathway <- rownames(pathabund.filter.adenoma2[[names]])
}
pathabund.filter.adenoma2.df <- pathabund.filter.adenoma2 %>% reduce(merge, by = "pathway", all = T)
pathabund.filter.adenoma2.df[is.na(pathabund.filter.adenoma2.df)] <- 0
pathabund.filter.adenoma2.df <- pathabund.filter.adenoma2.df[, names(pathabund.filter.adenoma2.df) %in% c("pathway", rownames(data.list$data.df.t.nohdc))]
pathabund.filter.adenoma2.df.t <- data.frame(t(pathabund.filter.adenoma2.df), stringsAsFactors = F)
pathabund.filter.adenoma2.df.t["Abb",] <- unlist(lapply(strsplit(as.character(pathabund.filter.adenoma2.df.t["pathway",]), split = ":", fixed = T), function(y){x <- y[[1]]}))
names(pathabund.filter.adenoma2.df.t) <- pathabund.filter.adenoma2.df.t["Abb",]

pathabund.filter.adenoma2.df.t <- pathabund.filter.adenoma2.df.t[!rownames(pathabund.filter.adenoma2.df.t) %in% c("pathway","Abb"),]
for (i in 1:ncol(pathabund.filter.adenoma2.df.t)) {
  pathabund.filter.adenoma2.df.t[,i] <- as.numeric(pathabund.filter.adenoma2.df.t[,i])
}
pathabund.filter.adenoma2.df.t$Group <- data.df.t.nohdc[rownames(pathabund.filter.adenoma2.df.t), "Group"]

## ---- remove adenoma data------------##
pathabund.filter2 <- pathabund.filter.adenoma2
for (pro in c("PRJEB6070","PRJEB7774" ,"PRJNA447983")) {
  pathabund.filter2[[pro]] <- pathabund.filter2[[pro]][, !names(pathabund.filter2[[pro]]) %in% adenoma.metadata.3$Sample_ID]
}

pathabund.filter.grp <- lapply(pathabund.filter2, function(x){
  x.t <- data.frame(t(x), stringsAsFactors = F)
  names(x.t) <- x.t["pathway",]
  x.t["Abb",] <- unlist(lapply(strsplit(as.character(x.t["pathway",]), split = ":", fixed = T), function(y){x <- y[[1]]}))
  names(x.t) <- x.t["Abb",]
  names(x.t) <- gsub(names(x.t), pattern = "-", replacement = "_", fixed = T)
  names(x.t) <- gsub(names(x.t), pattern = "+", replacement = "_", fixed = T)
  names(x.t) <- paste0("pathways_",names(x.t))
  x.t <- x.t[!rownames(x.t) %in% c("pathway", "Abb"), ]
  return(x.t)
})
## !!!!!!!!!!!!!----------change name --------------------##
names(pathabund.filter.grp)[5] <- "PRJNA447983_cohort1"

for(pro in names(pathabund.filter.grp)){
  for (i in seq_len(ncol(pathabund.filter.grp[[pro]]))) {
    pathabund.filter.grp[[pro]][,i] <- as.numeric(pathabund.filter.grp[[pro]][,i])
  }
}
for (pro in names(pathabund.filter.grp)) {
  pathabund.filter.grp[[pro]]$Group <- metadata[rownames(pathabund.filter.grp[[pro]]), "Group"]
}

## ------Sep 07, 2019 -------------##
## ------Differential pathways ----##
pathabund.filter.grp$SRP057027 <- pathabund.filter.grp$SRP057027[rownames(pathabund.filter.grp$SRP057027) %in% rownames(data.df.t.nohdc),]


CD.study <- c("PRJNA389280", "SRP057027", "PRJNA400072")
UC.study <- c("PRJNA389280", "PRJEB1220", "PRJNA400072")

kw.function <- function(data, case){
  x <- subset(data, Group %in% c("Control", case))
  x$Group <- factor(x$Group, levels = c("Control",case))
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
  x.result <- subset(x.result, !features %in% "Group")
  return(x.result)
}

## ----- pathways with differential abundance in CRC

path.crc.kw.list <- list()
for (pro in study) {
  path.crc.kw.list[[pro]] <- kw.function(pathabund.filter.grp[[pro]], case = "CRC")
  path.crc.kw.list[[pro]]$Project <- pro
  path.crc.kw.list[[pro]] <- path.crc.kw.list[[pro]] %>% arrange(p.value) %>% 
    mutate(fdr = p.adjust(p.value, method = "fdr"))
}
path.crc.kw.fdr <- path.crc.kw.list %>% reduce(rbind)
path.crc.kw.fdr$trend <- "Control"
path.crc.kw.fdr[which(path.crc.kw.fdr$less<0.05), "trend"] <- "Case"
path.crc.kw.sum <- path.crc.kw.fdr %>% filter(fdr <0.05) %>% group_by(features, trend) %>% summarise(count = n())
path.crc.kw.CTR <- subset(path.crc.kw.sum, trend %in% "Control")
path.crc.kw.CRC <- subset(path.crc.kw.sum, trend %in% "Case")
delete.spe <- intersect(path.crc.kw.CTR$features, path.crc.kw.CRC$features)
path.crc.kw.sum.2 <- subset(path.crc.kw.sum, !features %in% delete.spe)
##rm(crc.kw.sum.2)
path.crc.kw.sum.3 <- subset(path.crc.kw.sum.2, count >1)
path.crc.kw.sum.3 <- data.frame(path.crc.kw.sum.3, stringsAsFactors = F)
path.crc.kw.sum.3$features <- as.character(path.crc.kw.sum.3$features)
path.crc.kw.sum.3$Disease <- "CRC"
crc.kw.path <- path.crc.kw.sum.3$features #####CRC  53个

## ----- pathways with differential abundance in CD
## -----correct in Sep 07, 2019 --------------##
path.cd.kw.list <- list()
for (pro in CD.study) {
  path.cd.kw.list[[pro]] <- kw.function(pathabund.filter.grp[[pro]], case = "CD")
  path.cd.kw.list[[pro]]$Project <- pro
  path.cd.kw.list[[pro]] <- path.cd.kw.list[[pro]] %>% arrange(p.value) %>% 
    mutate(fdr = p.adjust(p.value, method = "fdr"))
}
path.cd.kw.fdr <- path.cd.kw.list %>% reduce(rbind)
path.cd.kw.fdr$trend <- "Control"
path.cd.kw.fdr[which(path.cd.kw.fdr$less<0.05), "trend"] <- "Case"
path.cd.kw.sum <- path.cd.kw.fdr %>% filter(fdr <0.05) %>% 
  group_by(features, trend) %>% summarise(count = n())
path.cd.kw.CTR <- subset(path.cd.kw.sum, trend %in% "Control")
path.cd.kw.CD <- subset(path.cd.kw.sum, trend %in% "Case")
delete.spe <- intersect(path.cd.kw.CTR$features, path.cd.kw.CD$features)
##rm(crc.kw.sum.2)
path.cd.kw.sum.3 <- subset(path.cd.kw.sum, count >1)
path.cd.kw.sum.3 <- data.frame(path.cd.kw.sum.3, stringsAsFactors = F)
path.cd.kw.sum.3$features <- as.character(path.cd.kw.sum.3$features)
path.cd.kw.sum.3$Disease <- "CD"
cd.kw.path <- path.cd.kw.sum.3$features #####CD 27个

## ----- pathways with differential abundance in UC
path.uc.kw.list <- list()
for (pro in UC.study) {
  path.uc.kw.list[[pro]] <- kw.function(pathabund.filter.grp[[pro]], case = "UC")
  path.uc.kw.list[[pro]]$Project <- pro
  path.uc.kw.list[[pro]] <- path.uc.kw.list[[pro]] %>% arrange(p.value) %>% 
    mutate(fdr = p.adjust(p.value, method = "fdr"))
}
path.uc.kw.fdr <- path.uc.kw.list %>% reduce(rbind)
path.uc.kw.fdr$trend <- "Control"
path.uc.kw.fdr[which(path.uc.kw.fdr$less<0.05), "trend"] <- "Case"
path.uc.kw.sum <- path.uc.kw.fdr %>% filter(fdr <0.05) %>% 
  group_by(features, trend) %>% summarise(count = n())
path.uc.kw.CTR <- subset(path.uc.kw.sum, trend %in% "Control")
path.uc.kw.UC <- subset(path.uc.kw.sum, trend %in% "Case")
delete.spe <- intersect(path.uc.kw.CTR$features, path.uc.kw.UC$features)
path.uc.kw.sum.2 <- subset(path.uc.kw.sum, !features %in% delete.spe)

path.uc.kw.sum.3 <- subset(path.uc.kw.sum.2, count >1)
path.uc.kw.sum.3 <- data.frame(path.uc.kw.sum.3, stringsAsFactors = F)
path.uc.kw.sum.3$features <- as.character(path.uc.kw.sum.3$features)
path.uc.kw.sum.3$Disease <- "UC"
uc.kw.path <- path.uc.kw.sum.3$features #####UC 2个

## -----------HDC correlated pathways -----------##
pathabund.filter.HDC <- pathabund.filter.grp
for (pro in names(pathabund.filter.HDC)) {
  pathabund.filter.HDC[[pro]] <- cbind(metadata[rownames(pathabund.filter.HDC[[pro]]), "HDC"],
                                       pathabund.filter.HDC[[pro]])
  pathabund.filter.HDC[[pro]] <- data.frame(pathabund.filter.HDC[[pro]] , stringsAsFactors = F)
  names(pathabund.filter.HDC[[pro]])[1] <- "HDC"
}
pathabund.filter.HDC[[5]] <- NULL

for (pro in c("PRJEB6070","PRJEB7774" ,"PRJNA447983")) {
  x.t <- data.frame(t(pathabund.filter.adenoma2[[pro]]), stringsAsFactors = F);
  names(x.t) <- x.t["pathway",]
  x.t["Abb",] <- unlist(lapply(strsplit(as.character(x.t["pathway",]), split = ":", fixed = T), function(y){x <- y[[1]]}))
  names(x.t) <- x.t["Abb",]
  names(x.t) <- gsub(names(x.t), pattern = "-", replacement = "_", fixed = T)
  names(x.t) <- gsub(names(x.t), pattern = "+", replacement = "_", fixed = T)
  names(x.t) <- paste0("pathways_",names(x.t))
  x.t <- x.t[!rownames(x.t) %in% c("pathway", "Abb"), ]
  pathabund.filter.HDC[[pro]] <- x.t
}

crc.meta.adeno[which(crc.meta.adeno$Group %in% "CTR"),"Group"] <- "Control"

for (pro in c("PRJEB6070","PRJEB7774" ,"PRJNA447983")) {
  pathabund.filter.HDC[[pro]] <- cbind(crc.meta.adeno[rownames(pathabund.filter.HDC[[pro]]), c("HDC","Group")],
                                       pathabund.filter.HDC[[pro]])
  pathabund.filter.HDC[[pro]] <- data.frame(pathabund.filter.HDC[[pro]] , stringsAsFactors = F)
  names(pathabund.filter.HDC[[pro]])[1:2] <- c("HDC","Group")
}
names(pathabund.filter.HDC)[11] <- "PRJNA447983_cohort1"

pathabund.filter.HDC.result <- list()
cor.func2 <- function(data, case, names){
  x <- subset(data, Group %in% c("Control", case))
  x.2 <- x[,! names(x) %in% "Group"]
  x.cor.result <- cor.func(data = x.2, type = case, project = names)
  return(x.cor.result)
}
for (pro in study) {
  pathabund.filter.HDC.result[[pro]] <- cor.func2(pathabund.filter.HDC[[pro]], 
                                                  case = "CRC", names = pro)
}
for (pro in CD.study) {
  pathabund.filter.HDC.result[[paste0(pro, "_CD")]] <- cor.func2(pathabund.filter.HDC[[pro]], 
                                                                 case = "CD", names = pro)
}
for (pro in UC.study) {
  pathabund.filter.HDC.result[[paste0(pro, "_UC")]] <- cor.func2(pathabund.filter.HDC[[pro]], 
                                                                 case = "UC", names = pro)
}
path.cor.data <- lapply(pathabund.filter.HDC.result, function(data){y <- data$data.cor.sum}) %>% reduce(rbind)
path.cor.data.sig <- subset(path.cor.data, P_value<0.05)
path.cor.data.posi <- subset(path.cor.data.sig, coef>0)
path.cor.data.nega <- subset(path.cor.data.sig, coef<0)
delete.spe <- intersect(path.cor.data.posi$Species, path.cor.data.nega$Species)
path.cor.data.sig <- subset(path.cor.data.sig, !Species %in% delete.spe)

path.cor.data.sig.sum <- path.cor.data.sig %>% group_by(Species, Type) %>% summarise(count= n())
path.cor.data.sig.sum.2 <- dcast(path.cor.data.sig.sum, Species ~ Type)
sig.path <- subset(path.cor.data.sig.sum.2, CRC>1 | UC>1 | CD>1)$Species
sig.path.df <- subset(path.cor.data.sig.sum.2, CRC>1 | UC>1 | CD>1)
