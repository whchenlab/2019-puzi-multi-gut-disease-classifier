kw.path.spearman.pval <- list()
kw.path.spearman.cor <- list()

all.path.spe.df <- cbind(log10(path.data.list$pathabund.df[,1:346]+1e-9),
                         log10(data.list$data.df.t.kw.nohdc[rownames(path.data.list$pathabund.df),kw.spe]+1e-6))

for (pro in setdiff(unique(metadata$Project), 
                    c("PRJNA389280","PRJNA400072"))) {
  
  metadata.pro <- subset(metadata, Project %in% pro) %>% 
    pull(Sample_ID)
  
  a <- rcorr(as.matrix(all.path.spe.df[metadata.pro,]), 
             type = "spearman")
  a.triu <- a$r[-c(1:346),c(1:346)]
  a.triu.df <- data.frame(a.triu, stringsAsFactors = F)
  a.triu.df$features <- rownames(a.triu.df) 
  a.triu.df <- melt(a.triu.df, id.vars = "features")
  kw.path.spearman.cor[[pro]] <- a.triu.df
  
  a.triu <- a$P[-c(1:346),c(1:346)]
  a.triu.df <- data.frame(a.triu, stringsAsFactors = F)
  a.triu.df$features <- rownames(a.triu.df) 
  a.triu.df <- melt(a.triu.df, id.vars = "features")
  kw.path.spearman.pval[[pro]] <- a.triu.df
}
for (pro in c("PRJNA389280","PRJNA400072")) {
  for (grp in c("UC","CD")) {
    metadata.pro <- subset(metadata, Project %in% pro & Group %in% c("Control", grp)) %>% 
      pull(Sample_ID)
    
    a <- rcorr(as.matrix(all.path.spe.df[metadata.pro,]), 
               type = "spearman")
    a.triu <- a$r[-c(1:346),c(1:346)]
    a.triu.df <- data.frame(a.triu, stringsAsFactors = F)
    a.triu.df$features <- rownames(a.triu.df) 
    a.triu.df <- melt(a.triu.df, id.vars = "features")
    kw.path.spearman.cor[[paste(pro, grp , sep = "_")]] <- a.triu.df
    
    a.triu <- a$P[-c(1:346),c(1:346)]
    a.triu.df <- data.frame(a.triu, stringsAsFactors = F)
    a.triu.df$features <- rownames(a.triu.df) 
    a.triu.df <- melt(a.triu.df, id.vars = "features")
    kw.path.spearman.pval[[paste(pro, grp , sep = "_")]] <- a.triu.df
  }
}


kw.path.spearman.cor.filter <- lapply(kw.path.spearman.cor, function(data){
  y <- data %>% filter(!is.na(value))
  return(y)
})
kw.path.spearman.pval.filter <- lapply(kw.path.spearman.pval, function(data){
  y <- data %>%  filter(!is.na(value))
  return(y)
})

for (pro in names(kw.path.spearman.cor.filter)) {
  kw.path.spearman.cor.filter[[pro]]$pro <- pro
  kw.path.spearman.cor.filter[[pro]]$combine <- paste(kw.path.spearman.cor.filter[[pro]]$features, 
                                                      kw.path.spearman.cor.filter[[pro]]$variable,
                                                             sep = "_")
  kw.path.spearman.pval.filter[[pro]]$pro <- pro
  kw.path.spearman.pval.filter[[pro]]$combine <- paste(kw.path.spearman.pval.filter[[pro]]$features, 
                                                       kw.path.spearman.pval.filter[[pro]]$variable,
                                                              sep = "_")
}

for (pro in names(kw.path.spearman.cor.filter)) {
  if (pro %in% study) {
    kw.path.spearman.cor.filter[[pro]]$disease <- "CRC"
    kw.path.spearman.pval.filter[[pro]]$disease <- "CRC"
  }else if(pro %in% c("PRJEB1220","PRJNA389280_UC","PRJNA400072_UC")){
    kw.path.spearman.cor.filter[[pro]]$disease <- "UC"
    kw.path.spearman.pval.filter[[pro]]$disease <- "UC"
  }else{
    kw.path.spearman.cor.filter[[pro]]$disease <- "CD"
    kw.path.spearman.pval.filter[[pro]]$disease <- "CD"
  }
}

kw.path.spearman.pval.filter.df <- kw.path.spearman.pval.filter %>% 
  reduce(rbind) %>% group_by(pro) %>% arrange(value) %>% 
  mutate(fdr = p.adjust(value, method = "fdr"))
kw.path.spearman.pval.final <- kw.path.spearman.pval.filter.df %>% 
  filter(fdr <0.05)
names(kw.path.spearman.pval.final)[3] <- "pvalue"

kw.path.spearman.cor.filter.df <- kw.path.spearman.cor.filter %>% 
  reduce(rbind)

kw.path.spearman.cor.final.merge <- merge(kw.path.spearman.pval.final, 
                                          kw.path.spearman.cor.filter.df,
                                          by = c("features", "variable",
                                                             "combine","disease",
                                                             "pro"),
                                          all.x = T)
kw.path.spearman.cor.final.merge.2 <- kw.path.spearman.cor.final.merge %>% 
  filter(!is.na(value))

save(kw.path.spearman.cor.final.merge.2,
     file = "permer_kw_cor_path.Rdata")

spe.cor.pathways.meta <- function(data, metadata){
  ## ----correlated species and pathways ------##
  ## ----meta-analysis of correlation ---------##
  ## ----data: contains features, targets, combine, disease, r^2, fdr, "project", -----##
  ## ----metadata: contains nSamples, project, 
  merge.data <- merge(data, metadata,
                      by = "project",
                      all.x = T)
  merge.data$nSamples <- as.numeric(merge.data$nSamples)
  
  data.cor.correct <- list()
  for (comp in unique(merge.data$combine)) {
    degra.sub <- merge.data %>% filter(combine %in% comp)
    if (nrow(degra.sub)>1) {
      degra.sub.zscore <- metacor(value, nSamples, project, 
                                  data = degra.sub, sm = "ZCOR")
      degra.sub.correct.r <- (exp(2* degra.sub.zscore$TE.random)-1)/(exp(2* degra.sub.zscore$TE.random)+1)
      
      data.cor.correct[[comp]]$meta_z <- degra.sub.zscore
      data.cor.correct[[comp]]$correct_r <- degra.sub.correct.r
    }
  }
  
  
  ## ------------corrected coefficient ----------##
  data.cor.correct.r.df <- lapply(data.cor.correct,
                                  function(data){
                                    y <- data$correct_r;
                                    return(y)})  %>% reduce(rbind)
  
  data.cor.correct.r.df <- data.frame(data.cor.correct.r.df,
                                                         stringsAsFactors = F)
  names(data.cor.correct.r.df)[1] <- "random.r"
  data.cor.correct.r.df$target <- names(data.cor.correct)
  
  data.merge.3 <- merge(data,
                        data.cor.correct.r.df,
                        by.x = c("combine"),
                        by.y= c("target"),
                        all.y= T) %>% distinct(combine, .keep_all = T)
  
  data.merge.dcast <- dcast(data.merge.3[,c("features","variable","random.r")],
                            features ~ variable)
  data.merge.dcast[is.na(data.merge.dcast)] <- "0"
  data.merge.dcast[,-1] <- apply(data.merge.dcast[,-1] ,2,as.numeric)
  rownames(data.merge.dcast) <- data.merge.dcast$features
  
  return(list(data.merge.dcast = data.merge.dcast,
    data.cor.correct.r.df = data.cor.correct.r.df,
    data.cor.correct = data.cor.correct,
    merge.data = merge.data))
}
metadata.3 <- data.frame(project = unique(kw.path.spearman.cor.final.merge.2$pro),
                         nSamples = c(104,107,110,80, 307,149,128,
                                      144, 132, 82, 216, 53, 157))
                                      
kw.path.spearman.cor.result <- spe.cor.pathways.meta(kw.path.spearman.cor.final.merge.2,
                                                     metadata.3)

kw.path.spearman.cor.final.dcast <- kw.path.spearman.cor.result$data.merge.dcast[,-1]
## ------¼ÆËãpvalue ------------------##
kw.path.spearman.cor.pval <- data.frame(unlist(lapply(kw.path.spearman.cor.result$data.cor.correct,
                                                      function(data){y <- data$meta_z$pval.random})), 
                                        stringsAsFactors = F)
kw.path.spearman.cor.pval$target <- names(kw.path.spearman.cor.result$data.cor.correct)
names(kw.path.spearman.cor.pval)[1] <- "correct.pval"
kw.path.spearman.cor.pval$features <- unlist(lapply(strsplit(kw.path.spearman.cor.pval$target,
                                                             split="_pathways",fixed = T),
                                                    function(data){
                                                      y<- data[[1]]}))
kw.path.spearman.cor.pval$variable <- NA
for (i in seq_len(nrow(kw.path.spearman.cor.pval))) {
  kw.path.spearman.cor.pval[i,"variable"] <- gsub(kw.path.spearman.cor.pval[i,"target"],
                                                  pattern = paste0(kw.path.spearman.cor.pval[i,"features"], "_"),
                                                  replacement = "",fixed = T)}

kw.path.spearman.cor.pval.dcast <- dcast(kw.path.spearman.cor.pval[,
                                                                   c("features","variable","correct.pval")],
                                         features ~ variable)
kw.path.spearman.cor.pval.dcast[is.na(kw.path.spearman.cor.pval.dcast)] <- 0
kw.path.spearman.cor.pval.dcast.2 <- apply(kw.path.spearman.cor.pval.dcast[,-1], 2, function(data){
  y <- cut(data, breaks = c(0,0.01,1),
           include.lowest = F,
           labels=c("*"," "));
  return(y)})
kw.path.spearman.cor.pval.dcast.2[is.na(kw.path.spearman.cor.pval.dcast.2)] <- " "
rownames(kw.path.spearman.cor.pval.dcast.2) <- kw.path.spearman.cor.pval.dcast$features

## -----------
## -----------只取起码在3个数据集中原来是显著关联的关联-----------------------##
filter.combine <- names(table(kw.path.spearman.cor.final.merge.2$combine)[table(kw.path.spearman.cor.final.merge.2$combine)>2])

kw.path.spearman.cor.filter.dcast <- merge(kw.path.spearman.cor.final.merge.2,
                                           kw.path.spearman.cor.result$data.cor.correct.r.df,
                                           by.x = c("combine"),
                                           by.y= c("target"),
                                           all.y= T) %>% distinct(combine, .keep_all = T)
kw.path.spearman.cor.filter.dcast <- kw.path.spearman.cor.filter.dcast %>%
  filter(combine %in% filter.combine)
kw.path.spearman.cor.filter.dcast.2 <- dcast(kw.path.spearman.cor.filter.dcast[,c("features","variable","random.r")],
                                             features ~ variable)
rownames(kw.path.spearman.cor.filter.dcast.2) <- kw.path.spearman.cor.filter.dcast.2$features
kw.path.spearman.cor.filter.dcast.2 <- kw.path.spearman.cor.filter.dcast.2[,-1]

kw.path.spearman.cor.pval.filter.dcast.2 <- kw.path.spearman.cor.pval.dcast.2[rownames(kw.path.spearman.cor.filter.dcast.2),
                                                                              names(kw.path.spearman.cor.filter.dcast.2)]
## -----
## ---------------------------------
## -----

for (i in 1:nrow(kw.path.spearman.cor.filter.dcast.2)) {
  for (j in 1:ncol(kw.path.spearman.cor.filter.dcast.2)) {
    if (is.na(kw.path.spearman.cor.filter.dcast.2[i,j])) {
      kw.path.spearman.cor.pval.filter.dcast.2[i,j] <- " "
    }
  }
}

kw.path.spearman.cor.filter.dcast.2[is.na(kw.path.spearman.cor.filter.dcast.2)] <- 0

## -----------½öÌáÈ¡degradationÊý¾Ý
degradation.df <- read.delim("degradation_path.txt",
                             header = T, sep = "\t", as.is = T)
rownames(degradation.df) <- degradation.df$Abb.in.df

kw.path.spearman.cor.degrad.dcast <- kw.path.spearman.cor.filter.dcast.2[, 
                                                                         names(kw.path.spearman.cor.filter.dcast.2) %in% rownames(degradation.df)]
kw.path.spearman.cor.degrad.pval.dcast <- kw.path.spearman.cor.pval.filter.dcast.2[, 
                                                                                   dimnames(kw.path.spearman.cor.pval.filter.dcast.2)[[2]] %in% rownames(degradation.df)]

## -----------
dimnames(kw.path.spearman.cor.degrad.dcast)[[2]] <- degradation.df[names(kw.path.spearman.cor.degrad.dcast),
                                                                         "description2"]
dimnames(kw.path.spearman.cor.degrad.pval.dcast)[[2]] <- degradation.df[dimnames(kw.path.spearman.cor.degrad.pval.dcast)[[2]],
                                                                              "description2"]

unqualified.spe <- names(apply(kw.path.spearman.cor.degrad.pval.dcast, 1, 
                               function(data){y <- sum(data =="*")})[apply(kw.path.spearman.cor.degrad.pval.dcast, 1, 
                                                                           function(data){y <- sum(data =="*")})==0])

kw.path.spearman.cor.degrad.dcast <- kw.path.spearman.cor.degrad.dcast[!rownames(kw.path.spearman.cor.degrad.dcast) %in% unqualified.spe, ]
kw.path.spearman.cor.degrad.pval.dcast <- kw.path.spearman.cor.degrad.pval.dcast[!rownames(kw.path.spearman.cor.degrad.pval.dcast) %in% unqualified.spe, ]

## cluster setting -------##
## ------Species setting ---------##
## ------First, control / case trend  -----##
kw.sum.df <- read.delim("kwspe_sum.txt",
                        header = T, sep = "\t", as.is = T)
kw.spe.annot <- kw.sum.df %>% 
  distinct(features, .keep_all = T) %>% 
  select(features, trend)
rownames(kw.spe.annot) <- kw.spe.annot$features

kw.spe.annot$type <- "CD"

kw.spe.annot[unique(subset(kw.sum.df, Disease %in% "CRC")$features),
                    "type"] <- "CRC"
kw.spe.annot[unique(subset(kw.sum.df, Disease %in% "UC")$features),
             "type"] <- "UC"
kw.spe.annot[names(table(kw.sum.df$features)[table(kw.sum.df$features) >1]),
             "type"] <- "common"
kw.spe.annot <- kw.spe.annot[,c("trend","type")]


## ---clusters setting --------##
## ---´ý¶¨
kw.spe.annot$clusters <- "Control"
kw.spe.annot[c("Escherichia_coli",
               "Escherichia_unclassified",
               "Klebsiella_pneumoniae"),
             "clusters"] <- "CD"
kw.spe.annot[c("Bacteroides_coprocola" ,
               "Paraprevotella_xylaniphila",
               "Paraprevotella_unclassified"),
             "clusters"] <- "UC"
kw.spe.annot[c("Porphyromonas_asaccharolytica",
               "Porphyromonas_somerae",
               "Porphyromonas_uenonis"),
             "clusters"] <- "CRC"
kw.spe.annot[c("Bifidobacterium_dentium","Clostridium_clostridioforme",
               "Clostridium_symbiosum", "Lachnospiraceae_bacterium_7_1_58FAA",
               "Clostridium_hathewayi", "Coprobacillus_unclassified"),
             "clusters"] <- "Mixed I"

kw.spe.annot[intersect(rownames(subset(kw.spe.annot, clusters %in% "Control")), subset(kw.sum.df, trend %in% "Case")$features),
             "clusters"] <- "Mixed II"
write.table(kw.spe.annot, 
            "kw_spe_annot.txt",
            col.names = T, row.names = T, 
            sep = "\t", quote = F)
## -----ÔÚexcelÖØÅÅÐò -----------##
kw.spe.annot <- read.delim("kw_spe_annot.txt",
                           header = T, sep = "\t", as.is = T)
rownames(kw.spe.annot) <- kw.spe.annot$Spe
kw.spe.annot.2 <- kw.spe.annot[,c("type","trend","clusters")]


## ------Pathways setting ---------##
## ------First, control / case trend  -----##
kw.path.sum.df <- read.delim("kwpath_sum.txt",
                             header = T, sep = "\t", as.is = T)
kw.path.annot <- kw.path.sum.df %>% 
  distinct(features, .keep_all = T) %>% 
  select(features, trend)
rownames(kw.path.annot) <- kw.path.annot$features

degrad.annot <- degradation.df[,c("majorclass", "type")]
degrad.annot$trend <- "Not markers"
degrad.annot[c("pathways_P162_PWY",
               "pathways_PWY0_1297",
               "pathways_FUC_RHAMCAT_PWY"), "trend"] <- "Case"
degrad.annot[subset(kw.path.annot, features %in% degradation.df$Abb.in.df & trend %in% "Control")$features, "trend"] <- "Control"
degrad.annot <- degrad.annot[,c("majorclass", "trend")]

rownames(degrad.annot) <- degradation.df[rownames(degrad.annot), "description2"]

## color setting -----##
degradation.annot.color <- list(majorclass = c(`Amine and Polyamine Degradation`="#999999",
                                               `Amino Acid Degradation`="#E69F00",
                                               `Sugar Degradation`="#56B4E9",
                                               `Sugar Derivative Degradation`="#009E73",
                                               `Polysaccharide Degradation`="#F0E442",
                                               `Alcohol Degradation`="#0072B2",
                                               `Aromatic Compound Degradation`="#D55E00",
                                               `Aldehyde Degradation`="#CC79A7",
                                               `Nucleoside and Nucleotide Degradation`="#e85a71",
                                               `Carboxylate Degradation`="#56A902",
                                               `Generation of Precursor Metabolite and Energy`= "lawngreen",
                                               `Inorganic Nutrient Metabolism`="#FFFF40FF"),
                                clusters = c(CD= "#869A00",
                                             UC= "#A9657A",
                                             CRC= "#A7C2D4",
                                             `Mixed I`= "#5F615B", 
                                             `Mixed II` = "#dadad8", 
                                             Control = "black"),
                                trend = c(Control ="#56B4E9",
                                          Case = "firebrick3",
                                          `Not markers` = "black"),
                                type = c(UC = "#EF6C9C",
                                         common = "#3D4756",
                                         CRC = "#A1ACBD",
                                         CD = "#00B053"))



bk.kw.degra <- c(seq(-1,-0.01,by=0.01),seq(0,1,by=0.01))

## -------½«Êý¾Ý°´ÕÕspecies orderÅÅÐò -----------##
kw.spe.order <- rownames(kw.spe.annot.2[rownames(kw.spe.annot.2) %in% rownames(kw.path.spearman.cor.degrad.dcast), ])

kw.path.spearman.cor.degrad.dcast <- kw.path.spearman.cor.degrad.dcast[kw.spe.order, ]
kw.path.spearman.cor.degrad.pval.dcast <- kw.path.spearman.cor.degrad.pval.dcast[kw.spe.order, ]

kw.degradation.spearman.plot <- pheatmap(t(as.matrix(kw.path.spearman.cor.degrad.dcast)),
                                         cluster_rows = T,
                                         cluster_cols = F, 
                                         annotation_colors = degradation.annot.color,
                                         annotation_col = kw.spe.annot.2,
                                         annotation_row = degrad.annot,
                                         clustering_method = "complete",
                                         color = c(colorRampPalette(colors = c("#016392","white"))(length(bk.kw.degra)/2),
                                                   colorRampPalette(colors = c("white","#c72e29"))(length(bk.kw.degra)/2)),
                                         legend_breaks=seq(-1,1,0.2),
                                         breaks=bk.kw.degra,
                                         display_numbers = t(kw.path.spearman.cor.degrad.pval.dcast),
                                         fontsize = 10, gaps_col = c(3,6,7,11, 44,45),
                                         cellwidth = 9, cellheight =9 , border=FALSE)

## -----------biosynthesis pathways -----------------------##
kw.path.spearman.cor.other.dcast <- kw.path.spearman.cor.filter.dcast.2[, 
                                                                        !names(kw.path.spearman.cor.filter.dcast.2) %in% rownames(degradation.df)]
kw.path.spearman.cor.other.pval.dcast <- kw.path.spearman.cor.pval.filter.dcast.2[, 
                                                                                  !dimnames(kw.path.spearman.cor.pval.filter.dcast.2)[[2]] %in% rownames(degradation.df)]

unqualified.spe.2 <- names(apply(kw.path.spearman.cor.other.pval.dcast, 1, function(data){y <- sum(data =="*")})[apply(kw.path.spearman.cor.other.pval.dcast, 1, function(data){y <- sum(data =="*")})==0])

kw.path.spearman.cor.other.dcast <- kw.path.spearman.cor.other.dcast[!rownames(kw.path.spearman.cor.other.dcast) %in% unqualified.spe.2, ]
kw.path.spearman.cor.other.pval.dcast <- kw.path.spearman.cor.other.pval.dcast[!rownames(kw.path.spearman.cor.other.pval.dcast) %in% unqualified.spe.2, ]

kw.spe.order.other <- rownames(kw.spe.annot.2[rownames(kw.spe.annot.2) %in% rownames(kw.path.spearman.cor.other.dcast), ])

kw.path.spearman.cor.other.dcast <- kw.path.spearman.cor.other.dcast[kw.spe.order.other, ]
kw.path.spearman.cor.other.pval.dcast <- kw.path.spearman.cor.other.pval.dcast[kw.spe.order.other, ]

biosyn.annot <- data.frame(Pathways = names(kw.path.spearman.cor.other.dcast), 
                           trend = "Not markers", 
                           stringsAsFactors = F)
rownames(biosyn.annot) <- biosyn.annot$Pathways
biosyn.annot[unique(subset(kw.path.annot, trend %in% "Control")$features), "trend"] <- "Control" 
biosyn.annot[unique(subset(kw.path.annot, trend %in% "Case")$features), "trend"] <- "Case" 

biosyn.annot.2 <- data.frame(trend = biosyn.annot$trend, 
                             stringsAsFactors = F)
rownames(biosyn.annot.2) <- rownames(biosyn.annot)
rownames(biosyn.annot.2) <- gsub(rownames(biosyn.annot.2), 
                                 pattern = "pathways_",
                                 replacement = "",
                                 fixed = T)
rownames(biosyn.annot.2) <- gsub(rownames(biosyn.annot.2), 
                                 pattern = "_",
                                 replacement = "-",
                                 fixed = T)

names(kw.path.spearman.cor.other.dcast) <- gsub(names(kw.path.spearman.cor.other.dcast), 
                                                pattern = "pathways_",
                                                replacement = "",
                                                fixed = T)
names(kw.path.spearman.cor.other.dcast) <- gsub(names(kw.path.spearman.cor.other.dcast), 
                                 pattern = "_",
                                 replacement = "-",
                                 fixed = T)

dimnames(kw.path.spearman.cor.other.pval.dcast)[[2]] <- gsub(dimnames(kw.path.spearman.cor.other.pval.dcast)[[2]], 
                                                pattern = "pathways_",
                                                replacement = "",
                                                fixed = T)
dimnames(kw.path.spearman.cor.other.pval.dcast)[[2]] <- gsub(dimnames(kw.path.spearman.cor.other.pval.dcast)[[2]], 
                                                pattern = "_",
                                                replacement = "-",
                                                fixed = T)


kw.biosyn.spearman.plot <- pheatmap(as.matrix(kw.path.spearman.cor.other.dcast),
         cluster_rows = F,
         cluster_cols = T, 
         annotation_colors = degradation.annot.color,
         annotation_row = kw.spe.annot.2,
         annotation_col = biosyn.annot.2,
         clustering_method = "complete",
         color = c(colorRampPalette(colors = c("#016392","white"))(length(bk.kw.degra)/2),
                   colorRampPalette(colors = c("white","#c72e29"))(length(bk.kw.degra)/2)),
         legend_breaks=seq(-1,1,0.2),
         breaks=bk.kw.degra,
         display_numbers = kw.path.spearman.cor.other.pval.dcast,
         gaps_row = c(3,6,9,13, 49,50),border=FALSE,
         fontsize = 8, cellwidth = 7, cellheight = 8)

save.image("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result/Figure4_data.RData")
