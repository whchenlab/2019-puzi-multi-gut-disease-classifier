## species-pathways correlation
install.packages("PResiduals", dependencies = T)
library(PResiduals)
load("Filtered_Pathways_abund.Rdata")

kw.path.spe.df <- cbind(log10(kw.spe.list$four.gp[, !(names(kw.spe.list$four.gp) %in% "Group")] + 1e-6),
                        log10(kw.path.list$four.gp[rownames(kw.spe.list$four.gp), !names(kw.path.list$four.gp) %in% "Group"] + 1e-9))
kw.path.spe.df[, c("Age", "Gender", "BMI")] <- metadata[rownames(kw.path.spe.df), c("Age", "Gender", "BMI")]

## ----------------------------------
adjusted.spearman.func <- function(metadata = metadata, 
                                   kw.path.spe.df = kw.path.spe.df,
                                   project = "PRJEB6070", 
                                   case = "CRC",
                                   covariates = c("Age")){
  library(dplyr)
  library(PResiduals) 
  samples <- metadata %>% filter(Project %in% project) %>%
    filter(Group %in% c("CTR", case)) %>%
    pull(Sample_ID)
  a <- kw.path.spe.df[samples, ]
  feature.name <- colnames(kw.path.spe.df)[1:(ncol(kw.path.spe.df)-3)]
  
  kw.path.spearman.cor <- matrix(data = NA,nrow = length(feature.name), 
                                 ncol = length(feature.name), 
                                 dimnames = list(feature.name, 
                                                 feature.name))
  
  kw.path.spearman.pval <- matrix(data = NA,nrow = length(feature.name), 
                                  ncol = length(feature.name), 
                                  dimnames = list(feature.name, 
                                                  feature.name))
  
  for (i in 1:length(feature.name)) {
    for (j in 1:length(feature.name)) {
      a.2 <- a[, c(feature.name[i], feature.name[j], covariates)]
      names(a.2) <- c("x", "y", 
                      "Age")
      if (length(unique(a.2[,1])) < 3 | length(unique(a.2[,2])) < 3 ) {
        kw.path.spearman.cor[i, j] <- NA
        kw.path.spearman.pval[i, j] <- NA
      }else if (feature.name[i] == feature.name[j]) {
        kw.path.spearman.cor[i, j] <- NA
        kw.path.spearman.pval[i, j] <- NA
      }else{
        a.2.result <- partial_Spearman(x|y ~ Age, a.2)
        kw.path.spearman.cor[i, j] <- a.2.result$TS$TB$ts
        kw.path.spearman.pval[i, j] <- a.2.result$TS$TB$pval
      }
    }
  }
  return(list(kw.path.spearman.cor = kw.path.spearman.cor,
              kw.path.spearman.pval = kw.path.spearman.pval))
}

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

## for CRC datasets
study <- c("PRJEB6070", "PRJEB27928", "PRJEB10878")

crc.adj.spear.list <- foreach(pro = study) %dopar% adjusted.spearman.func(project = pro, 
                                                                              metadata = metadata, 
                                                                              kw.path.spe.df = kw.path.spe.df,
                                                                              case = "CRC", 
                                                                              covariates = c("Age"))
save(crc.adj.spear.list, file = "partial_crc_spe_pathways.RData")

## for PRJNA400072 datasets
case_type <- c("CD", "UC")

ibd.adj.spear.list <- foreach(grp = case_type) %dopar% adjusted.spearman.func(project = "PRJNA400072", 
                                                                              metadata = metadata, 
                                                                              kw.path.spe.df = kw.path.spe.df,
                                                                              case = grp, 
                                                                              covariates = c("Age"))
save(ibd.adj.spear.list, file = "partial_prjna400072_spe_pathways.RData")

## for PRJEB1220 dataset
##
adjusted.spearman.func.2 <- function(metadata = metadata, 
                                   kw.path.spe.df = kw.path.spe.df,
                                   project, 
                                   case = "CRC",
                                   covariates = c("Age", "BMI")){
  library(dplyr)
  library(PResiduals)
  kw.path.spearman.cor <- list()
  kw.path.spearman.pval <- list()
  samples <- metadata %>% filter(Project %in% project) %>%
    filter(Group %in% c("CTR", case)) %>%
    pull(Sample_ID)
  a <- kw.path.spe.df[samples, ]
  feature.name <- colnames(kw.path.spe.df)[1:(ncol(kw.path.spe.df)-3)]
  
  kw.path.spearman.cor <- matrix(data = NA,nrow = length(feature.name), 
                             ncol = length(feature.name), 
                             dimnames = list(feature.name, 
                                             feature.name))
  
  kw.path.spearman.pval <- matrix(data = NA,nrow = length(feature.name), 
                                  ncol = length(feature.name), 
                                  dimnames = list(feature.name, 
                                                  feature.name))
  
  for (i in 1:length(feature.name)) {
    for (j in 1:length(feature.name)) {
      a.2 <- a[, c(feature.name[i], feature.name[j], covariates)]
      names(a.2) <- c("x", "y", 
                      "Age", "BMI")
      if (length(unique(a.2[,1])) < 3 | length(unique(a.2[,2])) < 3 ) {
        kw.path.spearman.cor[i, j] <- NA
        kw.path.spearman.pval[i, j] <- NA
      }else if (feature.name[i] == feature.name[j]) {
        kw.path.spearman.cor[i, j] <- NA
        kw.path.spearman.pval[i, j] <- NA
      }else{
        a.2.result <- partial_Spearman(x|y ~ Age + BMI, a.2)
        kw.path.spearman.cor[i, j] <- a.2.result$TS$TB$ts
        kw.path.spearman.pval[i, j] <- a.2.result$TS$TB$pval
      }
    }
  }
  return(list(kw.path.spearman.cor = kw.path.spearman.cor,
              kw.path.spearman.pval = kw.path.spearman.pval))
}
uc.adj.spear.list <- list()
uc.adj.spear.list[["PRJEB1220"]] <-  adjusted.spearman.func.2(project = "PRJEB1220", 
                                                            metadata = metadata, 
                                                            kw.path.spe.df = kw.path.spe.df,
                                                            case = "UC", 
                                                            covariates = c("Age", "BMI"))
save(uc.adj.spear.list, file = "partial_prjeb1220_spe_pathways.RData")
## ---------------------------------------------------------------------------------------------

## 
load("../RData/partial_prjeb1220_spe_pathways.RData")
load("../RData/partial_crc_spe_pathways.RData")
load("../RData/partial_prjna400072_spe_pathways.RData")

names(uc.adj.spear.list)
names(ibd.adj.spear.list)
names(crc.adj.spear.list)

cd.adj.spear.list <- list()
cd.adj.spear.list[["PRJNA400072"]] <- ibd.adj.spear.list$PRJNA400072_CD
uc.adj.spear.list[["PRJNA400072"]] <- ibd.adj.spear.list$PRJNA400072_UC

## CRC projects
for (pro in setdiff(crc.project, 
                    c("PRJEB6070","PRJEB27928", "PRJEB10878"))) {
  
  metadata.pro <- subset(metadata, Project %in% pro) %>% 
    pull(Sample_ID)
  
  a <- rcorr(as.matrix(kw.path.spe.df[metadata.pro, !names(kw.path.spe.df) %in% c("Age","Gender" ,"BMI")]), 
             type = "spearman")
  crc.adj.spear.list[[pro]]$kw.path.spearman.cor <- a$r
  crc.adj.spear.list[[pro]]$kw.path.spearman.pval <- a$P
  
}

## CD projects

for (pro in setdiff(cd.project, 
                    "PRJNA400072")) {
  
  metadata.pro <- subset(metadata, Project %in% pro) %>% 
    filter(Group %in% c("CTR", "CD")) %>%
    pull(Sample_ID)
  
  a <- rcorr(as.matrix(kw.path.spe.df[metadata.pro, !names(kw.path.spe.df) %in% c("Age","Gender" ,"BMI")]), 
             type = "spearman")
  cd.adj.spear.list[[pro]]$kw.path.spearman.cor <- a$r
  cd.adj.spear.list[[pro]]$kw.path.spearman.pval <- a$P
  
}

## UC projects

for (pro in "PRJNA389280") {
  
  metadata.pro <- subset(metadata, Project %in% pro) %>% 
    filter(Group %in% c("CTR", "UC")) %>%
    pull(Sample_ID)
  
  a <- rcorr(as.matrix(kw.path.spe.df[metadata.pro, !names(kw.path.spe.df) %in% c("Age","Gender" ,"BMI")]), 
             type = "spearman")
  uc.adj.spear.list[[pro]]$kw.path.spearman.cor <- a$r
  uc.adj.spear.list[[pro]]$kw.path.spearman.pval <- a$P
  
}

##
names(uc.adj.spear.list) <- paste0("UC_", names(uc.adj.spear.list))
names(cd.adj.spear.list) <- paste0("CD_", names(cd.adj.spear.list))
names(crc.adj.spear.list) <- paste0("CRC_", names(crc.adj.spear.list))

all.adj.spear.list <- c(crc.adj.spear.list, 
                        cd.adj.spear.list, 
                        uc.adj.spear.list)

kw.path.spearman.cor <- list()
kw.path.spearman.pval <- list()

for (i in names(all.adj.spear.list)) {
  a.triu <- all.adj.spear.list[[i]]$kw.path.spearman.cor
  a.triu <- a.triu[-c(1:87), c(1:87)]
  a.triu.df <- data.frame(a.triu, stringsAsFactors = F)
  a.triu.df$features <- rownames(a.triu.df) 
  a.triu.df <- melt(a.triu.df, id.vars = "features")
  a.triu.df$data <- i
  a.triu.df <- a.triu.df %>%
    filter(!is.na(value))
  kw.path.spearman.cor[[i]] <- a.triu.df
  
  
  a.triu <- all.adj.spear.list[[i]]$kw.path.spearman.pval
  a.triu <- a.triu[-c(1:87), c(1:87)]
  a.triu.df <- data.frame(a.triu, stringsAsFactors = F)
  a.triu.df$features <- rownames(a.triu.df) 
  a.triu.df <- melt(a.triu.df, id.vars = "features")
  a.triu.df$data <- i
  a.triu.df <- a.triu.df %>%
    filter(!is.na(value))
  kw.path.spearman.pval[[i]] <- a.triu.df
}

## first
## filter those correlation with P-value < 0.01
kw.path.spearman.pval.filter <- kw.path.spearman.pval %>% 
  reduce(rbind)

kw.path.spearman.pval.final <- kw.path.spearman.pval.filter %>% 
  filter(value <0.01)
names(kw.path.spearman.pval.final)[3] <- "p_value"

## second
## 
kw.path.spearman.cor.df <- kw.path.spearman.cor %>% 
  reduce(rbind)
names(kw.path.spearman.cor.df)[3] <- "coeff"

kw.path.spearman.cor.merge <- merge(kw.path.spearman.pval.final, 
              kw.path.spearman.cor.df,
              by = c("features", "variable",
                     "data"),
              all.x = T)
## 9798
sum(is.na(kw.path.spearman.cor.merge$coeff))

## metadata n_samples
kw.path.spearman.metadata <- data.frame(Project = unique(kw.path.spearman.cor.merge$data), 
                                        N_samples = c(237, 144, 157, 128, 132, 82, 80,
                                                      110, 58, 104, 53, 83, 107))
save(all.adj.spear.list, 
     kw.path.spearman.cor.merge,
     kw.path.spearman.metadata, 
     file = "/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/RData/kw_cor_path.Rdata")

## meta-analysis
load("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/RData/kw_cor_path.Rdata")
kw.path.spearman.cor.merge$combined <- paste(kw.path.spearman.cor.merge$features, 
                                             kw.path.spearman.cor.merge$variable, 
                                             sep = "_")
kw.path.spearman.cor.merge.2 <- merge(kw.path.spearman.cor.merge, 
                                      kw.path.spearman.metadata, 
                                      by.x = "data", 
                                      by.y = "Project",
                                      all.x = T)

kw.path.spearman.meta <- list()

for (comp in unique(kw.path.spearman.cor.merge.2$combined)) {
  degra.sub <- kw.path.spearman.cor.merge.2 %>% 
    filter(combined %in% comp)
  if (nrow(degra.sub)>1) {
    degra.sub.zscore <- metacor(coeff, N_samples, data, 
                                data = degra.sub, sm = "ZCOR")
    degra.sub.correct.r <- (exp(2* degra.sub.zscore$TE.random)-1)/(exp(2* degra.sub.zscore$TE.random)+1)
    
    kw.path.spearman.meta[[comp]]$meta_z <- degra.sub.zscore
    kw.path.spearman.meta[[comp]]$correct_r <- degra.sub.correct.r
  }
}
data.cor.correct.r.df <- lapply(kw.path.spearman.meta,
                                function(data){
                                  y <- data$correct_r;
                                  return(y)})  %>% reduce(rbind)

data.cor.correct.r.df <- data.frame(data.cor.correct.r.df,
                                    stringsAsFactors = F)
names(data.cor.correct.r.df)[1] <- "random.r"
data.cor.correct.r.df$target <- names(kw.path.spearman.meta)

kw.path.spearman.cor.merge.3 <- merge(kw.path.spearman.cor.merge.2,
                      data.cor.correct.r.df,
                      by.x = c("combined"),
                      by.y= c("target"),
                      all.y= T) %>% distinct(combined, .keep_all = T)

kw.path.spearman.meta.dcast <- dcast(kw.path.spearman.cor.merge.3[,c("features","variable","random.r")],
                          features ~ variable)
kw.path.spearman.meta.dcast[is.na(kw.path.spearman.meta.dcast)] <- "0"
kw.path.spearman.meta.dcast[,-1] <- apply(kw.path.spearman.meta.dcast[,-1] ,2, as.numeric)
rownames(kw.path.spearman.meta.dcast) <- kw.path.spearman.meta.dcast$features


data.cor.correct.pval.df <- data.frame(unlist(lapply(kw.path.spearman.meta,
                                                      function(data){y <- data$meta_z$pval.random})), 
                                        stringsAsFactors = F)
data.cor.correct.pval.df$target <- names(kw.path.spearman.meta)
names(data.cor.correct.pval.df)[1] <- "random.pvalue"

data.cor.correct.pval.df <- data.cor.correct.pval.df %>%
  arrange(random.pvalue) %>% 
  mutate(fdr = p.adjust(random.pvalue, method = "fdr")) 

data.cor.correct.pval.merge <- merge(kw.path.spearman.cor.merge.2,
                                     data.cor.correct.pval.df,
                                      by.x = c("combined"),
                                      by.y= c("target"),
                                      all.y= T) %>% distinct(combined, .keep_all = T)

kw.path.spearman.meta.pval.dcast <- dcast(data.cor.correct.pval.merge[,c("features","variable","fdr")],
                                         features ~ variable)

kw.path.spearman.meta.pval.dcast[is.na(kw.path.spearman.meta.pval.dcast)] <- 0
kw.path.spearman.meta.pval.dcast.2 <- apply(kw.path.spearman.meta.pval.dcast[,-1], 2, function(data){
  y <- cut(data, breaks = c(0,0.05,1),
           include.lowest = F,
           labels=c("*"," "));
  return(y)})
kw.path.spearman.meta.pval.dcast.2[is.na(kw.path.spearman.meta.pval.dcast.2)] <- " "

rownames(kw.path.spearman.meta.pval.dcast.2) <- kw.path.spearman.meta.pval.dcast$features

save(kw.path.spearman.meta.dcast, 
     kw.path.spearman.meta.pval.dcast.2, 
     file = "/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/RData/kw_cor_path_metaresults.Rdata")


## load results
load("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/RData/kw_cor_path_metaresults.Rdata")
kw.path.spearman.meta.pval.dcast.2[1:5,1:5]
kw.path.spearman.meta.dcast[1:5,1:5]

kw.path.meta.coef.final <- kw.path.spearman.meta.dcast[, -1]
kw.path.meta.pval.final <- kw.path.spearman.meta.pval.dcast.2

## extract the full names of the pathways
sig.path.df.merge <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SuppleData/HDC_correlated_DifPathways.txt",
                                header = T, sep = "\t", as.is = T)

kw.pathways.map <- sig.path.df.merge %>% select(Abb, feature, Superclass) %>%
  distinct(Abb, feature, .keep_all = T)
rownames(kw.pathways.map) <- kw.pathways.map$Abb

## testing if matching
length(intersect(rownames(kw.path.meta.coef.final), 
                 rownames(kw.pathways.map)))
kw.pathways.map <- kw.pathways.map[rownames(kw.path.meta.coef.final), ]
rownames(kw.path.meta.coef.final) <- kw.pathways.map$feature

## testing if matching
length(intersect(dimnames(kw.path.meta.pval.final)[[1]], 
                 rownames(kw.pathways.map)))
kw.pathways.map <- kw.pathways.map[dimnames(kw.path.meta.pval.final)[[1]], ]
dimnames(kw.path.meta.pval.final)[[1]] <- kw.pathways.map$feature

## split species
kw.spe.df.2 <- data.frame(Count = table(kw.spe.df$feature), 
                          stringsAsFactors = F)
names(kw.spe.df.2) <- c("Species", "Count")

kw.spe.df.sum <- kw.spe.df %>% 
  select(feature, Trend, exposure) %>%
  distinct(.keep_all = T)
kw.spe.df.2.merge <- merge(kw.spe.df.2, 
                           kw.spe.df.sum, 
                           by.x = "Species", 
                           by.y = "feature", 
                           all = T)
kw.spe.df.2.merge[which(kw.spe.df.2.merge$Count >1), "exposure"] <- "Common"
kw.spe.df.2.merge[which(kw.spe.df.2.merge$Species %in% c("Alistipes_onderdonkii","Ruminococcus_torques")), "Trend"] <- "Conflict"

kw.spe.df.2.merge <- kw.spe.df.2.merge %>%
  distinct(.keep_all = T) %>% arrange(Trend, exposure)

kw.spe.df.2.merge <- kw.spe.df.2.merge %>%
  mutate(Trend = factor(Trend, levels = c("Case-enriched", "Control-enriched", "Conflict"))) %>%
  mutate(exposure = factor(exposure, levels = c("UC", "CD", "CRC", "Common")))
write.table(kw.spe.df.2.merge, 
            "/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SuppleData/Dif_species_annotation.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)


## annotation based on dagradation pathways
kw.spe.df <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SuppleData/Dif_species_annotation.txt" ,
                                header = T, sep = "\t", as.is = T)
rownames(kw.spe.df) <- kw.spe.df$Species
kw.spe.df$final.spe <- gsub(kw.spe.df$Species, pattern = "_", replacement = " ", fixed = T)
kw.spe.df$final.spe <- gsub(kw.spe.df$final.spe, pattern = "unclassified", replacement = "spp.", fixed = T)

length(intersect(rownames(kw.spe.df), names(kw.path.meta.coef.final)))
length(intersect(rownames(kw.spe.df), dimnames(kw.path.meta.pval.final)[[2]]))

names(kw.path.meta.coef.final) == dimnames(kw.path.meta.pval.final)[[2]]

names(kw.path.meta.coef.final) <- kw.spe.df[names(kw.path.meta.coef.final), "final.spe"]
dimnames(kw.path.meta.pval.final)[[2]] <- kw.spe.df[dimnames(kw.path.meta.pval.final)[[2]], "final.spe"]

rownames(kw.spe.df) <- kw.spe.df$final.spe

kw.spe.df.annot <- kw.spe.df[, 3:4]

## degradation pathways
degradation.annot <- kw.pathways.map %>%
  filter(Superclass %in% c("Amino Acid Degradation", 
                           "Carbohydrate Degradation", 
                           "Carboxylate Degradation", 
                           "Nucleoside and Nucleotide Degradation", 
                           "Secondary Metabolite Degradation"))
degradation.annot.2 <- data.frame(Superclass = degradation.annot$Superclass,
                                  stringsAsFactors = F)
rownames(degradation.annot.2) <- degradation.annot$feature

kw.degradation.path.sum <- kw.path.df.merge %>% 
  filter(feature %in% rownames(degradation.annot)) %>% 
  group_by(feature, Trend) %>% summarise(count = n()) %>% data.frame

rownames(kw.degradation.path.sum) <- kw.degradation.path.sum$feature

degradation.annot.2$Trend <- kw.degradation.path.sum[rownames(degradation.annot.2), "Trend"]


pathways.annot.color <- list(Superclass = c(`Amino Acid Degradation`="#8DD3C7",
                    `Carbohydrate Degradation`="#FFED6F",
                    `Carboxylate Degradation`="#BEBADA",
                    `Nucleoside and Nucleotide Degradation`="#FDB462",
                    `Secondary Metabolite Degradation`="#B3DE69"),
     exposure = c(UC = "#d0a727",
     CRC = "#9CC5C9",
     Common = "#A08689",
     CD = "#87a44f"),
     Trend = c(`Control-enriched` ="#56B4E9",
               `Case-enriched` = "firebrick3",
               Conflict = "gray"))

## other pathways
biosynthesis.annot <- kw.pathways.map %>%
  filter(!Superclass %in% c("Amino Acid Degradation", 
                           "Carbohydrate Degradation", 
                           "Carboxylate Degradation", 
                           "Nucleoside and Nucleotide Degradation", 
                           "Secondary Metabolite Degradation"))
biosynthesis.annot.2 <- data.frame(Superclass = biosynthesis.annot$Superclass,
                                  stringsAsFactors = F)
rownames(biosynthesis.annot.2) <- biosynthesis.annot$feature

kw.biosyn.path.sum <- kw.path.df.merge %>% 
  filter(feature %in% rownames(biosynthesis.annot.2)) %>% 
  group_by(feature, Trend) %>% summarise(count = n()) %>% data.frame

rownames(kw.biosyn.path.sum) <- kw.biosyn.path.sum$feature

biosynthesis.annot.2$Trend <- kw.biosyn.path.sum[rownames(biosynthesis.annot.2), "Trend"]


pathways.annot.color.other <- list(Superclass = c(`Fermentation` = "#8DD3C7",
                                            `Generation of Precursor Metabolites and Energy` = "#FFFFB3",
                                            `Inorganic Nutrient Metabolism` = "#BEBADA",
                                            `Aminoacyl-tRNA Charging` = "#FB8072",
                                            `Amino Acid Biosynthesis` = "#80B1D3", 
                                            `Glycolysis` = "#FFED6F",
                                            `Aromatic Compound Biosynthesis` = "#FDB462", 
                                            `Carbohydrate Biosynthesis` = "#B3DE69",
                                            `Secondary Metabolite Biosynthesis` = "#FCCDE5",
                                            `Cofactor, Carrier, and Vitamin Biosynthesis` = "#D9D9D9" ,
                                            `Fatty Acid and Lipid Biosynthesis` = "#BC80BD",
                                            `Nucleoside and Nucleotide Biosynthesis` = "#CCEBC5"),
                             exposure = c(UC = "#d0a727",
                                          CRC = "#9CC5C9",
                                          Common = "#A08689",
                                          CD = "#87a44f"),
                             Trend = c(`Control-enriched` ="#56B4E9",
                                       `Case-enriched` = "firebrick3",
                                       Conflict = "gray"))

## degradation pheatmap
bk <- c(seq(-0.8, -0.01, by = 0.01), 
        seq(0, 0.8, by = 0.01))
kw.degradation.pheatmap <- pheatmap(as.matrix(kw.path.meta.coef.final[rownames(degradation.annot), kw.spe.df$final.spe]),
         cluster_rows = T,
         cluster_cols = F, 
         annotation_colors = pathways.annot.color,
         annotation_col = kw.spe.df.annot,
         annotation_row = degradation.annot.2,
         clustering_method = "mcquitty",
         legend_breaks = seq(-0.8,0.8,0.4), breaks = bk, 
         color = c(colorRampPalette(colors = c("#016392","white"))(80),
                   colorRampPalette(colors = c("white","#c72e29"))(80)), 
         display_numbers = kw.path.meta.pval.final[rownames(degradation.annot), kw.spe.df$final.spe],
         fontsize = 11,  
         cellwidth = 12, cellheight = 12, border=FALSE, 
         gaps_col = c(6, 16, 45, 51, 53, 
                      55, 77, 79, 85),
         filename = NA)

pdf("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/pdf/Species_degradation_pheatmap.pdf",
    width = 30,
    height = 16)
kw.degradation.pheatmap
dev.off()

kw.biosyn.pheatmap <- pheatmap(as.matrix(kw.path.meta.coef.final[rownames(biosynthesis.annot.2), kw.spe.df$final.spe]),
                               cluster_rows = T,
                               cluster_cols = F, 
                               annotation_colors = pathways.annot.color.other,
                               annotation_col = kw.spe.df.annot,
                               annotation_row = biosynthesis.annot.2,
                               clustering_method = "mcquitty",
                               legend_breaks = seq(-0.8,0.8,0.4), breaks = bk, 
                               color = c(colorRampPalette(colors = c("#016392","white"))(80),
                                         colorRampPalette(colors = c("white","#c72e29"))(80)), 
                               display_numbers = kw.path.meta.pval.final[rownames(biosynthesis.annot.2), kw.spe.df$final.spe],
                               fontsize = 11,  
                               cellwidth = 12, cellheight = 12, border=FALSE, 
                               gaps_col = c(6, 16, 45, 51, 53, 
                                            55, 77, 79, 85),
                               filename = NA)
                               

pdf("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/pdf/Species_biosynthesis_pheatmap.pdf",
    width = 30,
    height = 16)
kw.biosyn.pheatmap
dev.off()

save(kw.path.meta.coef.final, 
     kw.path.meta.pval.final, 
     pathways.annot.color,
     pathways.annot.color.other, 
     kw.spe.df.annot, 
     biosynthesis.annot.2, 
     degradation.annot,
     kw.spe.df, 
     file = "../RData/Figure3_species_pathways.RData")