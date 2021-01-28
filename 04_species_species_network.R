## 
library(SpiecEasi)
## import counts

## IBD datasets 
metaphlan.counts <- list()
metaphlan.counts[["PRJEB1220"]] <- read.delim("/mnt/raid5/puzi/CDdata/PRJEB1220/02_bowtie_reads.txt",
                                              header = F, sep = " ", stringsAsFactors = F)
metaphlan.counts[["PRJNA400072"]] <- read.delim("/mnt/raid5/puzi/CDdata/PRJNA400072/02_bowtie_reads.txt",
                                                header = F, sep = " ", stringsAsFactors = F)
metaphlan.counts[["PRJNA389280"]] <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/IBD/PRJNA389280/02_bowtie_reads.txt",
                                                header = F, sep = " ", stringsAsFactors = F)
metaphlan.counts[["SRP057027"]] <- read.delim("/mnt/raid5/puzi/CDdata/SRP057027/02_bowtie_reads.txt",
                                                header = F, sep = " ", stringsAsFactors = F)

metaphlan.counts <- lapply(metaphlan.counts, 
       function(data){
         names(data) <- c("Sample_ID", "counts");
         data$Sample_ID <- lapply(strsplit(data$Sample_ID,
                                    split = ".bowtie2",
                                    fixed = T), function(x){ y <- x[1]}) %>% unlist;
         return(data)
       })


metaphlan.counts.cd <- metaphlan.counts %>% reduce(rbind)
metaphlan.counts.cd$Sample_ID <- gsub(metaphlan.counts.cd$Sample_ID, pattern = "-", 
                                      replacement = ".", fixed = T)
metaphlan.counts.cd <- metaphlan.counts.cd %>%
  filter(Sample_ID %in% metadata$Sample_ID)

rownames(metaphlan.counts.cd) <- metaphlan.counts.cd$Sample_ID


## CRC datasets 
metaphlan.counts[["PRJEB10878"]] <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CRC/PRJEB10878/02_bowtie_reads.txt",
                                              header = F, sep = " ", stringsAsFactors = F)
metaphlan.counts[["PRJEB12449"]] <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CRC/PRJEB12449/02_bowtie_reads.txt",
                                                header = F, sep = " ", stringsAsFactors = F)
metaphlan.counts[["PRJEB27928"]] <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CRC/PRJEB27928/02_bowtie_reads.txt",
                                                header = F, sep = " ", stringsAsFactors = F)
metaphlan.counts[["PRJDB4176"]] <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CRC/PRJDB4176/02_bowtie_reads.txt",
                                              header = F, sep = " ", stringsAsFactors = F)
metaphlan.counts[["PRJEB6070"]] <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CRC/PRJEB6070/02_bowtie_reads.txt",
                                              header = F, sep = " ", stringsAsFactors = F)
metaphlan.counts[["PRJEB7774"]] <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CRC/PRJEB7774/02_bowtie_reads.txt",
                                              header = F, sep = " ", stringsAsFactors = F)
metaphlan.counts[["PRJNA447983"]] <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CRC/PRJNA447983/02_bowtie_reads.txt",
                                              header = F, sep = " ", stringsAsFactors = F)

metaphlan.counts.2 <- lapply(metaphlan.counts[5:11], 
                           function(data){
                             names(data) <- c("Sample_ID", "counts");
                             data$Sample_ID <- lapply(strsplit(data$Sample_ID,
                                                               split = ".bowtie2",
                                                               fixed = T), function(x){ y <- x[1]}) %>% unlist;
                             return(data)
                           })

## CRC datasets
PRJNA447983.run <- read.delim("/mnt/raid5/puzi/IBD/otherproject/metadata/prjna447983_runinfo.txt",
                              header = T, sep = "\t", as.is = T)
rownames(PRJNA447983.run) <- PRJNA447983.run$Run
rownames(metaphlan.counts.2[["PRJNA447983"]]) <- metaphlan.counts.2[["PRJNA447983"]]$Sample_ID
metaphlan.counts.2[["PRJNA447983"]]$Sample_ID <- PRJNA447983.run[rownames(metaphlan.counts.2[["PRJNA447983"]]), "Sample_Name"]
rm(PRJNA447983.run)

PRJDB4176.run <- read.delim("/mnt/raid5/puzi/IBD/otherproject/metadata/prjdb4176_runinfo.txt",
                            header = T, sep = "\t", as.is = T)
rownames(PRJDB4176.run) <- PRJDB4176.run$Run
rownames(metaphlan.counts.2[["PRJDB4176"]]) <- metaphlan.counts.2[["PRJDB4176"]]$Sample_ID
metaphlan.counts.2[["PRJDB4176"]]$Sample_ID <- PRJDB4176.run[rownames(metaphlan.counts.2[["PRJDB4176"]]), "BioSample"]
rm(PRJDB4176.run)

metaphlan.counts.crc <- metaphlan.counts.2 %>% reduce(rbind)
length(intersect(metaphlan.counts.crc$Sample_ID, metadata$Sample_ID))
rownames(metaphlan.counts.crc) <- metaphlan.counts.crc$Sample_ID


## separation
project <- list(CRC = crc.project, 
                CD = cd.project, 
                UC = uc.project)
data.sepa.list <- list()

for (disease in c("CRC", "CD", "UC")) {
  data <- abund.s.list.filter[[disease]]
  data.t <- t(data[, !names(data) %in% "features"]);
  data.t <- data.frame(data.t, 
                       stringsAsFactors = F);
  study <- project[[disease]]
  for (i in study) {
    meta.ctr <- metadata %>% filter(Project %in% i) %>%
      filter(Group %in% "CTR") %>% 
      pull(Sample_ID)
    
    meta.case <- metadata %>% filter(Project %in% i) %>%
      filter(Group %in% disease)  %>% 
      pull(Sample_ID)
    
    data.sepa.list[[paste(i, "CTR", sep = "_")]] <- data.t[meta.ctr , ]
    data.sepa.list[[paste(i, disease, sep = "_")]] <- data.t[meta.case, ]
  }
}

## IBD dataframes
data.sepa.list.ibd <- data.sepa.list[15:24]
for (i in 1:10) {
  data.sepa.list.ibd[[i]] <- data.sepa.list.ibd[[i]][, -which(colSums(data.sepa.list.ibd[[i]]) == 0)]
  data.sepa.list.ibd[[i]] <- data.sepa.list.ibd[[i]] * metaphlan.counts.cd[rownames(data.sepa.list.ibd[[i]]), "counts"]
  
  data.sepa.list.ibd[[i]] <- data.frame(features = names(data.sepa.list.ibd[[i]]), 
                                        t(data.sepa.list.ibd[[i]]),
                                        stringsAsFactors = F)
  write.table(data.sepa.list.ibd[[i]], 
              file = paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SparCC/", 
                           names(data.sepa.list.ibd)[i], "_abund"),
              col.names = T, row.names = F, sep = "\t", quote = F)
}

## CRC dataframes
data.sepa.list.crc <- data.sepa.list[1:14]
for (i in 1:14) {
  data.sepa.list.crc[[i]] <- data.sepa.list.crc[[i]][, -which(colSums(data.sepa.list.crc[[i]]) == 0)]
  data.sepa.list.crc[[i]] <- data.sepa.list.crc[[i]] * metaphlan.counts.crc[rownames(data.sepa.list.crc[[i]]), "counts"]
  
  data.sepa.list.crc[[i]] <- data.frame(features = names(data.sepa.list.crc[[i]]), 
                                        t(data.sepa.list.crc[[i]]),
                                        stringsAsFactors = F)
  write.table(data.sepa.list.crc[[i]], 
              file = paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SparCC/", 
                            names(data.sepa.list.crc)[i], "_abund"),
              col.names = T, row.names = F, sep = "\t", quote = F)
}

## compute prevalence
abund.s.filter.prevalence <- list()
kw.spe.prevalence <- list()

for (grp in c("CRC", "CD", "UC")) {
  # for disease samples
  meta.sub<- metadata %>%
    filter(Group %in% grp) %>%
    filter(Project %in% project[[grp]]) %>%
    pull(Sample_ID)
  abund.sub <- abund.s.list.filter.df.t[meta.sub, ]
  abund.sub.prev <- apply(abund.sub, 2, 
                          function(data){ y <- length(data[data>0])/length(data); return(y)})
  abund.sub.prev <- data.frame(Prevalence = abund.sub.prev, 
                               stringsAsFactors = F)
  abund.sub.prev$Species <- rownames(abund.sub.prev)
  abund.sub.prev$Group <- grp
  abund.s.filter.prevalence[[grp]] <- abund.sub.prev
  
  # for control samples
  meta.sub<- metadata %>%
    filter(Group %in% "CTR") %>%
    filter(Project %in% project[[grp]]) %>%
    pull(Sample_ID)
  abund.sub <- abund.s.list.filter.df.t[meta.sub, ]
  abund.sub.prev <- apply(abund.sub, 2, 
                          function(data){ y <- length(data[data>0])/length(data); return(y)})
  abund.sub.prev <- data.frame(Prevalence = abund.sub.prev, 
                               stringsAsFactors = F)
  abund.sub.prev$Species <- rownames(abund.sub.prev)
  abund.sub.prev$Group <- "CTR"
  abund.s.filter.prevalence[[paste0(grp, "_CTR")]] <- abund.sub.prev
  
  ## filter the kw species
  a <- kw.spe.df %>%
    filter(exposure %in% grp) 
  rownames(a) <- a$feature
  
  kw.spe.prevalence[[grp]] <- abund.s.filter.prevalence[[grp]] %>%
    filter(Species %in% a$feature)
  kw.spe.prevalence[[grp]]$Type <- a[rownames(kw.spe.prevalence[[grp]]), "Trend"]
  
  kw.spe.prevalence[[paste0(grp, "_CTR")]] <- abund.s.filter.prevalence[[paste0(grp, "_CTR")]] %>%
    filter(Species %in% a$feature)
  kw.spe.prevalence[[paste0(grp, "_CTR")]]$Type <- a[rownames(kw.spe.prevalence[[paste0(grp, "_CTR")]]), "Trend"]
}



for(i in names(kw.spe.prevalence)){
  write.table(kw.spe.prevalence[[i]],
              file = paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/Network/",
                            i, "_prevalence.txt"),
              row.names = F,
              col.names = T, sep = "\t",
              quote = F)
}


## sh /mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SparCC/submit.sh

## IBD datasets
## load P-value 

sparcc.ibd.cor <- list()
for (i in names(data.sepa.list.ibd)) {
  sparcc.ibd.cor[[i]] <- read.delim(paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SparCC/", 
                                            i, "_abund_cor.out"),
                                     header = T, sep = "\t", as.is = T)
  rownames(sparcc.ibd.cor[[i]]) <- sparcc.ibd.cor[[i]]$features
  sparcc.ibd.cor[[i]] <- sparcc.ibd.cor[[i]][,-1]
  sparcc.ibd.cor[[i]] <- as.matrix(sparcc.ibd.cor[[i]])
  sparcc.ibd.cor[[i]][upper.tri(sparcc.ibd.cor[[i]])] <- NA
  
  sparcc.ibd.cor[[i]] <- data.frame(sparcc.ibd.cor[[i]], 
                                    stringsAsFactors = F)
  sparcc.ibd.cor[[i]]$features <- rownames(sparcc.ibd.cor[[i]]) 

  sparcc.ibd.cor[[i]] <- sparcc.ibd.cor[[i]] %>%
    melt(id.vars = "features") %>%
    filter(!features == variable) %>%
    filter(!is.na(value))
  sparcc.ibd.cor[[i]]$file <- i
}

sparcc.ibd.pval <- list()
for (i in names(data.sepa.list.ibd)) {
  sparcc.ibd.pval[[i]] <- read.delim(paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SparCC/04_pvalue/", 
                                     i, "_abund_final_pval.txt"),
                                     header = T, sep = "\t", as.is = T)
  
  rownames(sparcc.ibd.pval[[i]]) <- sparcc.ibd.pval[[i]]$features
  sparcc.ibd.pval[[i]] <- sparcc.ibd.pval[[i]][,-1]
  sparcc.ibd.pval[[i]] <- as.matrix(sparcc.ibd.pval[[i]])
  sparcc.ibd.pval[[i]][upper.tri(sparcc.ibd.pval[[i]])] <- NA
  
  sparcc.ibd.pval[[i]] <- data.frame(sparcc.ibd.pval[[i]], 
                                    stringsAsFactors = F)
  sparcc.ibd.pval[[i]]$features <- rownames(sparcc.ibd.pval[[i]]) 
  
  sparcc.ibd.pval[[i]] <- sparcc.ibd.pval[[i]] %>%
    melt(id.vars = "features") %>%
    filter(!features == variable) %>%
    filter(!is.na(value))
  sparcc.ibd.pval[[i]]$file <- i
}

sparcc.ibd.cor <- lapply(sparcc.ibd.cor, 
       function(data){
         data$Project <- lapply(strsplit(data$file, 
                                         split = "_",fixed = T),
                                function(y){ y.2 <- y[1]}) %>% unlist
         
         data$Group <- lapply(strsplit(data$file, 
                                         split = "_",fixed = T),
                                function(y){ y.2 <- y[2]}) %>% unlist
         
         names(data) <- c("features", "variable", "r2", 
                          "file", "Project", "Group")
         
         return(data)
       })
sparcc.ibd.cor.df <- sparcc.ibd.cor %>%
  reduce(rbind)


sparcc.ibd.pval <- lapply(sparcc.ibd.pval, 
                         function(data){
                           data$Project <- lapply(strsplit(data$file, 
                                                           split = "_",fixed = T),
                                                  function(y){ y.2 <- y[1]}) %>% unlist
                           
                           data$Group <- lapply(strsplit(data$file, 
                                                         split = "_",fixed = T),
                                                function(y){ y.2 <- y[2]}) %>% unlist
                           
                           names(data) <- c("features", "variable", "P_value", 
                                            "file", "Project", "Group")
                           
                           return(data)
                         })
sparcc.ibd.pval.df <- sparcc.ibd.pval %>% reduce(rbind)

sparcc.ibd.merge <- merge(sparcc.ibd.cor.df, 
                          sparcc.ibd.pval.df, 
                          by = c("features", "variable", 
                                 "file", "Project", "Group"),
                          all = T)
sparcc.ibd.merge.sig <- sparcc.ibd.merge %>%
  filter(P_value < 0.05) 


data.sepa.list.ibd.nsamples <- lapply(data.sepa.list.ibd, function(y){
  y.2 <- ncol(y)-1;
  return(y.2)
}) %>% reduce(c) %>% data.frame
data.sepa.list.ibd.nsamples$file <- names(data.sepa.list.ibd)
names(data.sepa.list.ibd.nsamples) <- c("nSamples", "file")

sparcc.ibd.merge.sig <- merge(sparcc.ibd.merge.sig, 
                              data.sepa.list.ibd.nsamples, 
                              by = "file", all = T)

save(sparcc.ibd.merge.sig, 
      project, 
     file = "../RData/species_sparcc_IBD.Rdata")

## --------------------------------------------------------
## CRC datasets
sparcc.crc.cor <- list()
for (i in names(data.sepa.list.crc)) {
  sparcc.crc.cor[[i]] <- read.delim(paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SparCC/01_cor/", 
                                           i, "_abund_cor.out"),
                                    header = T, sep = "\t", as.is = T)
  rownames(sparcc.crc.cor[[i]]) <- sparcc.crc.cor[[i]]$features
  sparcc.crc.cor[[i]] <- sparcc.crc.cor[[i]][,-1]
  sparcc.crc.cor[[i]] <- as.matrix(sparcc.crc.cor[[i]])
  sparcc.crc.cor[[i]][upper.tri(sparcc.crc.cor[[i]])] <- NA
  
  sparcc.crc.cor[[i]] <- data.frame(sparcc.crc.cor[[i]], 
                                    stringsAsFactors = F)
  sparcc.crc.cor[[i]]$features <- rownames(sparcc.crc.cor[[i]]) 
  
  sparcc.crc.cor[[i]] <- sparcc.crc.cor[[i]] %>%
    melt(id.vars = "features") %>%
    filter(!features == variable) %>%
    filter(!is.na(value))
  sparcc.crc.cor[[i]]$file <- i
}

sparcc.crc.pval <- list()
for (i in names(data.sepa.list.crc)) {
  sparcc.crc.pval[[i]] <- read.delim(paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/SparCC/04_pvalue/", 
                                            i, "_abund_final_pval.txt"),
                                     header = T, sep = "\t", as.is = T)
  
  rownames(sparcc.crc.pval[[i]]) <- sparcc.crc.pval[[i]]$features
  sparcc.crc.pval[[i]] <- sparcc.crc.pval[[i]][,-1]
  sparcc.crc.pval[[i]] <- as.matrix(sparcc.crc.pval[[i]])
  sparcc.crc.pval[[i]][upper.tri(sparcc.crc.pval[[i]])] <- NA
  
  sparcc.crc.pval[[i]] <- data.frame(sparcc.crc.pval[[i]], 
                                     stringsAsFactors = F)
  sparcc.crc.pval[[i]]$features <- rownames(sparcc.crc.pval[[i]]) 
  
  sparcc.crc.pval[[i]] <- sparcc.crc.pval[[i]] %>%
    melt(id.vars = "features") %>%
    filter(!features == variable) %>%
    filter(!is.na(value))
  sparcc.crc.pval[[i]]$file <- i
}

sparcc.crc.cor <- lapply(sparcc.crc.cor, 
                         function(data){
                           data$Project <- lapply(strsplit(data$file, 
                                                           split = "_",fixed = T),
                                                  function(y){ y.2 <- y[1]}) %>% unlist
                           
                           data$Group <- lapply(strsplit(data$file, 
                                                         split = "_",fixed = T),
                                                function(y){ y.2 <- y[2]}) %>% unlist
                           
                           names(data) <- c("features", "variable", "r2", 
                                            "file", "Project", "Group")
                           
                           return(data)
                         })
sparcc.crc.cor.df <- sparcc.crc.cor %>%
  reduce(rbind)


sparcc.crc.pval <- lapply(sparcc.crc.pval, 
                          function(data){
                            data$Project <- lapply(strsplit(data$file, 
                                                            split = "_",fixed = T),
                                                   function(y){ y.2 <- y[1]}) %>% unlist
                            
                            data$Group <- lapply(strsplit(data$file, 
                                                          split = "_",fixed = T),
                                                 function(y){ y.2 <- y[2]}) %>% unlist
                            
                            names(data) <- c("features", "variable", "P_value", 
                                             "file", "Project", "Group")
                            
                            return(data)
                          })
sparcc.crc.pval.df <- sparcc.crc.pval %>% reduce(rbind)

sparcc.crc.merge <- merge(sparcc.crc.cor.df, 
                          sparcc.crc.pval.df, 
                          by = c("features", "variable", 
                                 "file", "Project", "Group"),
                          all = T)
sparcc.crc.merge.sig <- sparcc.crc.merge %>%
  filter(P_value < 0.05) 


data.sepa.list.crc.nsamples <- lapply(data.sepa.list.crc, function(y){
  y.2 <- ncol(y)-1;
  return(y.2)
}) %>% reduce(c) %>% data.frame
data.sepa.list.crc.nsamples$file <- names(data.sepa.list.crc)
names(data.sepa.list.crc.nsamples) <- c("nSamples", "file")

sparcc.crc.merge.sig <- merge(sparcc.crc.merge.sig, 
                              data.sepa.list.crc.nsamples, 
                              by = "file", all = T)

save(sparcc.crc.merge.sig, 
     project, 
     file = "../RData/species_sparcc_CRC.Rdata")


## in totoro
rm(list = ls())
library(metafor)
library(meta)

load("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/RData/species_sparcc_IBD.Rdata")

sparcc.ibd.merge.sig$target <- paste(sparcc.ibd.merge.sig$features,
                                     sparcc.ibd.merge.sig$variable,
                                     sep = "_")

sparcc.ibd.metacor.case <- list()
sparcc.ibd.metacor.ctr <- list()

for (grp in c("CD","UC")) {
  
  ## for disease samples
  targetspe <- sparcc.ibd.merge.sig %>% 
    filter(Group %in% grp) %>%
    filter(Project %in% project[[grp]]) %>%
    pull(target)
  
  for (comp in targetspe) {
    sparcc.sub <- sparcc.ibd.merge.sig %>% 
      filter(Group %in% grp & target %in% comp) %>% 
      filter(Project %in% project[[grp]]) 
    
    if (nrow(sparcc.sub)>1) {
      sub.zscore <- metacor(r2, nSamples, Project, 
                            data = sparcc.sub, 
                            sm = "ZCOR")
      sub.correct.r <- (exp(2* sub.zscore$TE.random)-1)/(exp(2* sub.zscore$TE.random)+1)
      
      sparcc.ibd.metacor.case[[paste(comp, grp, sep = "_")]]$meta_z <- sub.zscore
      sparcc.ibd.metacor.case[[paste(comp, grp, sep = "_")]]$correct_r <- sub.correct.r
    }
  }
  
  ## for control samples 
  targetspe <- sparcc.ibd.merge.sig %>% 
    filter(Group %in% "CTR") %>% 
    filter(Project %in% project[[grp]]) %>%
    pull(target)
  
  for (comp in targetspe) {
    sparcc.sub <- sparcc.ibd.merge.sig %>% 
      filter(Project %in% project[[grp]]) %>%
      filter(Group %in% "CTR" & target %in% comp) 
    
    if (nrow(sparcc.sub)>1) {
      sub.zscore <- metacor(r2, nSamples, Project, 
                            data = sparcc.sub, 
                            sm = "ZCOR")
      sub.correct.r <- (exp(2* sub.zscore$TE.random)-1)/(exp(2* sub.zscore$TE.random)+1)
      
      sparcc.ibd.metacor.ctr[[paste(comp, grp, sep = "_")]]$meta_z <- sub.zscore
      sparcc.ibd.metacor.ctr[[paste(comp, grp, sep = "_")]]$correct_r <- sub.correct.r
    }
  }
}

sparcc.ibd.metacor.case.result <- meta.correct.func(sparcc.ibd.metacor.case)
sparcc.ibd.metacor.ctr.result <- meta.correct.func(sparcc.ibd.metacor.ctr)

save(sparcc.ibd.metacor.case.result, 
     sparcc.ibd.metacor.ctr.result,
     file = "/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/RData/species_sparcc_IBD_meta.Rdata")

## in crick
load("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/RData/species_sparcc_IBD_meta.Rdata")

## correlation strength quantile
cut_edge <- function(r) {
  out <- cut(abs(r), breaks = c(0, 0.4, 0.6,1), include.lowest = F, 
             labels = c("1", "2", "3"))
  return(out)
}

## -------------IBD datasets---------------- ##

sparcc.ibd.metacor.case.merge <- merge(subset(sparcc.ibd.metacor.case.result$cor.merge, fdr <0.05),
                                       sparcc.ibd.merge.sig[,c("features","variable","target")] %>% 
                                         distinct(target, .keep_all = T),
                                       by.x = "target_real", by.y = "target",
                                       all.x = T)
sparcc.ibd.metacor.case.merge$edgestrength <- cut_edge(sparcc.ibd.metacor.case.merge$random.r)
sparcc.ibd.metacor.case.merge$edgetrend <- 1
sparcc.ibd.metacor.case.merge[which(sparcc.ibd.metacor.case.merge$random.r <0), "edgetrend"] <- -1

sparcc.ibd.metacor.case.merge$variable <- as.character(sparcc.ibd.metacor.case.merge$variable)
sparcc.ibd.metacor.case.merge$edgestrength <- as.character(sparcc.ibd.metacor.case.merge$edgestrength)


sparcc.ibd.metacor.case.kw <- list()
sparcc.ibd.metacor.case.kw[["CD"]] <- sparcc.ibd.metacor.case.merge %>%
  filter(disease %in% "CD") %>%
  filter(features %in% cd.kw.spe) %>% 
  filter(variable %in% cd.kw.spe)

sparcc.ibd.metacor.case.kw[["UC"]] <- sparcc.ibd.metacor.case.merge  %>%
  filter(disease %in% "UC") %>%
  filter(features %in% uc.kw.spe) %>%
  filter(variable %in% uc.kw.spe)

add.spe <- cd.kw.spe[!cd.kw.spe %in% unique(c(sparcc.ibd.metacor.case.kw[["CD"]]$features, 
                                              sparcc.ibd.metacor.case.kw[["CD"]]$variable))]

for (add in add.spe) {
  sparcc.ibd.metacor.case.kw[["CD"]][nrow(sparcc.ibd.metacor.case.kw[["CD"]])+1, 
                                    ] <- c(NA, NA, "CD",
                                           0, 0, 0, 
                                           add, add, 0, 0)
}

add.spe <- uc.kw.spe[!uc.kw.spe %in% unique(c(sparcc.ibd.metacor.case.kw[["UC"]]$features, 
                                              sparcc.ibd.metacor.case.kw[["UC"]]$variable))]

for (add in add.spe) {
  sparcc.ibd.metacor.case.kw[["UC"]][nrow(sparcc.ibd.metacor.case.kw[["UC"]])+1, 
                                    ] <- c(NA, NA, "UC",
                                           0, 0, 0, 
                                           add, add, 0, 0)
}

sparcc.ibd.metacor.ctr.merge <- merge(subset(sparcc.ibd.metacor.ctr.result$cor.merge, fdr <0.05),
                                       sparcc.ibd.merge.sig[,c("features","variable","target")] %>% 
                                         distinct(target, .keep_all = T),
                                       by.x = "target_real", by.y = "target",
                                       all.x = T)

sparcc.ibd.metacor.ctr.merge$edgestrength <- cut_edge(sparcc.ibd.metacor.ctr.merge$random.r)
sparcc.ibd.metacor.ctr.merge$edgetrend <- 1
sparcc.ibd.metacor.ctr.merge[which(sparcc.ibd.metacor.ctr.merge$random.r <0), "edgetrend"] <- -1

sparcc.ibd.metacor.ctr.merge$variable <- as.character(sparcc.ibd.metacor.ctr.merge$variable)
sparcc.ibd.metacor.ctr.merge$edgestrength <- as.character(sparcc.ibd.metacor.ctr.merge$edgestrength)



sparcc.ibd.metacor.ctr.kw <- list()
sparcc.ibd.metacor.ctr.kw[["CD"]] <- sparcc.ibd.metacor.ctr.merge  %>%
  filter(disease %in% "CD") %>%
  filter(features %in% cd.kw.spe) %>%
  filter(variable %in% cd.kw.spe)

sparcc.ibd.metacor.ctr.kw[["UC"]] <- sparcc.ibd.metacor.ctr.merge %>%
  filter(disease %in% "UC") %>%
  filter(features %in% uc.kw.spe)%>%
  filter(variable %in% uc.kw.spe)

add.spe <- cd.kw.spe[!cd.kw.spe %in% unique(c(sparcc.ibd.metacor.ctr.kw[["CD"]]$features, 
                                              sparcc.ibd.metacor.ctr.kw[["CD"]]$variable))]

for (add in add.spe) {
  sparcc.ibd.metacor.ctr.kw[["CD"]][nrow(sparcc.ibd.metacor.ctr.kw[["CD"]])+1, 
                            ] <- c(NA, NA, "CD",
                                   0, 0, 0, 
                                   add, add, 0, 0)
}

add.spe <- uc.kw.spe[!uc.kw.spe %in% unique(c(sparcc.ibd.metacor.ctr.kw[["UC"]]$features, 
                                              sparcc.ibd.metacor.ctr.kw[["UC"]]$variable))]

for (add in add.spe) {
  sparcc.ibd.metacor.ctr.kw[["UC"]][nrow(sparcc.ibd.metacor.ctr.kw[["UC"]])+1, 
                                    ] <- c(NA, NA, "UC",
                                           0, 0, 0, 
                                           add, add, 0, 0)
}


for(i in names(sparcc.ibd.metacor.ctr.kw)){
  write.table(sparcc.ibd.metacor.ctr.kw[[i]],
              file = paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/Network/",
                            i, "_CTR_network.txt"),
              row.names = F,
              col.names = T, sep = "\t",
              quote = F)
}

for(i in names(sparcc.ibd.metacor.case.kw)){
  write.table(sparcc.ibd.metacor.case.kw[[i]],
              file = paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/Network/",
                            i, "_Case_network.txt"),
              row.names = F,
              col.names = T, sep = "\t",
              quote = F)
}

## -------------CRC datasets---------------- ##
load("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/RData/species_sparcc_CRC_meta.Rdata")

sparcc.crc.merge.sig$target <- paste(sparcc.crc.merge.sig$features,
                                     sparcc.crc.merge.sig$variable, 
                                     sep = "_")

sparcc.crc.metacor.case.merge <- merge(subset(sparcc.crc.metacor.case.result$cor.merge, fdr <0.05),
                                       sparcc.crc.merge.sig[,c("features","variable","target")] %>% 
                                         distinct(target, .keep_all = T),
                                       by.x = "target_real", by.y = "target",
                                       all.x = T)
sparcc.crc.metacor.case.merge$edgestrength <- cut_edge(sparcc.crc.metacor.case.merge$random.r)
sparcc.crc.metacor.case.merge$edgetrend <- 1
sparcc.crc.metacor.case.merge[which(sparcc.crc.metacor.case.merge$random.r <0), "edgetrend"] <- -1

sparcc.crc.metacor.case.merge$variable <- as.character(sparcc.crc.metacor.case.merge$variable)
sparcc.crc.metacor.case.merge$edgestrength <- as.character(sparcc.crc.metacor.case.merge$edgestrength)

sparcc.crc.metacor.case.merge.pval <- merge(subset(sparcc.crc.metacor.case.result$cor.merge, 
                                                   random.pvalue <0.05),
                                            sparcc.crc.merge.sig[,c("features","variable","target")] %>% 
                                              distinct(target, .keep_all = T),
                                            by.x = "target_real", by.y = "target",
                                            all.x = T)

sparcc.crc.metacor.case.kw <- sparcc.crc.metacor.case.merge %>%
  filter(disease %in% "CRC") %>%
  filter(features %in% crc.kw.spe) %>% 
  filter(variable %in% crc.kw.spe)

add.spe <- crc.kw.spe[!crc.kw.spe %in% unique(c(sparcc.crc.metacor.case.kw$features, 
                                              sparcc.crc.metacor.case.kw$variable))]

## control data
sparcc.crc.metacor.ctr.merge <- merge(subset(sparcc.crc.metacor.ctr.result$cor.merge, fdr <0.05),
                                      sparcc.crc.merge.sig[,c("features","variable","target")] %>% 
                                        distinct(target, .keep_all = T),
                                      by.x = "target_real", by.y = "target",
                                      all.x = T)

sparcc.crc.metacor.ctr.merge$edgestrength <- cut_edge(sparcc.crc.metacor.ctr.merge$random.r)
sparcc.crc.metacor.ctr.merge$edgetrend <- 1
sparcc.crc.metacor.ctr.merge[which(sparcc.crc.metacor.ctr.merge$random.r <0), "edgetrend"] <- -1

sparcc.crc.metacor.ctr.merge$variable <- as.character(sparcc.crc.metacor.ctr.merge$variable)
sparcc.crc.metacor.ctr.merge$edgestrength <- as.character(sparcc.crc.metacor.ctr.merge$edgestrength)


sparcc.crc.metacor.ctr.merge.pval <- merge(subset(sparcc.crc.metacor.ctr.result$cor.merge, 
                                                  random.pvalue <0.05),
                                           sparcc.crc.merge.sig[,c("features","variable","target")] %>% 
                                             distinct(target, .keep_all = T),
                                           by.x = "target_real", by.y = "target",
                                           all.x = T)

sparcc.crc.metacor.ctr.kw <-  sparcc.crc.metacor.ctr.merge %>%
  filter(disease %in% "CRC") %>%
  filter(features %in% crc.kw.spe) %>% 
  filter(variable %in% crc.kw.spe)

add.spe <- crc.kw.spe[!crc.kw.spe %in% unique(c(sparcc.crc.metacor.ctr.kw$features, 
                                                sparcc.crc.metacor.ctr.kw$variable))]


for (add in add.spe) {
  sparcc.crc.metacor.ctr.kw[nrow(sparcc.crc.metacor.ctr.kw)+1, 
                                    ] <- c(NA, NA, "CRC",
                                           0, 0, 0, 
                                           add, add, 0, 0)
}

write.table(sparcc.crc.metacor.ctr.kw,
              file = paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/Network/",
                            "CRC_CTR_network.txt"),
              row.names = F,
              col.names = T, sep = "\t",
              quote = F)
write.table(sparcc.crc.metacor.case.kw,
            file = paste0("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/Network/",
                          "CRC_CRC_network.txt"),
            row.names = F,
            col.names = T, sep = "\t",
            quote = F)

## merge the data
sparcc.metacor.ctr.kw <- list(CRC = sparcc.crc.metacor.ctr.kw, 
                              CD = sparcc.ibd.metacor.ctr.kw$CD, 
                              UC = sparcc.ibd.metacor.ctr.kw$UC)
sparcc.metacor.case.kw <- list(CRC = sparcc.crc.metacor.case.kw, 
                               CD = sparcc.ibd.metacor.case.kw$CD, 
                               UC = sparcc.ibd.metacor.case.kw$UC)



