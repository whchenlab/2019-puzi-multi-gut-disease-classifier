## Preparing abundance data from MetaPhlAn2
library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(cowplot)

data.phylum.func <- function(data){
  ##-----data为metaphlan2生成的结果，i为sample名字长度
  ##colnames(data) <- substr(colnames(data),1,i)
  clade <- data$ID
  ##-----将ID替换成具体的分类名函数--------
  parse_taxonomy <- function (clade, derep = TRUE) {
    # The clade string has the form
    # "k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanosphaera|s__Methanosphaera_stadtmanae|t__GCF_000012545"
    # truncated at whatever the smallest rank is. We want to be able to parse
    # the clade string for any possible smallest rank, from Kingdom to Strain.
    rank_letters <- c("k", "p", "c", "o", "f", "g", "s", "t")
    tax_pattern  <- paste0("(?:", rank_letters, "__(\\w+))?") %>%
      paste(collapse = "\\|?")
    if (derep) {
      tax <- clade %>%
        unique %>%
        stringr::str_match(tax_pattern)
    } else {
      tax <- clade %>%
        stringr::str_match(tax_pattern)
    }
    colnames(tax) <- c("Clade", "Kingdom", "Phylum", "Class", "Order", "Family",
                       "Genus", "Species", "Strain")
    tax %>% as_tibble
  }
  ##-----
  clade.result <- parse_taxonomy(clade)
  data.merge <- cbind(clade.result, data[,!colnames(data) %in% "ID"])
  ##-----生成phylum.to.genus为后面画热图做准备
  
  phylum.to.genus <- subset(clade.result,is.na(Species) & !(is.na(Genus)) & Kingdom %in% "Bacteria" , select = c(Phylum, Genus))
  phylum.to.genus <- phylum.to.genus %>% distinct()
  annot.genus <- data.frame(phylum = phylum.to.genus$Phylum, stringsAsFactors = F)
  rownames(annot.genus) <- phylum.to.genus$Genus
  
  ##-----生成phylum.to.species为后面画热图做准备
  
  phylum.to.species <- subset(clade.result, is.na(Strain) & !(is.na(Species)) & Kingdom %in% "Bacteria" , select = c(Phylum, Species))
  phylum.to.species <- phylum.to.species %>% distinct()
  annot.species <- data.frame(phylum = phylum.to.species$Phylum, 
                              stringsAsFactors = F)
  rownames(annot.species) <- phylum.to.species$Species
  
  ##-----phylum水平与genus水平相对丰度
  data.abund.p <- subset(data.merge, is.na(Class) & (!is.na(Phylum)) & Kingdom %in% "Bacteria" )[,-c(1:2,4:9)]
  data.abund.g <- subset(data.merge, is.na(Species) & !(is.na(Genus)) & Kingdom %in% "Bacteria" )[,-c(1:6,8:9)]
  data.abund.c <- subset(data.merge, is.na(Order) & (!is.na(Class)) & Kingdom %in% "Bacteria" )[,-c(1:4,6:9)]
  data.abund.o <- subset(data.merge, is.na(Family) & (!is.na(Order)) & Kingdom %in% "Bacteria" )[,-c(1:5,7:9)]
  data.abund.f <- subset(data.merge, is.na(Genus) & (!is.na(Family)) & Kingdom %in% "Bacteria" )[,-c(1:6,8:9)]
  data.abund.s <- subset(data.merge, (!is.na(Species)) & is.na(Strain) & Kingdom %in% "Bacteria" )[,-c(1:7,9)]
  
  ##-----relative abundant plot for phylum above 0.001%
  ##-----selecte mean phyla above .1% in all samples
  data.abund.p.melt <- melt(data.abund.p, id.vars = "Phylum")
  data.abund.p.melt.stats <- data.abund.p.melt %>% group_by(Phylum) %>% 
    summarise(PercentPhylum = mean(value), 
              N = length(value), sd = sd(value), 
              se = sd/sqrt(N))
  ##-----筛选相对丰度大于0.001的门
  data.abund.p.melt.stats.1 <- subset(data.abund.p.melt.stats, PercentPhylum > 0.1)
  data.abund.p.melt.stats.1 <- data.abund.p.melt.stats.1[order(data.abund.p.melt.stats.1$PercentPhylum, decreasing = T),]
  annot.genus[-which(annot.genus$phylum %in% data.abund.p.melt.stats.1$Phylum),"phylum"] <- "low.phylum"
  annot.species[-which(annot.species$phylum %in% data.abund.p.melt.stats.1$Phylum),"phylum"] <- "low.phylum"
  
  phylum.mean.relaabund.plot <- ggplot(data.abund.p.melt.stats.1, aes(x = Phylum, y = PercentPhylum)) + 
    geom_bar(stat = 'identity', position = position_dodge(), fill = "cornflowerblue", color = "black") + 
    theme_bw() + xlab("Phylum") + ylab("Mean Relative Abundance (%)") + 
    geom_errorbar(aes(ymin = PercentPhylum - se, ymax = PercentPhylum + se), width = 0.25) + coord_flip() + 
    theme(axis.text = element_text(size = 10, face = "bold"), 
          axis.title = element_text(face = "bold", size = 12), 
          plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
  return(list(data.abund.p = data.abund.p, data.abund.g = data.abund.g,
              data.abund.c = data.abund.c, data.abund.o = data.abund.o,
              data.abund.f = data.abund.f, data.abund.s = data.abund.s,
              phylum.to.genus = annot.genus, data.merge = data.merge,
              phylum.to.species = annot.species,
              phylum.mean.relaabund.plot = phylum.mean.relaabund.plot
  ))
}

## import taxonomic abundance data
setwd("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/abund_data")

crc.project <- c("PRJEB6070", "PRJDB4176", "PRJEB10878",
                 "PRJEB12449", "PRJEB27928",
                 "PRJEB7774","PRJNA447983")
cd.project <- c("PRJNA389280", "PRJNA400072", "SRP057027")
uc.project <- c("PRJNA389280", "PRJNA400072", "PRJEB1220")

##
abund.list <- list()

for(i in c("PRJEB6070", "PRJDB4176", "PRJEB10878",
           "PRJEB12449", "PRJEB27928",
           "PRJEB7774","PRJNA447983",
           "PRJNA389280", "PRJNA400072", "PRJEB1220")){
  abund.list[[i]] <- read.delim(file = paste0(i, "_abund.txt"),
                                header = T, sep = "\t", as.is = T)
}
abund.list[["SRP057027"]] <- read.csv("SRP057027_abund.csv",
                                      header = T)
names(abund.list[["SRP057027"]])[1] <- "ID"

names(abund.list[["PRJNA389280"]]) <- unlist(lapply(strsplit(names(abund.list[["PRJNA389280"]]),
                                               split = ".",
                                               fixed = T),
                                       function(data){y <- data[[1]]}))

for (i in c("PRJEB6070", "PRJEB12449", 
            "PRJEB27928", "PRJEB7774",
            "PRJNA400072", "PRJEB1220")) {
  names(abund.list[[i]]) <- unlist(lapply(strsplit(names(abund.list[[i]]),
                                                               split = "_meta",
                                                               fixed = T),
                                                      function(data){y <- data[[1]]}))
  
}

for (i in c("PRJNA447983", "PRJDB4176", 
            "PRJEB10878")) {
  names(abund.list[[i]]) <- unlist(lapply(strsplit(names(abund.list[[i]]),
                                                   split = ".meta",
                                                   fixed = T),
                                          function(data){y <- data[[1]]}))
  
}

## keep the consistence of the abundance data and metadata
names(abund.list$PRJEB6070) <- gsub(names(abund.list$PRJEB6070), 
                                      pattern = ".", replacement = "-",
                                      fixed = T)

names(abund.list$PRJEB12449) <- gsub(names(abund.list$PRJEB12449), 
                                       pattern = ".", replacement = "-",
                                       fixed = T)

names(abund.list$PRJEB27928) <- gsub(names(abund.list$PRJEB27928), 
                                       pattern = ".", replacement = "-",
                                       fixed = T)
## 
PRJNA447983.run <- read.delim("/mnt/raid5/puzi/IBD/otherproject/metadata/prjna447983_runinfo.txt",
                              header = T, sep = "\t", as.is = T)
rownames(PRJNA447983.run) <- PRJNA447983.run$Run
names(abund.list$PRJNA447983) <- c("ID", 
                                   PRJNA447983.run[(names(abund.list$PRJNA447983))[-1],]$Sample_Name)
rm(PRJNA447983.run)

##
PRJDB4176.run <- read.delim("/mnt/raid5/puzi/IBD/otherproject/metadata/prjdb4176_runinfo.txt",
                            header = T, sep = "\t", as.is = T)
rownames(PRJDB4176.run) <- PRJDB4176.run$Run
names(abund.list$PRJDB4176) <- c("ID", 
                                 PRJDB4176.run[(names(abund.list$PRJDB4176))[-1],]$Sample_Name)
rm(PRJDB4176.run)

## check names
## abund.names <- lapply(abund.list, names) %>% reduce(c)
## length(intersect(abund.names, crc.meta$Sample_ID))
## length(intersect(abund.names, ibd.meta$Sample_ID))
## rm(abund.names)

abund.list.2 <- lapply(abund.list, data.phylum.func)
abund.s.list <- lapply(abund.list.2, function(data){
  y <- data$data.abund.s
})

##
metadata <- rbind(crc.meta, ibd.meta)
## filter samples
for (i in names(abund.s.list)) {
  abund.s.list[[i]] <- abund.s.list[[i]][, names(abund.s.list[[i]]) %in% c("Species", metadata$Sample_ID)]
}
lapply(abund.s.list, function(data){y <- apply(data[,-1], 2, sum)[apply(data[,-1], 2, sum)<=50]})

abund.s.list[["PRJEB1220"]] <- abund.s.list[["PRJEB1220"]][,!names(abund.s.list[["PRJEB1220"]]) %in% c("MH0034","MH0036","MH0039","MH0042","MH0046", "O2.UC16.0", "V1.UC14.0")]
metadata <- metadata %>%
  filter(!Sample_ID %in% c("MH0034","MH0036","MH0039","MH0042","MH0046", "O2.UC16.0", "V1.UC14.0"))
ibd.meta <- ibd.meta %>%
  filter(!Sample_ID %in% c("MH0034","MH0036","MH0039","MH0042","MH0046", "O2.UC16.0", "V1.UC14.0"))

cd.meta <- ibd.meta %>%
  filter(Project %in% cd.project & Group %in% c("CTR", "CD"))
uc.meta <- ibd.meta %>%
  filter(Project %in% uc.project & Group %in% c("CTR", "UC"))

abund.names <- lapply(abund.s.list, names) %>% 
  reduce(c) %>%
  unique
length(intersect(abund.names, metadata$Sample_ID))

## filter the species with low abundance
for (i in names(abund.s.list)) {
  abund.s.list[[i]][,-1] <- abund.s.list[[i]][,-1]/100
}

abund.s.df <- abund.s.list %>% reduce(merge, by = "Species", all = T)
abund.s.df[is.na(abund.s.df)] <- 0
rownames(abund.s.df) <- abund.s.df$Species


abund.s.df.crc <- abund.s.df[, crc.meta$Sample_ID]
abund.s.df.cd <- abund.s.df[, ibd.meta[which(ibd.meta$Project %in% cd.project & ibd.meta$Group %in% c("CTR", "CD")), "Sample_ID"]]
abund.s.df.uc <- abund.s.df[, ibd.meta[which(ibd.meta$Project %in% uc.project & ibd.meta$Group %in% c("CTR", "UC")), "Sample_ID"]]

## filter species with low relative abundance
## filter those species with relative abundance <0.001 in at least two datasets of a specified disease
## filter crc
spe.max.crc <- t(sapply(row.names(abund.s.df.crc),
                       FUN=function(marker){sapply(unique(crc.meta$Project),
                                                   FUN=function(study, marker){
                                                     max.ab = max(abund.s.df.crc[marker, which(crc.meta$Project == study)])
                                                   },
                                                   marker=marker)}))

f.idx <- rowSums(spe.max.crc >= 0.001) >= 4 &
  row.names(spe.max.crc) != '-1'
abund.s.crc.filter <- abund.s.df.crc[f.idx, ]

## filter cd
spe.max.cd <- t(sapply(row.names(abund.s.df.cd),
                        FUN=function(marker){sapply(unique(cd.meta$Project),
                                                    FUN=function(study, marker){
                                                      max.ab = max(abund.s.df.cd[marker, which(cd.meta$Project == study)])
                                                    },
                                                    marker=marker)}))

f.idx <- rowSums(spe.max.cd >= 0.001) >= 2 &
  row.names(spe.max.cd) != '-1'
abund.s.cd.filter <- abund.s.df.cd[f.idx, ]

## filter uc
spe.max.uc <- t(sapply(row.names(abund.s.df.uc),
                       FUN=function(marker){sapply(unique(uc.meta$Project),
                                                   FUN=function(study, marker){
                                                     max.ab = max(abund.s.df.uc[marker, which(uc.meta$Project == study)])
                                                   },
                                                   marker=marker)}))

f.idx <- rowSums(spe.max.uc >= 0.001) >= 2 &
  row.names(spe.max.uc) != '-1'
abund.s.uc.filter <- abund.s.df.uc[f.idx, ]

abund.s.list.filter <- list(CRC = abund.s.crc.filter, 
                            CD = abund.s.cd.filter, 
                            UC = abund.s.uc.filter)
abund.s.list.filter <- lapply(abund.s.list.filter, 
                              function(data){data$features <- rownames(data);
                               return(data)})

## for phylum to species mapping file 
phylum.to.species <- lapply(abund.list.2, function(data){
  y <- data$data.merge %>% 
    filter(Kingdom %in% "Bacteria") %>%
    select(Phylum, Species) %>%
    filter(!is.na(Species));
  return(y)
}) %>% reduce(rbind)
phylum.to.species <- phylum.to.species[!duplicated(phylum.to.species$Species), ]


save(abund.s.list.filter, 
     phylum.to.species, 
     metadata, 
     crc.project, 
     cd.project, 
     uc.project, 
     file = "Filtered_Spe_abund.Rdata")

## import pathways abundance data
## pathways abundance
path.abund.list <- list()

 for(i in c("PRJEB6070", "PRJDB4176", "PRJEB10878",
           "PRJEB12449", "PRJEB27928",
           "PRJEB7774","PRJNA447983",
           "PRJNA389280", "PRJNA400072", 
           "PRJEB1220", "SRP057027")){
  path.abund.list[[i]] <- read.delim(file = paste0(i, "_rela_pathabund.txt"),
                                header = T, sep = "\t", as.is = T)
}
lapply(path.abund.list,function(data){y <- data[1:5,1:5]})
lapply(path.abund.list, names)

path.abund.list[["PRJEB10878"]]$Pathway <- rownames(path.abund.list[["PRJEB10878"]])
path.abund.list[["PRJNA447983"]]$Pathway <- rownames(path.abund.list[["PRJNA447983"]])

for(i in c("PRJEB6070", "PRJDB4176", 
           "PRJEB12449", "PRJEB27928",
           "PRJEB7774", 
           "PRJNA389280", "PRJNA400072", 
           "PRJEB1220", "SRP057027")){
  names(path.abund.list[[i]])[1] <- "Pathway"
}


for(i in c("PRJEB6070", "PRJDB4176", "PRJEB10878",
           "PRJEB12449", "PRJEB27928",
           "PRJEB7774","PRJNA447983",
           "PRJNA389280", "PRJNA400072", 
           "PRJEB1220", "SRP057027")){
  names(path.abund.list[[i]]) <- lapply(strsplit(names(path.abund.list[[i]]), 
                                                 split = "_humann2", 
                                                 fixed = T), function(data){y <- data[1]}) %>% unlist
  names(path.abund.list[[i]]) <- gsub(names(path.abund.list[[i]]), 
                                      pattern = ".", replacement = "-", 
                                      fixed = T)
}

## rename the column names of the file
names(path.abund.list[["PRJEB1220"]]) <- gsub(names(path.abund.list[["PRJEB1220"]]), 
                                             pattern = "-", replacement = ".", 
                                             fixed = T)
names(path.abund.list[["PRJNA447983"]])[1:80] <- lapply(strsplit(names(path.abund.list[["PRJNA447983"]])[1:80], 
                                                           split = "_", 
                                                           fixed = T), function(data){y <- data[2]}) %>% unlist


PRJEB10878.run <- read.delim("/mnt/raid5/puzi/IBD/otherproject/metadata/prjeb10878_runinfo.txt",
                             header = T, sep = "\t", as.is = T)
rownames(PRJEB10878.run) <- PRJEB10878.run$BioSample
names(path.abund.list[["PRJEB10878"]])[1:128] <- PRJEB10878.run[names(path.abund.list[["PRJEB10878"]])[1:128], "Run"]
rm(PRJEB10878.run)

PRJDB4176.run <- read.delim("/mnt/raid5/puzi/IBD/otherproject/metadata/prjdb4176_runinfo.txt",
                            header = T, sep = "\t", as.is = T)
rownames(PRJDB4176.run) <- PRJDB4176.run$Run
names(path.abund.list[["PRJDB4176"]])[2:81] <- PRJDB4176.run[names(path.abund.list[["PRJDB4176"]])[2:81], "BioSample"]
rm(PRJDB4176.run)

metadata$Sample_ID[ !metadata$Sample_ID %in% (lapply(path.abund.list, names) %>% unlist)]
length((lapply(path.abund.list, names) %>% unlist))

## filter matched samples
path.abund.list <- lapply(path.abund.list, 
                          function(data){
                            y <- data[, names(data) %in% c("Pathway", metadata$Sample_ID)]; 
                            return(y)
                          })
path.abund.list <- lapply(path.abund.list, 
                          function(data){
                            rownames(data) <- data$Pathway; 
                            data <- data[, !names(data) %in% "Pathway"];
                            return(data)
                          })

path.abund.filter <- lapply(path.abund.list, 
                            function(x){
  y <- x[-grep(rownames(x), pattern = "|", fixed = T),]
})

lapply(path.abund.filter, dim)

## filter pathways with low abundance
noise.removal.path <- function(data, top=NULL){
  Matrix <- data
  zero.num <- apply(Matrix, 1, function(data){y <- length(data[data ==0])})
  col.num <- ncol(data)
  bigones <- zero.num<=0.15*col.num
  Matrix_1 <- Matrix[bigones,]
  return(Matrix_1)
}

path.abund.filter.2 <- lapply(path.abund.filter, function(x){
  y <- noise.removal.path(x)
  y <- y[!rownames(y) %in% c("UNMAPPED","UNINTEGRATED"),];
  y$Pathway <- rownames(y)
  return(y)
})

path.abund.fil.df <- path.abund.filter.2 %>% 
  reduce(merge, by = "Pathway", all = T)
path.abund.fil.df[is.na(path.abund.fil.df)] <- 0
rownames(path.abund.fil.df) <- path.abund.fil.df$Pathway


path.abund.fil.crc <- path.abund.fil.df[, crc.meta$Sample_ID]
path.abund.fil.cd <- path.abund.fil.df[, ibd.meta[which(ibd.meta$Project %in% cd.project & ibd.meta$Group %in% c("CTR", "CD")), "Sample_ID"]]
path.abund.fil.uc <- path.abund.fil.df[, ibd.meta[which(ibd.meta$Project %in% uc.project & ibd.meta$Group %in% c("CTR", "UC")), "Sample_ID"]]

## filter out pathways with low abundance
## filter crc
path.max.crc <- t(sapply(row.names(path.abund.fil.crc),
                        FUN=function(marker){
                          sapply(unique(crc.meta$Project),
                                 FUN=function(study, marker){
                                   max.ab = max(path.abund.fil.crc[marker, which(crc.meta$Project == study)])
                                   },
                                 marker=marker)}))

f.idx <- rowSums(path.max.crc >= 0.000001) >= 4 &
  row.names(path.max.crc) != '-1'
path.abund.fil.crc.2 <- path.abund.fil.crc[f.idx, ]

## filter cd
path.max.cd <- t(sapply(row.names(path.abund.fil.cd),
                       FUN=function(marker){sapply(unique(cd.meta$Project),
                                                   FUN=function(study, marker){
                                                     max.ab = max(path.abund.fil.cd[marker, which(cd.meta$Project == study)])
                                                   },
                                                   marker=marker)}))

f.idx <- rowSums(path.max.cd >= 0.000001) >= 2 &
  row.names(path.max.cd) != '-1'
path.abund.fil.cd.2 <- path.abund.fil.cd[f.idx, ]

## filter uc
path.max.uc <- t(sapply(row.names(path.abund.fil.uc),
                       FUN=function(marker){sapply(unique(uc.meta$Project),
                                                   FUN=function(study, marker){
                                                     max.ab = max(path.abund.fil.uc[marker, which(uc.meta$Project == study)])
                                                   },
                                                   marker=marker)}))

f.idx <- rowSums(path.max.uc >= 0.000001) >= 2 &
  row.names(path.max.uc) != '-1'
path.abund.fil.uc.2 <- path.abund.fil.uc[f.idx, ]

## 
path.abund.fil.list <- list(CRC = path.abund.fil.crc.2, 
                            CD = path.abund.fil.cd.2, 
                            UC = path.abund.fil.uc.2)
path.abund.fil.list <- lapply(path.abund.fil.list, 
                              function(data){ data$features <- rownames(data);
                              return(data)})
save(path.abund.fil.list, 
     metadata, 
     crc.project, 
     cd.project, 
     uc.project, 
     file = "Filtered_Pathways_abund.Rdata")
