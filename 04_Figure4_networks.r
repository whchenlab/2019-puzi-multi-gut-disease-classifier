## igraph
library(igraph)

setwd("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/Network/mcode")
file.list.1 <- dir()[-grep(dir(), pattern = "edge")]
file.list.2 <- dir()[grep(dir(), pattern = "edge")]

mcode.list <- list()
for (file in file.list.1) {
  mcode.list[[gsub(file, pattern = ".csv", 
                   replacement = "", fixed = T)]] <- read.csv(file = file,
                                                              header = T, as.is = T)
}

mcode.edge.list <- list()
for (file in file.list.2) {
  mcode.edge.list[[gsub(file, pattern = "_edge.csv", 
                        replacement = "", fixed = T)]] <- read.csv(file = file,
                                                                   header = T, as.is = T)
}

cut_preva <- function(p) {
  out <- cut(p, 
             breaks = c(0, 1/3, 2/3, 1), 
             include.lowest = T,
             labels = c("1.5", "2.5", "3.5"))
  return(out)
}

for (file in names(mcode.edge.list)) {
  mcode.edge.list[[file]]$features <- lapply(strsplit(mcode.edge.list[[file]]$name, 
                                                      split = " (interacts with) ", fixed = T),
                                             function(data){ y <- data[1]}) %>% unlist
  mcode.edge.list[[file]]$variable <- lapply(strsplit(mcode.edge.list[[file]]$name, 
                                                      split = " (interacts with) ", fixed = T),
                                             function(data){ y <- data[2]}) %>% unlist
  mcode.edge.list[[file]]$features <- gsub(mcode.edge.list[[file]]$features, 
                                           pattern = "_", replacement = " ", fixed = T)
  mcode.edge.list[[file]]$variable <- gsub(mcode.edge.list[[file]]$variable, 
                                           pattern = "_", replacement = " ", fixed = T)
  
  mcode.edge.list[[file]]$features.2 <- mcode.edge.list[[file]]$features 
  mcode.edge.list[[file]]$variable.2 <- mcode.edge.list[[file]]$variable 
  
  
  for(i in 1:nrow(mcode.edge.list[[file]])){
    ## from
    feature <- mcode.edge.list[[file]][i, "features"]
    a <- gregexpr(feature, 
                  pattern = " ")[[1]][1] +1
    b <- nchar(feature)
    mcode.edge.list[[file]][i, "features.2"] <- paste(substr(feature, 1, 1), 
                                                      substr(feature, a, b), 
                                                      sep = ". ")
    ## to 
    feature <- mcode.edge.list[[file]][i, "variable"]
    a <- gregexpr(feature, 
                  pattern = " ")[[1]][1] +1
    b <- nchar(feature)
    mcode.edge.list[[file]][i, "variable.2"] <- paste(substr(feature, 1, 1), 
                                                      substr(feature, a, b), 
                                                      sep = ". ")
  }
  
  mcode.list[[file]]$name.2 <- mcode.list[[file]]$name
  mcode.list[[file]]$name.2 <- gsub(mcode.list[[file]]$name, 
                                    pattern = "_", replacement = " ", 
                                    fixed = T)
  for(i in 1:nrow(mcode.list[[file]])){
    feature <- mcode.list[[file]][i, "name.2"]
    a <- gregexpr(feature, 
                  pattern = " ")[[1]][1] + 1
    b <- nchar(feature)
    mcode.list[[file]][i, "name.2"] <- paste(substr(feature, 1, 1), 
                                             substr(feature, a, b), 
                                             sep = ". ")
  }
 
  rownames(mcode.list[[file]]) <- mcode.list[[file]]$name.2
  mcode.list[[file]]$Preva_group <- cut_preva(mcode.list[[file]]$Prevalence)
}

## convert dataframe to igraph format
data.frame <- mcode.edge.list$CD_case_cluster1

convert.igraph <- function(file = file, 
                           mcode.edge.list = mcode.edge.list, 
                           mcode.list = mcode.list){
  data.frame <- mcode.edge.list[[file]]
  a <- data.frame %>%
    filter(!features == variable) %>%
    select(features.2, variable.2, edgestrength, edgetrend)
  
  a$edgestrength <- as.numeric(a$edgestrength)
  #a$random.r <- abs(a$random.r)
  names(a) <- c("from", "to", "weight", "edgetrend")
  
  ## convert to igraph format
  data <- graph.data.frame(a[, 1:3], directed = FALSE)
  ## define the edge strength as the weight
  E(data)$weight <- a$weight
  ## define the trend of the edge
  E(data)$trend <- a$edgetrend
  ## define the type of the vectors
  V(data)$Trend <- mcode.list[[file]][names(V(data)), "Type"]
  return(list(igraph.data = data))
}

mcode.igraph.list <- list()
for (i in names(mcode.list)) {
  mcode.igraph.list[[i]] <- convert.igraph(file = i, 
                                           mcode.edge.list = mcode.edge.list, 
                                           mcode.list = mcode.list)
}


## -----see layouts
layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1] 
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
  
 
## gold:#ffd700 CRC
## green:#B3DE69 CD
## blue:#50C1E9 UC

## plot module
for (file in names(mcode.list)) {
  
  a <- mcode.igraph.list[[file]]$igraph.data;
  V(a)$Trend[which(V(a)$Trend=="Case-enriched")] <- 2
  V(a)$Trend[which(V(a)$Trend=="Control-enriched")] <- 1
  E(a)$trend[E(a)$trend == "-1"] <- 2
  
  ordering <- as.numeric(V(a))
  names(ordering) <- names(degree(a))
  
  ordering2 <- mcode.list[[file]] %>% 
    filter(name.2 %in% names(ordering)) %>%
    arrange(Type, name.2) %>% pull(name.2)
  
  ordering <- ordering[ordering2]
  
  
  if (grepl(file, pattern = "CD")) {
    ## CD
    vcolor <- "#B3DE69";
    l <- layout_in_circle(a, order = as.numeric(ordering));
  }else if (grepl(file, pattern = "UC")) { 
    ## UC
    vcolor <- "#50C1E9";
    l <- layout_in_circle(a, order = as.numeric(ordering));
  }else {
    ## CRC
    vcolor <- "#ffd700";
    l <- layout_in_circle(a, order = as.numeric(ordering));
  }
  pdf(paste0(file, ".pdf"), width = 5, height = 5)
  plot(a, 
       vertex.label.color = "black",
       vertex.frame.color = "white",
       vertex.label.dist=1.5, 
       vertex.label.size = 9,
       vertex.color=c( "slategrey", vcolor)[as.numeric(V(a)$Trend)] , 
       vertex.size = as.numeric(as.character(mcode.list[[file]][names(V(a)), "Preva_group"])) *3, 
       edge.width = E(a)$weight,
       edge.color = c( "gray", "#FF5F5F")[as.numeric(E(a)$trend)],
       #edge.curved = 0.1, 
       layout = l*0.4,
       vertex.label.family = "Times", 
       vertex.label.font = 3,
       rescale = FALSE)
  dev.off()
}

## plot module
for (file in names(mcode.list)) {
  
  a <- mcode.igraph.list[[file]]$igraph.data;
  V(a)$Trend[which(V(a)$Trend=="Case-enriched")] <- 2
  V(a)$Trend[which(V(a)$Trend=="Control-enriched")] <- 1
  E(a)$trend[E(a)$trend == "-1"] <- 2
  
  ordering <- as.numeric(V(a))
  names(ordering) <- names(degree(a))
  
  ordering2 <- mcode.list[[file]] %>% 
    filter(name.2 %in% names(ordering)) %>%
    arrange(Type, name.2) %>% pull(name.2)
  
  ordering <- ordering[ordering2]
  
  
  if (grepl(file, pattern = "CD")) {
    ## CD
    vcolor <- "#B3DE69";
    l <- layout_in_circle(a, order = as.numeric(ordering));
  }else if (grepl(file, pattern = "UC")) { 
    ## UC
    vcolor <- "#50C1E9";
    l <- layout_in_circle(a, order = as.numeric(ordering));
  }else {
    ## CRC
    vcolor <- "#ffd700";
    l <- layout_in_circle(a, order = as.numeric(ordering));
  }
  pdf(paste0(file, ".pdf"), width = 5, height = 5)
  plot(a, 
       vertex.label.color = "black",
       vertex.frame.color = "white",
       vertex.label.dist=1.5, 
       vertex.label.size = 9,
       vertex.color=c( "slategrey", vcolor)[as.numeric(V(a)$Trend)] , 
       vertex.size = as.numeric(as.character(mcode.list[[file]][names(V(a)), "Preva_group"])) *3, 
       edge.width = E(a)$weight,
       edge.color = c( "gray", "#FF5F5F")[as.numeric(E(a)$trend)],
       #edge.curved = 0.1, 
       layout = l*0.4,
       vertex.label.family = "Times", 
       vertex.label.font = 3,
       rescale = FALSE)
  dev.off()
}

setwd("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/Network")
network.file.list <- dir()[intersect(grep(dir(), pattern = "^UC_C|^CD_C|^CRC_C"), 
                             grep(dir(), pattern = "network"))]

network.list <- list()
for (file in network.file.list) {
  network.list[[gsub(file, pattern = "_network.txt", 
                        replacement = "", fixed = T)]] <- read.delim(file = file, 
                                                                     sep= "\t", header = T, as.is = T)
}

network.prev.file <- dir()[grep(dir(), pattern = "prevalence")]

network.prev.list <- list()
for (file in network.prev.file) {
  network.prev.list[[gsub(file, pattern = "_prevalence.txt", 
                     replacement = "", fixed = T)]] <- read.delim(file = file, 
                                                                  sep= "\t", header = T, as.is = T)
}
names(network.prev.list) <- c("CD_CTR", "CD_Case", "CRC_CTR", "CRC_CRC", "UC_CTR", "UC_Case")


for (file in names(network.list)) {
  ## prevalence 
  network.prev.list[[file]]$Preva_group <- cut_preva(network.prev.list[[file]]$Prevalence)
  network.prev.list[[file]]$name.2 <- network.prev.list[[file]]$Species
  network.prev.list[[file]]$name.2 <- gsub(network.prev.list[[file]]$Species, 
                                           pattern = "_", replacement = " ", 
                                           fixed = T)
  for(i in 1:nrow(network.prev.list[[file]])){
    feature <- network.prev.list[[file]][i, "name.2"]
    a <- gregexpr(feature, 
                  pattern = " ")[[1]][1] + 1
    b <- nchar(feature)
    network.prev.list[[file]][i, "name.2"] <- paste(substr(feature, 1, 1), 
                                             substr(feature, a, b), 
                                             sep = ". ")
  } 
 
  #rownames(network.prev.list[[file]]) <- network.prev.list[[file]]$name.2
  
  ## network 
  network.list[[file]]$features <- gsub(network.list[[file]]$features, 
                                        pattern = "_", replacement = " ", fixed = T)
  network.list[[file]]$variable <- gsub(network.list[[file]]$variable, 
                                        pattern = "_", replacement = " ", fixed = T)
  
  network.list[[file]]$features.2 <- network.list[[file]]$features 
  network.list[[file]]$variable.2 <- network.list[[file]]$variable 
  
  
  for(i in 1:nrow(network.list[[file]])){
    ## from
    feature <- network.list[[file]][i, "features"]
    a <- gregexpr(feature, 
                  pattern = " ")[[1]][1] +1
    b <- nchar(feature)
    network.list[[file]][i, "features.2"] <- paste(substr(feature, 1, 1), 
                                                      substr(feature, a, b), 
                                                      sep = ". ")
    ## to 
    feature <- network.list[[file]][i, "variable"]
    a <- gregexpr(feature, 
                  pattern = " ")[[1]][1] +1
    b <- nchar(feature)
    network.list[[file]][i, "variable.2"] <- paste(substr(feature, 1, 1), 
                                                      substr(feature, a, b), 
                                                      sep = ". ")
  }
}

lapply(network.prev.list, function(data){ length(unique(data$name.2)) == nrow(data)})

## CRC_CTR
a <- network.prev.list$CRC_CTR$name.2[duplicated(network.prev.list$CRC_CTR$name.2)]

n <- 1
for (i in which(network.prev.list$CRC_CTR$name.2 == a)) {
  network.prev.list$CRC_CTR[i, "name.2"] <- paste0(a, "_", n)
  n <- n+1
}
network.list$CRC_CTR[which(network.list$CRC_CTR$features %in% "Parvimonas unclassified"), "features.2"] <- "P. unclassified_1"
network.list$CRC_CTR[which(network.list$CRC_CTR$variable %in% "Parvimonas unclassified"), "variable.2"] <- "P. unclassified_1"

network.list$CRC_CTR[which(network.list$CRC_CTR$features %in% "Peptostreptococcus unclassified"), "features.2"] <- "P. unclassified_2"
network.list$CRC_CTR[which(network.list$CRC_CTR$variable %in% "Peptostreptococcus unclassified"), "variable.2"] <- "P. unclassified_2"

## CRC_CRC
a <- network.prev.list$CRC_CRC$name.2[duplicated(network.prev.list$CRC_CRC$name.2)]

n <- 1
for (i in which(network.prev.list$CRC_CRC$name.2 == a)) {
  network.prev.list$CRC_CRC[i, "name.2"] <- paste0(a, "_", n)
  n <- n+1
}

network.list$CRC_CRC[which(network.list$CRC_CRC$features %in% "Parvimonas unclassified"), "features.2"] <- "P. unclassified_1"
network.list$CRC_CRC[which(network.list$CRC_CRC$variable %in% "Parvimonas unclassified"), "variable.2"] <- "P. unclassified_1"

network.list$CRC_CRC[which(network.list$CRC_CRC$features %in% "Peptostreptococcus unclassified"), "features.2"] <- "P. unclassified_2"
network.list$CRC_CRC[which(network.list$CRC_CRC$variable %in% "Peptostreptococcus unclassified"), "variable.2"] <- "P. unclassified_2"

## rownames
network.prev.list <- lapply(network.prev.list, 
                            function(data){ rownames(data) <- data$name.2; return(data)})

## network igraph
convert.igraph <- function(file = file, 
                           mcode.edge.list = mcode.edge.list, 
                           mcode.list = mcode.list){
  data.frame <- mcode.edge.list[[file]]
  a <- data.frame %>%
    filter(!features == variable) %>%
    select(features.2, variable.2, edgestrength, edgetrend)
  
  a$edgestrength <- as.numeric(a$edgestrength)
  #a$random.r <- abs(a$random.r)
  names(a) <- c("from", "to", "weight", "edgetrend")
  
  ## convert to igraph format
  data <- graph.data.frame(a[, 1:3], directed = FALSE)
  ## define the edge strength as the weight
  E(data)$weight <- a$weight
  ## define the trend of the edge
  E(data)$trend <- a$edgetrend
  ## add vectors
  add.v <- mcode.list[[file]]$name.2[!mcode.list[[file]]$name.2 %in% names(V(data))]
  data <- data + vertices(add.v)
  
  ## define the type of the vectors
  
  V(data)$Trend <- mcode.list[[file]][names(V(data)), "Type"]
  ## 
  V(data)$Trend[which(V(data)$Trend == "Case-enriched")] <- 2
  V(data)$Trend[which(V(data)$Trend == "Control-enriched")] <- 1
  ## trend == -1 means negative correlation
  E(data)$trend[E(data)$trend == "-1"] <- 2 
  return(list(igraph.data = data))
}

network.igraph.list <- list()
for (i in names(network.list)) {
  network.igraph.list[[i]] <- convert.igraph(file = i, 
                                           mcode.edge.list = network.list, 
                                           mcode.list = network.prev.list)
}


for (file in names(network.igraph.list)) {
  
  a <- network.igraph.list[[file]]$igraph.data;
 
  ordering <- as.numeric(V(a))
  names(ordering) <- names(degree(a))
  
  ordering2 <- network.prev.list[[file]] %>% 
    filter(name.2 %in% names(ordering)) %>%
    arrange(Type, name.2) %>% pull(name.2)
  
  ordering <- ordering[ordering2]
  
  if (grepl(file, pattern = "CD")) {
    vcolor <- "#B3DE69";
    
    l <- layout_in_circle(a, order = as.numeric(ordering));
    pdf(paste0(file, ".pdf"), width = 5, height = 5)
    plot(a, 
         vertex.label.color = "black",
         vertex.frame.color = "white",
         vertex.label.dist=1.5, 
         vertex.label.cex = 0.6,
         vertex.color=c( "slategrey", vcolor)[as.numeric(V(a)$Trend)] , 
         vertex.size = as.numeric(as.character(network.prev.list[[file]][names(V(a)), "Preva_group"])) *3, 
         edge.width = E(a)$weight ,
         edge.color = c( "gray", "#FF5F5F")[as.numeric(E(a)$trend)],
         #edge.curved = 0.1, 
         layout = l,
         vertex.label.family = "Times", 
         vertex.label.font = 3)
    dev.off()
  }else if (grepl(file, pattern = "UC")) {
    
    vcolor <- "#50C1E9";
    l <- layout_in_circle(a, order = as.numeric(ordering));
    pdf(paste0(file, ".pdf"), width = 5, height = 5)
    plot(a, 
         vertex.label.color = "black",
         vertex.frame.color = "white",
         vertex.label.dist=1.5, 
         vertex.label.cex = 0.7,
         vertex.color=c( "slategrey", vcolor)[as.numeric(V(a)$Trend)] , 
         vertex.size = as.numeric(as.character(network.prev.list[[file]][names(V(a)), "Preva_group"])) *3, 
         edge.width = E(a)$weight ,
         edge.color = c( "gray", "#FF5F5F")[as.numeric(E(a)$trend)],
         #edge.curved = 0.1, 
         layout = l,
         vertex.label.family = "Times", 
         vertex.label.font = 3)
    dev.off()
  }else {
    ## CRC
    
    vcolor <- "#ffd700";
    
    l <- layout_in_circle(a, order = as.numeric(ordering));
    pdf(paste0(file, ".pdf"), width = 5, height = 5)
    plot(a, 
         vertex.label.color = "black",
         vertex.frame.color = "white",
         vertex.label.dist=1.5, 
         vertex.label.cex = 0.6,
         vertex.color=c( "slategrey", vcolor)[as.numeric(V(a)$Trend)] , 
         vertex.size = as.numeric(as.character(network.prev.list[[file]][names(V(a)), "Preva_group"])) *3, 
         edge.width = E(a)$weight ,
         edge.color = c( "gray", "#FF5F5F")[as.numeric(E(a)$trend)],
         #edge.curved = 0.1, 
         layout = l,
         vertex.label.family = "Times", 
         vertex.label.font = 3)
    dev.off()
  }
}

save(network.igraph.list, 
     network.prev.list, 
     network.list, 
     mcode.list, 
     mcode.igraph.list, 
     mcode.edge.list, 
     file = "Figure_network.RData")
     
## centrality
## IGRAPH
## for betweenness 
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
    theme_bw() +
    labs(y = "Eigenvector Centrality Scores") + 
    scale_fill_manual(values = c("#a1d4e2", "#e62739")) + 
  theme(legend.position = "none",
        axis.text.x = element_text(vjust = .5, hjust = 1, angle = 90))
  return(list(spear.net.eigen.order = spear.net.eigen.order,
              spear.net.eigen.plot = spear.net.eigen.plot))
}

kw.01.spearman.cor.ctr.ig <- lapply(network.list, network.func) 

spear.net.betweenness <- lapply(kw.01.spearman.cor.ctr.ig, function(data){
                                    y <- data.frame(data$betweenness)
                                    names(y) <- "Betweenness"
                                    y$Species <- rownames(y)
                                    return(y)
                                  })
spear.net.eigen <- lapply(kw.01.spearman.cor.ctr.ig, function(data){
                              y <- data.frame(data$eigen.data$vector)
                              names(y) <- "Eigen"
                              y$Species <- rownames(y)
                              return(y)
                            })

names(spear.net.eigen)[3] <- "CRC_Case"
names(spear.net.betweenness)[3] <- "CRC_Case"

for (grp in names(spear.net.eigen)) {
  spear.net.eigen[[grp]]$type <- grp
  spear.net.betweenness[[grp]]$type <- grp
}

## Eigenvalue
spear.net.eigen.df <- spear.net.eigen %>% reduce(rbind)
spear.net.eigen.df$type2 <- factor(spear.net.eigen.df$type,
                                   levels = c( "CD_CTR", "CD_Case",
                                               "CRC_CTR", "CRC_Case",
                                               "UC_CTR","UC_Case"))

spear.net.eigen.plot <- list()
spear.net.eigen.plot[["CRC"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.eigen.df, 
                                                                                type2 %in% c("CRC_CTR","CRC_Case")))
spear.net.eigen.plot[["CRC"]]$spear.net.eigen.plot <- spear.net.eigen.plot[["CRC"]]$spear.net.eigen.plot +
  labs(title = "CRC", x = "") 

spear.net.eigen.plot[["CD"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.eigen.df, 
                                                                               type2 %in% c("CD_CTR","CD_Case")))
spear.net.eigen.plot[["CD"]]$spear.net.eigen.plot <- spear.net.eigen.plot[["CD"]]$spear.net.eigen.plot +
  labs(title = "CD", x = "") 

spear.net.eigen.plot[["UC"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.eigen.df, 
                                                                               type2 %in% c("UC_CTR","UC_Case")))
spear.net.eigen.plot[["UC"]]$spear.net.eigen.plot <- spear.net.eigen.plot[["UC"]]$spear.net.eigen.plot +
  labs(title = "UC", x = "") 


## Betweenness
spear.net.betweenness.df <- spear.net.betweenness %>% reduce(rbind)
spear.net.betweenness.df$type2 <- factor(spear.net.betweenness.df$type,
                                         levels = c( "CD_CTR", "CD_Case",
                                                     "CRC_CTR", "CRC_Case",
                                                     "UC_CTR","UC_Case"))

names(spear.net.betweenness.df)[1] <- "Eigen"
spear.net.bet.plot <- list() 
spear.net.bet.plot[["CRC"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.betweenness.df, 
                                                                              type2 %in% c("CRC_CTR","CRC_Case")))
spear.net.bet.plot[["CRC"]]$spear.net.eigen.plot <- spear.net.bet.plot[["CRC"]]$spear.net.eigen.plot +
  labs(y = "Betweenness Centrality Score", 
       title = "CRC", x = "")

spear.net.bet.plot[["CD"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.betweenness.df, 
                                                                             type2 %in% c("CD_CTR","CD_Case")))
spear.net.bet.plot[["CD"]]$spear.net.eigen.plot <- spear.net.bet.plot[["CD"]]$spear.net.eigen.plot +
  labs(y = "Betweenness Centrality Score", 
       title = "CD", x = "")

spear.net.bet.plot[["UC"]] <- network.eigen.plot(spear.net.eigen.df = subset(spear.net.betweenness.df, 
                                                                             type2 %in% c("UC_CTR","UC_Case")))
spear.net.bet.plot[["UC"]]$spear.net.eigen.plot <- spear.net.bet.plot[["UC"]]$spear.net.eigen.plot +
  labs(y = "Betweenness Centrality Score", 
       title = "UC", x = "") 

## final figures
centrality.plot <- list()
centrality.plot[["CRC"]] <- plot_grid(spear.net.eigen.plot[["CRC"]]$spear.net.eigen.plot, 
                                      spear.net.bet.plot[["CRC"]]$spear.net.eigen.plot, 
                                      rel_widths = c(1/2, 1/2), nrow = 1, 
                                      labels = c('A', 'B'),
                                      align = "v")
centrality.plot[["CD"]] <- plot_grid(spear.net.eigen.plot[["CD"]]$spear.net.eigen.plot, 
                                     spear.net.bet.plot[["CD"]]$spear.net.eigen.plot, 
                                     rel_widths = c(1/2, 1/2), nrow = 1, 
                                     labels = c('C', 'D'),
                                     align = "v")
centrality.plot[["UC"]] <- plot_grid(spear.net.eigen.plot[["UC"]]$spear.net.eigen.plot, 
                                     spear.net.bet.plot[["UC"]]$spear.net.eigen.plot, 
                                     rel_widths = c(1/2, 1/2), nrow = 1, 
                                     labels = c('E', 'F'),
                                     align = "v")
centrality.plot.all <- plot_grid(centrality.plot[["CRC"]], 
                                 centrality.plot[["CD"]], 
                                 centrality.plot[["UC"]], 
                                 rel_heights = c(1.1/3,1.1/3,0.8/3), ncol = 1, 
                                 align = "h")

pdf("../pdf/Network_centrality_plot.pdf", 
    width = 14, 
    height = 16)
centrality.plot.all
dev.off()

## -----numbers of correlations ------------##
spear.net.cornum <- lapply(network.list, 
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

spear.net.cornum.df$type <- c("Case","CTR", "Case","CTR", "Case","CTR")
spear.net.cornum.df$data <- c("CD","CD", 
                              "CRC" ,"CRC" ,
                              "UC","UC")
names(spear.net.cornum.df)[1:2] <- c("posi.num", "nega.num")

spear.net.cornum.df$posi.frac <- spear.net.cornum.df$posi.num/(spear.net.cornum.df$posi.num + spear.net.cornum.df$nega.num)
spear.net.cornum.df$posi.ratio <- c(76/22, 1, 172/138, 1, 9/1, 1)

spear.net.cornum.df$nega.frac <- spear.net.cornum.df$nega.num/(spear.net.cornum.df$posi.num + spear.net.cornum.df$nega.num)
spear.net.cornum.df$nega.ratio <- c(50/15, 1, 45/36, 1, 5/2, 1)

## 
## position correlation numbers ratio
spear.net.cornum.plot.1.ratio <- spear.net.cornum.df %>%  
  mutate(data = factor(data, levels = c("UC","CD","CRC"))) %>%
  mutate(type = factor(type, levels = c("CTR", "Case"))) %>%
  ggplot(aes(x = data, y = posi.ratio, fill = type)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#a1d4e2", "#e62739")) +
  labs(x = "Disease", 
       y = "Case-to-control ratio of positive correlation numbers") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) + 
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4)) +
  geom_text(aes(label = posi.num, y = posi.ratio -0.15), 
            position = position_dodge(0.9), size = 5,
            color = "white")

spear.net.cornum.plot.2.ratio <- spear.net.cornum.df %>%  
  mutate(data = factor(data, levels = c("UC","CD","CRC"))) %>%
  mutate(type = factor(type, levels = c("CTR", "Case"))) %>%
  ggplot(aes(x = data, y = posi.ratio, fill = type)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#a1d4e2", "#e62739")) +
  labs(x = "Disease", 
       y = "Network density") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_cartesian(ylim = c(8, 9)) +
  scale_y_continuous(breaks = c(8, 9))  +
  geom_text(aes(label = posi.num, y = posi.ratio -0.15), 
            position = position_dodge(0.9), size = 5,
            color = "white")

spear.net.cornum.plot.all <- ggarrange(spear.net.cornum.plot.2.ratio,
                                       spear.net.cornum.plot.1.ratio,
                                       heights=c(1/5, 4/5),ncol = 1, 
                                       nrow = 2, common.legend = TRUE,
                                       legend="right",align = "v") 


## negative correlation numbers ratio
spear.net.neganum.plot.1.ratio <- spear.net.cornum.df %>%  
  mutate(data = factor(data, levels = c("UC","CD","CRC"))) %>%
  mutate(type = factor(type, levels = c("CTR", "Case"))) %>%
  ggplot(aes(x = data, y = nega.ratio, fill = type)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  theme_bw() +
  scale_fill_manual(values = c("#a1d4e2", "#e62739")) +
  labs(x = "Disease", 
       y = "Case-to-control ratio of negative correlation numbers") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)) + 
  coord_cartesian(ylim = c(0, 4)) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4)) +
  geom_text(aes(label = nega.num, y = nega.ratio -0.15), 
            position = position_dodge(0.9), size = 5,
            color = "white")


spear.net.cornum.plot.combined <- plot_grid(spear.net.cornum.plot.all,
                                            spear.net.neganum.plot.1.ratio,
                                       widths = c(1/2, 1/2),ncol = 2,  
                                       legend="right", align = "v", 
                                       labels = c("F", "G")) 
pdf("../pdf/Network_ratioPlot.pdf",
    width = 10, height = 10)
spear.net.cornum.plot.combined
dev.off()
