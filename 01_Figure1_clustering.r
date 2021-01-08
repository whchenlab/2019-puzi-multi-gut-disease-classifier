load("../RData/combined_list_rf.Rdata")
load("../RData/species_list_rf.Rdata")
dim(all.spe.list$four.gp)

library(vegan)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyr)
library(pheatmap)
library(gdata)

kw.spe.df <- read.delim("../SuppleData/Dif_species_annotation.txt",
                        header = T, sep = "\t", as.is = T)
kw.spe.case <- kw.spe.df %>% 
  filter(Trend %in% "Case-enriched") %>%
  pull(Species)
kw.spe.df.case <- kw.spe.df %>% 
  filter(Trend %in% "Case-enriched")

data.df.dif.case <- all.spe.list$four.gp[, kw.spe.case]


data.df.dif.case.matrix <- data.df.dif.case
for (i in 1:ncol(data.df.dif.case.matrix)) {
  y <- data.df.dif.case.matrix[metadata %>% filter(Group %in% "CTR") %>% pull(Sample_ID), i]
  y.quantile <- quantile(y, probs = 0.95)
  data.df.dif.case.matrix[,i] <- ifelse(data.df.dif.case.matrix[,i] > y.quantile, 1,0)
}
rm(y, y.quantile)

#metadata.case <- metadata %>% filter(Sample_ID %in% rownames(data.df.t.nohdc))
data.df.noctr <- all.spe.list$four.gp %>% filter(!Group %in% "CTR") 

jaccard.dist <- vegdist(data.frame(t(data.df.dif.case.matrix[rownames(data.df.noctr),]), 
                                   stringsAsFactors = F),
                        method = "jaccard", binary = T)
jacc.dend <- hclust(jaccard.dist, method='ward.D2')
jacc.dend.clustering <- cutree(jacc.dend, k=5)

jacc.dend.clustering <- kw.spe.df.case$exposure
names(jacc.dend.clustering) <- kw.spe.df.case$Species
##
## ------------clusters similarity -------------##
hm.mat <- as.matrix(jaccard.dist)
hm.mat <- 1- hm.mat

hm.mat.2 <- hm.mat
diag(hm.mat.2) <- NA
upperTriangle(hm.mat.2) <- NA

hm.mat.melt <- as_tibble(hm.mat.2) %>% 
  mutate(species=rownames(hm.mat.2)) %>% 
  mutate(cluster=jacc.dend.clustering[rownames(hm.mat.2)]) %>% 
  melt(id.vars = c("species", "cluster")) 

names(hm.mat.melt)[3:4] <- c("species.2", "jaccard")
comparisons <- list(c("bg", "Cluster UC"),
                    c("bg", "Cluster CD"),
                    c("bg", "Cluster CRC"),
                    c("bg", "Cluster Common") )

g.inset <-  hm.mat.melt %>% 
    mutate(cluster.2=jacc.dend.clustering[species.2]) %>% 
    mutate(background=
               case_when(cluster!=cluster.2~'bg',
                         TRUE~paste0('Cluster ', as.character(cluster)))) %>% 
    filter(complete.cases(.)) %>% 
    mutate(background = factor(background, levels = c("bg", "Cluster UC",
                                                      "Cluster CD", "Cluster CRC", "Cluster Common"))) %>%
    ggplot(aes(x=background, y=jaccard, color  = background)) +
    geom_boxplot(outlier.shape=NA, ) +   
    theme_classic() + 
    theme(axis.text = element_text(size = 12)) +
    xlab('') + ylab('Jaccard Similarity') +
    scale_color_manual(values = c("#74828F","#80B1D3", "#74c239","firebrick3" ,"#996699")) +
    theme(axis.text.x = element_text(vjust= .5, hjust = .5, angle = 90, size = 11),
          legend.position = "none") +
    stat_compare_means(comparisons = comparisons , label.y = c(0.4, 0.5, 0.6, 0.7 ),
                       label = "p.signif") + 
    stat_compare_means(label.y = 0.8,
                       label = "p.signif", 
                       label.x.npc = "middle" )  + 
    geom_segment(aes(x = 1, y = 0.8, xend = 5, yend = 0.8), size = .4, color = "black")


features.red <- t(data.df.dif.case.matrix)[, metadata$Sample_ID]
features.red <- t(sapply(unique(jacc.dend.clustering), FUN=function(x){
  as.numeric(colSums(features.red[which(jacc.dend.clustering == x),]) >= 1)
}))

## ------positive fraction in each disease ------## 
cluster.to.group.sum <-as_tibble(t(features.red)) %>% 
  mutate(conf = metadata[["Group"]], 
         project = metadata[["Project"]], 
         Sample_ID = metadata[["Sample_ID"]]) %>%
  melt(id.vars = c("conf", "project", "Sample_ID" ))

names(cluster.to.group.sum)[4:5] <- c("key", "group")
cluster.to.group.sum.2 <- cluster.to.group.sum %>%
  group_by(key, conf) %>%
  summarise(value = sum(group)/n())

cluster.to.group.pvals <- list()

a <- cluster.to.group.sum
a <- data.frame(a)
a$conf <- as.factor(a$conf)
a$project <- as.factor(a$project)
a$group <- as.factor(a$group)

study <- metadata %>%
  filter(Group %in% "CRC") %>% 
  pull(Project) %>% unique

a.test <- subset(a, project %in% study)
for (keyvalue in unique(a.test$key)) {
  a.test.key <- subset(a.test, key %in% keyvalue)
  a.pval <- cmh_test(group ~ conf |project, a.test.key)
  cluster.to.group.pvals[[paste("CRC",keyvalue, sep = "_")]] <- pvalue(a.pval)
}
 

a.test <- subset(a, project %in% cd.project & !conf %in% "UC")
for (keyvalue in unique(a.test$key)) {
  a.test.key <- subset(a.test, key %in% keyvalue)
  a.pval <- cmh_test(group ~ conf |project, a.test.key)
  cluster.to.group.pvals[[paste("CD",keyvalue, sep = "_")]] <- pvalue(a.pval)
}

a.test <- subset(a, project %in% uc.project & !conf %in% "CD")
for (keyvalue in unique(a.test$key)) {
  a.test.key <- subset(a.test, key %in% keyvalue)
  a.pval <- cmh_test(group ~ conf |project, a.test.key)
  cluster.to.group.pvals[[paste("UC",keyvalue, sep = "_")]] <- pvalue(a.pval)
}

cluster.to.group.pvals.df <- cluster.to.group.pvals %>% reduce(rbind)
cluster.to.group.pvals.df <- data.frame(cluster.to.group.pvals.df)
cluster.to.group.pvals.df$type <- names(cluster.to.group.pvals)
cluster.to.group.pvals.df$disease <- c("CRC","CRC","CRC","CRC",
                                       "CD","CD","CD","CD",
                                       "UC","UC","UC","UC")
cluster.to.group.pvals.df$type <- rep(c("UC","CD","CRC","Common"),3)

names(cluster.to.group.pvals.df)[1] <- "pval"
cluster.to.group.pvals.df$pval_sig <- " "
cluster.to.group.pvals.df[which(cluster.to.group.pvals.df$pval < 0.05), "pval_sig"] <- "*"
cluster.to.group.pvals.df[which(cluster.to.group.pvals.df$pval < 0.01), "pval_sig"] <- "**"

cluster.to.group.sum.2.merge <- merge(cluster.to.group.sum.2,
                                      cluster.to.group.pvals.df,
                                      by.x = c("key","conf"),
                                      by.y = c("type","disease"), 
                                      all.y = T)
posi.frac.plot <- cluster.to.group.sum.2.merge %>%
  mutate(key = paste("Cluster", key)) %>%
  mutate(conf = factor(conf, levels =  c("UC","CD","CRC"))) %>%
  mutate(key = factor(key, levels =  c("Cluster UC","Cluster CD","Cluster CRC","Cluster Common"))) %>%
  ggplot(aes(x = key, y = value, fill =  conf)) +
  geom_bar(stat = "identity",
           position = "dodge", 
           alpha = 0.8,
           color = "white" ) +
  geom_text(aes(label = pval_sig, y = value +0.03, group = conf), 
            position = position_dodge(0.9), size = 5) + 
  scale_fill_manual(values = c("#D9D9D9","#6BAED6","firebrick3")) +
  guides(fill = guide_legend(reverse=FALSE)) +
  ylim(0, 1)+
  labs(x = "clusters",
       y = "positive fraction") +
  theme_classic() +
  theme(legend.position = "top", 
        axis.text.x = element_text(vjust= .5, hjust = .5, angle = 90, size = 11),
        axis.text.y = element_text(size = 12)) 

cluster.plot <- plot_grid(g.inset, 
          posi.frac.plot, 
          rel_widths = c(2/5, 3/5),
          axis = "b",
          align = "h")

pdf("../pdf/figure1_cluster.pdf", 
    width = 8,
    height = 5)
cluster.plot
dev.off()
