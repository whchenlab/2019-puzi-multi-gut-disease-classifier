
## ------case-enriched species clustering ------------##

library(RColorBrewer)
library(tidyr)
library(pheatmap)
library(dplyr)
library(coin)

kw.spe.case <- kw.sum.df %>% filter(trend %in% "Case") %>% pull(features)
kw.spe.case <- unique(kw.spe.case)
data.df.dif.case <- data.df.t.nohdc[, kw.spe.case]


data.df.dif.case.matrix <- data.df.dif.case
for (i in 1:ncol(data.df.dif.case.matrix)) {
  y <- data.df.dif.case.matrix[metadata %>% filter(Group %in% "Control") %>% pull(Sample_ID), i]
  y.quantile <- quantile(y, probs = 0.95)
  data.df.dif.case.matrix[,i] <- ifelse(data.df.dif.case.matrix[,i] > y.quantile, 1,0)
}
rm(y, y.quantile)

## ---only use cases data to create species clustering

jaccard.dist <- vegdist(data.frame(t(data.df.dif.case.matrix[rownames(data.nocontrol.list$data.df.t.nohdc),]), 
                                   stringsAsFactors = F),
                        method = "jaccard", binary = T)
jacc.dend <- hclust(jaccard.dist, method='ward.D2')
jacc.dend.clustering <- cutree(jacc.dend, k=5)

hm.mat <- as.matrix(jaccard.dist)
hm.mat <- 1- hm.mat
hm.mat.heatmap <- pheatmap(hm.mat, cluster_rows = jacc.dend,
                   cluster_cols = jacc.dend,
                   cellwidth = 12, cellheight =12,
                   fontsize_col = 10,
                   fontsize_row = 10,
                   color = colorRampPalette(c("white", "navy", "firebrick3"))(300))

metadata.case <- metadata %>% filter(Sample_ID %in% rownames(data.df.t.nohdc))

features.red <- t(data.df.dif.case.matrix)[, metadata.case$Sample_ID]
features.red <- t(sapply(unique(jacc.dend.clustering), FUN=function(x){
  as.numeric(colSums(features.red[which(jacc.dend.clustering == x),]) >= 1)
}))
rownames(features.red) <- c('UC markers', 'mixed markers 1', 
                            'mixed markers 2', 'CD markers', 'CRC markers')

## ------positive fraction in each disease ------##
## ------每个项目的binary情况
cluster.to.group.sum <- as_tibble(t(features.red)) %>% 
  mutate(conf=metadata.case[["Group"]], project = metadata.case[["Project"]]) %>% 
  gather(key=key, value=group, -conf, -project) 

## -----每个疾病的positive fraction的情况--------
cluster.to.group.sum.2 <- cluster.to.group.sum %>%
  group_by(key, conf) %>%
  summarise(value = sum(group)/n())

## -------计算卡方block检验值
cluster.to.group.pvals <- list()

a <- cluster.to.group.sum
a <- data.frame(a)
a$conf <- as.factor(a$conf)
a$project <- as.factor(a$project)
a$group <- as.factor(a$group)

a.test <- subset(a, project %in% study)
for (keyvalue in unique(a.test$key)) {
  a.test.key <- subset(a.test, key %in% keyvalue)
  a.pval <- cmh_test(group ~ conf |project, a.test.key)
  cluster.to.group.pvals[[paste("CRC",keyvalue, sep = "_")]] <- pvalue(a.pval)
}

a.test <- subset(a, project %in% CD.study & !conf %in% "UC")
for (keyvalue in unique(a.test$key)) {
  a.test.key <- subset(a.test, key %in% keyvalue)
  a.pval <- cmh_test(group ~ conf |project, a.test.key)
  cluster.to.group.pvals[[paste("CD",keyvalue, sep = "_")]] <- pvalue(a.pval)
}

a.test <- subset(a, project %in% UC.study & !conf %in% "CD")
for (keyvalue in unique(a.test$key)) {
  a.test.key <- subset(a.test, key %in% keyvalue)
  a.pval <- cmh_test(group ~ conf |project, a.test.key)
  cluster.to.group.pvals[[paste("UC",keyvalue, sep = "_")]] <- pvalue(a.pval)
}

cluster.to.group.pvals.df <- cluster.to.group.pvals %>% reduce(rbind)
cluster.to.group.pvals.df <- data.frame(cluster.to.group.pvals.df)
cluster.to.group.pvals.df$type <- names(cluster.to.group.pvals)
cluster.to.group.pvals.df$disease <- c("CRC","CRC","CRC","CRC","CRC",
                                       "CD","CD","CD","CD","CD",
                                       "UC","UC","UC","UC","UC")
cluster.to.group.pvals.df$type <- c(unique(a.test$key), unique(a.test$key), unique(a.test$key))
names(cluster.to.group.pvals.df)[1] <- "pval"
cluster.to.group.pvals.df$pval_sig <- " "
cluster.to.group.pvals.df[which(cluster.to.group.pvals.df$pval < 0.05), "pval_sig"] <- "*"
cluster.to.group.pvals.df[which(cluster.to.group.pvals.df$pval < 0.01), "pval_sig"] <- "**"

cluster.to.group.sum.2.merge <- merge(cluster.to.group.sum.2,
                                      cluster.to.group.pvals.df,
                                      by.x = c("key","conf"),
                                      by.y = c("type","disease"), 
                                      all.y = T)
cluster.to.group.plot <- cluster.to.group.sum.2.merge %>%
  mutate(conf = factor(conf, levels = rev(c("UC","CD","CRC")))) %>%
  mutate(key = factor(key, levels = rev(c("UC markers",
                                          "CD markers",
                                          "CRC markers",
                                          "mixed markers 1",
                                          "mixed markers 2")))) %>%
  ggplot(aes(x = key, y = value, fill =  conf)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  geom_text(aes(label = pval_sig, y = value +0.03, group = conf), 
            position = position_dodge(0.9), size = 5) +
  coord_flip() + 
  scale_fill_manual(values = c("#8265CC","#74B347","#F2CC30")) +
  guides(fill = guide_legend(reverse=TRUE)) +
  ylim(0,0.8)+
  labs(x = "clusters",
       y = "positive fraction")

