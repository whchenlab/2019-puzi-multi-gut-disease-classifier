library(Hmisc)

cor.func <- function(abund.data, metadata = metadata, 
                     type = "CD", project = "PRJNA389280"){
  ##data, first column is humanper 
  cut_sig <- function(p) {
    out <- cut(p, breaks = c(0, 0.01,0.05,1), include.lowest = T, 
               labels = c("**", "*", ""))
    return(out)
  }
  ##-----require data, first column is humanper----
  metadata.subset <- metadata %>% filter(Project %in% project) %>%
    filter(Group %in% c("CTR", type))
  rownames(metadata.subset) <- metadata.subset$Sample_ID
  
  data <- data.frame(HDC = metadata.subset$humanper, 
                  abund.data[rownames(metadata.subset), !names(abund.data) %in% "Group"],
                  stringsAsFactors = F)

  data.cor <- rcorr(as.matrix(data), type = "spearman")
  data.coef <- data.cor$r[1,-1]
  data.P <- data.cor$P[1,-1]
  data.P.nona <- na.omit(data.P)
  data.coef.nona <- na.omit(data.coef)
  data.sig <- data.frame(apply(data.frame(data.P.nona), 2, cut_sig))
  rownames(data.sig) <- names(data.P.nona)
  ##samples in week 0
  ##coef results bewtween humanDNA% and taxa abundance
  data.cor.sum <- cbind(data.coef.nona, data.P.nona, data.sig)
  names(data.cor.sum) <- c("coef","P_value","SigLevel")
  data.cor.sum[,"Species"] <- rownames(data.cor.sum)
  data.cor.sum[,"Type"] <- type
  data.cor.sum[,"Project"] <- project
  sig.taxa <- subset(data.cor.sum,! SigLevel %in% "")
  return(list(data.cor.sum = data.cor.sum,
              sig.taxa = sig.taxa))
}
## Species level
## CRC
crc.HDC.cor.list <- list()

for (pro in crc.project) {
  crc.HDC.cor.list[[pro]] <- cor.func(abund.data = crc.spe.list$all, 
                                 metadata = metadata, 
                                 type = "CRC", project = pro)
}
## CD 
cd.HDC.cor.list <- list()
for(pro in cd.project){
  cd.HDC.cor.list[[pro]] <- cor.func(abund.data = cd.spe.list$all, 
                                     metadata = metadata, 
                                     type = "CD", project = pro)
}
## UC
uc.HDC.cor.list <- list()
for(pro in uc.project){
  uc.HDC.cor.list[[pro]] <- cor.func(abund.data = uc.spe.list$all, 
                                     metadata = metadata, 
                                     type = "UC", project = pro)
}
## combined
cor.data.list <- lapply(c(crc.HDC.cor.list,
                          cd.HDC.cor.list,
                          uc.HDC.cor.list), 
                        function(data){y <- data$data.cor.sum})
cor.data <- reduce(cor.data.list, rbind)

## ----filter results with P-value < 0.05
cor.data.sig <- subset(cor.data, P_value<0.05)
cor.data.posi <- subset(cor.data.sig, coef>0)
cor.data.nega <- subset(cor.data.sig, coef<0)
delete.spe <- intersect(cor.data.nega$Species, cor.data.posi$Species)
cor.data.sig <- subset(cor.data.sig, !Species %in% delete.spe)

cor.data.sig.sum <- cor.data.sig %>% 
  group_by(Species, Type) %>% 
  summarise(count= n())
cor.data.sig.sum.2 <- dcast(cor.data.sig.sum, 
                            Species ~Type)

sig.spe <- subset(cor.data.sig.sum.2, 
                  CRC>1 | UC>1 | CD>1)$Species
sig.spe.df <- subset(cor.data.sig.sum.2, 
                     CRC>1 | UC>1 | CD>1)
sig.sum.df <- melt(sig.spe.df, 
                   id.vars = "Species")
sig.sum.df <- sig.sum.df %>% filter(!is.na(value)) %>%
  filter(value > 1)


sig.spe.df.merge <- merge(kw.spe.df, 
                          sig.sum.df, 
                          by.x = c("feature", "exposure"), 
                          by.y = c("Species", "variable"),
                          all.x = T)
names(sig.spe.df.merge)[9] <- "count_sig"
sig.spe.df.merge$HDC_related <- "Yes"
sig.spe.df.merge[is.na(sig.spe.df.merge$count_sig), "HDC_related"] <- "NO"

sig.spe.df.merge <- sig.spe.df.merge %>% arrange(Phylum)

sig.spe.df.bottom <- sig.spe.df.merge %>% 
  mutate(final.spe = factor(final.spe, levels = unique(kw.spe.df$final.spe))) %>%
  mutate(exposure = factor(exposure, levels = c("CRC","CD","UC"))) %>%
  mutate(HDC_related = factor(HDC_related, levels = c("Yes","NO"))) %>%
  ggplot(aes(x = final.spe, y = exposure)) +
  geom_tile(aes(fill = HDC_related),
            color="#D9D9D9") + 
  scale_fill_manual(values = c("firebrick3","#6BAED6")) +
  theme_classic()+
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(),  
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text = element_text(size = 8, angle = 90, hjust = 1),
        #axis.text = element_text(vjust = .5, hjust = .5),
        axis.title = element_text(size = 13),
        #panel.border=element_blank(),
        legend.position = "bottom") 
sig.spe.df.top <- sig.spe.df.merge %>%
  mutate(final.spe = factor(final.spe, levels = unique(kw.spe.df$final.spe))) %>%
  ggplot(aes(x = final.spe, y = 1)) + geom_tile(aes(fill = Phylum)) +
  theme_classic()+
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(),  
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.border=element_blank(),
        legend.position = "top")+
  scale_fill_manual(values =  c("#B882BC","#E8C3DB", "#C5C5C5",
                                "#A8BC68","#8EC8B6","#4A9869"))

sig.spe.df.plot <-plot_grid(sig.spe.df.top, sig.spe.df.bottom, 
                             ncol = 1, align = 'v', rel_heights = c(0.3,1))

pdf("../pdf/HDC_related_species.pdf", 
    height = 5, width = 8)
sig.spe.df.plot
dev.off()

## Pathways level
## CRC
crc.HDC.cor.path.list <- list()

for (pro in crc.project) {
  crc.HDC.cor.path.list[[pro]] <- cor.func(abund.data = crc.path.list$all, 
                                      metadata = metadata, 
                                      type = "CRC", project = pro)
}
## CD 
cd.HDC.cor.path.list <- list()
for(pro in cd.project){
  cd.HDC.cor.path.list[[pro]] <- cor.func(abund.data = cd.path.list$all, 
                                     metadata = metadata, 
                                     type = "CD", project = pro)
}
## UC
uc.HDC.cor.path.list <- list()
for(pro in uc.project){
  uc.HDC.cor.path.list[[pro]] <- cor.func(abund.data = uc.path.list$all, 
                                     metadata = metadata, 
                                     type = "UC", project = pro)
}

## 
cor.path.list <- lapply(c(crc.HDC.cor.path.list,
                          cd.HDC.cor.path.list,
                          uc.HDC.cor.path.list), 
                        function(data){y <- data$data.cor.sum})
cor.path <- reduce(cor.path.list, rbind)

## ----filter results with P-value < 0.05
cor.path.sig <- subset(cor.path, P_value<0.05)
cor.path.posi <- subset(cor.path.sig, coef>0)
cor.path.nega <- subset(cor.path.sig, coef<0)
delete.spe <- intersect(cor.path.nega$Species, cor.path.posi$Species)
cor.path.sig <- subset(cor.path.sig, !Species %in% delete.spe)

cor.path.sig.sum <- cor.path.sig %>% 
  group_by(Species, Type) %>% 
  summarise(count= n())
cor.path.sig.sum.2 <- dcast(cor.path.sig.sum, 
                            Species ~Type)

sig.path <- subset(cor.path.sig.sum.2, 
                  CRC>1 | UC>1 | CD>1)$Species
sig.path.df <- subset(cor.path.sig.sum.2, 
                     CRC>1 | UC>1 | CD>1)

sig.sum.path.df <- melt(sig.path.df, 
                   id.vars = "Species")
sig.sum.path.df <- sig.sum.path.df %>% filter(!is.na(value)) %>%
  filter(value > 1)


sig.path.df.merge <- merge(kw.path.df.merge, 
                           sig.sum.path.df, 
                          by.x = c("Abb", "exposure"), 
                          by.y = c("Species", "variable"),
                          all.x = T)
names(sig.path.df.merge)[8] <- "count_sig"
sig.path.df.merge$HDC_related <- "Yes"
sig.path.df.merge[is.na(sig.path.df.merge$count_sig), "HDC_related"] <- "NO"

sig.path.df.bottom <- sig.path.df.merge  %>%
  mutate(exposure = factor(exposure, levels = c("CRC","CD","UC"))) %>%
  ggplot(aes(x = feature, y = exposure)) +
  geom_tile(aes(fill = HDC_related), 
            color="#D9D9D9") + 
  scale_fill_manual(values = c("firebrick3","#6BAED6")) +
  theme_classic()+
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(),  
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text = element_text(size = 8, angle = 90, hjust = 1),
        #axis.text = element_text(vjust = .5, hjust = .5),
        axis.title = element_text(size = 13),
        #panel.border=element_blank(),
        legend.position = "bottom") 

pdf("../pdf/HDC_related_pathways.pdf", 
    height = 7, width = 8)
sig.path.df.bottom
dev.off()

write.table(sig.path.df.merge, file = "../SuppleData/HDC_correlated_DifPathways.txt", 
            col.names = T, row.names = F, sep = "\t", quote = F)