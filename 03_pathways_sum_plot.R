## ----- pathways heatmap ------##
## ----- generalized fold change --------##
kw.path <- names(kw.path.list$four.gp)
kw.path <- kw.path[!kw.path %in% "Group"]


## summarise the sum of dagradation
degra.sum.df.plot <- kw.path.list$four.gp[,c(rownames(degradation.annot),"Group")]
degra.sum.df.plot$Sample_ID <- rownames(degra.sum.df.plot)
degra.sum.df.plot$Project <- metadata[rownames(degra.sum.df.plot),"Project"]

degra.sum.df.2 <- tibble(
  `Carbohydrate degradation` = 
    rowMeans(log10(degra.sum.df.plot[,c("Pathway_PWY_621", "Pathway_PWY_6317", "Pathway_PWY_6527",
                                        "Pathway_PWY_6737", "Pathway_PWY_6901", "Pathway_PWY66_422")] + 1e-09), na.rm=TRUE),
  `Amino Acid Degradation` = 
    rowMeans(log10(degra.sum.df.plot[,c("Pathway_P162_PWY", "Pathway_HISDEG_PWY")] + 1e-09), na.rm=TRUE),
  `Carboxylate Degradation` = 
    rowMeans(log10(degra.sum.df.plot[,c("Pathway_P441_PWY", "Pathway_P161_PWY")] + 1e-09), na.rm=TRUE),
  `Nucleoside and Nucleotide Degradation` = 
    rowMeans(log10(degra.sum.df.plot[,c("Pathway_PWY0_1296", "Pathway_PWY0_1297")] + 1e-09), na.rm=TRUE),
  `Secondary Metabolite Degradation` = 
    log10(degra.sum.df.plot[,c("Pathway_PWY_4702")] + 1e-09),
  
  Group=degra.sum.df.plot$Group,
  Study=degra.sum.df.plot$Project,
  Sample_ID = degra.sum.df.plot$Sample_ID)

degra.sum.df.2.melt <- melt(degra.sum.df.2, 
                            id.vars = c("Group","Study","Sample_ID"))

degra.sum.df.2.melt.list <- list()
cd.project <- c("SRP057027", "PRJNA389280", "PRJNA400072")
uc.project <- c("PRJEB1220", "PRJNA389280", "PRJNA400072")
crc.project <- metadata %>% filter(Group %in% "CRC") %>% pull(Project) %>% unique

degra.sum.df.2.melt.list[["CD"]] <- degra.sum.df.2.melt %>% 
  filter(Study %in% cd.project) %>%
  filter(!Group %in% "UC")
degra.sum.df.2.melt.list[["CD"]]$Group <- factor(degra.sum.df.2.melt.list[["CD"]]$Group)
degra.sum.df.2.melt.list[["CD"]]$Study <- factor(degra.sum.df.2.melt.list[["CD"]]$Study)

degra.sum.df.2.melt.list[["UC"]] <- degra.sum.df.2.melt %>% 
  filter(Study %in% uc.project) %>%
  filter(!Group %in% "CD")
degra.sum.df.2.melt.list[["UC"]]$Group <- factor(degra.sum.df.2.melt.list[["UC"]]$Group)
degra.sum.df.2.melt.list[["UC"]]$Study <- factor(degra.sum.df.2.melt.list[["UC"]]$Study)

degra.sum.df.2.melt.list[["CRC"]] <- degra.sum.df.2.melt %>% 
  filter(Study %in% crc.project)
degra.sum.df.2.melt.list[["CRC"]]$Group <- factor(degra.sum.df.2.melt.list[["CRC"]]$Group)
degra.sum.df.2.melt.list[["CRC"]]$Study <- factor(degra.sum.df.2.melt.list[["CRC"]]$Study)

degra.sum.pval <- list()

for (i in names(degra.sum.df.2.melt.list)) {
  for (j in unique(degra.sum.df.2.melt.list[[i]]$variable)) {
    a <- degra.sum.df.2.melt.list[[i]] %>%
      filter(variable %in% j)
    
    a.wilcox <- wilcox_test(value ~ Group | Study, a)
    degra.sum.pval[[i]][[j]] <- pvalue(a.wilcox)
  }
}

## despite of these degradation pathways
## summarise the sum of biosynthesis
biosyn.sum.df.plot <- kw.path.list$four.gp[, !names(kw.path.list$four.gp) %in% rownames(degradation.annot)]
biosyn.sum.df.plot$Sample_ID <- rownames(biosyn.sum.df.plot)
biosyn.sum.df.plot$Project <- metadata[rownames(biosyn.sum.df.plot),"Project"]

biosyn.sum.df.2 <- tibble(
  `Amino Acid Biosynthesis` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Amino Acid Biosynthesis")$Abb] + 1e-09), na.rm=TRUE),
  `Aromatic Compound Biosynthesis` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Aromatic Compound Biosynthesis")$Abb] + 1e-09), na.rm=TRUE),
  `Carbohydrate Biosynthesis` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Carbohydrate Biosynthesis")$Abb] + 1e-09), na.rm=TRUE),
  `Cofactor, Carrier, and Vitamin Biosynthesis` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Cofactor, Carrier, and Vitamin Biosynthesis")$Abb] + 1e-09), na.rm=TRUE),
  `Fatty Acid and Lipid Biosynthesis` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Fatty Acid and Lipid Biosynthesis")$Abb] + 1e-09), na.rm=TRUE),
  `Nucleoside and Nucleotide Biosynthesis` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Nucleoside and Nucleotide Biosynthesis")$Abb] + 1e-09), na.rm=TRUE),
  `Secondary Metabolite Biosynthesis` = 
    log10(biosyn.sum.df.plot[, "Pathway_NONMEVIPP_PWY"] + 1e-09), 
  
  Group=biosyn.sum.df.plot$Group,
  Study=biosyn.sum.df.plot$Project,
  Sample_ID = biosyn.sum.df.plot$Sample_ID)

biosyn.sum.df.2.melt <- melt(biosyn.sum.df.2, 
                             id.vars = c("Group","Study","Sample_ID"))

biosyn.sum.df.2.melt.list <- list()

biosyn.sum.df.2.melt.list[["CD"]] <- biosyn.sum.df.2.melt %>% 
  filter(Study %in% cd.project) %>%
  filter(!Group %in% "UC")
biosyn.sum.df.2.melt.list[["CD"]]$Group <- factor(biosyn.sum.df.2.melt.list[["CD"]]$Group)
biosyn.sum.df.2.melt.list[["CD"]]$Study <- factor(biosyn.sum.df.2.melt.list[["CD"]]$Study)

biosyn.sum.df.2.melt.list[["UC"]] <- biosyn.sum.df.2.melt %>% 
  filter(Study %in% uc.project) %>%
  filter(!Group %in% "CD")
biosyn.sum.df.2.melt.list[["UC"]]$Group <- factor(biosyn.sum.df.2.melt.list[["UC"]]$Group)
biosyn.sum.df.2.melt.list[["UC"]]$Study <- factor(biosyn.sum.df.2.melt.list[["UC"]]$Study)

biosyn.sum.df.2.melt.list[["CRC"]] <- biosyn.sum.df.2.melt %>% 
  filter(Study %in% crc.project)
biosyn.sum.df.2.melt.list[["CRC"]]$Group <- factor(biosyn.sum.df.2.melt.list[["CRC"]]$Group)
biosyn.sum.df.2.melt.list[["CRC"]]$Study <- factor(biosyn.sum.df.2.melt.list[["CRC"]]$Study)

biosyn.sum.pval <- list()

for (i in names(biosyn.sum.df.2.melt.list)) {
  for (j in unique(biosyn.sum.df.2.melt.list[[i]]$variable)) {
    a <- biosyn.sum.df.2.melt.list[[i]] %>%
      filter(variable %in% j)
    
    a.wilcox <- wilcox_test(value ~ Group | Study, a)
    biosyn.sum.pval[[i]][[j]] <- pvalue(a.wilcox)
  }
}

## despite of these degradation & biosynthesis pathways
## summarise the sum of distinct pathways

otherpath.sum.df.2 <- tibble(
  `Fermentation` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Fermentation")$Abb] + 1e-09), na.rm=TRUE),
  `Glycolysis` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Glycolysis")$Abb] + 1e-09), na.rm=TRUE),
  `Generation of Precursor Metabolites and Energy` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Generation of Precursor Metabolites and Energy")$Abb] + 1e-09), na.rm=TRUE),
  `Inorganic Nutrient Metabolism` = 
    rowMeans(log10(biosyn.sum.df.plot[,subset(kw.pathways.map, Superclass %in% "Inorganic Nutrient Metabolism")$Abb] + 1e-09), na.rm=TRUE),
  `Aminoacyl-tRNA Charging` = 
    log10(biosyn.sum.df.plot[,"Pathway_TRNA_CHARGING_PWY"] + 1e-09),
  
  Group=biosyn.sum.df.plot$Group,
  Study=biosyn.sum.df.plot$Project,
  Sample_ID = biosyn.sum.df.plot$Sample_ID)

otherpath.sum.df.2.melt <- melt(otherpath.sum.df.2, 
                             id.vars = c("Group","Study","Sample_ID"))

otherpath.sum.df.2.melt.list <- list()

otherpath.sum.df.2.melt.list[["CD"]] <- otherpath.sum.df.2.melt %>% 
  filter(Study %in% cd.project) %>%
  filter(!Group %in% "UC")
otherpath.sum.df.2.melt.list[["CD"]]$Group <- factor(otherpath.sum.df.2.melt.list[["CD"]]$Group)
otherpath.sum.df.2.melt.list[["CD"]]$Study <- factor(otherpath.sum.df.2.melt.list[["CD"]]$Study)

otherpath.sum.df.2.melt.list[["UC"]] <- otherpath.sum.df.2.melt %>% 
  filter(Study %in% uc.project) %>%
  filter(!Group %in% "CD")
otherpath.sum.df.2.melt.list[["UC"]]$Group <- factor(otherpath.sum.df.2.melt.list[["UC"]]$Group)
otherpath.sum.df.2.melt.list[["UC"]]$Study <- factor(otherpath.sum.df.2.melt.list[["UC"]]$Study)

otherpath.sum.df.2.melt.list[["CRC"]] <- otherpath.sum.df.2.melt %>% 
  filter(Study %in% crc.project)
otherpath.sum.df.2.melt.list[["CRC"]]$Group <- factor(otherpath.sum.df.2.melt.list[["CRC"]]$Group)
otherpath.sum.df.2.melt.list[["CRC"]]$Study <- factor(otherpath.sum.df.2.melt.list[["CRC"]]$Study)

otherpath.sum.pval <- list()

for (i in names(otherpath.sum.df.2.melt.list)) {
  for (j in unique(otherpath.sum.df.2.melt.list[[i]]$variable)) {
    a <- otherpath.sum.df.2.melt.list[[i]] %>%
      filter(variable %in% j)
    
    a.wilcox <- wilcox_test(value ~ Group | Study, a)
    otherpath.sum.pval[[i]][[j]] <- pvalue(a.wilcox)
  }
}


pathways.sum.pval <- cbind(degra.sum.pval %>% reduce(rbind), 
                           biosyn.sum.pval %>% reduce(rbind),
                           otherpath.sum.pval %>% reduce(rbind))
rownames(pathways.sum.pval) <- c("CD", "UC", "CRC")
write.table(pathways.sum.pval, 
            file = "../SuppleData/Pathways_sum_Pvalue.txt",
            row.names = T, col.names = T, sep = "\t", quote = F)

a <- merge(degra.sum.df.2, biosyn.sum.df.2, 
           by = c("Group", "Study", "Sample_ID"), 
           all = T)
a.2 <- merge(a, otherpath.sum.df.2, 
             by = c("Group", "Study", "Sample_ID"), 
             all = T)
rownames(a.2) <- a.2$Sample_ID

a.3 <- a.2[, -c(1:3)]
a.3$Group <- metadata[rownames(a.3), "Group"]

write.table(a.3, file = "../SuppleData/Pathways_sum_df.txt",
            col.names = T, row.names = T,
            sep = "\t", quote = F)

a.3.df <- a.2  %>% melt(id.vars = c("Sample_ID", "Group", "Study")) 
a.3.df.cd <- a.3.df %>% filter(Study %in% cd.project) %>%
  filter(!Group %in% "UC") %>%
  group_by(variable, Group) 

a.3.df.uc <- a.3.df %>% filter(Study %in% uc.project) %>%
  filter(!Group %in% "CD") %>%
  group_by(variable, Group)

fc <- matrix(NA, nrow=length(unique(a.3.df$variable)), ncol=3, 
             dimnames=list(unique(a.3.df$variable), c("UC","CD", "CRC")))

for (pathways in unique(a.3.df.uc$variable)) { 
    
    metadata.pro.ctr <- subset(a.3.df.uc,  Group %in% "CTR" & variable %in% pathways) 
    metadata.pro.case <- subset(a.3.df.uc,   variable %in% pathways & !Group %in% "CTR")
    q.p <- quantile(metadata.pro.ctr$value, probs=seq(.1, .9, .1), na.rm=TRUE)
    q.n <- quantile(metadata.pro.case$value, probs=seq(.1, .9, .1), na.rm=TRUE)
   fc[pathways, "UC"] <- sum(q.n - q.p)/length(q.p)

} 

fc.cd <- list()
for (pathways in unique(a.3.df.cd$variable)) { 
    
    metadata.pro.ctr <- subset(a.3.df.cd,  Group %in% "CTR" & variable %in% pathways) 
    metadata.pro.case <- subset(a.3.df.cd,  variable %in% pathways & !Group %in% "CTR")
    q.p <- quantile(metadata.pro.ctr$value, probs=seq(.1, .9, .1), na.rm=TRUE)
    q.n <- quantile(metadata.pro.case$value, probs=seq(.1, .9, .1), na.rm=TRUE)
    fc[pathways, "CD"] <- sum(q.n - q.p)/length(q.p) 
}

a.3.df.crc <- a.3.df %>% filter(Study %in% crc.project)  %>%
  group_by(variable, Group)

fc.crc <- list()
for (pathways in unique(a.3.df.crc$variable)) {
    
    metadata.pro.ctr <- subset(a.3.df.crc,  Group %in% "CTR" & variable %in% pathways) 
    metadata.pro.case <- subset(a.3.df.crc,  variable %in% pathways & !Group %in% "CTR")
    q.p <- quantile(metadata.pro.ctr$value, probs=seq(.1, .9, .1), na.rm=TRUE)
    q.n <- quantile(metadata.pro.case$value, probs=seq(.1, .9, .1), na.rm=TRUE)
    fc[pathways, "CRC"] <- sum(q.n - q.p)/length(q.p) 
}

fc.df <- data.frame(fc, stringsAsFactors = F)
fc.df$Modules <- rownames(fc.df)
names(fc.df) <- c("UC","CD", "CRC",
                  "Modules")

fc.df.merge <- fc.df %>% melt(id.vars = "Modules") 

fc.df.merge.2 <- merge(fc.df.merge, 
                      pathways.sum.pval.df, 
                         by.x = c("Modules", "variable"),
                         by.y = c( "variable", "Dataset"), 
                         all = T)
fc.df.merge.2[which(fc.df.merge.2$signif %in% "Not signif"), "signif"] <- " "

col.lodo.heatmap <- c( "#4F86C6",  "white", "firebrick3" )
pathways.sum.fc.plot <- fc.df.merge.2 %>%  
  mutate(Modules = factor(Modules, levels = c(names(degra.sum.df.2)[1:5], 
                                              names(biosyn.sum.df.2)[1:7],
                                              names(otherpath.sum.df.2)[1:5]))) %>%
  mutate(variable = factor(variable, levels = rev(c("UC", "CD", "CRC")))) %>%
  ggplot(aes(x = Modules, y = variable, fill = value.x)) +
  geom_tile()  + theme_bw()  +
  geom_text(aes(label= signif), col='black', size=6) +
  scale_fill_gradientn(colours = col.lodo.heatmap, limits=c(-1.2, 1.2), 
                       guide = "colourbar", breaks = c(-1.2, -0.6, 0, 0.6, 1.2)) + 
  theme(axis.line=element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid=element_blank(), 
        panel.border=element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        axis.text.x = element_text(size = 11,
                                   angle = 90,
                                   hjust = .5, vjust = .5),
        axis.text.y = element_text(size = 12)) + 
  ylab('') + xlab('')

pdf("../pdf/pathways_sum_fc.pdf",
    width = 12, height = 6)
pathways.sum.fc.plot
dev.off()
