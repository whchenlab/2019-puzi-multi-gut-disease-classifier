## based on R 4.0.0

setwd("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/abund_data")

.libPaths("/mnt/raid1/puzi/Rlib/bork")

library(MMUPHin, lib.loc = "/mnt/raid1/puzi/Rlib/bork")
library(rlang)
library(dplyr)
library(ggplot2)
library(purrr)
library(cowplot)

## from MMUPHin
rename_maaslin <- function(old_names, prefix) {
  if(is.null(old_names) | length(old_names) == 0) return(NULL)
  new_names <- paste0(prefix, seq_along(old_names))
  names(new_names) <- old_names
  return(new_names)
}

create_table_maaslin <- function(features, exposure, lvl_exposure) {
  if(is.null(lvl_exposure))
    values_exposure <- exposure
  else
    values_exposure <- lvl_exposure[-1]
  names(features) <- NULL
  table_maaslin <- expand.grid(features, exposure, values_exposure, 
                               stringsAsFactors = FALSE)
  names(table_maaslin) <- c("feature", "metadata", "value")
  return(table_maaslin)
}
catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), 
                  error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 

## for adjusting covariates and identifying features with differential abundance using Maaslin2
maaslin2_refunc <- function(metadata = metadata, 
                            pro = pro, case = "CRC", 
                            exposure = "Group",
                            covariate = "Age",
                            abund.s.list.filter.df = data,
                            output = getwd(),
                            normalization = "TSS",
                            transform = "AST",
                            analysis_method = "LM"){
  a <- metadata %>% filter(Project %in% pro) %>% 
    filter(Group %in% c("CTR", case))
  a$Group <- factor(a$Group, levels = c("CTR", case))
  rownames(a) <- a$Sample_ID
  abd <- abund.s.list.filter.df[, rownames(a)]
  
  data_rename <- a[, c(exposure, covariate), 
                   drop = FALSE]
  feature_abd_rename <- abd
  
  features_rename <- rename_maaslin(rownames(feature_abd_rename), prefix = "T")
  samples_rename <- rename_maaslin(colnames(feature_abd_rename), prefix = "S")
  exposure_rename <- rename_maaslin(exposure, prefix = "E")
  covariates_rename <- rename_maaslin(covariate, prefix = "X")
  #covariates_random_rename <- rename_maaslin(covariates_random, prefix = "RX")
  dimnames(feature_abd_rename) <- list(features_rename, samples_rename)
  
  dimnames(data_rename) <- list(samples_rename,
                                c(exposure_rename,
                                  covariates_rename))
  ind_features <- apply(feature_abd_rename> 0, 1, any)
  
  a.maaslin <- Maaslin2::Maaslin2(
    input_data = feature_abd_rename[ind_features, , drop = TRUE],
    input_metadata = data_rename,
    output = output,
    min_abundance = 0,
    min_prevalence = 0,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    max_significance = 1,
    fixed_effects = c(exposure_rename, covariates_rename),
    standardize = FALSE,
    plot_heatmap = FALSE,
    plot_scatter = FALSE)
  logtomaaslin <- catchToList(a.maaslin$results)
  
  res_rename <- logtomaaslin$value
  
  # Read Maaslin results
  lvl_exposure <- NULL
  if(is.factor(a[[exposure]]))
    lvl_exposure <- levels(a[[exposure]])
  
  table_maaslin <- dplyr::left_join(
    data.frame(feature = names(features_rename),
               feature_rename = features_rename,
               stringsAsFactors = FALSE),
    create_table_maaslin(features_rename,
                         exposure_rename,
                         lvl_exposure),
    by = c("feature_rename" = "feature"))
  
  res <- dplyr::left_join(table_maaslin, res_rename,
                          by = c("feature_rename" = "feature",
                                 "metadata",
                                 "value"))
  res <- dplyr::select(res, -feature_rename, -name)
  
  res$metadata <- exposure
  if(all(res$value == exposure_rename)) res$value <- exposure
  # Maaslin adjust p-values for all coefficients, modify to be for only the 
  # exposure
  #res <- subset(res, !is.na(coef))
  res <- res %>% arrange(pval) %>%
    mutate(qval = p.adjust(pval, method = "fdr"))
  
  return(res = res)
}

## meta-analysis for identifying features from multiple datasets 
rma_wrapper <- function(maaslin_fits, 
                        method = "REML",
                        output = getwd(),
                        forest_plot = NULL, 
                        rma_conv = 1e-6,
                        rma_maxit = 1000,
                        verbose = TRUE) {
  lvl_batch <- names(maaslin_fits)
  n_batch <- length(lvl_batch)
  exposure <- unique(maaslin_fits[[1]]$metadata)
  values_exposure <- unique(maaslin_fits[[1]]$value)
  features <- unique(maaslin_fits[[1]]$feature)
  
  for (i in 1:length(maaslin_fits)) {
    rownames(maaslin_fits[[i]]) <- maaslin_fits[[i]]$feature
    maaslin_fits[[i]] <- maaslin_fits[[i]][features, ]
  }
  
  l_results <- list()
  for(value_exposure in values_exposure) {
    i_result <- data.frame(matrix(NA,
                                  nrow = length(features),
                                  ncol = 11 + length(lvl_batch)))
    colnames(i_result) <- c("feature",
                            "exposure",
                            "coef",
                            "stderr",
                            "pval",
                            "k",
                            "tau2",
                            "stderr.tau2",
                            "pval.tau2",
                            "I2",
                            "H2",
                            paste0("weight_", lvl_batch))
    i_result$feature <- features
    i_result$exposure <- value_exposure
    rownames(i_result) <- i_result$feature
    if(!is.null(forest_plot)) 
      pdf(paste0(output, "/",
                 exposure, "_", value_exposure, "_",
                 forest_plot),
          width = 6,
          height = 4 + ifelse(n_batch > 4,
                              (n_batch - 4) * 0.5,
                              0))
    # sanity check
    if(any(features != maaslin_fits[[2]][
      maaslin_fits[[2]]$value == value_exposure, 
      "feature"]))
      stop("Feature names don't match between maaslin_fits components!")
    betas <- vapply(
      maaslin_fits, 
      function(i_maaslin_fit)
        i_maaslin_fit[i_maaslin_fit$value == value_exposure, 
                      "coef", drop = TRUE],
      rep_len(0.0, length(features))
    )
    sds <- vapply(
      maaslin_fits, 
      function(i_maaslin_fit)
        i_maaslin_fit[i_maaslin_fit$value == value_exposure, 
                      "stderr", drop = TRUE],
      rep_len(0.0, length(features))
    )
    pvals <- vapply(
      maaslin_fits, 
      function(i_maaslin_fit)
        i_maaslin_fit[i_maaslin_fit$value == value_exposure, 
                      "pval", drop = TRUE],
      rep_len(0.0, length(features))
    )
    rownames(betas) <- rownames(sds) <- rownames(pvals) <- features
    ind_features <- !is.na(betas) & !is.na(sds) & (sds != 0)
    count_feature <- apply(ind_features, 1, sum)
    for(feature in features) {
      if(count_feature[feature] >= 2) {
        i_log <- catchToList(
          metafor::rma.uni(yi = betas[feature, ind_features[feature, ]],
                           sei = sds[feature, ind_features[feature, ]],
                           slab = lvl_batch[ind_features[feature, ]],
                           method = method,
                           control = list(threshold = rma_conv,
                                          maxiter = rma_maxit))
        )
        if(!is.null(i_log$error)) {
          warning("Fitting rma on feature ", feature, ";\n",
                  i_log$error)
          next
        }
        if(!is.null(i_log$warnings))
          warning("Fitting rma on feature ", feature, ";\n",
                  i_log$warnings)
        i_rma_fit <- i_log$value
        wts <- metafor::weights.rma.uni(i_rma_fit)
        i_result[feature, c("coef",
                            "stderr",
                            "pval",
                            "k",
                            "tau2",
                            "stderr.tau2",
                            "pval.tau2",
                            "I2",
                            "H2",
                            paste0("weight_",
                                   names(wts))
        )] <- c(unlist(i_rma_fit[c("beta", 
                                   "se",
                                   "pval",
                                   "k",
                                   "tau2",
                                   "se.tau2",
                                   "QEp",
                                   "I2",
                                   "H2")]),
                wts)
        if(i_rma_fit$pval < 0.05 & !is.null(forest_plot))
          metafor::forest(
            i_rma_fit,
            xlab = shorten_name(feature, cutoff = 5),
            slab = shorten_name(lvl_batch[ind_features[feature, ]], cutoff = 3))
      }
      if(count_feature[feature] == 1) {
        i_ind_features <- ind_features[feature, ]
        tmp_batch <- lvl_batch[i_ind_features]
        i_result[feature, c("coef",
                            "stderr",
                            "pval",
                            "k",
                            paste0("weight_",
                                   tmp_batch)
        )] <- c(betas[feature, i_ind_features],
                sds[feature, i_ind_features],
                pvals[feature, i_ind_features],
                1,
                100)
      }
    }
    if(!is.null(forest_plot)) dev.off()
    i_result <- i_result %>% arrange(desc(pval)) %>%
      mutate(qval.bonf = p.adjust(pval, method = "bonf"))
    i_result <- i_result %>% arrange(desc(pval)) %>%
      mutate(qval.fdr = p.adjust(pval, method = "fdr"))
    
    l_results[[value_exposure]] <- i_result
  }
  results <- Reduce("rbind", l_results)
  return(results)
}

## identifying differential species
## confounding analysis
load("Filtered_Spe_abund.Rdata")
rownames(metadata) <- metadata$Sample_ID
metadata[which(metadata$Project %in% "PRJNA447983_cohort1"), "Project"] <- "PRJNA447983"

confounder.sum <- list()

for (pro in unique(metadata$Project)) {
  a <- metadata %>% filter(Project %in% pro)
  a.summ <- a %>% group_by(Group) %>% 
    summarise(age.median = median(Age, na.rm = T),
              bmi.median = median(BMI, na.rm = T),
              male = sum(Gender == "male", na.rm = T),
              female = sum(Gender == "female", na.rm = T), 
              HDC.median = median(humanper, na.rm = T))
  ## testing age
  if (!(is.na(a.summ$age.median[1]))){
    a.test <- kruskal.test(a$Age ~ a$Group)
    a.summ$age.pval <- a.test$p.value
  }else {
    a.summ$age.pval <- NA
  }
  ## testing bmi
  if (!(is.na(a.summ$bmi.median[1]))){
    a.test <- kruskal.test(a$BMI ~ a$Group)
    a.summ$BMI.pval <- a.test$p.value
  }else {
    a.summ$BMI.pval <- NA
  }
  ## testing Gender
  if (a.summ$male != 0 & a.summ$female != 0) {
    a.test <- a.summ[,c("male","female")]
    a.chis <- chisq.test(a.test)
    a.summ$Gender.pval <- a.chis$p.value
  }else {
    a.summ$Gender.pval <- NA
  }
  ## testing HDC
  a.test <- kruskal.test(a$humanper ~ a$Group)
  a.summ$HDC.pval <- a.test$p.value
  
  a.summ$Project <- pro
  confounder.sum[[pro]] <- a.summ
  
}
confounder.sum.df <- confounder.sum %>% reduce(rbind)

confounder.sum.df %>%
  filter(age.pval < 0.05)
confounder.sum.df %>%
  filter(BMI.pval < 0.05)

## identifying differential species
uc.maaslin <- list()
cd.maaslin <- list()
crc.maaslin <- list()

## for UC
uc.maaslin[["PRJEB1220"]] <- maaslin2_refunc(metadata = metadata, 
                                             pro = "PRJEB1220", case = "UC", 
                                             covariate = c("Age", "BMI"), 
                                             abund.s.list.filter.df = abund.s.list.filter$UC,
                                             transform = "LOG")

uc.maaslin[["PRJNA400072"]] <- maaslin2_refunc(metadata = metadata,
                                               pro = "PRJNA400072", case = "UC", 
                                               covariate = c("Age"), 
                                               abund.s.list.filter.df = abund.s.list.filter$UC,
                                               transform = "LOG")

uc.maaslin[["PRJNA389280"]] <- maaslin2_refunc(metadata = metadata, 
                                               pro = "PRJNA389280", case = "UC", 
                                                 covariate = NULL, 
                                                 abund.s.list.filter.df = abund.s.list.filter$UC,
                                               transform = "LOG")

uc.maaslin.meta <- rma_wrapper(uc.maaslin)
uc.kw.spe <- uc.maaslin.meta %>% filter(qval.fdr < 0.05) %>%
  pull(feature)

## for CD
cd.maaslin[["SRP057027"]] <- maaslin2_refunc(metadata = metadata, 
                                             pro = "SRP057027", case = "CD", 
                                             covariate = NULL, 
                                             abund.s.list.filter.df = abund.s.list.filter$CD,
                                             transform = "LOG")

cd.maaslin[["PRJNA400072"]] <- maaslin2_refunc(metadata = metadata,
                                               pro = "PRJNA400072", case = "CD", 
                                               covariate = c("Age"), 
                                               abund.s.list.filter.df = abund.s.list.filter$CD,
                                               transform = "LOG")
cd.maaslin[["PRJNA389280"]] <- maaslin2_refunc(metadata = metadata, 
                                               pro = "PRJNA389280", case = "CD", 
                                               covariate = NULL, 
                                               abund.s.list.filter.df = abund.s.list.filter$CD,
                                               transform = "LOG")

cd.maaslin.meta <- rma_wrapper(cd.maaslin)
cd.kw.spe <- cd.maaslin.meta %>% filter(qval.fdr < 0.05) %>%
  pull(feature)


## for CRC
for (i in c("PRJEB27928", "PRJEB6070", "PRJEB10878")) {
  crc.maaslin[[i]] <- maaslin2_refunc(metadata = metadata, 
                                      pro = i, case = "CRC", 
                                      covariate = c("Age"), 
                                      abund.s.list.filter.df = abund.s.list.filter$CRC,
                                      transform = "LOG")
}

for (i in setdiff(crc.project, c("PRJEB27928", "PRJEB6070", "PRJEB10878"))) {
  crc.maaslin[[i]] <- maaslin2_refunc(metadata = metadata, 
                                      pro = i, case = "CRC", 
                                      covariate = NULL, 
                                      abund.s.list.filter.df = abund.s.list.filter$CRC,
                                      transform = "LOG")
}

crc.maaslin.meta <- rma_wrapper(crc.maaslin)
crc.kw.spe <- crc.maaslin.meta %>% filter(qval.fdr < 0.05) %>%
  pull(feature)


## save results 
names(cd.maaslin) <- paste("CD", names(cd.maaslin), sep = "_")
names(uc.maaslin) <- paste("UC", names(uc.maaslin), sep = "_")

maaslin <- c(crc.maaslin, cd.maaslin, uc.maaslin)
maaslin.meta <- list(crc.maaslin.meta = crc.maaslin.meta, 
                     cd.maaslin.meta = cd.maaslin.meta,
                     uc.maaslin.meta = uc.maaslin.meta)


fit.lm.results.df <- rbind(crc.maaslin.meta[,c("feature","exposure","coef","k","qval.fdr")], 
                           cd.maaslin.meta[,c("feature","exposure","coef","k","qval.fdr")], 
                           uc.maaslin.meta[,c("feature","exposure","coef","k","qval.fdr")])
fit.lm.results.df.sig <- fit.lm.results.df %>% filter(qval.fdr < 0.05)

fit.lm.results.df.sig.merge <- merge(fit.lm.results.df.sig, 
                                     phylum.to.species, 
                                     by.x = "feature", 
                                     by.y = "Species", 
                                     all.x = T)

fit.lm.results.df.sig.merge$final.spe <- gsub(fit.lm.results.df.sig.merge$feature,
                                               pattern = "_",
                                               replacement = " ",
                                               fixed = T)
fit.lm.results.df.sig.merge$final.spe <- gsub(fit.lm.results.df.sig.merge$final.spe,
                                        pattern = "unclassified",
                                        replacement = "spp.",
                                        fixed = T)
fit.lm.results.df.sig.merge$Trend <- "Control-enriched"
fit.lm.results.df.sig.merge[which(fit.lm.results.df.sig.merge$coef >0),"Trend"] <- "Case-enriched"

fit.lm.results.df.sig.merge <- fit.lm.results.df.sig.merge %>%
  arrange(Phylum)

kw.spe.df.bottom <- fit.lm.results.df.sig.merge %>% 
  mutate(final.spe = factor(final.spe, 
                            levels = unique(fit.lm.results.df.sig.merge$final.spe))) %>%
  mutate(exposure = factor(exposure, levels = c("CRC","CD","UC"))) %>%
  ggplot(aes(x = final.spe, y = exposure)) +
  geom_tile(aes(fill = Trend), 
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

kw.spe.df.top <- fit.lm.results.df.sig.merge %>%
  mutate(final.spe = factor(final.spe, 
                            levels = unique(fit.lm.results.df.sig.merge$final.spe))) %>%
  ggplot(aes(x = final.spe, y = 1)) + geom_tile(aes(fill = Phylum)) +
  theme_classic() +
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

kw.spe.df.plot <- plot_grid(kw.spe.df.top, kw.spe.df.bottom, 
                            ncol = 1, align = 'v', rel_heights = c(0.3,1))

pdf("../pdf/differntial_species.pdf", 
    height = 5, width = 8)
kw.spe.df.plot
dev.off()


save(confounder.sum.df, maaslin, maaslin.meta, 
     fit.lm.results.df, fit.lm.results.df.sig,
     file = "../RData/Maaslin_confounders_bork.Rdata")

write.table(confounder.sum.df,  
            file = "../SuppleData/confounder.txt",
            col.names = T, row.names = F, 
            sep = "\t", quote = F)
write.table(fit.lm.results.df,  
            file = "../SuppleData/metafor_species_coef.txt",
            col.names = T, row.names = F, 
            sep = "\t", quote = F)
write.table(fit.lm.results.df.sig.merge,  
            file = "../SuppleData/metafor_species_coef_signif.txt",
            col.names = T, row.names = F, 
            sep = "\t", quote = F)

## identifying differential pathways
## confounding analysis
load("Filtered_Pathways_abund.Rdata")
rownames(metadata) <- metadata$Sample_ID
metadata[which(metadata$Project %in% "PRJNA447983_cohort1"), "Project"] <- "PRJNA447983"

## identifying differential pathways
uc.path.maaslin <- list()
cd.path.maaslin <- list()
crc.path.maaslin <- list()

## for UC
uc.path.maaslin[["PRJEB1220"]] <- maaslin2_refunc(metadata = metadata, 
                                             pro = "PRJEB1220", case = "UC", 
                                             covariate = c("Age", "BMI"), 
                                             abund.s.list.filter.df = path.abund.fil.list$UC,
                                             transform = "LOG")

uc.path.maaslin[["PRJNA400072"]] <- maaslin2_refunc(metadata = metadata,
                                               pro = "PRJNA400072", case = "UC", 
                                               covariate = c("Age"), 
                                               abund.s.list.filter.df = path.abund.fil.list$UC,
                                               transform = "LOG")

uc.path.maaslin[["PRJNA389280"]] <- maaslin2_refunc(metadata = metadata, 
                                               pro = "PRJNA389280", case = "UC", 
                                               covariate = NULL, 
                                               abund.s.list.filter.df = path.abund.fil.list$UC,
                                               transform = "LOG")

uc.path.maaslin.meta <- rma_wrapper(uc.path.maaslin)
uc.kw.path <- uc.path.maaslin.meta %>% filter(qval.fdr < 0.05) %>%
  pull(feature)

## for CD
cd.path.maaslin[["SRP057027"]] <- maaslin2_refunc(metadata = metadata, 
                                             pro = "SRP057027", case = "CD", 
                                             covariate = NULL, 
                                             abund.s.list.filter.df = path.abund.fil.list$CD,
                                             transform = "LOG")

cd.path.maaslin[["PRJNA400072"]] <- maaslin2_refunc(metadata = metadata,
                                               pro = "PRJNA400072", case = "CD", 
                                               covariate = c("Age"), 
                                               abund.s.list.filter.df = path.abund.fil.list$CD,
                                               transform = "LOG")
cd.path.maaslin[["PRJNA389280"]] <- maaslin2_refunc(metadata = metadata, 
                                               pro = "PRJNA389280", case = "CD", 
                                               covariate = NULL, 
                                               abund.s.list.filter.df = path.abund.fil.list$CD,
                                               transform = "LOG")

cd.path.maaslin.meta <- rma_wrapper(cd.path.maaslin)
cd.kw.path <- cd.path.maaslin.meta %>% filter(qval.fdr < 0.05) %>%
  pull(feature)


## for CRC
for (i in c("PRJEB27928", "PRJEB6070", "PRJEB10878")) {
  crc.path.maaslin[[i]] <- maaslin2_refunc(metadata = metadata, 
                                      pro = i, case = "CRC", 
                                      covariate = c("Age"), 
                                      abund.s.list.filter.df = path.abund.fil.list$CRC,
                                      transform = "LOG")
}

for (i in setdiff(crc.project, c("PRJEB27928", "PRJEB6070", "PRJEB10878"))) {
  crc.path.maaslin[[i]] <- maaslin2_refunc(metadata = metadata, 
                                      pro = i, case = "CRC", 
                                      covariate = NULL, 
                                      abund.s.list.filter.df = path.abund.fil.list$CRC,
                                      transform = "LOG")
}

crc.path.maaslin.meta <- rma_wrapper(crc.path.maaslin)
crc.kw.path <- crc.path.maaslin.meta %>% filter(qval.fdr < 0.05) %>%
  pull(feature)


## save results 
names(cd.path.maaslin) <- paste("CD", names(cd.path.maaslin), sep = "_")
names(uc.path.maaslin) <- paste("UC", names(uc.path.maaslin), sep = "_")

path.maaslin <- c(crc.path.maaslin, cd.path.maaslin, uc.path.maaslin)
path.maaslin.meta <- list(crc.path.maaslin.meta = crc.path.maaslin.meta, 
                     cd.path.maaslin.meta = cd.path.maaslin.meta,
                     uc.path.maaslin.meta = uc.path.maaslin.meta)


path.fit.lm.results.df <- rbind(crc.path.maaslin.meta[,c("feature","exposure","coef","k","qval.fdr")], 
                                cd.path.maaslin.meta[,c("feature","exposure","coef","k","qval.fdr")], 
                                uc.path.maaslin.meta[,c("feature","exposure","coef","k","qval.fdr")])
path.fit.lm.results.df.sig <- path.fit.lm.results.df %>% filter(qval.fdr < 0.05)

path.fit.lm.results.df.sig$Trend <- "Control-enriched"
path.fit.lm.results.df.sig[which(path.fit.lm.results.df.sig$coef >0),"Trend"] <- "Case-enriched"

kw.path.df.bottom <- path.fit.lm.results.df.sig  %>%
  mutate(exposure = factor(exposure, levels = c("CRC","CD","UC"))) %>%
  ggplot(aes(x = feature, y = exposure)) +
  geom_tile(aes(fill = Trend), 
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

pdf("../pdf/differntial_pathways.pdf", 
    height = 7, width = 8)
kw.path.df.bottom
dev.off()

save(confounder.sum.df, path.maaslin, path.maaslin.meta, 
     path.fit.lm.results.df, path.fit.lm.results.df.sig,
     file = "../RData/Maaslin_confounders_pathways_bork.Rdata")

write.table(path.fit.lm.results.df,  
            file = "../SuppleData/metafor_pathways_coef.txt",
            col.names = T, row.names = F, 
            sep = "\t", quote = F)
write.table(path.fit.lm.results.df.sig,  
            file = "../SuppleData/metafor_pathways_coef_signif.txt",
            col.names = T, row.names = F, 
            sep = "\t", quote = F)