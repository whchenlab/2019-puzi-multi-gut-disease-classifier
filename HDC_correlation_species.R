## -----load data------##
load("data_prepare.RData")

## -----taxonomic level ----##
## -----crc.abund.merge.list: HDC & abundances of species as columns;----
## -----crc.abund.merge.list: sample names as rows; ---------------------
## -----cd.abund.merge.list & uc.abund.merge.list : same as crc.abund.merge.list ########
## -----CRC list -------##
crc.abund.cor.list <- list()
for(pro in names(crc.abund.merge.list)){
  crc.abund.cor.list[[pro]] <- cor.func(crc.abund.merge.list[[pro]], 
                                        type = "CRC", project = pro)
}

## -----CD list --------##
cd.abund.cor.list <- list()
for(pro in names(cd.abund.merge.list)){
  cd.abund.cor.list[[pro]] <- cor.func(cd.abund.merge.list[[pro]], 
                                        type = "CD", project = pro)
}

## -----UC list --------##
uc.abund.cor.list <- list()
for(pro in names(uc.abund.merge.list)){
  uc.abund.cor.list[[pro]] <- cor.func(uc.abund.merge.list[[pro]], 
                                          type = "UC", project = pro)
}

cor.data.list <- lapply(c(crc.abund.cor.list,
                          uc.abund.cor.list,
                          cd.abund.cor.list), 
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

save.image("HDC_related_species.Rdata")
