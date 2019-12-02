## -----load data------##
load("data_prepare.RData")

## -----taxonomic level ----##
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


## -----functional level ----##
crc.path.cor.list <- list()
for(pro in names(crc.abund.path.merge.list)){
  crc.path.cor.list[[pro]] <- cor.func(crc.abund.path.merge.list[[pro]], 
                                        type = "CRC", project = pro)
}

## -----CD list --------##
cd.path.cor.list <- list()
for(pro in names(cd.abund.path.merge.list)){
  cd.path.cor.list[[pro]] <- cor.func(cd.abund.path.merge.list[[pro]], 
                                        type = "CD", project = pro)
}

## -----UC list --------##
uc.path.cor.list <- list()
for(pro in names(uc.abund.path.merge.list)){
  uc.path.cor.list[[pro]] <- cor.func(uc.abund.path.merge.list[[pro]], 
                                          type = "UC", project = pro)
}

cor.path.data.list <- lapply(c(crc.path.cor.list,
                          uc.path.cor.list,
                          cd.path.cor.list), 
                        function(data){y <- data$data.cor.sum})
cor.path.data <- reduce(cor.path.data.list, rbind)

## ----filter results with P-value < 0.05
cor.path.data.sig <- subset(cor.path.data, P_value<0.05)
cor.path.data.posi <- subset(cor.path.data.sig, coef>0)
cor.path.data.nega <- subset(cor.path.data.sig, coef<0)
delete.spe <- intersect(cor.path.data.posi$Species, cor.path.data.nega$Species)
cor.path.data.sig <- subset(cor.path.data.sig, !Species %in% delete.spe)

cor.path.sig.sum <- cor.path.data.sig %>% 
  group_by(Species, Type) %>% 
  summarise(count= n())
cor.path.sig.sum.2 <- dcast(cor.path.sig.sum, 
                            Species ~Type)

sig.path <- subset(cor.path.sig.sum.2, 
                  CRC>1 | UC>1 | CD>1)$Species
sig.path.df <- subset(cor.path.sig.sum.2, 
                     CRC>1 | UC>1 | CD>1)


save.image("HDC_related.Rdata")