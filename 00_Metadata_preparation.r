## load HDC reads
crc.reads <- read.delim("/mnt/raid5/puzi/IBD/otherproject/result/all_reads.txt",
                        header = T, sep = "\t", as.is = T)
crc.reads.statistics <- crc.reads %>%
  group_by(Project, Sample_ID, Group, type) %>%
  summarise(reads_sum = sum(reads))
crc.reads.statistics.2 <- dcast(crc.reads.statistics, 
                                Project + Sample_ID + Group ~ type, value.var = "reads_sum")
crc.reads.statistics.2$HDC <- crc.reads.statistics.2$conta/(crc.reads.statistics.2$clean + crc.reads.statistics.2$conta)
crc.reads.statistics.2$sum <- crc.reads.statistics.2$clean + crc.reads.statistics.2$conta

## load metadata
crc.meta <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/metadata/allCRC_metadata.txt",
                       header = T, sep = "\t", as.is = T)
crc.meta.merge <- merge(crc.meta, 
                        crc.reads.statistics.2[,c("Sample_ID", "clean", "conta", 
                                                  "sum")], 
                        by = "Sample_ID",
                        all = T)
write.table(crc.meta.merge, 
            file = "/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/metadata/allCRC_metadata.txt",
           col.names = T, row.names = F, sep = "\t", quote = F)

# ********
crc.meta <- crc.meta.merge
rm(crc.meta.merge, crc.reads.statistics, crc.reads.statistics.2,
   crc.reads)

crc.meta[which(crc.meta$Country %in% "Germany"), "Country"] <- "GER"
crc.meta[which(crc.meta$Country %in% "ITA_1"), "Country"] <- "ITA"
crc.meta[which(crc.meta$Country %in% "France"), "Country"] <- "FRA"
crc.meta[which(crc.meta$Country %in% "Japan"), "Country"] <- "JAP"

## filter the first collections of PRJNA389280

PRJNA389280.meta <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/metadata/PRJNA389280_meta.txt",
                               header = T, sep = "\t", as.is = T)
PRJNA389280.sample.1 <- PRJNA389280.meta %>% 
  filter(Participant.ID.x %in% names(table(PRJNA389280.meta$Participant.ID.x)[table(PRJNA389280.meta$Participant.ID.x)==1])) %>%
  pull(Run)

PRJNA389280.sample.2 <- names(table(PRJNA389280.meta$Participant.ID.x)[table(PRJNA389280.meta$Participant.ID.x)>1])

PRJNA389280.sample.3 <- c()

for (i in PRJNA389280.sample.2) {
  PRJNA389280.meta.2 <- PRJNA389280.meta %>%
    filter(Participant.ID.x %in% i) 
  a <- PRJNA389280.meta.2[which.min(PRJNA389280.meta.2$week_num), "Run"]
  PRJNA389280.sample.3 <- c(PRJNA389280.sample.3, a)
}

PRJNA389280.sample.final <- c(PRJNA389280.sample.1, PRJNA389280.sample.3)
PRJNA389280.meta.final <- subset(PRJNA389280.meta, 
                                 Run %in% PRJNA389280.sample.final)

PRJNA389280.csv <- read.csv("/mnt/raid5/puzi/IBD/389280/hmp2_metadata.csv",
                            header = T, sep = ",", as.is = T)
PRJNA389280.csv.filter <- PRJNA389280.csv %>%
  filter(site_sub_coll %in% PRJNA389280.meta.final$site_sub_coll) %>%
  filter(data_type %in% "metagenomics")
PRJNA389280.csv.filter.meta <- PRJNA389280.csv.filter %>%
  select(Participant.ID, site_sub_coll, 
         consent_age, sex, BMI, diagnosis) %>%
  distinct(.keep_all = T)

PRJNA389280.meta.final.merge <- merge(PRJNA389280.meta.final, 
                                      PRJNA389280.csv.filter.meta, 
                                      by.x = c("site_sub_coll", 
                                             "Participant.ID.x"),
                                      by.y = c("site_sub_coll", 
                                               "Participant.ID"))
PRJNA389280.meta.final.merge <- PRJNA389280.meta.final.merge[,-c(38, 42)]
write.table(PRJNA389280.meta.final.merge,
            file = "/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/metadata/PRJNA389280_meta_final_complete.txt",
            col.names = T, row.names = F, sep = "\t",
            quote = F)
rm(PRJNA389280.meta, PRJNA389280.sample.1, PRJNA389280.sample.2, PRJNA389280.sample.3, 
   PRJNA389280.csv, PRJNA389280.meta.final)

PRJNA389280.meta <- PRJNA389280.meta.final.merge[,c("Run", "site_sub_coll", "Participant.ID.x",
                                                    "BMI", "consent_age", "sex",
                                                    "diagnosis", "clean",
                                                    "conta", "sum", "human_per")]
names(PRJNA389280.meta) <- c("Sample_ID", "site_sub_coll", "Participant.ID.HMP",
                             "BMI", "Age", "Gender",
                             "Group", "clean",
                             "conta", "sum", "humanper")
PRJNA389280.meta$Country <- "USA"
PRJNA389280.meta$Project <- "PRJNA389280"
PRJNA389280.meta <- PRJNA389280.meta[,c("Sample_ID","Age","Gender","BMI",
                                        "Country","Project","Group","humanper",
                                        "site_sub_coll", "Participant.ID.HMP",
                                        "clean", "conta", "sum")]

rm(PRJNA389280.csv.filter, PRJNA389280.csv.filter.meta,
   PRJNA389280.meta.2, PRJNA389280.meta.week0, PRJNA389280.meta.final.merge)

# ********
PRJNA389280.meta[which(PRJNA389280.meta$Group %in% "nonIBD"), "Group"] <- "CTR"

## filter the first collections of PRJEB1220
PRJEB1220.meta <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/metadata/PRJEB1220_meta.txt",
                             header = T, sep = "\t", as.is = T)
PRJEB1220.meta$Time.Sample <- unlist(lapply(strsplit(PRJEB1220.meta$Sample.ID, 
                                                     split = ".",
                                                     fixed = T),
                                            function(data){ y <- data[3]}))
PRJEB1220.meta[is.na(PRJEB1220.meta$Time.Sample),"Time.Sample"] <- 0
PRJEB1220.meta.2 <- subset(PRJEB1220.meta, 
                           Health.Status %in% c("Healthy", "Ulcerative colitis"))

## ------------------------calculating HUMAN DNA--------------------------------
PRJEB1220.reads <- read.delim("/mnt/raid5/puzi/CDdata/PRJEB1220/PRJEB1220_reads.txt",
                                header = F, sep = " ", as.is = T)
PRJEB1220.run <- read.delim("/mnt/raid5/puzi/CDdata/PRJEB1220/runinfo",
                            header = T, sep = "\t", as.is = T)
names(PRJEB1220.reads) <- c("filenames","reads")
PRJEB1220.reads$sample <- substr(PRJEB1220.reads$filenames, 1,9)
print(nrow(PRJEB1220.reads))
PRJEB1220.reads$reads <- as.numeric(PRJEB1220.reads$reads)
PRJEB1220.reads$reads <- PRJEB1220.reads$reads/4
PRJEB1220.reads$type <- substr(PRJEB1220.reads$filenames, 9+14, 9+18)

print(length(intersect(unique(PRJEB1220.reads$sample), 
                       PRJEB1220.run$Run)))

PRJEB1220.reads <- subset(PRJEB1220.reads, 
                          sample %in% PRJEB1220.run$Run)
##merge metadata and reads.data
PRJEB1220.reads.merge <- merge(PRJEB1220.reads, 
                               PRJEB1220.run[,c("Alias", "Run")], 
                               by.x = "sample", 
                               by.y = "Run", all.x = T)
PRJEB1220.reads.all <- PRJEB1220.reads.merge %>% 
  group_by(Alias, type) %>% 
  summarise(all = sum(reads))
PRJEB1220.reads.all.dcast <- dcast(PRJEB1220.reads.all, 
                                   Alias ~ type)
PRJEB1220.reads.all.dcast$all <- PRJEB1220.reads.all.dcast$clean+PRJEB1220.reads.all.dcast$conta
PRJEB1220.reads.all.dcast$humanper <- PRJEB1220.reads.all.dcast$conta/PRJEB1220.reads.all.dcast$all
PRJEB1220.reads.all.dcast <- subset(PRJEB1220.reads.all.dcast, !all == 0)
PRJEB1220.reads.all.dcast$Alias <- gsub(PRJEB1220.reads.all.dcast$Alias, 
                                        pattern = "-", replacement = ".",
                                        fixed = T)

rm(PRJEB1220.reads, PRJEB1220.run, 
   PRJEB1220.reads.merge, PRJEB1220.reads.all)

write.table(PRJEB1220.reads.all.dcast, 
            file = "/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/metadata/PRJEB1220_HDC.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)

PRJEB1220.meta.2.merge <- merge(PRJEB1220.meta.2, 
                                PRJEB1220.reads.all.dcast, 
                                by.x = "Sample.ID", 
                                by.y = "Alias", 
                                all.x = T)
PRJEB1220.meta.2.merge <- subset(PRJEB1220.meta.2.merge, 
                                 !is.na(humanper))
repeated.sample <- names(table(PRJEB1220.meta.2.merge$Individual.ID)[table(PRJEB1220.meta.2.merge$Individual.ID)>1])

## ------------------------calculating HUMAN DNA--------------------------------

## filter the final samples of PRJEB1220
table(subset(PRJEB1220.meta.2.merge, Individual.ID %in% repeated.sample)$Time.Sample)
PRJEB1220.meta.3 <- subset(PRJEB1220.meta.2.merge, 
                           !(Individual.ID %in% repeated.sample & !(Time.Sample %in% 0)))

PRJEB1220.meta.final <- PRJEB1220.meta.3[,c("Sample.ID","Age","Gender","BMI",
                                            "Nationality", "Health.Status", "humanper",
                                            "clean", "conta","all")]
names(PRJEB1220.meta.final) <- c("Sample_ID","Age","Gender",
                                 "BMI", "Country", "Group",
                                 "humanper", "clean", "conta","sum")
PRJEB1220.meta.final[which(PRJEB1220.meta.final$Group %in% "Healthy"), "Group"] <- "CTR"
PRJEB1220.meta.final[which(PRJEB1220.meta.final$Group %in% "Ulcerative colitis"), "Group"] <- "UC"

PRJEB1220.meta.final$Project <- "PRJEB1220"

# ********
PRJEB1220.meta.final <- PRJEB1220.meta.final[,c("Sample_ID","Age","Gender","BMI",
                                                "Country", "Project","Group",
                                                "humanper", "clean", "conta","sum")]
PRJEB1220.meta <- PRJEB1220.meta.final

rm(PRJEB1220.meta.final, PRJEB1220.meta.2, PRJEB1220.meta.2.merge, 
   PRJEB1220.meta.3, PRJEB1220.meta.3.merge, 
   PRJEB1220.reads.all.dcast)


## filter the first collections of SRP057027
SRP057027.meta <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/metadata/SRP057027_meta.txt",
                             header = T, sep = "\t", as.is = T)
SRP057027.meta.final <- subset(SRP057027.meta, 
                               Sample.group %in% c("Baseline","Control"))

SRP057027.meta.final[,c("Age","Gender","BMI","Country")] <- NA

##
SRP057027.reads <- read.delim("/mnt/raid5/puzi/IBD/otherproject/CD_reads_HDC.txt",
                              header = T, sep = "\t", as.is = T)
SRP057027.meta.final.merge <- merge(SRP057027.meta.final, 
                                    SRP057027.reads[,c("sample", "allclean",
                                                       "allcontam", "all")],
                                    by.x = "run_id", 
                                    by.y = "sample",
                                    all.x = T)

SRP057027.meta <- SRP057027.meta.final.merge[,c("run_id", "Age","Gender","BMI",
                                          "Country","ProjectID","Sample.group",
                                          "our_HDC", "allclean",
                                          "allcontam", "all")]
names(SRP057027.meta) <- c("Sample_ID", "Age","Gender","BMI",
                           "Country", "Project", "Group",
                           "humanper", "clean","conta","sum")
rm(SRP057027.meta.final, SRP057027.reads, SRP057027.meta.final.merge)

# ********
# lack basic info, but source literature showed it is a case-control study
SRP057027.meta[which(SRP057027.meta$Group %in% "Baseline"), "Group"] <- "CD"
SRP057027.meta[which(SRP057027.meta$Group %in% "Control"), "Group"] <- "CTR"


##
PRJNA400072.meta <- read.delim("/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/metadata/PRJNA400072_meta.txt",
                               header = T, sep = "\t", as.is = T)
PRJNA400072.meta[,c("Gender","BMI")] <- NA


## ----------------------calculating human DNA ----------------------
PRJNA400072.reads <- read.delim("/mnt/raid5/puzi/CDdata/PRJNA400072/prjna400072_reads.txt",
                                header = F, sep = " ", as.is = T)
names(PRJNA400072.reads) <- c("filenames","reads")
PRJNA400072.reads$sample <- substr(PRJNA400072.reads$filenames, 1,10)
print(nrow(PRJNA400072.reads))
PRJNA400072.reads$reads <- as.numeric(PRJNA400072.reads$reads)
PRJNA400072.reads$reads <- PRJNA400072.reads$reads/4
PRJNA400072.reads$type <- substr(PRJNA400072.reads$filenames, 10+14, 10+18)

print(length(intersect(unique(PRJNA400072.reads$sample), 
                       PRJNA400072.meta$Run)))


PRJNA400072.reads.all <- PRJNA400072.reads %>% 
  group_by(sample, type) %>% 
  summarise(all = sum(reads))
PRJNA400072.reads.all.dcast <- dcast(PRJNA400072.reads.all, 
                                   sample ~ type)
PRJNA400072.reads.all.dcast$all <- PRJNA400072.reads.all.dcast$clean+PRJNA400072.reads.all.dcast$conta
PRJNA400072.reads.all.dcast$humanper <- PRJNA400072.reads.all.dcast$conta/PRJNA400072.reads.all.dcast$all
## ----------------------calculating human DNA ----------------------

PRJNA400072.meta.merge <- merge(PRJNA400072.meta, 
                                PRJNA400072.reads.all.dcast,
                                by.x = "Run", 
                                by.y = "sample", 
                                all = T)
PRJNA400072.meta.2 <- PRJNA400072.meta.merge[,c("Run", "Sample_Data.Age",
                                          "Gender","BMI", "geo_loc_name",
                                          "BioProject", "Diagnosis", "humanper",
                                          "clean","conta","all")]
write.table(PRJNA400072.meta.merge, 
            file = "/mnt/raid6/puzi/2019intestine/IBD_CRC/CDdata/result_ArticleII/metadata/PRJNA400072_meta.txt",
            col.names = T, row.names = F, sep = "\t", quote = F)
rm(PRJNA400072.reads.all, 
   PRJNA400072.reads,
   PRJNA400072.reads.all.dcast)

names(PRJNA400072.meta.2) <- c("Sample_ID", "Age",
                               "Gender","BMI", "Country",
                               "Project", "Group", "humanper",
                               "clean","conta","sum")
# ********
PRJNA400072.meta.2[which(PRJNA400072.meta.2$Group %in% "Control"), "Group"] <- "CTR"
PRJNA400072.meta <- PRJNA400072.meta.2

rm(PRJNA400072.meta.2, PRJNA400072.meta.merge)


# ********final Metadata for IBD******** #
ibd.meta <- rbind(PRJEB1220.meta, 
                  PRJNA389280.meta[,-c(9:10)],
                  PRJNA400072.meta, 
                  SRP057027.meta)
ibd.meta %>% group_by(Project, Group) %>% 
  summarise(count = n())

## statisticals of reads
crc.meta.reads.stat <- crc.meta %>%
  group_by(Project, Group) %>%
  summarise(median_reads = median(clean, na.rm = T), 
            min_reads = min(clean, na.rm = T),
            max_reads = max(clean, na.rm = T))

ibd.meta.reads.stat <- ibd.meta %>%
  group_by(Project, Group) %>%
  summarise(median_reads = median(clean, na.rm = T), 
            min_reads = min(clean, na.rm = T),
            max_reads = max(clean, na.rm = T))