setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep1-replace-by0829data/")
samples <- list.dirs(recursive = F)
total_sample_stat <- list()
for(sample in samples){
  setwd(sample)
    mapping_rate <- read.table("CRISPResso_mapping_statistics.txt", 
                               sep="\t", header = T)
    total_sample_stat[[sample]] <- c(unlist(mapping_rate[1,c(1,3)]))
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
total_sample_stat <- data.frame(do.call(rbind, total_sample_stat))
total_sample_stat$unmapped <- total_sample_stat[,1] - total_sample_stat[,2]
total_sample_stat$Sample <- rownames(total_sample_stat)
total_sample_stat$Sample <- str_remove(total_sample_stat$Sample, "./CRISPResso_on_")
rownames(total_sample_stat) <- total_sample_stat$Sample
rep1_map_stat <- total_sample_stat



setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep2-replace-by0829data/")
samples <- list.dirs(recursive = F)
total_sample_stat <- list()
for(sample in samples){
  setwd(sample)
  mapping_rate <- read.table("CRISPResso_mapping_statistics.txt", 
                             sep="\t", header = T)
  total_sample_stat[[sample]] <- c(unlist(mapping_rate[1,c(1,3)]))
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
total_sample_stat <- data.frame(do.call(rbind, total_sample_stat))
total_sample_stat$unmapped <- total_sample_stat[,1] - total_sample_stat[,2]
total_sample_stat$Sample <- rownames(total_sample_stat)
total_sample_stat$Sample <- str_remove(total_sample_stat$Sample, "./CRISPResso_on_")
rownames(total_sample_stat) <- total_sample_stat$Sample
rep2_map_stat <- total_sample_stat


setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep3-replace-by0829data/")
samples <- list.dirs(recursive = F)
total_sample_stat <- list()
for(sample in samples){
  setwd(sample)
  mapping_rate <- read.table("CRISPResso_mapping_statistics.txt", 
                             sep="\t", header = T)
  total_sample_stat[[sample]] <- c(unlist(mapping_rate[1,c(1,3)]))
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
total_sample_stat <- data.frame(do.call(rbind, total_sample_stat))
total_sample_stat$unmapped <- total_sample_stat[,1] - total_sample_stat[,2]
total_sample_stat$Sample <- rownames(total_sample_stat)
total_sample_stat$Sample <- str_remove(total_sample_stat$Sample, "./CRISPResso_on_")
rownames(total_sample_stat) <- total_sample_stat$Sample
rep3_map_stat <- total_sample_stat

load("../../batch1/Result/group_info.rda")
total_info <- data.frame(rbind(cell_info, tissue_info))

rep1_map_stat <- merge(rep1_map_stat, total_info, by.x = "Sample", by.y = "id")
rep1_map_stat$a <- paste0(rep1_map_stat$a, rep1_map_stat$duplicateID)
rep1_map_stat <- rep1_map_stat[,c(5, 1 : 4)]
colnames(rep1_map_stat)[1] <- "ID"
rep1_map_stat$Rep <- "Rep1"

rep2_map_stat <- merge(rep2_map_stat, total_info, by.x = "Sample", by.y = "id")
rep2_map_stat$a <- paste0(rep2_map_stat$a, rep2_map_stat$duplicateID)
rep2_map_stat <- rep2_map_stat[,c(5, 1 : 4)]
colnames(rep2_map_stat)[1] <- "ID"
rep2_map_stat$Rep <- "Rep2"

rep3_map_stat <- merge(rep3_map_stat, total_info, by.x = "Sample", by.y = "id")
rep3_map_stat$a <- paste0(rep3_map_stat$a, rep3_map_stat$duplicateID)
rep3_map_stat <- rep3_map_stat[,c(5, 1 : 4)]
colnames(rep3_map_stat)[1] <- "ID"
rep3_map_stat$Rep <- "Rep3"

total_map_stat <- data.frame(rbind(rep1_map_stat, rep2_map_stat, rep3_map_stat))
total_map_stat <- total_map_stat[,c(1,2,6, 3 : 5)]
total_map_stat <- total_map_stat[order(total_map_stat$ID,total_map_stat$Sample,total_map_stat$Rep),]
total_map_stat$unmapped_ratio <- total_map_stat$unmapped / total_map_stat$READS.IN.INPUTS
openxlsx::write.xlsx(total_map_stat, file="../Previous_data_map_stat.xlsx", colNames=T, rowNames=F)
