setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep1-replace-by0829data/")
samples <- list.dirs(recursive = F)
total_sample_stat <- list()
for(sample in samples){
  setwd(sample)

  mapping_rate <- read.table("CRISPResso_quantification_of_editing_frequency.txt", 
                             sep="\t", header = T)
  total_sample_stat[[sample]] <- c(unlist(mapping_rate[1,c(5 : 11)]))

  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
output_res<- data.frame(do.call(rbind, total_sample_stat))
output_res$Sample <- rownames(output_res)
output_res$Sample <- str_remove(output_res$Sample, "./CRISPResso_on_")
rownames(output_res) <- output_res$Sample
rep1_edit_stat <- output_res
#####Rep2

setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep2-replace-by0829data/")
samples <- list.dirs(recursive = F)
total_sample_stat <- list()
for(sample in samples){
  setwd(sample)
  
  mapping_rate <- read.table("CRISPResso_quantification_of_editing_frequency.txt", 
                             sep="\t", header = T)
  total_sample_stat[[sample]] <- c(unlist(mapping_rate[1,c(5 : 11)]))
  
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
output_res<- data.frame(do.call(rbind, total_sample_stat))
output_res$Sample <- rownames(output_res)
output_res$Sample <- str_remove(output_res$Sample, "./CRISPResso_on_")
rownames(output_res) <- output_res$Sample
rep2_edit_stat <- output_res

####Rep3
setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep3-replace-by0829data/")
samples <- list.dirs(recursive = F)
total_sample_stat <- list()
for(sample in samples){
  setwd(sample)
  
  mapping_rate <- read.table("CRISPResso_quantification_of_editing_frequency.txt", 
                             sep="\t", header = T)
  total_sample_stat[[sample]] <- c(unlist(mapping_rate[1,c(5 : 11)]))
  
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
output_res<- data.frame(do.call(rbind, total_sample_stat))
output_res$Sample <- rownames(output_res)
output_res$Sample <- str_remove(output_res$Sample, "./CRISPResso_on_")
rownames(output_res) <- output_res$Sample
rep3_edit_stat <- output_res


load("../../batch1/Result/group_info.rda")
total_info <- data.frame(rbind(cell_info, tissue_info))

rep1_edit_stat <- merge(rep1_edit_stat, total_info, by.x = "Sample", by.y = "id")
rep1_edit_stat$a <- paste0(rep1_edit_stat$a, rep1_edit_stat$duplicateID)
rep1_edit_stat <- rep1_edit_stat[,c(9, 1 : 8)]
colnames(rep1_edit_stat)[1] <- "ID"
rep1_edit_stat$Rep <- "Rep1"

rep2_edit_stat <- merge(rep2_edit_stat, total_info, by.x = "Sample", by.y = "id")
rep2_edit_stat$a <- paste0(rep2_edit_stat$a, rep2_edit_stat$duplicateID)
rep2_edit_stat <- rep2_edit_stat[,c(9, 1 : 8)]
colnames(rep2_edit_stat)[1] <- "ID"
rep2_edit_stat$Rep <- "Rep2"


rep3_edit_stat <- merge(rep3_edit_stat, total_info, by.x = "Sample", by.y = "id")
rep3_edit_stat$a <- paste0(rep3_edit_stat$a, rep3_edit_stat$duplicateID)
rep3_edit_stat <- rep3_edit_stat[,c(9, 1 : 8)]
colnames(rep3_edit_stat)[1] <- "ID"
rep3_edit_stat$Rep <- "Rep3"


total_edit_stat <- data.frame(rbind(rep1_edit_stat, rep2_edit_stat, rep3_edit_stat))
total_edit_stat <- total_edit_stat[,c(1,2,10, 3 : 9)]
total_edit_stat <- total_edit_stat[order(total_edit_stat$ID,total_edit_stat$Sample,total_edit_stat$Rep),]
openxlsx::write.xlsx(total_edit_stat, file="../Previous_data_edit_stat.xlsx", colNames=T, rowNames=F)
