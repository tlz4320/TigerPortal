#在2.0_total_First_pop_peak_indel1.R
####Read batch1_1 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/Second/Batch1/240318-A00599B/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
sgRNA_tmp <- read.table("../sg_info.txt")
sgRNA_tmp <- sgRNA_tmp[c(19:30,32,33,35,36),]
total_indel_table <- list()
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    setwd(paste0("CRISPResso_on_",sg))
    edit_info <- read.table("Indel_histogram.txt", 
                            sep="\t", header = T)
    total_indel_table[[str_remove(sample, "40")]][[sg]] <- edit_info
    setwd("..")
    
  }
  setwd("..")
}
###Read batch1_2 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/Second/Batch2/240426-A00599B/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
sgRNA_tmp <- read.table("../sg_info.txt")
sgRNA_tmp <- sgRNA_tmp[c(1 : 18, 37 : 160), ]
rm_list <- c("169","170","177","178","207","208","209","210","221","224","225","240","250","256","258","268","279","281","284","285","299","311")
sgRNA_tmp <- sgRNA_tmp[!sgRNA_tmp$V1 %in% paste0("Sg_24_", rm_list),]
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    setwd(paste0("CRISPResso_on_",sg))
    edit_info <- read.table("Indel_histogram.txt", 
                            sep="\t", header = T)
    total_indel_table[[sample]][[sg]] <- edit_info
    setwd("..")
  }
  setwd("..")
}
rm(sample, samples, sg, sgRNA_tmp, edit_info)


tmp <- unique(unlist(lapply(total_indel_table, names)))
total_indel_table_rev <- lapply(tmp, function(x){
  res <- list()
  for(name in names(total_indel_table)){
    if(x %in% names(total_indel_table[[name]])){
      res[[name]] <- total_indel_table[[name]][[x]]
    }
  }
  res
})

names(total_indel_table_rev) <- tmp
total_indel_table_second <- total_indel_table
total_indel_table_second_rev <- total_indel_table_rev
save(total_indel_table_second, total_indel_table_second_rev, 
     file="~/data/project/ear_project/gene_therapy_ll/Result/total_indel_table_second.rda")
