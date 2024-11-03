###接着2.1 绘制version6之后，需要增加新的OTOF数据，在该脚本中增加统计和绘图
#load cell line data
setwd("~/data/project/ear_project/gene_therapy_lyh/batch13/clean_data/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
total_insert_stat <- list()
total_restore_stat <- list()
for(sample in samples){
  setwd(sample)
  setwd("CRISPResso_on_nhej")
  mapping_rate <- read.table("Indel_histogram.txt", 
                             sep="\t", header = T)
  filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
  edit_table <- read.table(filename, 
                           sep="\t", header = T, comment.char = "")
  total_restore_stat[[sample]] <- edit_table
  total_insert_stat[[sample]] <- mapping_rate
  setwd("../..")
  rm(mapping_rate)
}
rm(sample, samples)
otof_cell_sg12_15 <- total_insert_stat
otof_cell_sg12_15_restore <- total_restore_stat
rm(total_insert_stat)
otof_cell_sg12_15_split <- paste0(str_remove(names(otof_cell_sg12_15), "[-].*"),"-Otof-1236dC")

#load tissue data2
setwd("~/data/project/ear_project/gene_therapy_qsw/batch2/fastq")

samples <- openxlsx::read.xlsx("OTOF_around_sg_info.xlsx", colNames=F)
total_insert_stat <- list()
total_restore_stat <- list()
for(sample in samples$X1){
  setwd(paste0(sample, "_res"))
  setwd(paste0("CRISPResso_on_", sample))
  mapping_rate <- read.table("Indel_histogram.txt", 
                             sep="\t", header = T)
  filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
  edit_table <- read.table(filename, 
                           sep="\t", header = T, comment.char = "")
  id <- samples[samples[,1] == sample, 2]
  total_restore_stat[[id]] <- edit_table
  total_insert_stat[[id]] <- mapping_rate
  setwd("../..")
  rm(mapping_rate)
}
rm(sample, samples)
otof_tissue_sg12_15 <- total_insert_stat
otof_tissue_sg12_15_resote <- total_restore_stat
rm(total_insert_stat, total_restore_stat)
otof_tissue_sg12_15_split <- paste0(str_remove(str_to_lower(names(otof_tissue_sg12_15))
                                              , "[-_].*"),"-Otof-1236dC")
###统计新的restore比例
otof_cell_sg12_15_edit <- lapply(otof_cell_sg12_15_restore, function(tmp){
  tmp[tmp$Unedited == "False",]
})
otof_tissue_sg12_15_edit <- lapply(otof_tissue_sg12_15_resote, function(tmp){
  tmp[tmp$Unedited == "False",]
})

otof_cell_sg12_15_framerestore <- lapply(names(otof_cell_sg12_15_edit), function(n){
  tmp <- otof_cell_sg12_15_edit[[n]]
  type <- -1
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_reads <- sum(tmp[restore, 7])
  unrestore_reads <- sum(tmp[!restore, 7])
  c(restore_reads, unrestore_reads, restore_reads / (unrestore_reads + restore_reads))
})
names(otof_cell_sg12_15_framerestore) <-  names(otof_cell_sg12_15_edit)
otof_tissue_sg12_15_framerestore <- lapply(names(otof_tissue_sg12_15_edit), function(n){
  tmp <- otof_tissue_sg12_15_edit[[n]]
  type <- -1
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_reads <- sum(tmp[restore, 7])
  unrestore_reads <- sum(tmp[!restore, 7])
  c(restore_reads, unrestore_reads, restore_reads / (unrestore_reads + restore_reads))
})
names(otof_tissue_sg12_15_framerestore) <-  names(otof_tissue_sg12_15_edit)

otof_cell_sg12_15_framerestore <- unlist(lapply(otof_cell_sg12_15_framerestore, function(x){
  x[3]
}))
otof_tissue_sg12_15_framerestore <- unlist(lapply(otof_tissue_sg12_15_framerestore, function(x){
  x[3]
}))





save(otof_cell_sg12_15, otof_tissue_sg12_15,
     otof_cell_sg12_15_edit, otof_cell_sg12_15_restore,
     otof_tissue_sg12_15_edit, otof_tissue_sg12_15_resote,
     otof_cell_sg12_15_split, otof_tissue_sg12_15_split,
     otof_tissue_sg12_15_framerestore, otof_cell_sg12_15_framerestore,
     file= "~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data_newer_newer.rda")


