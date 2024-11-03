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

#load tissue data
setwd("~/data/project/ear_project/gene_therapy_wgg/batch5_24_2_1/clean_data/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
samples <- samples[samples != "w17GJB6_res"]
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
otof_tissue_sg7_14 <- total_insert_stat
otof_tissue_sg7_14_resote <- total_restore_stat
rm(total_insert_stat, total_restore_stat)
names(otof_tissue_sg7_14) <- str_remove(names(otof_tissue_sg7_14), "WT-OTOF-")
names(otof_tissue_sg7_14_resote) <- str_remove(names(otof_tissue_sg7_14_resote), "WT-OTOF-")
otof_tissue_sg7_14_split <- paste0(str_remove(str_to_lower(names(otof_tissue_sg7_14))
                                              , "[-_].*"),"-Otof-1236dC")

#load tissue data2
setwd("~/data/project/ear_project/gene_therapy_lyh/batch15/clean_data/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
samples <- samples[grep("OTOF", samples)]
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
otof_tissue_sg7_15 <- total_insert_stat
otof_tissue_sg7_15_resote <- total_restore_stat
rm(total_insert_stat, total_restore_stat)
names(otof_tissue_sg7_15) <- paste0(str_remove(str_remove(names(otof_tissue_sg7_15), 
                                                   "OTOF1236delC-"), "-WGG_res"), "-2")
names(otof_tissue_sg7_15_resote) <- paste0(str_remove(str_remove(names(otof_tissue_sg7_15_resote), 
                                                          "OTOF1236delC-"), "-WGG_res"), "-2")
otof_tissue_sg7_15_split <- paste0(str_remove(str_to_lower(names(otof_tissue_sg7_15))
                                              , "[-_].*"),"-Otof-1236dC")
for(name in names(otof_tissue_sg7_14)){
  otof_tissue_sg7_15[[name]] <- otof_tissue_sg7_14[[name]]
  otof_tissue_sg7_15_resote[[name]] <- otof_tissue_sg7_14_resote[[name]]
}
otof_tissue_sg7_15_split <- paste0(str_remove(str_to_lower(names(otof_tissue_sg7_15))
                                              , "[-_].*"),"-Otof-1236dC")
###统计新的restore比例
otof_cell_sg12_15_edit <- lapply(otof_cell_sg12_15_restore, function(tmp){
  tmp[tmp$Unedited == "False",]
})
otof_tissue_sg7_15_edit <- lapply(otof_tissue_sg7_14_resote, function(tmp){
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
otof_tissue_sg7_15_framerestore <- lapply(names(otof_tissue_sg7_15_edit), function(n){
  tmp <- otof_tissue_sg7_15_edit[[n]]
  type <- -1
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_reads <- sum(tmp[restore, 7])
  unrestore_reads <- sum(tmp[!restore, 7])
  c(restore_reads, unrestore_reads, restore_reads / (unrestore_reads + restore_reads))
})
names(otof_tissue_sg7_15_framerestore) <-  names(otof_tissue_sg7_15_edit)

otof_cell_sg12_15_framerestore <- unlist(lapply(otof_cell_sg12_15_framerestore, function(x){
  x[3]
}))
otof_tissue_sg7_15_framerestore <- unlist(lapply(otof_tissue_sg7_15_framerestore, function(x){
  x[3]
}))





# save(otof_cell_sg12_15, otof_tissue_sg7_15,
#      otof_cell_sg12_15_edit, otof_cell_sg12_15_restore,
#      otof_tissue_sg7_15_edit, otof_tissue_sg7_15_resote,
#      otof_cell_sg12_15_split, otof_tissue_sg7_15_split,
#      otof_tissue_sg7_15_framerestore, otof_cell_sg12_15_framerestore,
#      file= "~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data_newer.rda")


###为每个样本单独画图
region <- c(-19: 19)

cell_line_data <- lapply(otof_cell_sg12_15, function(y){
  sel_sg <- y
  if(is.null(sel_sg))
    return(sel_sg)
  long_del <- sum(sel_sg[sel_sg[,1] < min(region), 2])
  long_ins <- sum(sel_sg[sel_sg[,1] > max(region), 2])
  long_res <- data.frame(indel_size = c("<=-20", ">=20"), fq = c(long_del, long_ins))
  tmp_region <- region[!region %in% sel_sg[,1]]
  sel_sg <- sel_sg[sel_sg[,1] %in% region,]
  if(length(tmp_region) != 0){
    sel_sg <- data.frame(rbind(sel_sg, data.frame(indel_size = tmp_region, fq = 0)))
    sel_sg <- sel_sg[order(sel_sg[,1]),]
  }
  sel_sg <- data.frame(rbind(sel_sg, long_res))
  sel_sg[sel_sg[,1] == 0, 2] <- 0 
  sel_sg
})

tissue_data <- lapply(otof_tissue_sg7_14, function(y){
  sel_sg <- y
  if(is.null(sel_sg))
    return(sel_sg)
  long_del <- sum(sel_sg[sel_sg[,1] < min(region), 2])
  long_ins <- sum(sel_sg[sel_sg[,1] > max(region), 2])
  long_res <- data.frame(indel_size = c("<=-20", ">=20"), fq = c(long_del, long_ins))
  tmp_region <- region[!region %in% sel_sg[,1]]
  sel_sg <- sel_sg[sel_sg[,1] %in% region,]
  if(length(tmp_region) != 0){
    sel_sg <- data.frame(rbind(sel_sg, data.frame(indel_size = tmp_region, fq = 0)))
    sel_sg <- sel_sg[order(sel_sg[,1]),]
  }
  sel_sg <- data.frame(rbind(sel_sg, long_res))
  sel_sg[sel_sg[,1] == 0, 2] <- 0 
  sel_sg
})


pdf("../../batch1/Result/cell_line_pattern_each_sample.pdf", width = 20, height = 8)
for(name in names(cell_line_data)){
  tmp_data <- cell_line_data[[name]]
  max_pct <- max(tmp_data[,2] / sum(tmp_data[,2]))
  x <- tmp_data
  x$indel_size <- factor(x$indel_size, levels = c("<=-20", -19 : 19, ">=20"))
  x$pct <- x[,2] / sum(x[,2])
  print(ggplot(x, aes(x = indel_size, y = pct * 100)) + geom_bar(stat="identity") + 
          scale_y_continuous(limits = c(0, max_pct * 100)) + 
          ylab("% in indel reads") + xlab("") + 
          theme_bw() + ggtitle(name))
}
dev.off()


pdf("../../batch1/Result/tissue_pattern_each_sample.pdf", width = 20, height = 8)
for(name in names(tissue_data)){
  tmp_data <- tissue_data[[name]]
  max_pct <- max(tmp_data[,2] / sum(tmp_data[,2]))
  x <- tmp_data
  x$indel_size <- factor(x$indel_size, levels = c("<=-20", -19 : 19, ">=20"))
  x$pct <- x[,2] / sum(x[,2])
  print(ggplot(x, aes(x = indel_size, y = pct * 100)) + geom_bar(stat="identity") + 
          scale_y_continuous(limits = c(0, max_pct * 100)) + 
          ylab("% in indel reads") + xlab("") + 
          theme_bw() + ggtitle(name))
}
dev.off()
