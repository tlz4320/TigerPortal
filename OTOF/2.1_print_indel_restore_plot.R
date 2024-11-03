setwd("~/data/project/ear_project/gene_therapy_ll/batch1/Result/")
sample_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)

load("indel_stat_total.rda")
names(old_tissue) <- str_remove(names(old_tissue), "CRISPResso_on_")
names(old_cell) <- str_remove(names(old_cell), "CRISPResso_on_")
old_cell_order <- old_cell[names(old_cell) %in% sample_info$`Cell-w-W`]
split_name <- data.frame(cellid = sample_info$`Cell-w-W`, 
                         tissueid = sample_info$`Tissue-w-W(Mix)`,
                         mutation = sample_info$GRNAa)
split_name$mutation <- unlist(lapply(split_name$mutation, function(x){
  x <- unlist(strsplit(x, "[-]"))
  paste(x[-length(x)], collapse = "-")
}))

split_name_cell <- split_name[split_name$cellid %in% names(old_cell),]
split_name_cell <- na.omit(split_name_cell)
old_cell_order <- old_cell[split_name_cell$cellid]


split_name_tissue <- split_name[split_name$tissueid %in% names(old_tissue),]
split_name_tissue <- na.omit(split_name_tissue)
old_tissue_order <- old_tissue[split_name_tissue$tissueid]

load("old_rep123_edit.rda")
old_rep1_edited <- lapply(old_rep1_edit, function(tmp){
  tmp[tmp$Unedited == "False",]
})
rm(old_rep1_edit)
old_rep2_edited <- lapply(old_rep2_edit, function(tmp){
  tmp[tmp$Unedited == "False",]
})
rm(old_rep2_edit)
old_rep3_edited <- lapply(old_rep3_edit, function(tmp){
  tmp[tmp$Unedited == "False",]
})
rm(old_rep3_edit)

rep1_framerestore <- lapply(names(old_rep1_edited), function(n){
  tmp <- old_rep1_edited[[n]]
  n <- str_remove(n, "CRISPResso_on_")
  type <- str_sub(split_name$mutation[split_name$cellid == n | split_name$tissueid == n], -2, -2)
  type <- ifelse(type == 'd', -1, 1)
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_reads <- sum(tmp[restore, 7])
  unrestore_reads <- sum(tmp[!restore, 7])
  c(restore_reads, unrestore_reads, restore_reads / (unrestore_reads + restore_reads))
})
names(rep1_framerestore) <-  names(old_rep1_edited)
rep2_framerestore <- lapply(names(old_rep2_edited), function(n){
  tmp <- old_rep2_edited[[n]]
  n <- str_remove(n, "CRISPResso_on_")
  type <- str_sub(split_name$mutation[split_name$cellid == n | split_name$tissueid == n], -2, -2)
  type <- ifelse(type == 'd', -1, 1)
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_reads <- sum(tmp[restore, 7])
  unrestore_reads <- sum(tmp[!restore, 7])
  c(restore_reads, unrestore_reads, restore_reads / (unrestore_reads + restore_reads))
})
names(rep2_framerestore) <-  names(old_rep2_edited)
rep3_framerestore <- lapply(names(old_rep3_edited), function(n){
  tmp <- old_rep3_edited[[n]]
  n <- str_remove(n, "CRISPResso_on_")
  type <- str_sub(split_name$mutation[split_name$cellid == n | split_name$tissueid == n], -2, -2)
  type <- ifelse(type == 'd', -1, 1)
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_reads <- sum(tmp[restore, 7])
  unrestore_reads <- sum(tmp[!restore, 7])
  c(restore_reads, unrestore_reads, restore_reads / (unrestore_reads + restore_reads))
})
names(rep3_framerestore) <-  names(old_rep3_edited)



rep123_framerestore <- list()

for(name in names(rep1_framerestore)){
  tmp_rep1 <- rep1_framerestore[[name]]
  if(!name %in% names(rep2_framerestore) & !name %in% names(rep3_framerestore)){
    rep123_framerestore[[name]] <- tmp_rep1
    next
  }
  restore_reads <- tmp_rep1[1]
  unrestore_reads <- tmp_rep1[2]
  has_rep2 <- name %in% names(rep2_framerestore)
  if(has_rep2){
    tmp_rep2 <- rep2_framerestore[[name]]
    restore_reads <- restore_reads + tmp_rep2[1]
    unrestore_reads <- unrestore_reads + tmp_rep2[2]
  }
  has_rep3 <- name %in% names(rep3_framerestore)
  if(has_rep3){
    tmp_rep3 <- rep3_framerestore[[name]]
    restore_reads <- restore_reads + tmp_rep3[1]
    unrestore_reads <- unrestore_reads + tmp_rep3[2]
  }
  
  rep123_framerestore[[name]] <- c(restore_reads, unrestore_reads, restore_reads / (restore_reads + unrestore_reads))
}
rm(tmp_rep1, tmp_rep2, tmp_rep3)
rep123_framerestore_list <- unlist(lapply(rep123_framerestore, function(x){
  x[3]
}))
names(rep123_framerestore_list) <- str_remove(names(rep123_framerestore_list), "CRISPResso_on_")


split_name_cell$mutation <- unlist(lapply(split_name_cell$mutation, function(x){
  x <- unlist(str_split(x, "[-]"))
  x[2] <- stringr::str_to_title(x[2])
  paste(x, collapse = "-")
}))
split_name_tissue$mutation <- unlist(lapply(split_name_tissue$mutation, function(x){
  x <- unlist(str_split(x, "[-]"))
  x[2] <- stringr::str_to_title(x[2])
  paste(x, collapse = "-")
}))

rep123_framerestore_cell <- rep123_framerestore_list[split_name_cell$cellid]
rep123_framerestore_tissue <- rep123_framerestore_list[split_name_tissue$tissueid]
maxRestore <- max(unlist(c(rep123_framerestore_cell, rep123_framerestore_tissue)))
pdf("indel_stat_heatmap_old_version9.pdf", width = 30, height = 7)
plotStat(old_cell_order, region = -18 : 5, title = "Cell", showname = T,
         split_col = split_name_cell$mutation,
         shownumber = T, restore = rep123_framerestore_cell, 
         maxRestore = maxRestore)
plotStat(old_tissue_order, region = -18 : 5, title = "Tissue", showname = T,
         split_col = split_name_tissue$mutation,
         shownumber = T, restore = rep123_framerestore_tissue, 
         maxRestore = maxRestore)
dev.off()





library(ggpubr)

pdf("indel_stat_point_old_2_20_cell.pdf", width = 16, height = 5)
plotStat2(rep12_framerestore_cell, split_name_cell$mutation)
dev.off()
pdf("indel_stat_point_old_2_20_tissue.pdf", width = 16, height = 5)

plotStat2(rep12_framerestore_tissue, split_name_tissue$mutation)
dev.off()

cell_restore_data <- data.frame(restore = rep12_framerestore_cell, 
                               sgRNA = split_name_cell$mutation,
                               name = names(rep12_framerestore_cell))
tissue_restore_data <- data.frame(restore = rep12_framerestore_tissue, 
                                sgRNA = split_name_tissue$mutation,
                                name = names(rep12_framerestore_tissue))
openxlsx::write.xlsx(list("Tissue" = tissue_restore_data, "Cell" = cell_restore_data), 
                     file="Restore_stat_result.xlsx", rowNames=F, colNames=T)
