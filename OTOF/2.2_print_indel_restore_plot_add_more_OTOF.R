###接着1.4 + 2.1脚本跑
####合并编辑pattern数据
print(load("~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data.rda"))

names(otof_cell_sg12_15)
otof_cell_sg12_15_rep1 <- otof_cell_sg12_15[seq(1, length(otof_cell_sg12_15), 2)]
otof_cell_sg12_15_rep2 <- otof_cell_sg12_15[seq(2, length(otof_cell_sg12_15), 2)]

otof_cell_sg12_15_ave <- list()
otof_cell_sg12_15_ave_list <- otof_cell_sg12_15_split[seq(1, 8, 2)]
for(i in 1 : length(otof_cell_sg12_15_rep1)){
  tmp_rep1 <- otof_cell_sg12_15_rep1[[i]]
  tmp_rep2 <- otof_cell_sg12_15_rep2[[i]]
  name <- otof_cell_sg12_15_ave_list[i]
  region <- min(c(tmp_rep1[,1], tmp_rep2[,1])) : 
    max(c(tmp_rep1[,1], tmp_rep2[,1]))
  tmp_region <- region[!region %in% tmp_rep1[,1]]
  if(length(tmp_region) != 0){
    tmp_rep1 <- data.frame(rbind(tmp_rep1, data.frame(indel_size = tmp_region, fq = 0)))
  }
  tmp_region <- region[!region %in% tmp_rep2[,1]]
  if(length(tmp_region) != 0){
    tmp_rep2 <- data.frame(rbind(tmp_rep2, data.frame(indel_size = tmp_region, fq = 0)))
  }
  tmp_rep1 <- tmp_rep1[order(tmp_rep1[,1]),]
  tmp_rep2 <- tmp_rep2[order(tmp_rep2[,1]),]
  tmp_rep1[,2] <- (tmp_rep1[,2] + tmp_rep2[,2]) / 2
  otof_cell_sg12_15_ave[[name]] <- tmp_rep1
  rm(tmp_rep1, tmp_rep2, tmp_region, region)
}

otof_cell_sg12_15_framerestore_rep1 <- otof_cell_sg12_15_framerestore[seq(1, 8, 2)]
otof_cell_sg12_15_framerestore_rep2 <- otof_cell_sg12_15_framerestore[seq(2, 8, 2)]

otof_cell_sg12_15_framerestore_ave <- list()
for(i in 1 : length(otof_cell_sg12_15_framerestore_rep1)){
  name <- otof_cell_sg12_15_ave_list[i]
  otof_cell_sg12_15_framerestore_ave[name] <- (otof_cell_sg12_15_framerestore_rep1[i] + 
    otof_cell_sg12_15_framerestore_rep2[i]) / 2
}
otof_cell_sg12_15_framerestore_ave <- unlist(otof_cell_sg12_15_framerestore_ave)


merged_cell_line_result <- old_cell_order
for(name in names(otof_cell_sg12_15_ave)){
  merged_cell_line_result[[name]] <- otof_cell_sg12_15_ave[[name]]
}
merged_cell_line_split <- c(split_name_cell$mutation, otof_cell_sg12_15_ave_list)
merged_cell_line_restore <- c(rep123_framerestore_cell, otof_cell_sg12_15_framerestore_ave)

merged_tissue_result <- old_tissue_order
names(otof_tissue_sg7_14) <- otof_tissue_sg7_14_split
names(otof_tissue_sg7_14_framerestore) <- otof_tissue_sg7_14_split
for(name in names(otof_tissue_sg7_14)){
  merged_tissue_result[[name]] <- otof_tissue_sg7_14[[name]]
}
merged_tissue_split <- c(split_name_tissue$mutation, otof_tissue_sg7_14_split)
merged_tissue_restore <- c(rep123_framerestore_tissue, otof_tissue_sg7_14_framerestore)

maxRestore <- max(c(merged_cell_line_restore, merged_tissue_restore))
setwd("~/data/project/ear_project/gene_therapy_ll/batch1/Result/")

merged_cell_line_split <- str_remove(merged_cell_line_split, "(w|wm|mw|m[0-9]+)[-]")
merged_tissue_split <- str_remove(merged_tissue_split, "(w|wm|mw|m[0-9]+)[-]")

###重新排序
library(gtools)
print(load("group_info.rda"))
cell_info <- cell_info[!duplicated(cell_info$id),]
merged_cell_line_result_reorder <- merged_cell_line_result[cell_info$id]
merged_cell_line_restore_reorder <- merged_cell_line_restore[cell_info$id]
merged_cell_line_split_reorder <- cell_info$a

tissue_info <- tissue_info[!duplicated(tissue_info$id),]
merged_tissue_result_reorder <- merged_tissue_result[tissue_info$id]
merged_tissue_restore_reorder <- merged_tissue_restore[tissue_info$id]
merged_tissue_split_reorder <- tissue_info$a




pdf("indel_stat_heatmap_old_version12_cell.pdf", width = 30, height = 7)
plotStat(merged_cell_line_result_reorder, region = -18 : 5, title = "Cell", showname = T,
         split_col = merged_cell_line_split_reorder,
         shownumber = T, restore = merged_cell_line_restore_reorder, 
         maxRestore = maxRestore, convertGene = function(x){
           x <- unlist(strsplit(x, "[-]"))
           x[1] <- paste0("*", x[1], "*")
           paste0(x, collapse = "-")
         })
dev.off()
pdf("indel_stat_heatmap_old_version12_tissue.pdf", width = 30, height = 7)

plotStat(merged_tissue_result_reorder, region = -18 : 5, title = "Tissue", showname = T,
         split_col = merged_tissue_split_reorder,
         shownumber = T, restore = merged_tissue_restore_reorder, 
         maxRestore = maxRestore, convertGene = function(x){
           x <- unlist(strsplit(x, "[-]"))
           x[1] <- paste0("*", x[1], "*")
           paste0(x, collapse = "-")
         })
dev.off()

####输出回码的结果

load("../../batch1/Result/group_info.rda")
total_info <- data.frame(rbind(cell_info, tissue_info))

rep1_framerestore_output <- data.frame(do.call(rbind, rep1_framerestore))
rep1_framerestore_output$sgRNA <- rownames(rep1_framerestore_output)
rep1_framerestore_output$Rep <- "Rep1"

rep2_framerestore_output <- data.frame(do.call(rbind, rep2_framerestore))
rep2_framerestore_output$sgRNA <- rownames(rep2_framerestore_output)
rep2_framerestore_output$Rep <- "Rep2"

rep3_framerestore_output <- data.frame(do.call(rbind, rep3_framerestore))
rep3_framerestore_output$sgRNA <- rownames(rep3_framerestore_output)
rep3_framerestore_output$Rep <- "Rep3"

rep123_framestore_output <- data.frame(rbind(rep1_framerestore_output, rep2_framerestore_output, rep3_framerestore_output))
colnames(rep123_framestore_output) <- c("Restore Reads", "Unrestore Reads", "Percent", "sample", "Rep")
rep123_framestore_output <- rep123_framestore_output[, c(4,5,1,2,3)]
rep123_framestore_output$sample <- str_remove(rep123_framestore_output$sample, "CRISPResso_on_")
rep123_framestore_output <- merge(total_info, rep123_framestore_output, by.x="id", by.y="sample")
rep123_framestore_output <- rep123_framestore_output[,c(2,1,5 : ncol(rep123_framestore_output))]
colnames(rep123_framestore_output) <- c("ID", "Sample", "Rep", "Restore Reads", "Unrestore Reads", "Percent")
rep123_framestore_output <- rep123_framestore_output[order(rep123_framestore_output$ID, rep123_framestore_output$Sample, rep123_framestore_output$Rep),]
openxlsx::write.xlsx(rep123_framestore_output,file="../Rep123_framerestore_stat.xlsx", rowNames=F, colNames=T)


