###接着1.4 + 2.1脚本跑
load("old_rep12_edit.rda")
old_rep1_edited <- lapply(old_rep1_edit, function(tmp){
  tmp[tmp$Unedited == "False",]
})
rm(old_rep1_edit)
old_rep2_edited <- lapply(old_rep2_edit, function(tmp){
  tmp[tmp$Unedited == "False",]
})
rm(old_rep2_edit)

rep1_aa_change <- lapply(names(old_rep1_edited), function(n){
  tmp <- old_rep1_edited[[n]]
  n <- str_remove(n, "CRISPResso_on_")
  type <- str_sub(split_name$mutation[split_name$cellid == n | split_name$tissueid == n], -2, -2)
  type <- ifelse(type == 'd', -1, 1)
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_tmp <- tmp[restore, ]
  aa_change <- abs(type + restore_tmp$n_inserted - restore_tmp$n_deleted) / 3
  aa_change <- lapply(split(restore_tmp[,7], aa_change), sum)
  aa_change <- data.frame(aa = names(aa_change), Counts = unlist(aa_change))
  aa_change
})
names(rep1_aa_change) <-  names(old_rep1_edited)
rep2_aa_change <- lapply(names(old_rep2_edited), function(n){
  tmp <- old_rep2_edited[[n]]
  n <- str_remove(n, "CRISPResso_on_")
  type <- str_sub(split_name$mutation[split_name$cellid == n | split_name$tissueid == n], -2, -2)
  type <- ifelse(type == 'd', -1, 1)
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_tmp <- tmp[restore, ]
  aa_change <- abs(type + restore_tmp$n_inserted - restore_tmp$n_deleted) / 3
  aa_change <- lapply(split(restore_tmp[,7], aa_change), sum)
  aa_change <- data.frame(aa = names(aa_change), Counts = unlist(aa_change))
  aa_change
})
names(rep2_aa_change) <-  names(old_rep2_edited)


load("~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data.rda")
otof_cell_sg12_15_edit <- lapply(otof_cell_sg12_15_restore, function(tmp){
  tmp[tmp$Unedited == "False",]
})
rm(otof_cell_sg12_15_restore)
otof_tissue_sg7_14_edit <- lapply(otof_tissue_sg7_14_resote, function(tmp){
  tmp[tmp$Unedited == "False",]
})
rm(otof_tissue_sg7_14_resote)

otof_cell_sg12_15_aa_change <- lapply(names(otof_cell_sg12_15_edit), function(n){
  tmp <- otof_cell_sg12_15_edit[[n]]
  type <- -1
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_tmp <- tmp[restore, ]
  aa_change <- abs(type + restore_tmp$n_inserted - restore_tmp$n_deleted) / 3
  aa_change <- lapply(split(restore_tmp[,7], aa_change), sum)
  aa_change <- data.frame(aa = names(aa_change), Counts = unlist(aa_change))
  aa_change
})
names(otof_cell_sg12_15_aa_change) <-  names(otof_cell_sg12_15_edit)
otof_tissue_sg7_14_aa_change <- lapply(names(otof_tissue_sg7_14_edit), function(n){
  tmp <- otof_tissue_sg7_14_edit[[n]]
  type <- -1
  tmp <- tmp[tmp$n_inserted != 0 | tmp$n_deleted != 0,]
  restore <- abs(type + tmp$n_inserted - tmp$n_deleted) %% 3 == 0
  restore_tmp <- tmp[restore, ]
  aa_change <- abs(type + restore_tmp$n_inserted - restore_tmp$n_deleted) / 3
  aa_change <- lapply(split(restore_tmp[,7], aa_change), sum)
  aa_change <- data.frame(aa = names(aa_change), Counts = unlist(aa_change))
  aa_change
})
names(otof_tissue_sg7_14_aa_change) <-  names(otof_tissue_sg7_14_edit)

# tmp <- ls(pattern="aa_change")
# save(list = tmp, file="AA_change_stat.rda")

rep12_aa_change <- list()

for(name in names(rep1_aa_change)){
  tmp_rep1 <- rep1_aa_change[[name]]
  if(!name %in% names(rep2_aa_change)){
    rep12_aa_change[[name]] <- tmp_rep1
    next
  }
  tmp_rep2 <- rep2_aa_change[[name]]
  tmp_rep12 <- data.frame(rbind(tmp_rep1, tmp_rep2))
  if(nrow(tmp_rep12) == 0){
    rep12_aa_change[[name]] <- tmp_rep12
    next
  }
  tmp_rep12 <- lapply(split(tmp_rep12$Counts, tmp_rep12$aa), sum)
  tmp_rep12 <- data.frame(aa = names(tmp_rep12), Counts = unlist(tmp_rep12))
  rep12_aa_change[[name]] <- tmp_rep12
}

rm(tmp_rep1, tmp_rep2, tmp_rep12)

names(rep12_aa_change) <- str_remove(names(rep12_aa_change), "CRISPResso_on_")

rep12_aa_change_cell <- rep12_aa_change[split_name_cell$cellid]
rep12_aa_change_tissue <- rep12_aa_change[split_name_tissue$tissueid]



####合并蛋白质长度数据
merged_cell_line_aa_change <- rep12_aa_change_cell
for(name in names(otof_cell_sg12_15_aa_change)){
  merged_cell_line_aa_change[[name]] <- otof_cell_sg12_15_aa_change[[name]]
}
merged_cell_line_id <- c(split_name_cell$cellid, otof_cell_sg12_15_split)
merged_cell_line_split <- c(split_name_cell$mutation, otof_cell_sg12_15_split)
merged_tissue_aa_change <- rep12_aa_change_tissue
for(name in names(otof_tissue_sg7_14_aa_change)){
  merged_tissue_aa_change[[name]] <- otof_tissue_sg7_14_aa_change[[name]]
}
merged_tissue_split <- c(split_name_tissue$mutation, otof_tissue_sg7_14_split)
merged_tissue_id <- c(split_name_tissue$tissueid, otof_tissue_sg7_14_split)
merged_cell_line_split <- str_remove(merged_cell_line_split, "(w|wm|mw|m[0-9]+)[-]")
merged_tissue_split <- str_remove(merged_tissue_split, "(w|wm|mw|m[0-9]+)[-]")

###重新排序
library(gtools)
merged_cell_line_aa_change_reorder <- merged_cell_line_aa_change[mixedorder(merged_cell_line_split)]
merged_cell_line_id_reorder <- merged_cell_line_id[mixedorder(merged_cell_line_split)]
merged_cell_line_split_reorder <- merged_cell_line_split[mixedorder(merged_cell_line_split)]

merged_tissue_aa_change_reorder <- merged_tissue_aa_change[mixedorder(merged_tissue_split)]
merged_tissue_id_reorder <- merged_tissue_id[mixedorder(merged_tissue_split)]
merged_tissue_split_reorder <- merged_tissue_split[mixedorder(merged_tissue_split)]

pdf("aa_change_heatmap_old_version1_cell.pdf", width = 30, height = 7)
plot_data <- plotStat_AA(merged_cell_line_aa_change_reorder, region = 5, title = "Cell", showname = T,
         split_col = merged_cell_line_split_reorder,
         shownumber = T, convertGene = function(x){
           x <- unlist(strsplit(x, "[-]"))
           x[1] <- paste0("*", x[1], "*")
           paste0(x, collapse = "-")
         })
dev.off()

library(dplyr)
df1 <- data.frame(a = merged_cell_line_split_reorder, b = merged_cell_line_split_reorder)
df1 <- df1 %>% 
  dplyr::group_by(a, b) %>%
  dplyr::mutate(duplicateID = dplyr::row_number()) %>%
  dplyr::ungroup()
table(df1$a == merged_cell_line_split_reorder)
plot_data$Pct <- data.frame(plot_data$Pct)
colnames(plot_data$Counts) <- paste0(df1$a, df1$duplicateID)
colnames(plot_data$Pct) <- paste0(df1$a, df1$duplicateID)

plot_data$Counts$Change <- rownames(plot_data$Counts)
plot_data$Pct$Change <- rownames(plot_data$Pct)

final_result <- data.frame(rbind(plot_data$Counts, plot_data$Pct))

final_result$Type <- rep(c("Counts", "Percent"), c(6,6))
final_result <- final_result[,c(ncol(final_result), (ncol(final_result) - 1), 1 : (ncol(final_result) - 2))]
colnames(final_result) <- c("Type", "Change", paste0(df1$a, df1$duplicateID))
final_result[13,] <- colnames(final_result)
final_result[14,] <- c("Type","Change",merged_cell_line_id_reorder)
final_result <- final_result[c(13, 14 ,1:12),]
openxlsx::write.xlsx(final_result, file="cell_AA_change_stat.xlsx", rowNames=F, colNames=F)


pdf("aa_change_heatmap_old_version1_tissue.pdf", width = 30, height = 7)

plot_data <- plotStat_AA(merged_tissue_aa_change_reorder, region = 5, title = "Tissue", showname = R,
         split_col = merged_tissue_split_reorder,
         shownumber = T, convertGene = function(x){
           x <- unlist(strsplit(x, "[-]"))
           x[1] <- paste0("*", x[1], "*")
           paste0(x, collapse = "-")
         })
dev.off()


df2 <- data.frame(a = merged_tissue_split_reorder, b = merged_tissue_split_reorder)
df2 <- df2 %>% 
  dplyr::group_by(a, b) %>%
  dplyr::mutate(duplicateID = dplyr::row_number()) %>%
  dplyr::ungroup()
# cell_info <- df1
# cell_info$id <- merged_cell_line_id_reorder
# tissue_info <- df2
# tissue_info$id <- merged_tissue_id_reorder
# save(cell_info, tissue_info, file="../../batch1/Result/group_info.rda")
table(df2$a == merged_tissue_split_reorder)
plot_data$Pct <- data.frame(plot_data$Pct)
colnames(plot_data$Counts) <- paste0(df2$a, df2$duplicateID)
colnames(plot_data$Pct) <- paste0(df2$a, df2$duplicateID)

plot_data$Counts$Change <- rownames(plot_data$Counts)
plot_data$Pct$Change <- rownames(plot_data$Pct)

final_result <- data.frame(rbind(plot_data$Counts, plot_data$Pct))

final_result$Type <- rep(c("Counts", "Percent"), c(6,6))
final_result <- final_result[,c(ncol(final_result), (ncol(final_result) - 1), 1 : (ncol(final_result) - 2))]
colnames(final_result) <- c("Type", "Change", paste0(df2$a, df2$duplicateID))
final_result[13,] <- colnames(final_result)
final_result[14,] <- c("Type","Change",merged_tissue_id_reorder)
final_result <- final_result[c(13, 14 ,1:12),]
openxlsx::write.xlsx(final_result, file="tissue_AA_change_stat.xlsx", rowNames=F, colNames=F)
