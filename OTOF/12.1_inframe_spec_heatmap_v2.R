###load 之前114个sgRNA的结果
library(stringr)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
load(file="~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda")
sp_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)
sp_info <- split_name <- data.frame(
  tissueid = c(sp_info$`Tissue-w-W(Mix)`,
               paste0("m", 12:15, "-tissue")),
  mutation = c(sp_info$GRNAa, paste0("OTOF-1236dC-g", 3:6)))
rownames(sp_info) <- sp_info$tissueid
inframe_result_processed$aa <- inframe_result_processed$indel / 3
inframe_result_processed$aa[abs(inframe_result_processed$aa) > 5] <- "Other"

cell_inframe_result_processed$aa <- cell_inframe_result_processed$indel / 3
cell_inframe_result_processed$aa[abs(cell_inframe_result_processed$aa) > 5] <- "Other"

inframe_tissue_stat <- read.xlsx("~/Nutstore Files/Tobin/Previous/previous_tissue_inframe_in_indel.xlsx")
inframe_cell_stat <- read.xlsx("~/Nutstore Files/Tobin/Previous/previous_cell_inframe_in_indel.xlsx")

rownames(inframe_tissue_stat) <- str_remove(inframe_tissue_stat$result, "[-].*")
rownames(inframe_cell_stat) <- str_remove(inframe_cell_stat$result, "[-].*")

tmp <- cell_inframe_result_processed[cell_inframe_result_processed$aa == "0",]
rownames(tmp) <- tmp$id
tmp <- tmp[inframe_cell_stat$result,]
inframe_cell_stat$noAA <- tmp$pct
inframe_cell_stat$noAAse <- tmp$se


tmp <- inframe_result_processed[inframe_result_processed$aa == "0",]
rownames(tmp) <- tmp$id
tmp <- tmp[inframe_tissue_stat$result,]
inframe_tissue_stat$noAA <- tmp$pct
inframe_tissue_stat$noAAse <- tmp$se
inframe_tissue_stat <- inframe_tissue_stat[order(inframe_tissue_stat$inframePct, decreasing = T),]
inframe_cell_stat <- inframe_cell_stat[rownames(inframe_tissue_stat),]
sgname <- sp_info[inframe_tissue_stat$result,]$mutation

paired_id <- intersect(unique(str_remove(inframe_result_processed$id, "[-].*")),
                       unique(str_remove(cell_inframe_result_processed$id, "[-].*")))

inframe_cor <- lapply(paired_id, function(id){
  cell_inframe <- cell_inframe_result_processed[str_remove(cell_inframe_result_processed$id, "[-].*") == id,]
  tissue_inframe <- inframe_result_processed[str_remove(inframe_result_processed$id, "[-].*") == id,]
  
  cell_inframe_array <- rep(0, 12)
  names(cell_inframe_array) <- as.character(c(-5: 5, "Other"))
  tissue_inframe_array <- rep(0, 12)
  names(tissue_inframe_array) <- as.character(c(-5: 5, "Other"))
  for(i in 1 : nrow(cell_inframe)){
    cell_inframe_array[cell_inframe$aa[i]]<-cell_inframe$pct[i]
  }
  for(i in 1 : nrow(tissue_inframe)){
    tissue_inframe_array[tissue_inframe$aa[i]]<-tissue_inframe$pct[i]
  }
  cor(cell_inframe_array, tissue_inframe_array)
})
inframe_cor <- data.frame(id = paired_id, cor = unlist(inframe_cor))
inframe_cor[is.na(inframe_cor)] <- 0
rownames(inframe_cor) <- inframe_cor$id
inframe_cor <- inframe_cor[rownames(inframe_tissue_stat),]

inframe_data <- matrix(0, ncol=12 * 2 + 1, nrow = length(paired_id))
for(index in 1 : nrow(inframe_cor)){
  id <- inframe_cor$id[index]
  cell_inframe <- cell_inframe_result_processed[str_remove(cell_inframe_result_processed$id, "[-].*") == id,]
  tissue_inframe <- inframe_result_processed[str_remove(inframe_result_processed$id, "[-].*") == id,]
  
  cell_inframe_array <- rep(0, 12)
  names(cell_inframe_array) <- as.character(c(5: -5, "Other"))
  tissue_inframe_array <- rep(0, 12)
  names(tissue_inframe_array) <- as.character(c("Other",-5: 5))
  for(i in 1 : nrow(cell_inframe)){
    cell_inframe_array[cell_inframe$aa[i]]<-cell_inframe$pct[i]
  }
  for(i in 1 : nrow(tissue_inframe)){
    tissue_inframe_array[tissue_inframe$aa[i]]<-tissue_inframe$pct[i]
  }
  inframe_data[index,] <- c(tissue_inframe_array, inframe_cor$cor[index], cell_inframe_array)
}
inframe_data[is.na(inframe_data)] <- 0


####分开画  先画delete突变
which_delete <- which(inframe_tissue_stat$mutType == "d")
inframe_tissue_stat_sel <- inframe_tissue_stat[which_delete,]
inframe_cell_stat_sel <- inframe_cell_stat[which_delete,]
inframe_data_sel <- inframe_data[which_delete,]
inframe_cor_sel <- inframe_cor[which_delete,]
sgname_sel <- sgname[which_delete]

max_0aa <- ceiling(max(c(max(inframe_tissue_stat_sel$noAA + 
                               inframe_tissue_stat_sel$noAAse, na.rm = T),
                         max(inframe_cell_stat_sel$noAA + 
                               inframe_cell_stat_sel$noAAse, na.rm = T))))
max_inframe <- max(c(max(inframe_tissue_stat_sel$inframePct + 
                           inframe_tissue_stat_sel$inframeSe, na.rm = T), 
                     max(
                       inframe_cell_stat_sel$inframePct + 
                         inframe_cell_stat_sel$inframeSe, na.rm = T)))
anno_data1 <- matrix(0, ncol = 2, nrow = nrow(inframe_tissue_stat_sel))
anno_data1[, 1]<- inframe_tissue_stat_sel$noAA
anno_data1[, 2] <- inframe_tissue_stat_sel$inframePct
anno_data1[,2] <- anno_data1[,2] - anno_data1[,1]
left_anno <- rowAnnotation(
  sgname1 = anno_empty(width = unit(4, "cm"), border = F),
  tissue_0aa = anno_barplot(anno_data1, 
                            border = T,
                            gp = gpar(fill = c("#494490","#FF8C00")), 
                            bar_width = 0.5,
                            axis_param = list(direction = "reverse", 
                                              "side"="top"),
                            ylim = c(0, max_0aa), 
                            width = unit(2, "cm")),
  show_annotation_name = F)
anno_data2 <- matrix(0, ncol = 2, nrow = nrow(inframe_cell_stat_sel))
anno_data2[, 1]<- inframe_cell_stat_sel$noAA
anno_data2[, 2] <- inframe_cell_stat_sel$inframePct
anno_data2[,2] <- anno_data2[,2] - anno_data2[,1]
right_anno <- rowAnnotation(
  cell_0aa = anno_barplot(anno_data2, 
                          border = T,
                          gp = gpar(fill = c("#494490","#FF8C00")), 
                          bar_width = 0.5,
                          axis_param = list("side"="top"),
                          ylim = c(0, max_0aa), 
                          width = unit(2, "cm")),
  
  show_annotation_name = F)
col_acc <- colSums(inframe_data_sel)
col_acc[13] <- NA
col_mean <- colMeans(inframe_data_sel)
col_se <- apply(inframe_data_sel, 2, function(x){
  plotrix::std.error(x)
})


top_anno <- HeatmapAnnotation(acc = anno_barplot(col_mean[1:12], 
                                                 border = T,
                                                 gp = gpar(fill = "#494490"), 
                                                 bar_width = 0.5,
                                                 width = unit(2, "cm"), 
                                                 ylim = c(0, max(col_mean + col_se, na.rm = T))),
                              pos = anno_empty(border = F, width = unit(2, "mm")))
plot_mat_left <- inframe_data_sel[,c(1:12)]
indel_color <- colorRamp2(c(0, 100), c("white", "blue"))
ht <- Heatmap(plot_mat_left, 
              name = "heatmap1",
              cluster_rows = F,
              cluster_columns = F, 
              top_annotation = top_anno,
              col = c("white", "white"), 
              left_annotation = left_anno,
              cell_fun = function(j, i, x, y, w, h, fill){
                v <- pindex(plot_mat_left, i, j)
                grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                                col = "#AAAAAA"))
              }, show_heatmap_legend = F)

plot_mat_right <- inframe_data_sel[,c(14 : 25)]
top_anno2 <- HeatmapAnnotation(acc2 = anno_barplot(col_mean[14: 25], 
                                                  border = T,
                                                  gp = gpar(fill = "#494490"), 
                                                  bar_width = 0.5,
                                                  width = unit(2, "cm"), 
                                                  ylim = c(0, max(col_mean + col_se, na.rm = T))),
                               pos2 = anno_empty(border = F, width = unit(2, "mm")))
left_anno2 <- rowAnnotation(
  Cor = anno_barplot(inframe_cor_sel$cor, 
                     border = T,
                     gp = gpar(fill = "red"), 
                     bar_width = 0.5,
                     axis_param = list(direction = "reverse", 
                                       "side"="top"), 
                     width = unit(2, "cm")), 
  show_annotation_name = T)
ht2 <- Heatmap(plot_mat_right, 
               name = "heatmap2",
               cluster_rows = F,
               cluster_columns = F, 
               top_annotation = top_anno2,
               col = c("white", "white"), 
               left_annotation = left_anno2,
               right_annotation =right_anno,
               cell_fun = function(j, i, x, y, w, h, fill){
                 v <- pindex(plot_mat_right, i, j)
                 
                 grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                                 col = "#AAAAAA"))
                 
               }, show_heatmap_legend = F)



pdf("~/Nutstore Files/Tobin/Previous/118_pair_inframe_table_del1mut_0aa_order_v3.pdf", 
    width = 12, height = 33.5)

draw(ht + ht2)


decorate_annotation("acc", slice = 1, {
  od = 1 : ncol(plot_mat_left)
  grid.segments(seq_along(od), 
                col_mean[od],
                seq_along(od), 
                col_mean[od] + col_se[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(seq_along(od) - 0.2, 
                col_mean[od] + col_se[od],
                seq_along(od) + 0.2, 
                col_mean[od] + col_se[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
})
decorate_annotation("acc2", slice = 1, {
  od = 1 : ncol(plot_mat_right)
  grid.segments(seq_along(od), col_mean[13 + od],
                seq_along(od), 
                col_mean[13 + od] + col_se[13 + od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(seq_along(od) - 0.2, 
                col_mean[13 + od] + col_se[13 + od],
                seq_along(od) + 0.2, 
                col_mean[13 + od] + col_se[13 + od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
})

decorate_annotation("tissue_0aa", slice = 1, {
  od = nrow(inframe_data_sel) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - inframe_tissue_stat_sel$noAA[od],seq_along(od), 
                xcale - (inframe_tissue_stat_sel$noAA[od] + inframe_tissue_stat_sel$noAAse[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (inframe_tissue_stat_sel$noAA[od] + inframe_tissue_stat_sel$noAAse[od]), 
                seq_along(od) - 0.2, 
                xcale - (inframe_tissue_stat_sel$noAA[od] + inframe_tissue_stat_sel$noAAse[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - inframe_tissue_stat_sel$inframePct[od],seq_along(od), 
                xcale - (inframe_tissue_stat_sel$inframePct[od] + inframe_tissue_stat_sel$inframeSe[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#FF8C00"))
  grid.segments(xcale - (inframe_tissue_stat_sel$inframePct[od] + inframe_tissue_stat_sel$inframeSe[od]), 
                seq_along(od) - 0.2, 
                xcale - (inframe_tissue_stat_sel$inframePct[od] + inframe_tissue_stat_sel$inframeSe[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#FF8C00"))
  
})

indel_text <- c("Other", -5 : 0, paste0("+", 1 : 5))
decorate_annotation("pos",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(indel_text)), 
                 rot = 0, 
                 x = unit(c(1 : length(indel_text)) / length(indel_text) -
                            0.5 / length(indel_text), "npc"),
                 y=unit(0.5, "npc"), hjust = 0.5, 
                 gp = gpar(fontsize = 14, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("pos2",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(rev(indel_text))), 
                 rot = 0, 
                 x = unit(c(1 : length(indel_text)) / length(indel_text) -
                            0.5 / length(indel_text), "npc"),
                 y=unit(0.5, "npc"), hjust = 0.5, 
                 gp = gpar(fontsize = 14, col = "black"))
  grid.draw(tg)
  invisible(tg)
})

decorate_annotation("cell_0aa", slice = 1, {
  od = nrow(inframe_data_sel) : 1
  grid.segments(inframe_cell_stat_sel$noAA[od],seq_along(od), 
                inframe_cell_stat_sel$noAA[od] + inframe_cell_stat_sel$noAAse[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(inframe_cell_stat_sel$noAA[od] + inframe_cell_stat_sel$noAAse[od], 
                seq_along(od) - 0.2, 
                inframe_cell_stat_sel$noAA[od] + inframe_cell_stat_sel$noAAse[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(inframe_cell_stat_sel$inframePc[od],seq_along(od), 
                inframe_cell_stat_sel$inframePc[od] + inframe_cell_stat_sel$inframeSe[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#FF8C00"))
  grid.segments(inframe_cell_stat_sel$inframePc[od] + inframe_cell_stat_sel$inframeSe[od], 
                seq_along(od) - 0.2, 
                inframe_cell_stat_sel$inframePc[od] + inframe_cell_stat_sel$inframeSe[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#FF8C00"))
})
decorate_annotation("sgname1",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(rev(sgname_sel))), 
                 rot = 0, 
                 x = unit(0, "npc"),
                 y=unit(c(1 : length(sgname_sel)) / length(sgname_sel) -
                          0.5 / length(sgname_sel), "npc"), hjust = 0, 
                 gp = gpar(fontsize = 10, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
dev.off()

####分开画  再画insert突变
which_insert <- which(inframe_tissue_stat$mutType == "i")
inframe_tissue_stat_sel <- inframe_tissue_stat[which_insert,]
inframe_cell_stat_sel <- inframe_cell_stat[which_insert,]
inframe_data_sel <- inframe_data[which_insert,]
inframe_cor_sel <- inframe_cor[which_insert,]
sgname_sel <- sgname[which_insert]


max_0aa <- ceiling(max(c(max(inframe_tissue_stat_sel$noAA + 
                               inframe_tissue_stat_sel$noAAse, na.rm = T),
                         max(inframe_cell_stat_sel$noAA + 
                               inframe_cell_stat_sel$noAAse, na.rm = T))))
max_inframe <- max(c(max(inframe_tissue_stat_sel$inframePct + 
                           inframe_tissue_stat_sel$inframeSe, na.rm = T), 
                     max(
                       inframe_cell_stat_sel$inframePct + 
                         inframe_cell_stat_sel$inframeSe, na.rm = T)))
anno_data1 <- matrix(0, ncol = 2, nrow = nrow(inframe_tissue_stat_sel))
anno_data1[, 1]<- inframe_tissue_stat_sel$noAA
anno_data1[, 2] <- inframe_tissue_stat_sel$inframePct
anno_data1[,2] <- anno_data1[,2] - anno_data1[,1]
left_anno <- rowAnnotation(
  sgname1 = anno_empty(width = unit(4, "cm"), border = F),
  tissue_0aa = anno_barplot(anno_data1, 
                            border = T,
                            gp = gpar(fill = c("#494490","#FF8C00")), 
                            bar_width = 0.5,
                            axis_param = list(direction = "reverse", 
                                              "side"="top"),
                            ylim = c(0, max_0aa), 
                            width = unit(2, "cm")),
  show_annotation_name = F)
anno_data2 <- matrix(0, ncol = 2, nrow = nrow(inframe_cell_stat_sel))
anno_data2[, 1]<- inframe_cell_stat_sel$noAA
anno_data2[, 2] <- inframe_cell_stat_sel$inframePct
anno_data2[,2] <- anno_data2[,2] - anno_data2[,1]
right_anno <- rowAnnotation(
  cell_0aa = anno_barplot(anno_data2, 
                          border = T,
                          gp = gpar(fill = c("#494490","#FF8C00")), 
                          bar_width = 0.5,
                          axis_param = list("side"="top"),
                          ylim = c(0, max_0aa), 
                          width = unit(2, "cm")),
  
  show_annotation_name = F)
col_acc <- colSums(inframe_data_sel)
col_acc[13] <- NA
col_mean <- colMeans(inframe_data_sel)
col_se <- apply(inframe_data_sel, 2, function(x){
  plotrix::std.error(x)
})


top_anno <- HeatmapAnnotation(acc = anno_barplot(col_mean[1:12], 
                                                 border = T,
                                                 gp = gpar(fill = "#494490"), 
                                                 bar_width = 0.5,
                                                 width = unit(2, "cm"), 
                                                 ylim = c(0, max(col_mean + col_se, na.rm = T))),
                              pos = anno_empty(border = F, width = unit(2, "mm")))
plot_mat_left <- inframe_data_sel[,c(1:12)]
indel_color <- colorRamp2(c(0, 100), c("white", "blue"))
ht <- Heatmap(plot_mat_left, 
              name = "heatmap1",
              cluster_rows = F,
              cluster_columns = F, 
              top_annotation = top_anno,
              col = c("white", "white"), 
              left_annotation = left_anno,
              cell_fun = function(j, i, x, y, w, h, fill){
                v <- pindex(plot_mat_left, i, j)
                grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                                col = "#AAAAAA"))
              }, show_heatmap_legend = F)

plot_mat_right <- inframe_data_sel[,c(14 : 25)]
top_anno2 <- HeatmapAnnotation(acc2 = anno_barplot(col_mean[14: 25], 
                                                  border = T,
                                                  gp = gpar(fill = "#494490"), 
                                                  bar_width = 0.5,
                                                  width = unit(2, "cm"), 
                                                  ylim = c(0, max(col_mean + col_se, na.rm = T))),
                               pos2 = anno_empty(border = F, width = unit(2, "mm")))
left_anno2 <- rowAnnotation(
  Cor = anno_barplot(inframe_cor_sel$cor, 
                     border = T,
                     gp = gpar(fill = "red"), 
                     bar_width = 0.5,
                     axis_param = list(direction = "reverse", 
                                       "side"="top"), 
                     width = unit(2, "cm")), 
  show_annotation_name = T)
ht2 <- Heatmap(plot_mat_right, 
               name = "heatmap2",
               cluster_rows = F,
               cluster_columns = F, 
               top_annotation = top_anno2,
               col = c("white", "white"), 
               left_annotation = left_anno2,
               right_annotation =right_anno,
               cell_fun = function(j, i, x, y, w, h, fill){
                 v <- pindex(plot_mat_right, i, j)
                 
                 grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                                 col = "#AAAAAA"))
                 
               }, show_heatmap_legend = F)



pdf("~/Nutstore Files/Tobin/Previous/118_pair_inframe_table_ins1mut_0aa_order_v3.pdf", 
    width = 12, height = 11.5)

draw(ht + ht2)

decorate_annotation("acc", slice = 1, {
  od = 1 : ncol(plot_mat_left)
  grid.segments(seq_along(od), 
                col_mean[od],
                seq_along(od), 
                col_mean[od] + col_se[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(seq_along(od) - 0.2, 
                col_mean[od] + col_se[od],
                seq_along(od) + 0.2, 
                col_mean[od] + col_se[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
})
decorate_annotation("acc2", slice = 1, {
  od = 1 : ncol(plot_mat_right)
  grid.segments(seq_along(od), col_mean[13 + od],
                seq_along(od), 
                col_mean[13 + od] + col_se[13 + od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(seq_along(od) - 0.2, 
                col_mean[13 + od] + col_se[13 + od],
                seq_along(od) + 0.2, 
                col_mean[13 + od] + col_se[13 + od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
})


decorate_annotation("tissue_0aa", slice = 1, {
  od = nrow(inframe_data_sel) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - inframe_tissue_stat_sel$noAA[od],seq_along(od), 
                xcale - (inframe_tissue_stat_sel$noAA[od] + inframe_tissue_stat_sel$noAAse[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (inframe_tissue_stat_sel$noAA[od] + inframe_tissue_stat_sel$noAAse[od]), 
                seq_along(od) - 0.2, 
                xcale - (inframe_tissue_stat_sel$noAA[od] + inframe_tissue_stat_sel$noAAse[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - inframe_tissue_stat_sel$inframePct[od],seq_along(od), 
                xcale - (inframe_tissue_stat_sel$inframePct[od] + inframe_tissue_stat_sel$inframeSe[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#FF8C00"))
  grid.segments(xcale - (inframe_tissue_stat_sel$inframePct[od] + inframe_tissue_stat_sel$inframeSe[od]), 
                seq_along(od) - 0.2, 
                xcale - (inframe_tissue_stat_sel$inframePct[od] + inframe_tissue_stat_sel$inframeSe[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#FF8C00"))
  
})

indel_text <- c("Other", -5 : 0, paste0("+", 1 : 5))
decorate_annotation("pos",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(indel_text)), 
                 rot = 0, 
                 x = unit(c(1 : length(indel_text)) / length(indel_text) -
                            0.5 / length(indel_text), "npc"),
                 y=unit(0.5, "npc"), hjust = 0.5, 
                 gp = gpar(fontsize = 14, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("pos2",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(rev(indel_text))), 
                 rot = 0, 
                 x = unit(c(1 : length(indel_text)) / length(indel_text) -
                            0.5 / length(indel_text), "npc"),
                 y=unit(0.5, "npc"), hjust = 0.5, 
                 gp = gpar(fontsize = 14, col = "black"))
  grid.draw(tg)
  invisible(tg)
})

decorate_annotation("cell_0aa", slice = 1, {
  od = nrow(inframe_data_sel) : 1
  grid.segments(inframe_cell_stat_sel$noAA[od],seq_along(od), 
                inframe_cell_stat_sel$noAA[od] + inframe_cell_stat_sel$noAAse[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(inframe_cell_stat_sel$noAA[od] + inframe_cell_stat_sel$noAAse[od], 
                seq_along(od) - 0.2, 
                inframe_cell_stat_sel$noAA[od] + inframe_cell_stat_sel$noAAse[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(inframe_cell_stat_sel$inframePc[od],seq_along(od), 
                inframe_cell_stat_sel$inframePc[od] + inframe_cell_stat_sel$inframeSe[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#FF8C00"))
  grid.segments(inframe_cell_stat_sel$inframePc[od] + inframe_cell_stat_sel$inframeSe[od], 
                seq_along(od) - 0.2, 
                inframe_cell_stat_sel$inframePc[od] + inframe_cell_stat_sel$inframeSe[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#FF8C00"))
})
decorate_annotation("sgname1",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(rev(sgname_sel))), 
                 rot = 0, 
                 x = unit(0, "npc"),
                 y=unit(c(1 : length(sgname_sel)) / length(sgname_sel) -
                          0.5 / length(sgname_sel), "npc"), hjust = 0, 
                 gp = gpar(fontsize = 10, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
dev.off()




