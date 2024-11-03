load("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_spec_result.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt", sep="\t")
tmp <- bulge_pos[,c(4,5)]
tmp <- merge(tmp, sgCmp, by.x="V4", by.y="id")
tmp <- tmp[,c("V4","id2", "V5")]
inframe_result_add_pos <- merge(inframe_result_processed, tmp, by.x="id", by.y="id2")
paired_id <-inframe_result_add_pos[!duplicated(inframe_result_add_pos$id),]
paired_id <- split(paired_id$id, paired_id$V4)
inframe_cor <- lapply(paired_id, function(id){
  id <- unlist(id)
  one_inframe <- inframe_result_add_pos[inframe_result_add_pos$id == id[1],]
  two_inframe <- inframe_result_add_pos[inframe_result_add_pos$id == id[2],]
  
  one_inframe_array <- rep(0, 12)
  names(one_inframe_array) <- as.character(c(-5: 5, "Other"))
  two_inframe_array <- rep(0, 12)
  names(two_inframe_array) <- as.character(c(-5: 5, "Other"))
  for(i in 1 : nrow(one_inframe)){
    one_inframe_array[one_inframe$aa[i]]<-one_inframe$pct[i]
  }
  for(i in 1 : nrow(two_inframe)){
    two_inframe_array[two_inframe$aa[i]]<-two_inframe$pct[i]
  }
  cor(one_inframe_array, two_inframe_array)
})

inframe_cor <- data.frame(pair = names(inframe_cor), cor = unlist(inframe_cor))
inframe_cor <- merge(inframe_cor, tmp, by.x="pair", by.y="V4")
inframe_cor <- inframe_cor[!duplicated(inframe_cor$pair),]
rownames(inframe_cor) <- inframe_cor$pair

noBulge_sg <- sgCmp$id2[sgCmp$isInsert == "Delete"]
withBulge_sg <- sgCmp$id2[sgCmp$isInsert != "Delete"]
inframe_ratio <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Merged1NT/Ins1_inframe_ratio_info.xlsx")
noBulge_inframe_ratio <- inframe_ratio[inframe_ratio$id %in% noBulge_sg,]
withBulge_inframe_ratio <- inframe_ratio[!inframe_ratio$id %in% noBulge_sg,]

noBulge_inframe_spec <- inframe_result_processed[inframe_result_processed$id %in% noBulge_sg,]
withBulge_inframe_spec <- inframe_result_processed[!inframe_result_processed$id %in% noBulge_sg,]


tmp <- inframe_result_processed[inframe_result_processed$aa == "0",]
tmp <- tmp[tmp$id %in% noBulge_sg,]
rownames(tmp) <- tmp$id
tmp <- tmp[noBulge_inframe_ratio$id,]
noBulge_inframe_ratio$noAA <- tmp$pct
noBulge_inframe_ratio$noAAse <- tmp$se
noBulge_inframe_ratio <- noBulge_inframe_ratio[order(noBulge_inframe_ratio$noAA, decreasing = T),]

tmp <- inframe_result_processed[inframe_result_processed$aa == "0",]
tmp <- tmp[tmp$id %in% withBulge_sg,]
rownames(tmp) <- tmp$id
tmp <- tmp[withBulge_inframe_ratio$id,]
withBulge_inframe_ratio$noAA <- tmp$pct
withBulge_inframe_ratio$noAAse <- tmp$se
rownames(withBulge_inframe_ratio) <- withBulge_inframe_ratio$pair
withBulge_inframe_ratio <- withBulge_inframe_ratio[noBulge_inframe_ratio$pair,]

inframe_cor <- inframe_cor[noBulge_inframe_ratio$pair,]

inframe_data <- matrix(0, ncol=12 * 2 + 1, nrow = length(paired_id))
for(index in 1 : nrow(inframe_cor)){
  id <- sgCmp$id2[sgCmp$id == inframe_cor$pair[index]]
  
  noBulge_inframe <- noBulge_inframe_spec[noBulge_inframe_spec$id %in% id,]
  withBulge_inframe <- withBulge_inframe_spec[withBulge_inframe_spec$id %in% id,]
  
  withBulge_inframe_array <- rep(0, 12)
  names(withBulge_inframe_array) <- as.character(c(5: -5, "Others"))
  noBulge_inframe_array <- rep(0, 12)
  names(noBulge_inframe_array) <- as.character(c("Others",-5: 5))
  for(i in 1 : nrow(withBulge_inframe)){
    withBulge_inframe_array[withBulge_inframe$aa[i]]<-withBulge_inframe$pct[i]
  }
  for(i in 1 : nrow(noBulge_inframe)){
    noBulge_inframe_array[noBulge_inframe$aa[i]]<-noBulge_inframe$pct[i]
  }
  inframe_data[index,] <- c(noBulge_inframe_array, inframe_cor$cor[index], withBulge_inframe_array)
}
inframe_data[is.na(inframe_data)] <- 0



max_0aa <- ceiling(max(c(max(noBulge_inframe_ratio$noAA + 
                               noBulge_inframe_ratio$noAAse, na.rm = T),
                         max(withBulge_inframe_ratio$noAA + 
                               withBulge_inframe_ratio$noAAse, na.rm = T))))
max_inframe <- max(c(max(noBulge_inframe_ratio$inframePct + 
                           noBulge_inframe_ratio$inframeSe, na.rm = T), 
                     max(
                       withBulge_inframe_ratio$inframePct + 
                         withBulge_inframe_ratio$inframeSe, na.rm = T)))
left_anno <- rowAnnotation(
  sgname1 = anno_empty(width = unit(4, "cm"), border = F),
  tissue_0aa = anno_barplot(noBulge_inframe_ratio$noAA, 
                            border = T,
                            gp = gpar(fill = "#494490"), 
                            bar_width = 0.5,
                            axis_param = list(direction = "reverse", 
                                              "side"="top"),
                            ylim = c(0, max_0aa), 
                            width = unit(2, "cm")),
  noBulge_inframe = anno_barplot(noBulge_inframe_ratio$inframePct, 
                                border = T,
                                gp = gpar(fill = "#494490"), 
                                bar_width = 0.5,
                                axis_param = list(direction = "reverse", 
                                                  "side"="top"),
                                ylim = c(0, max_inframe), 
                                width = unit(2, "cm")), 
  show_annotation_name = F)

right_anno <- rowAnnotation(
  withBulge_inframe = anno_barplot(withBulge_inframe_ratio$inframePct, 
                              border = T,
                              gp = gpar(fill = "#494490"), 
                              bar_width = 0.5,
                              axis_param = list("side"="top"),
                              ylim = c(0, max_inframe), 
                              width = unit(2, "cm")), 
  cell_0aa = anno_barplot(withBulge_inframe_ratio$noAA, 
                          border = T,
                          gp = gpar(fill = "#494490"), 
                          bar_width = 0.5,
                          axis_param = list("side"="top"),
                          ylim = c(0, max_0aa), 
                          width = unit(2, "cm")),
  
  show_annotation_name = F)
col_acc <- colSums(inframe_data)
col_acc[13] <- NA



top_anno <- HeatmapAnnotation(acc = anno_barplot(col_acc[1:12] / 100, 
                                                 border = T,
                                                 gp = gpar(fill = "#494490"), 
                                                 bar_width = 0.5,
                                                 width = unit(2, "cm"), 
                                                 ylim = c(0, max(col_acc, na.rm = T) / 100)),
                              pos = anno_empty(border = F, width = unit(2, "mm")))
plot_mat_left <- inframe_data[,c(1:12)]
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

plot_mat_right <- inframe_data[,c(14 : 25)]
top_anno2 <- HeatmapAnnotation(acc = anno_barplot(col_acc[14: 25] /100, 
                                                  border = T,
                                                  gp = gpar(fill = "#494490"), 
                                                  bar_width = 0.5,
                                                  width = unit(2, "cm"), 
                                                  ylim = c(0, max(col_acc, na.rm = T) / 100)),
                               pos2 = anno_empty(border = F, width = unit(2, "mm")))
left_anno2 <- rowAnnotation(
  Cor = anno_barplot(inframe_cor$cor, 
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



pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_table_del1mut_0aa_order.pdf", 
    width = 14, height = 36)

draw(ht + ht2)
decorate_annotation("tissue_0aa", slice = 1, {
  od = nrow(inframe_data) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - noBulge_inframe_ratio$noAA[od],seq_along(od), 
                xcale - (noBulge_inframe_ratio$noAA[od] + noBulge_inframe_ratio$noAAse[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (noBulge_inframe_ratio$noAA[od] + noBulge_inframe_ratio$noAAse[od]), 
                seq_along(od) - 0.2, 
                xcale - (noBulge_inframe_ratio$noAA[od] + noBulge_inframe_ratio$noAAse[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  
})
decorate_annotation("noBulge_inframe", slice = 1, {
  od = nrow(inframe_data) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - noBulge_inframe_ratio$inframePct[od],seq_along(od), 
                xcale - (noBulge_inframe_ratio$inframePct[od] + noBulge_inframe_ratio$inframeSe[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (noBulge_inframe_ratio$inframePct[od] + noBulge_inframe_ratio$inframeSe[od]), 
                seq_along(od) - 0.2, 
                xcale - (noBulge_inframe_ratio$inframePct[od] + noBulge_inframe_ratio$inframeSe[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  
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

decorate_annotation("withBulge_inframe", slice = 1, {
  od = nrow(inframe_data) : 1
  grid.segments(withBulge_inframe_ratio$inframePc[od],seq_along(od), 
                withBulge_inframe_ratio$inframePc[od] + withBulge_inframe_ratio$inframeSe[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(withBulge_inframe_ratio$inframePc[od] + withBulge_inframe_ratio$inframeSe[od], 
                seq_along(od) - 0.2, 
                withBulge_inframe_ratio$inframePc[od] + withBulge_inframe_ratio$inframeSe[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
})
decorate_annotation("cell_0aa", slice = 1, {
  od = nrow(inframe_data) : 1
  grid.segments(withBulge_inframe_ratio$noAA[od],seq_along(od), 
                withBulge_inframe_ratio$noAA[od] + withBulge_inframe_ratio$noAAse[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(withBulge_inframe_ratio$noAA[od] + withBulge_inframe_ratio$noAAse[od], 
                seq_along(od) - 0.2, 
                withBulge_inframe_ratio$noAA[od] + withBulge_inframe_ratio$noAAse[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
})
sgname_sel <- noBulge_inframe_ratio$id
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



####画一个过滤0.6的结果
load("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_spec_result.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt", sep="\t")
tmp <- bulge_pos[,c(4,5)]
tmp <- merge(tmp, sgCmp, by.x="V4", by.y="id")
tmp <- tmp[,c("V4","id2", "V5")]
inframe_result_add_pos <- merge(inframe_result_processed, tmp, by.x="id", by.y="id2")
paired_id <-inframe_result_add_pos[!duplicated(inframe_result_add_pos$id),]
paired_id <- split(paired_id$id, paired_id$V4)
inframe_cor <- lapply(paired_id, function(id){
  id <- unlist(id)
  one_inframe <- inframe_result_add_pos[inframe_result_add_pos$id == id[1],]
  two_inframe <- inframe_result_add_pos[inframe_result_add_pos$id == id[2],]
  
  one_inframe_array <- rep(0, 12)
  names(one_inframe_array) <- as.character(c(-5: 5, "Other"))
  two_inframe_array <- rep(0, 12)
  names(two_inframe_array) <- as.character(c(-5: 5, "Other"))
  for(i in 1 : nrow(one_inframe)){
    one_inframe_array[one_inframe$aa[i]]<-one_inframe$pct[i]
  }
  for(i in 1 : nrow(two_inframe)){
    two_inframe_array[two_inframe$aa[i]]<-two_inframe$pct[i]
  }
  cor(one_inframe_array, two_inframe_array)
})

inframe_cor <- data.frame(pair = names(inframe_cor), cor = unlist(inframe_cor))
inframe_cor <- merge(inframe_cor, tmp, by.x="pair", by.y="V4")
inframe_cor <- inframe_cor[!duplicated(inframe_cor$pair),]
rownames(inframe_cor) <- inframe_cor$pair
inframe_cor <- inframe_cor[inframe_cor$cor > 0.6,]

noBulge_sg <- sgCmp$id2[sgCmp$isInsert == "Delete"]
withBulge_sg <- sgCmp$id2[sgCmp$isInsert != "Delete"]
inframe_ratio <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Merged1NT/Ins1_inframe_ratio_info.xlsx")
noBulge_inframe_ratio <- inframe_ratio[inframe_ratio$id %in% noBulge_sg,]
withBulge_inframe_ratio <- inframe_ratio[!inframe_ratio$id %in% noBulge_sg,]

noBulge_inframe_spec <- inframe_result_processed[inframe_result_processed$id %in% noBulge_sg,]
withBulge_inframe_spec <- inframe_result_processed[!inframe_result_processed$id %in% noBulge_sg,]


tmp <- inframe_result_processed[inframe_result_processed$aa == "0",]
tmp <- tmp[tmp$id %in% noBulge_sg,]
rownames(tmp) <- tmp$id
tmp <- tmp[noBulge_inframe_ratio$id,]
noBulge_inframe_ratio$noAA <- tmp$pct
noBulge_inframe_ratio$noAAse <- tmp$se
noBulge_inframe_ratio <- noBulge_inframe_ratio[order(noBulge_inframe_ratio$noAA, decreasing = T),]

tmp <- inframe_result_processed[inframe_result_processed$aa == "0",]
tmp <- tmp[tmp$id %in% withBulge_sg,]
rownames(tmp) <- tmp$id
tmp <- tmp[withBulge_inframe_ratio$id,]
withBulge_inframe_ratio$noAA <- tmp$pct
withBulge_inframe_ratio$noAAse <- tmp$se
rownames(withBulge_inframe_ratio) <- withBulge_inframe_ratio$pair
withBulge_inframe_ratio <- withBulge_inframe_ratio[noBulge_inframe_ratio$pair,]

noBulge_inframe_ratio <- noBulge_inframe_ratio[noBulge_inframe_ratio$pair %in% inframe_cor$pair,]
withBulge_inframe_ratio <- withBulge_inframe_ratio[withBulge_inframe_ratio$pair %in% inframe_cor$pair,]

inframe_cor <- inframe_cor[noBulge_inframe_ratio$pair,]

inframe_data <- matrix(0, ncol=12 * 2 + 1, nrow = nrow(inframe_cor))
for(index in 1 : nrow(inframe_cor)){
  id <- sgCmp$id2[sgCmp$id == inframe_cor$pair[index]]
  
  noBulge_inframe <- noBulge_inframe_spec[noBulge_inframe_spec$id %in% id,]
  withBulge_inframe <- withBulge_inframe_spec[withBulge_inframe_spec$id %in% id,]
  
  withBulge_inframe_array <- rep(0, 12)
  names(withBulge_inframe_array) <- as.character(c(5: -5, "Others"))
  noBulge_inframe_array <- rep(0, 12)
  names(noBulge_inframe_array) <- as.character(c("Others",-5: 5))
  for(i in 1 : nrow(withBulge_inframe)){
    withBulge_inframe_array[withBulge_inframe$aa[i]]<-withBulge_inframe$pct[i]
  }
  for(i in 1 : nrow(noBulge_inframe)){
    noBulge_inframe_array[noBulge_inframe$aa[i]]<-noBulge_inframe$pct[i]
  }
  inframe_data[index,] <- c(noBulge_inframe_array, inframe_cor$cor[index], withBulge_inframe_array)
}
inframe_data[is.na(inframe_data)] <- 0



max_0aa <- ceiling(max(c(max(noBulge_inframe_ratio$noAA + 
                               noBulge_inframe_ratio$noAAse, na.rm = T),
                         max(withBulge_inframe_ratio$noAA + 
                               withBulge_inframe_ratio$noAAse, na.rm = T))))
max_inframe <- max(c(max(noBulge_inframe_ratio$inframePct + 
                           noBulge_inframe_ratio$inframeSe, na.rm = T), 
                     max(
                       withBulge_inframe_ratio$inframePct + 
                         withBulge_inframe_ratio$inframeSe, na.rm = T)))
left_anno <- rowAnnotation(
  sgname1 = anno_empty(width = unit(4, "cm"), border = F),
  tissue_0aa = anno_barplot(noBulge_inframe_ratio$noAA, 
                            border = T,
                            gp = gpar(fill = "#494490"), 
                            bar_width = 0.5,
                            axis_param = list(direction = "reverse", 
                                              "side"="top"),
                            ylim = c(0, max_0aa), 
                            width = unit(2, "cm")),
  noBulge_inframe = anno_barplot(noBulge_inframe_ratio$inframePct, 
                                 border = T,
                                 gp = gpar(fill = "#494490"), 
                                 bar_width = 0.5,
                                 axis_param = list(direction = "reverse", 
                                                   "side"="top"),
                                 ylim = c(0, max_inframe), 
                                 width = unit(2, "cm")), 
  show_annotation_name = F)

right_anno <- rowAnnotation(
  withBulge_inframe = anno_barplot(withBulge_inframe_ratio$inframePct, 
                                   border = T,
                                   gp = gpar(fill = "#494490"), 
                                   bar_width = 0.5,
                                   axis_param = list("side"="top"),
                                   ylim = c(0, max_inframe), 
                                   width = unit(2, "cm")), 
  cell_0aa = anno_barplot(withBulge_inframe_ratio$noAA, 
                          border = T,
                          gp = gpar(fill = "#494490"), 
                          bar_width = 0.5,
                          axis_param = list("side"="top"),
                          ylim = c(0, max_0aa), 
                          width = unit(2, "cm")),
  
  show_annotation_name = F)
col_acc <- colSums(inframe_data)
col_acc[13] <- NA



top_anno <- HeatmapAnnotation(acc = anno_barplot(col_acc[1:12] / 100, 
                                                 border = T,
                                                 gp = gpar(fill = "#494490"), 
                                                 bar_width = 0.5,
                                                 width = unit(2, "cm"), 
                                                 ylim = c(0, max(col_acc, na.rm = T) / 100)),
                              pos = anno_empty(border = F, width = unit(2, "mm")))
plot_mat_left <- inframe_data[,c(1:12)]
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

plot_mat_right <- inframe_data[,c(14 : 25)]
top_anno2 <- HeatmapAnnotation(acc = anno_barplot(col_acc[14: 25] /100, 
                                                  border = T,
                                                  gp = gpar(fill = "#494490"), 
                                                  bar_width = 0.5,
                                                  width = unit(2, "cm"), 
                                                  ylim = c(0, max(col_acc, na.rm = T) / 100)),
                               pos2 = anno_empty(border = F, width = unit(2, "mm")))
left_anno2 <- rowAnnotation(
  Cor = anno_barplot(inframe_cor$cor, 
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



pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_table_del1mut_0aa_order_large60cor.pdf", 
    width = 14, height = 36)

draw(ht + ht2)
decorate_annotation("tissue_0aa", slice = 1, {
  od = nrow(inframe_data) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - noBulge_inframe_ratio$noAA[od],seq_along(od), 
                xcale - (noBulge_inframe_ratio$noAA[od] + noBulge_inframe_ratio$noAAse[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (noBulge_inframe_ratio$noAA[od] + noBulge_inframe_ratio$noAAse[od]), 
                seq_along(od) - 0.2, 
                xcale - (noBulge_inframe_ratio$noAA[od] + noBulge_inframe_ratio$noAAse[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  
})
decorate_annotation("noBulge_inframe", slice = 1, {
  od = nrow(inframe_data) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - noBulge_inframe_ratio$inframePct[od],seq_along(od), 
                xcale - (noBulge_inframe_ratio$inframePct[od] + noBulge_inframe_ratio$inframeSe[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (noBulge_inframe_ratio$inframePct[od] + noBulge_inframe_ratio$inframeSe[od]), 
                seq_along(od) - 0.2, 
                xcale - (noBulge_inframe_ratio$inframePct[od] + noBulge_inframe_ratio$inframeSe[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  
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

decorate_annotation("withBulge_inframe", slice = 1, {
  od = nrow(inframe_data) : 1
  grid.segments(withBulge_inframe_ratio$inframePc[od],seq_along(od), 
                withBulge_inframe_ratio$inframePc[od] + withBulge_inframe_ratio$inframeSe[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(withBulge_inframe_ratio$inframePc[od] + withBulge_inframe_ratio$inframeSe[od], 
                seq_along(od) - 0.2, 
                withBulge_inframe_ratio$inframePc[od] + withBulge_inframe_ratio$inframeSe[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
})
decorate_annotation("cell_0aa", slice = 1, {
  od = nrow(inframe_data) : 1
  grid.segments(withBulge_inframe_ratio$noAA[od],seq_along(od), 
                withBulge_inframe_ratio$noAA[od] + withBulge_inframe_ratio$noAAse[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(withBulge_inframe_ratio$noAA[od] + withBulge_inframe_ratio$noAAse[od], 
                seq_along(od) - 0.2, 
                withBulge_inframe_ratio$noAA[od] + withBulge_inframe_ratio$noAAse[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
})
sgname_sel <- noBulge_inframe_ratio$id
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




###重新画一个按照相关性排序的结果

load("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_spec_result.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt", sep="\t")
tmp <- bulge_pos[,c(4,5)]
tmp <- merge(tmp, sgCmp, by.x="V4", by.y="id")
tmp <- tmp[,c("V4","id2", "V5")]
inframe_result_add_pos <- merge(inframe_result_processed, tmp, by.x="id", by.y="id2")
paired_id <-inframe_result_add_pos[!duplicated(inframe_result_add_pos$id),]
paired_id <- split(paired_id$id, paired_id$V4)
inframe_cor <- lapply(paired_id, function(id){
  id <- unlist(id)
  one_inframe <- inframe_result_add_pos[inframe_result_add_pos$id == id[1],]
  two_inframe <- inframe_result_add_pos[inframe_result_add_pos$id == id[2],]
  
  one_inframe_array <- rep(0, 12)
  names(one_inframe_array) <- as.character(c(-5: 5, "Other"))
  two_inframe_array <- rep(0, 12)
  names(two_inframe_array) <- as.character(c(-5: 5, "Other"))
  for(i in 1 : nrow(one_inframe)){
    one_inframe_array[one_inframe$aa[i]]<-one_inframe$pct[i]
  }
  for(i in 1 : nrow(two_inframe)){
    two_inframe_array[two_inframe$aa[i]]<-two_inframe$pct[i]
  }
  cor(one_inframe_array, two_inframe_array)
})

inframe_cor <- data.frame(pair = names(inframe_cor), cor = unlist(inframe_cor))
inframe_cor <- merge(inframe_cor, tmp, by.x="pair", by.y="V4")
inframe_cor <- inframe_cor[!duplicated(inframe_cor$pair),]
rownames(inframe_cor) <- inframe_cor$pair
inframe_cor <- inframe_cor[order(inframe_cor$cor,decreasing = T),]
noBulge_sg <- sgCmp$id2[sgCmp$isInsert == "Delete"]
withBulge_sg <- sgCmp$id2[sgCmp$isInsert != "Delete"]
inframe_ratio <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Merged1NT/Ins1_inframe_ratio_info.xlsx")
noBulge_inframe_ratio <- inframe_ratio[inframe_ratio$id %in% noBulge_sg,]
withBulge_inframe_ratio <- inframe_ratio[!inframe_ratio$id %in% noBulge_sg,]

noBulge_inframe_spec <- inframe_result_processed[inframe_result_processed$id %in% noBulge_sg,]
withBulge_inframe_spec <- inframe_result_processed[!inframe_result_processed$id %in% noBulge_sg,]


tmp <- inframe_result_processed[inframe_result_processed$aa == "0",]
tmp <- tmp[tmp$id %in% noBulge_sg,]
rownames(tmp) <- tmp$id
tmp <- tmp[noBulge_inframe_ratio$id,]
noBulge_inframe_ratio$noAA <- tmp$pct
noBulge_inframe_ratio$noAAse <- tmp$se
rownames(noBulge_inframe_ratio) <- noBulge_inframe_ratio$pair
noBulge_inframe_ratio <- noBulge_inframe_ratio[inframe_cor$pair,]

tmp <- inframe_result_processed[inframe_result_processed$aa == "0",]
tmp <- tmp[tmp$id %in% withBulge_sg,]
rownames(tmp) <- tmp$id
tmp <- tmp[withBulge_inframe_ratio$id,]
withBulge_inframe_ratio$noAA <- tmp$pct
withBulge_inframe_ratio$noAAse <- tmp$se
rownames(withBulge_inframe_ratio) <- withBulge_inframe_ratio$pair
withBulge_inframe_ratio <- withBulge_inframe_ratio[noBulge_inframe_ratio$pair,]


inframe_data <- matrix(0, ncol=12 * 2 + 1, nrow = length(paired_id))
for(index in 1 : nrow(inframe_cor)){
  id <- sgCmp$id2[sgCmp$id == inframe_cor$pair[index]]
  
  noBulge_inframe <- noBulge_inframe_spec[noBulge_inframe_spec$id %in% id,]
  withBulge_inframe <- withBulge_inframe_spec[withBulge_inframe_spec$id %in% id,]
  
  withBulge_inframe_array <- rep(0, 12)
  names(withBulge_inframe_array) <- as.character(c(5: -5, "Others"))
  noBulge_inframe_array <- rep(0, 12)
  names(noBulge_inframe_array) <- as.character(c("Others",-5: 5))
  for(i in 1 : nrow(withBulge_inframe)){
    withBulge_inframe_array[withBulge_inframe$aa[i]]<-withBulge_inframe$pct[i]
  }
  for(i in 1 : nrow(noBulge_inframe)){
    noBulge_inframe_array[noBulge_inframe$aa[i]]<-noBulge_inframe$pct[i]
  }
  inframe_data[index,] <- c(noBulge_inframe_array, inframe_cor$cor[index], withBulge_inframe_array)
}
inframe_data[is.na(inframe_data)] <- 0


####分开画  先画delete突变

max_0aa <- ceiling(max(c(max(noBulge_inframe_ratio$noAA + 
                               noBulge_inframe_ratio$noAAse, na.rm = T),
                         max(withBulge_inframe_ratio$noAA + 
                               withBulge_inframe_ratio$noAAse, na.rm = T))))
max_inframe <- max(c(max(noBulge_inframe_ratio$inframePct + 
                           noBulge_inframe_ratio$inframeSe, na.rm = T), 
                     max(
                       withBulge_inframe_ratio$inframePct + 
                         withBulge_inframe_ratio$inframeSe, na.rm = T)))
left_anno <- rowAnnotation(
  sgname1 = anno_empty(width = unit(4, "cm"), border = F),
  tissue_0aa = anno_barplot(noBulge_inframe_ratio$noAA, 
                            border = T,
                            gp = gpar(fill = "#494490"), 
                            bar_width = 0.5,
                            axis_param = list(direction = "reverse", 
                                              "side"="top"),
                            ylim = c(0, max_0aa), 
                            width = unit(2, "cm")),
  noBulge_inframe = anno_barplot(noBulge_inframe_ratio$inframePct, 
                                 border = T,
                                 gp = gpar(fill = "#494490"), 
                                 bar_width = 0.5,
                                 axis_param = list(direction = "reverse", 
                                                   "side"="top"),
                                 ylim = c(0, max_inframe), 
                                 width = unit(2, "cm")), 
  show_annotation_name = F)

right_anno <- rowAnnotation(
  withBulge_inframe = anno_barplot(withBulge_inframe_ratio$inframePct, 
                                   border = T,
                                   gp = gpar(fill = "#494490"), 
                                   bar_width = 0.5,
                                   axis_param = list("side"="top"),
                                   ylim = c(0, max_inframe), 
                                   width = unit(2, "cm")), 
  cell_0aa = anno_barplot(withBulge_inframe_ratio$noAA, 
                          border = T,
                          gp = gpar(fill = "#494490"), 
                          bar_width = 0.5,
                          axis_param = list("side"="top"),
                          ylim = c(0, max_0aa), 
                          width = unit(2, "cm")),
  
  show_annotation_name = F)
col_acc <- colSums(inframe_data)
col_acc[13] <- NA



top_anno <- HeatmapAnnotation(acc = anno_barplot(col_acc[1:12] / 100, 
                                                 border = T,
                                                 gp = gpar(fill = "#494490"), 
                                                 bar_width = 0.5,
                                                 width = unit(2, "cm"), 
                                                 ylim = c(0, max(col_acc, na.rm = T) / 100)),
                              pos = anno_empty(border = F, width = unit(2, "mm")))
plot_mat_left <- inframe_data[,c(1:12)]
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

plot_mat_right <- inframe_data[,c(14 : 25)]
top_anno2 <- HeatmapAnnotation(acc = anno_barplot(col_acc[14: 25] /100, 
                                                  border = T,
                                                  gp = gpar(fill = "#494490"), 
                                                  bar_width = 0.5,
                                                  width = unit(2, "cm"), 
                                                  ylim = c(0, max(col_acc, na.rm = T) / 100)),
                               pos2 = anno_empty(border = F, width = unit(2, "mm")))
left_anno2 <- rowAnnotation(
  Cor = anno_barplot(inframe_cor$cor, 
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



pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_table_del1mut_cor_order.pdf", 
    width = 14, height = 36)

draw(ht + ht2)
decorate_annotation("tissue_0aa", slice = 1, {
  od = nrow(inframe_data) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - noBulge_inframe_ratio$noAA[od],seq_along(od), 
                xcale - (noBulge_inframe_ratio$noAA[od] + noBulge_inframe_ratio$noAAse[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (noBulge_inframe_ratio$noAA[od] + noBulge_inframe_ratio$noAAse[od]), 
                seq_along(od) - 0.2, 
                xcale - (noBulge_inframe_ratio$noAA[od] + noBulge_inframe_ratio$noAAse[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  
})
decorate_annotation("noBulge_inframe", slice = 1, {
  od = nrow(inframe_data) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - noBulge_inframe_ratio$inframePct[od],seq_along(od), 
                xcale - (noBulge_inframe_ratio$inframePct[od] + noBulge_inframe_ratio$inframeSe[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (noBulge_inframe_ratio$inframePct[od] + noBulge_inframe_ratio$inframeSe[od]), 
                seq_along(od) - 0.2, 
                xcale - (noBulge_inframe_ratio$inframePct[od] + noBulge_inframe_ratio$inframeSe[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  
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

decorate_annotation("withBulge_inframe", slice = 1, {
  od = nrow(inframe_data) : 1
  grid.segments(withBulge_inframe_ratio$inframePc[od],seq_along(od), 
                withBulge_inframe_ratio$inframePc[od] + withBulge_inframe_ratio$inframeSe[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(withBulge_inframe_ratio$inframePc[od] + withBulge_inframe_ratio$inframeSe[od], 
                seq_along(od) - 0.2, 
                withBulge_inframe_ratio$inframePc[od] + withBulge_inframe_ratio$inframeSe[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
})
decorate_annotation("cell_0aa", slice = 1, {
  od = nrow(inframe_data) : 1
  grid.segments(withBulge_inframe_ratio$noAA[od],seq_along(od), 
                withBulge_inframe_ratio$noAA[od] + withBulge_inframe_ratio$noAAse[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(withBulge_inframe_ratio$noAA[od] + withBulge_inframe_ratio$noAAse[od], 
                seq_along(od) - 0.2, 
                withBulge_inframe_ratio$noAA[od] + withBulge_inframe_ratio$noAAse[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
})
sgname_sel <- noBulge_inframe_ratio$id
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

