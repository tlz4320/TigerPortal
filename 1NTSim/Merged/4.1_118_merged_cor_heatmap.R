###load 之前114个sgRNA的结果
print(load("~/data/project/ear_project/gene_therapy_ll/Previews/Result/Rep123_indel_stat.rda"))
sample_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)

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

split_name_tissue <- split_name[split_name$tissueid %in% names(old_tissue),]
split_name_tissue <- na.omit(split_name_tissue)



split_name_cell$mutation <- unlist(lapply(split_name_cell$mutation, function(x){
  x <- unlist(str_split(x, "[-]"))
  x[2] <- stringr::str_to_title(x[2])
  paste(x, collapse = "-")
}))
split_name_cell$type <- unlist(lapply(split_name_cell$mutation, function(x){
  str_sub(x, str_length(x) - 1,  str_length(x) - 1)
}))


split_name_tissue$mutation <- unlist(lapply(split_name_tissue$mutation, function(x){
  x <- unlist(str_split(x, "[-]"))
  x[2] <- stringr::str_to_title(x[2])
  paste(x, collapse = "-")
}))
split_name_tissue$type <- unlist(lapply(split_name_tissue$mutation, function(x){
  str_sub(x, str_length(x) - 1,  str_length(x) - 1)
}))

###读取新测的4个sgRNA结果 然后合并到之前的结果里面去
print(load("~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data_newer.rda"))
rep1_new_otof <- otof_cell_sg12_15[seq(1,8, 2)]
rep2_new_otof <- otof_cell_sg12_15[seq(2,8, 2)]
for(name in names(rep1_new_otof)){
  name2 <- str_remove(name, "[-].*")
  old_rep1[[name2]] <- rep1_new_otof[[name]]
}
for(name in names(rep2_new_otof)){
  name2 <- str_remove(name, "[-].*")
  old_rep2[[name2]] <- rep2_new_otof[[name]]
}
split_name_cell <- data.frame(rbind(split_name_cell, 
                                    data.frame("cellid" = paste0("m", 12 : 15),
                                               "tissueid" = paste0("m", 12 : 15),
                                               "mutation" = 'w-Otof-1236dC', 
                                               "type" = 'd')))
otof_tissue_sg7_15_rep1 <- otof_tissue_sg7_15[c(9,10,11)]
otof_tissue_sg7_15_rep2 <- otof_tissue_sg7_15[c(2:5)]
names(otof_tissue_sg7_15_rep1) <- c("m12-tissue", "m13-tissue", "m14-tissue")
names(otof_tissue_sg7_15_rep2) <- c("m12-tissue", "m13-tissue", "m14-tissue", "m15-tissue")


for(name in names(otof_tissue_sg7_15_rep1)){
  old_rep1[[name]] <- otof_tissue_sg7_15_rep1[[name]]
}
for(name in names(otof_tissue_sg7_15_rep2)){
  old_rep2[[name]] <- otof_tissue_sg7_15_rep2[[name]]
}
split_name_tissue <- data.frame(rbind(split_name_tissue, 
                                      data.frame("cellid" = paste0("m", 12 : 15, "-tissue"),
                                                 "tissueid" = paste0("m", 12 : 15, "-tissue"),
                                                 "mutation" = 'w-Otof-1236dC', 
                                                 "type" = 'd')))

###因为加入了几个新的样本 所以重新开始计算均值
old_rep123 <- list()
total_samples <- unique(c(names(old_rep1), names(old_rep2), names(old_rep3)))
old_result <- list(Rep1 = old_rep1, Rep2 = old_rep2, Rep3 = old_rep3)

for(name in total_samples){
  
  sample_count <- sum(unlist(lapply(old_result, function(x){
    name %in% names(x)
  })))
  tmp_rep123 <- list()
  for(sp in names(old_result)){
    if(name %in% names(old_result[[sp]])){
      tmp_rep123[[sp]] <- old_result[[sp]][[name]]
    }
  }
  tmp_rep123 <- data.frame(do.call(rbind, tmp_rep123))
  tmp_rep123 <- lapply(split(tmp_rep123$fq, tmp_rep123$indel_size), sum)
  tmp_rep123 <- data.frame(indel_size = as.numeric(names(tmp_rep123)), fq = unlist(tmp_rep123) / sample_count)
  tmp_rep123 <- tmp_rep123[order(tmp_rep123$indel_size),]
  old_rep123[[name]] <- tmp_rep123
}
old_tissue <- old_rep123[str_remove(names(old_rep123), "CRISPResso_on_") %in% split_name_tissue$tissueid]
old_cell <- old_rep123[str_remove(names(old_rep123), "CRISPResso_on_") %in% split_name_cell$cellid]

names(old_tissue) <- str_remove(names(old_tissue), "CRISPResso_on_")
names(old_cell) <- str_remove(names(old_cell), "CRISPResso_on_")


split_name_tissue$id <- str_remove(split_name_tissue$tissueid, "[-].*")
split_name_cell$id <- str_remove(split_name_cell$cellid, "[-].*")
split_name <- merge(split_name_tissue[,-1], split_name_cell[,-2], by="id")
split_name <- split_name[,c(-1, -3,-4)]
colnames(split_name) <- str_remove(colnames(split_name), "[.].*")
###接下来根据是Insert突变还是Delete突变来计算Del1与Ins1的编辑比例

###计算均值的相关性，因为本次数据每次样本是分开测得 不好AA BB CC一一对比了
max_region <- c(-200, 30)
total_pair_sg_cor <- lapply(1 : nrow(split_name), function(i){
  tissue_id <- split_name$tissueid[i]
  cell_id <- split_name$cellid[i]
  sel_sg1 <- old_tissue[[tissue_id]]
  sel_sg2 <- old_cell[[cell_id]]
  max_del <- min(c(sel_sg1[,1], sel_sg2[,1]))
  max_in <- max(c(sel_sg1[,1], sel_sg2[,1]))
  region <- max_del : max_in
  tmp_region <- region[!region %in% sel_sg1[,1]]
  if(length(tmp_region) != 0){
    sel_sg1 <- data.frame(rbind(sel_sg1, data.frame(indel_size = tmp_region, fq = 0)))
  }
  tmp_region <- region[!region %in% sel_sg2[,1]]
  if(length(tmp_region) != 0){
    sel_sg2 <- data.frame(rbind(sel_sg2, data.frame(indel_size = tmp_region, fq = 0)))
  }
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  common_indel <- intersect(sel_sg1[,1], sel_sg2[,1])
  sel_sg1 <- sel_sg1[sel_sg1[,1] %in% common_indel,]
  sel_sg2 <- sel_sg2[sel_sg2[,1] %in% common_indel,]
  sel_sg1 <- sel_sg1[order(sel_sg1[,1]),]
  sel_sg2 <- sel_sg2[order(sel_sg2[,1]),]
  # rm0 <- sel_sg1[,2] == 0 & sel_sg2[,2] == 0
  #cor.test(sel_sg1[!rm0,2], sel_sg2[!rm0,2])
  cor.test(sel_sg1[,2], sel_sg2[,2])
})

total_pair_sg_cor <- data.frame(tissueid = split_name$tissueid,
                                cellid = split_name$cellid,
                                mutation = split_name$mutation,
                                type = split_name$type,
                                cor = unlist(lapply(total_pair_sg_cor, function(x){x$estimate})),
                                pval = unlist(lapply(total_pair_sg_cor, function(x){x$p.value})))
total_pair_sg_cor <- total_pair_sg_cor[order(total_pair_sg_cor$cor, decreasing = T),]


rownames(total_pair_sg_cor) <- total_pair_sg_cor$tissueid
tmp <- sample_info[sample_info$`Tissue-w-W(Mix)` %in% total_pair_sg_cor$tissueid,]
total_pair_sg_cor$name <- total_pair_sg_cor$mutation
total_pair_sg_cor[tmp$`Tissue-w-W(Mix)`, "name"] <- tmp$GRNAa
total_pair_sg_cor[total_pair_sg_cor$mutation == "w-Otof-1236dC" & 
                    total_pair_sg_cor$cellid != "105-001", "name"] <- 
  paste0("w-Otof-1236dC-", total_pair_sg_cor[total_pair_sg_cor$mutation == "w-Otof-1236dC" & 
                      total_pair_sg_cor$cellid != "105-001", "cellid"])
  

plot_region <- c(-8: 5)

indel_table_tissue <- lapply(total_pair_sg_cor$tissueid, function(x){
  res <- old_tissue[[x]]
  res[res[,1] == "0",2] <- 0
  res[,3] <- res[,2] / sum(res[,2]) * 100

  res[res[,1] %in% plot_region,]
})
names(indel_table_tissue) <- total_pair_sg_cor$tissueid
indel_table_cell <- lapply(total_pair_sg_cor$cellid, function(x){
  res <- old_cell[[x]]
  res[res[,1] == "0",2] <- 0
  res[,3] <- res[,2] / sum(res[,2]) * 100
  res[res[,1] %in% plot_region,]
})
names(indel_table_cell) <- total_pair_sg_cor$cellid


plot_mat <- matrix(0, ncol = length(plot_region) * 2 - 1, nrow = nrow(total_pair_sg_cor))
for(i in 1 : nrow(total_pair_sg_cor)){
  tmp1 <- indel_table_tissue[[i]]
  tmp1 <- tmp1[tmp1[,1] != 0,]
  tmp2 <- indel_table_cell[[i]]
  tmp2 <- tmp2[tmp2[,1] != 0,]
  plot_mat[i, 1 : (length(plot_region) - 1)] <- rev(tmp1[,3]) / 100
  plot_mat[i, length(plot_region)] <- total_pair_sg_cor$cor[i]
  plot_mat[i, (length(plot_region) + 1) : ncol(plot_mat)] <- tmp2[,3] / 100
}

library(ComplexHeatmap)
max(unlist(plot_mat[,-length(plot_region)]))
indel_color <- colorRamp2(c(0, 1), c("white", "blue"))
cor_color <- colorRamp2(c(0, 1), c("white", "red"))
indel_reads_tissue <- lapply(total_pair_sg_cor$tissueid, function(x){
  tmp <- c()
  for(i in 1 : 3){
    names <- str_remove(names(old_result[[i]]), "CRISPResso_on_")
    if(x %in% names){
      y <- old_result[[i]][[which(names == x)]]
      tmp <- c(sum(y[y[,1] %in% plot_region & y[,1] != 0,2]), tmp)
    }
  }
  
  data.frame(count = log(mean(tmp)), se = plotrix::std.error(log(tmp)))
})
indel_reads_tissue <- data.frame(do.call(rbind, indel_reads_tissue))
indel_reads_tissue$se[is.na(indel_reads_tissue$se)] <- 0
indel_reads_cell <- lapply(total_pair_sg_cor$cellid, function(x){
  tmp <- c()
  for(i in 1 : 3){
    names <- str_remove(names(old_result[[i]]), "CRISPResso_on_")
    if(x %in% names){
      y <- old_result[[i]][[which(names == x)]]
      tmp <- c(sum(y[y[,1] %in% plot_region & y[,1] != 0,2]), tmp)
    }
  }
  
  data.frame(count = log(mean(tmp)), se = plotrix::std.error(log(tmp)))
})
indel_reads_cell <- data.frame(do.call(rbind, indel_reads_cell))
indel_reads_cell$se[is.na(indel_reads_cell$se)] <- 0

left_anno <- rowAnnotation(
  sgname1 = anno_empty(width = unit(2, "cm"), border = F),
  left = anno_barplot(indel_reads_tissue$count, 
                                               border = T,
                                               gp = gpar(fill = "#494490"), 
                                               bar_width = 0.5,
                                               axis_param = list(direction = "reverse", 
                                                                 "side"="top"),
                                               ylim = c(0,
                                                        ceiling(
                                                          max(
                                                            indel_reads_tissue$count + 
                                                              indel_reads_tissue$se))), 
                                               width = unit(2, "cm")), 
                           show_annotation_name = F)

right_anno <- rowAnnotation(right = anno_barplot(indel_reads_cell$count, 
                                                 border = T,
                                                 gp = gpar(fill = "#494490"), 
                                                 axis_param = c("side"="top"), 
                                                 bar_width = 0.5,
                                                 ylim = c(0,
                                                          ceiling(
                                                            max(
                                                              indel_reads_cell$count + 
                                                                indel_reads_cell$se))), 
                                                 width = unit(2, "cm")), 
                            sgname2 = anno_empty(width = unit(2, "cm"), border = F),
                            show_annotation_name = F)
col_acc <- colSums(plot_mat)
col_acc[length(plot_region)] <- NA
top_anno <- HeatmapAnnotation(acc = anno_barplot(col_acc, 
                                                 border = T,
                                                 gp = gpar(fill = "#494490"), 
                                                 bar_width = 0.5,
                                                 width = unit(2, "cm")),
                              pos = anno_empty(border = F, width = unit(2, "mm")))

pdf("~/data/project/ear_project/gene_therapy_ll/Result/118_pair_indel_table_v2.pdf", 
    width = 11, height = 44)
Heatmap(plot_mat, 
        name = "heatmap",
        cluster_rows = F,
        cluster_columns = F, 
        top_annotation = top_anno,
        column_split = c(rep("",length(plot_region) - 1), " ", rep("  ", length(plot_region) - 1)),
        col = c("white", "white"), 
        left_annotation = left_anno,
        right_annotation =right_anno,
        cell_fun = function(j, i, x, y, w, h, fill){
          v <- pindex(plot_mat, i, j)
          if(j == length(plot_region)){
            grid.rect(x, y, w, h, gp = gpar(fill = cor_color(v), 
                                            col = "#AAAAAA"))
          }
          else{
            grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                            col = "#AAAAAA"))
          }
        }, show_heatmap_legend = F)

decorate_annotation("left", slice = 1, {
  od = nrow(plot_mat) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - indel_reads_tissue$count[od],seq_along(od), 
                xcale - (indel_reads_tissue$count[od] + indel_reads_tissue$se[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (indel_reads_tissue$count[od] + indel_reads_tissue$se[od]), 
                seq_along(od) - 0.2, 
                xcale - (indel_reads_tissue$count[od] + indel_reads_tissue$se[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  
})
decorate_annotation("right",slice = 1, {
  od = nrow(plot_mat) : 1
  grid.segments(indel_reads_cell$count[od],seq_along(od), 
                indel_reads_cell$count[od] + indel_reads_cell$se[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(indel_reads_cell$count[od] + indel_reads_cell$se[od], 
                seq_along(od) - 0.2, 
                indel_reads_cell$count[od] + indel_reads_cell$se[od], 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
})

indel_text <- c(paste0("+", plot_region[plot_region > 0]), as.character(plot_region[plot_region < 0]))
indel_text <- gtools::mixedsort(indel_text, decreasing = T)
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
decorate_annotation("pos",slice = 3, {
  
  tg <- textGrob(gt_render(as.character(rev(indel_text))), 
                 rot = 0, 
                 x = unit(c(1 : length(indel_text)) / length(indel_text) -
                            0.5 / length(indel_text), "npc"),
                 y=unit(0.5, "npc"), hjust = 0.5, 
                 gp = gpar(fontsize = 14, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("pos",slice = 2, {
  
  tg <- textGrob(gt_render("Cor"), 
                 rot = 0, 
                 x = unit(0.5, "npc"),
                 y=unit(0.5, "npc"), hjust = 0.5, 
                 gp = gpar(fontsize = 14, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("sgname1",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(rev(total_pair_sg_cor$name))), 
                 rot = 0, 
                 x = unit(0, "npc"),
                 y=unit(c(1 : nrow(total_pair_sg_cor)) / nrow(total_pair_sg_cor) -
                          0.5 / nrow(total_pair_sg_cor), "npc"), hjust = 0, 
                 gp = gpar(fontsize = 10, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("sgname2",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(rev(total_pair_sg_cor$name))), 
                 rot = 0, 
                 x = unit(0, "npc"),
                 y=unit(c(1 : nrow(total_pair_sg_cor)) / nrow(total_pair_sg_cor) -
                          0.5 / nrow(total_pair_sg_cor), "npc"), hjust = 0, 
                 gp = gpar(fontsize = 10, col = "black"))
  grid.draw(tg)
  invisible(tg)
})

dev.off()

