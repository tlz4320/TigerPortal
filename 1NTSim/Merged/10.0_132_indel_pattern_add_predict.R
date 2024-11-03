load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/merged_indel_table_first_second.Rda")
###去掉那些不匹配的样本
sgCmp_sel <- sgCmp[!sgCmp$id2 %in%c("Sg_21_144","Sg_6_37","Sg_6_38","Sg_7_45","Sg_17_115"), ]
sgRNA_pair <- split(sgCmp_sel$id2, sgCmp_sel$id)

sgRNA_pair_rm <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(total_mean_first_second)) != 2)
}))]
sgRNA_pair_remain <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(total_mean_first_second)) == 2)
}))]
# sgCmp132 <- sgCmp[sgCmp$id2 %in% unlist(sgRNA_pair_remain),]
# save(sgCmp132 , file = "~/data/project/ear_project/gene_therapy_ll/Result/sgCmp132.rda")
c("Sg_21_144","Sg_6_37","Sg_6_38","Sg_7_45","Sg_17_115") %in% unlist(sgRNA_pair_remain)
###计算均值的相关性，因为本次数据每次样本是分开测得 不好AA BB CC一一对比了
max_region <- c(-200, 30)
total_pair_sg_cor <- lapply(sgRNA_pair_remain, function(y){
  ids <- unlist(y)
  sel_sg1 <- total_mean_first_second[[ids[1]]]
  sel_sg2 <- total_mean_first_second[[ids[2]]]
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
names(total_pair_sg_cor) <- unlist(lapply(sgRNA_pair_remain, function(x){
  paste(x, collapse = "-")
}))
total_pair_sg_cor <- data.frame(pair = names(total_pair_sg_cor), 
                                cor = unlist(lapply(total_pair_sg_cor, function(x){x$estimate})),
                                pval = unlist(lapply(total_pair_sg_cor, function(x){x$p.value})))

sg_name <- data.frame(id = names(sgRNA_pair_remain), 
                      pair = unlist(lapply(sgRNA_pair_remain, function(x){
                        paste(x, collapse = "-")
                      })))
sg_name$pos <- unlist(lapply(sg_name$id, function(x){
  if(!is.na(str_match(x, "last"))){
    return(1)
  }
  unlist(strsplit(x, "[_-]"))[2]
}))
total_pair_sg_cor <- merge(total_pair_sg_cor, sg_name, by = "pair")

total_pair_sg_cor$pos <- as.integer(total_pair_sg_cor$pos)
plot_data <- total_pair_sg_cor[order(total_pair_sg_cor$pos, total_pair_sg_cor$cor),]
plot_data$pair <- factor(plot_data$pair, levels = plot_data$pair)
plot_data$pos <- factor(as.character(plot_data$pos), levels = as.character(unique(plot_data$pos)))
sample_counts <- data.frame(table(plot_data$pos))
table(plot_data$cor >= 0.6)


plot_region <- c(-8: 5)
with_bulge_sg <- unlist(lapply(sgRNA_pair_remain, function(x){
  x <- unlist(x)
  x <- sgCmp[sgCmp$id2 %in% x,]
  x$id2[x$isInsert == "Insert"]
}))
without_bulge_sg <- unlist(lapply(sgRNA_pair_remain, function(x){
  x <- unlist(x)
  x <- sgCmp[sgCmp$id2 %in% x,]
  x$id2[x$isInsert != "Insert"]
}))

mean_indel_table_sel_with_bulge <- lapply(with_bulge_sg, function(x){
  res <- total_mean_first_second[[x]]
  res[,3] <- res[,2] / sum(res[,2]) * 100
  res[res[,1] %in% plot_region,]
})
names(mean_indel_table_sel_with_bulge) <- with_bulge_sg
mean_indel_table_sel_without_bulge <- lapply(without_bulge_sg, function(x){
  res <- total_mean_first_second[[x]]
  res[,3] <- res[,2] / sum(res[,2]) * 100
  res[res[,1] %in% plot_region,]
})
names(mean_indel_table_sel_without_bulge) <- without_bulge_sg
total_pair_sg_cor$bulgeOne <- unlist(lapply(total_pair_sg_cor$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))
  x <- sgCmp[sgCmp$id2 %in% x,]
  x$id2[x$isInsert == "Insert"]
}))
total_pair_sg_cor$notBulgeOne <- unlist(lapply(total_pair_sg_cor$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))
  x <- sgCmp[sgCmp$id2 %in% x,]
  x$id2[x$isInsert != "Insert"]
}))
total_pair_sg_cor <- total_pair_sg_cor[order(total_pair_sg_cor$cor, decreasing = T),]

plot_mat <- matrix(0, ncol = length(plot_region) * 2 - 1, nrow = nrow(total_pair_sg_cor))
for(i in 1 : nrow(total_pair_sg_cor)){
  tmp1 <- mean_indel_table_sel_with_bulge[[total_pair_sg_cor$bulgeOne[i]]]
  tmp1 <- tmp1[tmp1[,1] != 0,]
  tmp2 <- mean_indel_table_sel_without_bulge[[total_pair_sg_cor$notBulgeOne[i]]]
  tmp2 <- tmp2[tmp2[,1] != 0,]
  plot_mat[i, 1 : (length(plot_region) - 1)] <- rev(tmp1[,3]) / 100
  plot_mat[i, length(plot_region)] <- total_pair_sg_cor$cor[i]
  plot_mat[i, (length(plot_region) + 1) : ncol(plot_mat)] <- tmp2[,3] / 100
}

library(ComplexHeatmap)
library(circlize)
max(unlist(plot_mat[,-length(plot_region)]))
indel_color <- colorRamp2(c(0, 1), c("white", "blue"))
cor_color <- colorRamp2(c(0, 1), c("white", "red"))

indel_reads_with_bulge <- lapply(with_bulge_sg, function(x){
  tmp <- total_indel_table_first_second_rev[[x]]
  tmp <- unlist(lapply(tmp, function(y){
    sum(y[y[,1] %in% plot_region & y[,1] != 0,2])
  }))
  data.frame(count = log(mean(tmp)), se = plotrix::std.error(log(tmp)))
})
indel_reads_with_bulge <- data.frame(do.call(rbind, indel_reads_with_bulge))
rownames(indel_reads_with_bulge) <- with_bulge_sg
indel_reads_with_bulge <- indel_reads_with_bulge[total_pair_sg_cor$bulgeOne,]

indel_reads_without_bulge <- lapply(without_bulge_sg, function(x){
  tmp <- total_indel_table_first_second_rev[[x]]
  tmp <- unlist(lapply(tmp, function(y){
    sum(y[y[,1] %in% plot_region & y[,1] != 0,2])
  }))
  data.frame(count = log(mean(tmp)), se = plotrix::std.error(log(tmp)))
})
indel_reads_without_bulge <- data.frame(do.call(rbind, indel_reads_without_bulge))
rownames(indel_reads_without_bulge) <- without_bulge_sg
indel_reads_without_bulge <- indel_reads_without_bulge[total_pair_sg_cor$notBulgeOne,]
left_anno <- rowAnnotation(
  sgname1 = anno_empty(width = unit(3, "cm"), border = F),
  left1 = anno_barplot(indel_reads_with_bulge$count, 
                      border = T,
                      gp = gpar(fill = "#494490"), 
                      bar_width = 0.5,
                      axis_param = list(direction = "reverse", 
                                        "side"="top"),
                      ylim = c(0,
                               ceiling(
                                 max(
                                   indel_reads_with_bulge$count + 
                                     indel_reads_with_bulge$se))), 
                      width = unit(2, "cm")), 
  show_annotation_name = F)

right_anno <- rowAnnotation(right1 = anno_barplot(indel_reads_without_bulge$count, 
                                                 border = T,
                                                 gp = gpar(fill = "#494490"), 
                                                 axis_param = c("side"="top"), 
                                                 bar_width = 0.5,
                                                 ylim = c(0,
                                                          ceiling(
                                                            max(
                                                              indel_reads_without_bulge$count + 
                                                                indel_reads_without_bulge$se))), 
                                                 width = unit(2, "cm")), 
                            sgname2 = anno_empty(width = unit(3, "cm"), border = F),
                            show_annotation_name = F)

left_name <- total_pair_sg_cor$bulgeOne
right_name <- total_pair_sg_cor$notBulgeOne

output_indel_table <- plot_mat * 100
output_indel_table <- data.frame(rbind(output_indel_table[, 1 : 13], 
                                       output_indel_table[,ncol(output_indel_table) : 15]))
colnames(output_indel_table) <- c(paste0("+", 5 : 1), paste0("-", 1 : 8))
output_indel_table$sgName <- c(left_name, right_name)
output_indel_table <- output_indel_table[,c(ncol(output_indel_table), 1 : (ncol(output_indel_table) -1))]
openxlsx::write.xlsx(output_indel_table, "~/Nutstore Files/Tobin/sTable/s2Table.xlsx", rowNames=F, colNames=T)
col_acc <- colSums(plot_mat)
col_acc[length(plot_region)] <- NA
col_mean <- colMeans(plot_mat) * 100
col_se <- apply(plot_mat, 2, function(x){
  plotrix::std.error(x * 100)
})
col_mean[length(plot_region)] <- NA
top_anno <- HeatmapAnnotation(acc = anno_barplot(col_mean, 
                                                 border = T,
                                                 gp = gpar(fill = "#494490"), 
                                                 bar_width = 0.5,
                                                 width = unit(2, "cm"), 
                                                 ylim = c(0, max(col_mean + col_se, na.rm = T))),
                              pos = anno_empty(border = F, width = unit(2, "mm")))
which60 <- min(which(total_pair_sg_cor$cor < 0.6))

h1 <- Heatmap(plot_mat, 
        name = "heatmap",
        column_gap = unit(2.5, "mm"),
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
            grid.rect(x, y, unit(1.4, "npc"), h, gp = gpar(fill = cor_color(v), 
                                                           col = "#AAAAAA"))
            grid.text(round(v, 2), x, y)
          }
          else{
            grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                            col = "#AAAAAA"))
          }
        }, show_heatmap_legend = F)



###indelphi
load("~/data/project/ear_project/gene_therapy_ll/Result/indelphi_132_result.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/132_forecast.rda")
max_region <- c(-200, 30)

total_pair_sg_cor$indelphi_bulge <- unlist(lapply(total_pair_sg_cor$bulgeOne, function(id){
  sel_sg1 <- total_mean_first_second[[id]]
  sel_sg2 <- indelphi_132_res[[id]]
  sel_sg2$indel_size <- as.numeric(sel_sg2$indel_size)
  sel_sg2 <- sel_sg2[,c(1,2)]
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
  cor(sel_sg1[,2], sel_sg2[,2])
}))


total_pair_sg_cor$indelphi_notbulge <- unlist(lapply(total_pair_sg_cor$notBulgeOne, function(id){
  sel_sg1 <- total_mean_first_second[[id]]
  sel_sg2 <- indelphi_132_res[[id]]
  sel_sg2$indel_size <- as.numeric(sel_sg2$indel_size)
  sel_sg2 <- sel_sg2[,c(1,2)]
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
  cor(sel_sg1[,2], sel_sg2[,2])
}))

total_pair_sg_cor$forecast_bulge <- unlist(lapply(total_pair_sg_cor$bulgeOne, function(id){
  sel_sg1 <- total_mean_first_second[[id]]
  sel_sg2 <- forecast_132_res[[id]]
  sel_sg2$indel_size <- as.numeric(sel_sg2$indel_size)
  sel_sg2 <- sel_sg2[,c(1,2)]
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
  cor(sel_sg1[,2], sel_sg2[,2])
}))

total_pair_sg_cor$forecast_notbulge <- unlist(lapply(total_pair_sg_cor$notBulgeOne, function(id){
  sel_sg1 <- total_mean_first_second[[id]]
  sel_sg2 <- forecast_132_res[[id]]
  sel_sg2$indel_size <- as.numeric(sel_sg2$indel_size)
  sel_sg2 <- sel_sg2[,c(1,2)]
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
  cor(sel_sg1[,2], sel_sg2[,2])
}))






plot_region <- c(-8: 5)


indel_table_indelphi_bulge <- lapply(total_pair_sg_cor$bulgeOne, function(x){
  res <- indelphi_132_res[[x]]
  res[res[,1] == "0",2] <- 0
  res[,3] <- res[,2] / sum(res[,2]) * 100
  res <- res[res[,1] %in% plot_region,]
  append_region <- plot_region[!plot_region %in% res[,1]]
  res <- data.frame(rbind(res, data.frame(indel_size = append_region, fq = 0, Pct = 0)))
  res <- res[gtools::mixedorder(res$indel_size, decreasing = T),]
})
names(indel_table_indelphi_bulge) <- total_pair_sg_cor$bulgeOne


plot_mat3 <- matrix(0, ncol = length(plot_region) - 1, nrow = nrow(total_pair_sg_cor))
for(i in 1 : nrow(total_pair_sg_cor)){
  tmp1 <- indel_table_indelphi_bulge[[i]]
  tmp1 <- tmp1[tmp1[,1] != 0,]
  plot_mat3[i, ] <- tmp1[,3] / 100
}


indel_table_indelphi_nobulge <- lapply(total_pair_sg_cor$notBulgeOne, function(x){
  res <- indelphi_132_res[[x]]
  res[res[,1] == "0",2] <- 0
  res[,3] <- res[,2] / sum(res[,2]) * 100
  res <- res[res[,1] %in% plot_region,]
  append_region <- plot_region[!plot_region %in% res[,1]]
  res <- data.frame(rbind(res, data.frame(indel_size = append_region, fq = 0, Pct = 0)))
  res <- res[gtools::mixedorder(res$indel_size),]
})
names(indel_table_indelphi_nobulge) <- total_pair_sg_cor$notBulgeOne


plot_mat4 <- matrix(0, ncol = length(plot_region) - 1, nrow = nrow(total_pair_sg_cor))
for(i in 1 : nrow(total_pair_sg_cor)){
  tmp1 <- indel_table_indelphi_nobulge[[i]]
  tmp1 <- tmp1[tmp1[,1] != 0,]
  plot_mat4[i, ] <- tmp1[,3] / 100
}


indel_table_forecast_bulge <- lapply(total_pair_sg_cor$bulgeOne, function(x){
  res <- forecast_132_res[[x]]
  res[res[,1] == "0",2] <- 0
  res[,3] <- res[,2] / sum(res[,2]) * 100
  res <- res[res[,1] %in% plot_region,]
  append_region <- plot_region[!plot_region %in% res[,1]]
  res <- data.frame(rbind(res, data.frame(indel_size = append_region, fq = 0, Pct = 0)))
  res <- res[gtools::mixedorder(res$indel_size, decreasing = T),]
})
names(indel_table_forecast_bulge) <- total_pair_sg_cor$bulgeOne


plot_mat5 <- matrix(0, ncol = length(plot_region) - 1, nrow = nrow(total_pair_sg_cor))
for(i in 1 : nrow(total_pair_sg_cor)){
  tmp1 <- indel_table_forecast_bulge[[i]]
  tmp1 <- tmp1[tmp1[,1] != 0,]
  plot_mat5[i, ] <- tmp1[,3] / 100
}


indel_table_forecast_nobulge <- lapply(total_pair_sg_cor$notBulgeOne, function(x){
  res <- forecast_132_res[[x]]
  res[res[,1] == "0",2] <- 0
  res[,3] <- res[,2] / sum(res[,2]) * 100
  res <- res[res[,1] %in% plot_region,]
  append_region <- plot_region[!plot_region %in% res[,1]]
  res <- data.frame(rbind(res, data.frame(indel_size = append_region, fq = 0, Pct = 0)))
  res <- res[gtools::mixedorder(res$indel_size),]
})
names(indel_table_forecast_nobulge) <- total_pair_sg_cor$notBulgeOne


plot_mat6 <- matrix(0, ncol = length(plot_region) - 1, nrow = nrow(total_pair_sg_cor))
for(i in 1 : nrow(total_pair_sg_cor)){
  tmp1 <- indel_table_forecast_nobulge[[i]]
  tmp1 <- tmp1[tmp1[,1] != 0,]
  plot_mat6[i, ] <- tmp1[,3] / 100
}

left_anno3 <- rowAnnotation(
  right = anno_barplot(total_pair_sg_cor$indelphi_bulge, 
                      border = T,
                      gp = gpar(fill = "red"), 
                      bar_width = 0.5,
                      axis_param = list("side"="top"), 
                      width = unit(2, "cm")), 
  show_annotation_name = F)

col_acc3 <- colSums(plot_mat3)
col_mean3 <- colMeans(plot_mat3) * 100
col_se3 <- apply(plot_mat3, 2, function(x){
  plotrix::std.error(x * 100)
})
top_anno3 <- HeatmapAnnotation(acc3 = anno_barplot(col_mean3, 
                                                  border = T,
                                                  gp = gpar(fill = "#494490"), 
                                                  bar_width = 0.5,
                                                  width = unit(2, "cm"), 
                                                  ylim = c(0, max(col_mean3 + col_se3, na.rm = T))),
                               pos3 = anno_empty(border = F, width = unit(2, "mm")))

ht3 <- Heatmap(plot_mat3, 
               name = "heatmap3",
               cluster_rows = F,
               cluster_columns = F, 
               top_annotation = top_anno3,
               col = c("white", "white"), 
               right_annotation = left_anno3,
               cell_fun = function(j, i, x, y, w, h, fill){
                 v <- pindex(plot_mat3, i, j)
                 
                 grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                                 col = "#AAAAAA"))
                 
               }, show_heatmap_legend = F)


left_anno4 <- rowAnnotation(
  right = anno_barplot(total_pair_sg_cor$indelphi_notbulge, 
                       border = T,
                       gp = gpar(fill = "red"), 
                       bar_width = 0.5,
                       axis_param = list(direction = "reverse", 
                         "side"="top"), 
                       width = unit(2, "cm")), 
  show_annotation_name = F)

col_acc4 <- colSums(plot_mat4)
col_mean4 <- colMeans(plot_mat4) * 100
col_se4 <- apply(plot_mat4, 2, function(x){
  plotrix::std.error(x * 100)
})
top_anno4 <- HeatmapAnnotation(acc4 = anno_barplot(col_mean4, 
                                                  border = T,
                                                  gp = gpar(fill = "#494490"), 
                                                  bar_width = 0.5,
                                                  width = unit(2, "cm"), 
                                                  ylim = c(0, max(col_mean4 + col_se4, na.rm = T))),
                               pos3 = anno_empty(border = F, width = unit(2, "mm")))

ht4 <- Heatmap(plot_mat4, 
               name = "heatmap4",
               cluster_rows = F,
               cluster_columns = F, 
               top_annotation = top_anno4,
               col = c("white", "white"), 
               left_annotation = left_anno4,
               cell_fun = function(j, i, x, y, w, h, fill){
                 v <- pindex(plot_mat4, i, j)
                 
                 grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                                 col = "#AAAAAA"))
                 
               }, show_heatmap_legend = F)




left_anno5 <- rowAnnotation(
  right = anno_barplot(total_pair_sg_cor$forecast_bulge, 
                       border = T,
                       gp = gpar(fill = "red"), 
                       bar_width = 0.5,
                       axis_param = list("side"="top"), 
                       width = unit(2, "cm")), 
  show_annotation_name = F)

col_acc5 <- colSums(plot_mat5)
col_mean5 <- colMeans(plot_mat5) * 100
col_se5 <- apply(plot_mat5, 2, function(x){
  plotrix::std.error(x * 100)
})
top_anno5 <- HeatmapAnnotation(acc5 = anno_barplot(col_mean5, 
                                                  border = T,
                                                  gp = gpar(fill = "#494490"), 
                                                  bar_width = 0.5,
                                                  width = unit(2, "cm"), 
                                                  ylim = c(0, max(col_mean5 + col_se5, na.rm = T))),
                               pos3 = anno_empty(border = F, width = unit(2, "mm")))

ht5 <- Heatmap(plot_mat5, 
               name = "heatmap5",
               cluster_rows = F,
               cluster_columns = F, 
               top_annotation = top_anno5,
               col = c("white", "white"), 
               right_annotation = left_anno5,
               cell_fun = function(j, i, x, y, w, h, fill){
                 v <- pindex(plot_mat5, i, j)
                 
                 grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                                 col = "#AAAAAA"))
                 
               }, show_heatmap_legend = F)


left_anno6 <- rowAnnotation(
  right = anno_barplot(total_pair_sg_cor$forecast_notbulge, 
                       border = T,
                       gp = gpar(fill = "red"), 
                       bar_width = 0.5,
                       axis_param = list(direction = "reverse", 
                                         "side"="top"), 
                       width = unit(2, "cm")), 
  show_annotation_name = F)

col_acc6 <- colSums(plot_mat6)
col_mean6 <- colMeans(plot_mat6) * 100
col_se6 <- apply(plot_mat6, 2, function(x){
  plotrix::std.error(x * 100)
})
top_anno6 <- HeatmapAnnotation(acc6 = anno_barplot(col_mean6, 
                                                  border = T,
                                                  gp = gpar(fill = "#494490"), 
                                                  bar_width = 0.5,
                                                  width = unit(2, "cm"), 
                                                  ylim = c(0, max(col_mean6 + col_se6, na.rm = T))),
                               pos3 = anno_empty(border = F, width = unit(2, "mm")))

ht6 <- Heatmap(plot_mat6, 
               name = "heatmap6",
               cluster_rows = F,
               cluster_columns = F, 
               top_annotation = top_anno6,
               col = c("white", "white"), 
               left_annotation = left_anno6,
               cell_fun = function(j, i, x, y, w, h, fill){
                 v <- pindex(plot_mat6, i, j)
                 
                 grid.rect(x, y, w, h, gp = gpar(fill = indel_color(v), 
                                                 col = "#AAAAAA"))
                 
               }, show_heatmap_legend = F)




pdf("~/Nutstore Files/Tobin/Merged1NT/132_pair_indel_table_v5.pdf", 
    width = 35, height = 42)
draw(ht5 + ht3 + h1 + ht4+ ht6)
decorate_annotation("left1", slice = 1, {
  od = nrow(plot_mat) : 1
  vp = current.viewport()
  xcale <- vp$xscale[2]
  grid.segments(xcale - indel_reads_with_bulge$count[od],seq_along(od), 
                xcale - (indel_reads_with_bulge$count[od] + indel_reads_with_bulge$se[od]),
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(xcale - (indel_reads_with_bulge$count[od] + indel_reads_with_bulge$se[od]), 
                seq_along(od) - 0.2, 
                xcale - (indel_reads_with_bulge$count[od] + indel_reads_with_bulge$se[od]), 
                seq_along(od) + 0.2, default.units = "native", 
                gp = gpar(col = "#494490"))
  
})
decorate_annotation("right1",slice = 1, {
  od = nrow(plot_mat) : 1
  
  grid.segments(indel_reads_without_bulge$count[od],seq_along(od), 
                indel_reads_without_bulge$count[od] + indel_reads_without_bulge$se[od],
                seq_along(od), default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(indel_reads_without_bulge$count[od] + indel_reads_without_bulge$se[od], 
                seq_along(od) - 0.2, 
                indel_reads_without_bulge$count[od] + indel_reads_without_bulge$se[od], 
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
                 gp = gpar(fontsize = 15, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("pos",slice = 3, {
  
  tg <- textGrob(gt_render(as.character(rev(indel_text))), 
                 rot = 0, 
                 x = unit(c(1 : length(indel_text)) / length(indel_text) -
                            0.5 / length(indel_text), "npc"),
                 y=unit(0.5, "npc"), hjust = 0.5, 
                 gp = gpar(fontsize = 15, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("pos",slice = 2, {
  
  tg <- textGrob(gt_render("Cor"), 
                 rot = 0, 
                 x = unit(0.5, "npc"),
                 y=unit(0.5, "npc"), hjust = 0.5, 
                 gp = gpar(fontsize = 15, col = "black"))
  grid.draw(tg)
  invisible(tg)
})

for(i in 1 : 3){
  decorate_heatmap_body("heatmap",slice = 1,column_slice = i, {
    grid.lines(y = unit(c(1 - (which60 - 1) / nrow(plot_data), 1 - (which60 - 1) / nrow(plot_data)), "npc"))
  })
}

decorate_annotation("sgname1",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(rev(left_name))), 
                 rot = 0, 
                 x = unit(0, "npc"),
                 y=unit(c(1 : length(left_name)) / length(left_name) -
                          0.5 / length(left_name), "npc"), hjust = 0, 
                 gp = gpar(fontsize = 13, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("sgname2",slice = 1, {
  
  tg <- textGrob(gt_render(as.character(rev(right_name))), 
                 rot = 0, 
                 x = unit(0, "npc"),
                 y=unit(c(1 : length(left_name)) / length(left_name) -
                          0.5 / length(left_name), "npc"), hjust = 0, 
                 gp = gpar(fontsize = 13, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
for(i in 1 : 3){
  decorate_annotation("acc", slice = i, {
    od = column_order(h1)[[i]]
    grid.segments(seq_along(od), col_mean[od],
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
}
decorate_annotation("acc3", slice = 1, {
  od = 1 : length(col_mean3)
  grid.segments(seq_along(od), col_mean3[od],
                seq_along(od), 
                col_mean3[od] + col_se3[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(seq_along(od) - 0.2, 
                col_mean3[od] + col_se3[od],
                seq_along(od) + 0.2, 
                col_mean3[od] + col_se3[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
})
decorate_annotation("acc4", slice = 1, {
  od = 1 : length(col_mean4)
  grid.segments(seq_along(od), col_mean4[od],
                seq_along(od), 
                col_mean4[od] + col_se4[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(seq_along(od) - 0.2, 
                col_mean4[od] + col_se4[od],
                seq_along(od) + 0.2, 
                col_mean4[od] + col_se4[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
})
decorate_annotation("acc5", slice = 1, {
  od = 1 : length(col_mean5)
  grid.segments(seq_along(od), col_mean5[od],
                seq_along(od), 
                col_mean5[od] + col_se5[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(seq_along(od) - 0.2, 
                col_mean5[od] + col_se5[od],
                seq_along(od) + 0.2, 
                col_mean5[od] + col_se5[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
})
decorate_annotation("acc6", slice = 1, {
  od = 1 : length(col_mean6)
  grid.segments(seq_along(od), col_mean6[od],
                seq_along(od), 
                col_mean6[od] + col_se6[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
  grid.segments(seq_along(od) - 0.2, 
                col_mean6[od] + col_se6[od],
                seq_along(od) + 0.2, 
                col_mean6[od] + col_se6[od],
                default.units = "native", 
                gp = gpar(col = "#494490"))
})

dev.off()

pdf("~/Nutstore Files/Tobin/Merged1NT/132_cor_compare_predict.pdf", width = 6, height = 4)
tmp <- melt(total_pair_sg_cor, id.vars = "pair", measure.vars = c("cor", "indelphi_bulge", "forecast_bulge", "indelphi_notbulge", "forecast_notbulge"))
ggplot(tmp) + geom_boxplot(aes(x = variable, y = value)) + ggtitle("Total 132")
tmp <- melt(total_pair_sg_cor[total_pair_sg_cor$cor > 0.6,], id.vars = "pair", measure.vars = c("cor", "indelphi_bulge", "forecast_bulge", "indelphi_notbulge", "forecast_notbulge"))
ggplot(tmp) + geom_boxplot(aes(x = variable, y = value)) + ggtitle("Cor > 0.6")
dev.off()

pdf("~/Nutstore Files/Tobin/Merged1NT/heatmap_legend.pdf", width = 3, height = 2)
pushViewport(viewport(width = 0.9, height = 0.9))
grid.rect()

indel_color <- colorRamp2(c(0, 1), c("white", "blue"))
cor_color <- colorRamp2(c(0, 1), c("white", "red"))
lgd = Legend(col_fun = indel_color, title = "Indel Color", direction = "horizontal")
lgd2 = Legend(col_fun = cor_color, title = "Cor Color", direction = "horizontal")

draw(lgd, x = unit(0.5, "npc"), y = unit(0.8, "npc"))
draw(lgd2, x = unit(0.5, "npc"), y = unit(0.3, "npc"))
popViewport()
dev.off()


