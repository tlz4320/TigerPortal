setwd("~/data/project/ear_project/gene_therapy_qsw/mOTOF/Result/")
samples <- list.files(pattern = "CRISPResso_on")
samples <- samples[-grep('[.]', samples)]
otof_result <- list()
for(sample in samples){
  setwd(sample)
  filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
  filename <- filename[which.max(str_length(filename))]
  edit_table <- read.table(filename, 
                           sep="\t", header = T, comment.char = "")
  otof_result[[str_remove(sample, "CRISPResso_on_")]] <- edit_table
  setwd("..")
}

otof_wt_seqs <- lapply(names(otof_result), function(id){
  edit_table <- otof_result[[id]]
  edit_table <- edit_table[edit_table$Unedited == "True",]
  getWtSeq(edit_table)
})
names(otof_wt_seqs) <- names(otof_result)

control_result <- otof_result$`m-OTOF-HO-BM`
control_result <- control_result[control_result$n_inserted + control_result$n_deleted != 0,]
control_result$id <- paste(control_result$Aligned_Sequence, control_result$Reference_Sequence,sep = "-")
control_result <- control_result[!duplicated(control_result$id),]
rownames(control_result) <- control_result$id
otof_result <- otof_result[-1]
otof_result_indel <- lapply(otof_result, function(x){
  indel <- x$n_inserted + x$n_deleted
  x <- x[indel != 0 & indel < 15,]
  id <- paste(x[,1],x[,2], sep = "-")
  for(i in 1 : nrow(x)){
    if(id[i] %in% control_result$id){
      x[i, 8] <- x[i, 8] - control_result[id[i],8]
    }
  }
  x <- x[x[,8] >= 0,]
  x$Pct <- x[,8]/ sum(x[,8])* 100
  x
})

plot_data <- lapply(names(otof_result_indel)[-1], function(id){
  x <- otof_result_indel[[id]]
  insert <- x[x$n_inserted != 0,]
  delete <- x[x$n_deleted != 0,]
  if(nrow(insert) != 0){
    insert <- lapply(split(insert$Pct, insert$n_inserted), sum)
    insert <- data.frame(indel = as.numeric(names(insert)), Pct = unlist(insert))
  } else{
    insert <- data.frame(indel = 0, Pct = 0)
  }
  if(nrow(delete) != 0){
    delete <- lapply(split(delete$Pct, delete$n_deleted), sum)
    delete <- data.frame(indel = -as.numeric(names(delete)), Pct = unlist(delete))
  } else{
    delete <- data.frame(indel = 0, Pct = 0)
  }
  result <- data.frame(rbind(insert, delete))
  result$sample <- id
  result[result$indel != 0,]
})

plot_data <- data.frame(do.call(rbind, plot_data))
plot_data <- lapply(split(plot_data,plot_data$indel), function(x){
  if(nrow(x) == 1){
    x$se <- 0
    return(x)
  }
  x$se <- plotrix::std.error(x$Pct)
  x$Pct <- mean(x$Pct)
  x[1,]
})
plot_data <- data.frame(do.call(rbind, plot_data))
range <- -18 : 18
range <- range[!range %in% plot_data$indel]
plot_data <- data.frame(rbind(plot_data, data.frame(indel = range, 
                                                    Pct = 0, 
                                                    sample = "",
                                                    se = 0)))
plot_data <- plot_data[order(plot_data$indel),]
plot_data <- plot_data[plot_data$indel %in% as.character(-18:18),]
plot_data$se2 <- plot_data$se
plot_data$se2[plot_data$se < 0.01] <- NA
plot_data$Pct2 <- plot_data$Pct
plot_data$Pct2[plot_data$se < 0.01] <- NA
plot_data$fill <- "#C9C9C9"
plot_data$fill[plot_data$indel == "1"] <- "#C92126"
plot_data$indel <- as.numeric(as.character(plot_data$indel))
x_label <- rep("", 37)
x_label[seq(1, length(x_label), 6)] <- as.character(seq(-18, 18, 6))
stick_len <- rep(0.1, 37)
stick_len[seq(1, length(x_label), 6)] <- 0.2
library(ggh4x)
# pdf("~/Nutstore Files/Tobin/Merged1NT/OTOF_virus_edit_pattern.pdf", width = 5, height = 5.4)
ggplot(plot_data,aes(x = indel, y = Pct)) + 
  geom_errorbar(aes(ymin=Pct2-se2, ymax=Pct2+se2), width=.5, na.rm = T, size = 0.3) + 
  geom_bar(stat="identity", fill =plot_data$fill, color = "black", size = 0.3, width = 0.7) + 
  scale_y_continuous(breaks = c(0, 10, 20,30,40,50,60,70), expand = c(0,0)) + 
  scale_x_continuous(minor_breaks = c(-18 : 18)[-seq(1, length(x_label), 6)],
                   breaks = seq(-18, 18, 6),
                   guide = "axis_minor", expand = c(0,0)) + 
  theme_classic2() + theme(panel.grid = element_blank(), 
                           axis.text.x = element_text(size = 14), 
                           axis.text.y = element_text(size = 14), 
                           axis.ticks.length.x=unit(.2, "cm"))
# dev.off()


#remove control
otof_result_indel <- lapply(otof_result, function(x){
  indel <- x$n_inserted + x$n_deleted
  x <- x[indel != 0 & indel < 15,]
  id <- paste(x[,1],x[,2], sep = "-")
  for(i in 1 : nrow(x)){
    if(id[i] %in% control_result$id){
      x[i, 8] <- x[i, 8] - control_result[id[i],8]
    }
  }
  x <- x[x[,8] >= 0,]
  x$Pct <- x[,8]/ sum(x[,8])* 100
  x
})

#dont remove control
otof_result_indel <- lapply(otof_result, function(x){
  indel <- x$n_inserted + x$n_deleted
  x <- x[indel != 0 & indel < 15,]
  id <- paste(x[,1],x[,2], sep = "-")
  x$Pct <- x[,8]/ sum(x[,8])* 100
  x
})

#ins_result
ins_result <- lapply(otof_result_indel, function(x){
  tmp <- x[x$n_inserted == 1,]
  reframe <- x[(x$n_inserted - x$n_deleted) %% 3== 1, ]
  tmp$Pct <- tmp[,8] / sum(reframe[,8]) * 100
  tmp
})
ins_result <- data.frame(do.call(rbind, ins_result))
ins_result_mean <- lapply(split(ins_result, ins_result$Aligned_Sequence), function(x){
  x$se <- plotrix::std.error(x$Pct)
  x$Pct <- mean(x$Pct)
  x[1,]
})
ins_result_mean <- data.frame(do.call(rbind, ins_result_mean))
ins_result_mean <- ins_result_mean[order(ins_result_mean$Pct, decreasing = T),]
ins_result_mean <- ins_result_mean[1:4,]
ins_result_mean$Reference_Sequence[4] <- "CTCCCCGAGGGCGTGCCCCCGAA-CGGCAGTGGGCACGGT"
seq_matrix <- matrix("", nrow = 4, ncol = 40)
for(i in 1 : nrow(ins_result_mean)){
  seq <- unlist(strsplit(ins_result_mean$Aligned_Sequence[i], "*"))
  seq_matrix[i,] <- seq
}
insert_matrix <- matrix(0, nrow = 4, ncol = 40)

for(i in 1 : nrow(ins_result_mean)){
  seq <- unlist(strsplit(ins_result_mean$Reference_Sequence[i], "*"))
  seq <- which(seq == '-')
  seq <- seq[which.min(abs(seq - 20))]
  insert_matrix[i,seq] <- 1
}

col_stat <-  lapply(otof_result_indel, function(x){
  tmp <- x[x[,1] %in% ins_result_mean$Aligned_Sequence,]
  tmp <- tmp[order(tmp[,1]),]
  tmp2 <- ins_result_mean[order(ins_result_mean[,1]),]
  tmp[,2] <- tmp2[,2]
  #tmp$Pct <- tmp[,8] / sum(tmp[,8]) * 100
  tmp$pos <- unlist(lapply(tmp$Reference_Sequence, function(x){
    x <- unlist(strsplit(x, "*"))
    x <- which(x == "-")
    x[which.min(abs(x - 20))]
  }))
  tmp <- lapply(split(tmp,tmp$pos), function(x){
    x$Pct <- sum(x$Pct)
    x[1,]
  })
  data.frame(do.call(rbind,tmp))
})


col_stat <- data.frame(do.call(rbind, col_stat))
col_stat <- lapply(split(col_stat,col_stat$pos),function(x){
  x$se <- plotrix::std.error(x$Pct)
  x$Pct <- mean(x$Pct)
  x[1,]
})
col_stat <- data.frame(do.call(rbind, col_stat))
col_stat <- col_stat[order(col_stat$Pct, decreasing = T),]
col_stat_data <- data.frame(Pct = rep(NA, 40), se = NA)
col_stat_data$Pct[col_stat$pos] <- col_stat$Pct
col_stat_data$se[col_stat$pos] <- col_stat$se
col_stat_data$pos <- rownames(col_stat_data)
openxlsx::write.xlsx(na.omit(col_stat_data), file= "~/data/project/ear_project/gene_therapy_qsw/mOTOF/Result/OTOF_ins1_col_stat_v2.xlsx", rowNames=F)
library(ComplexHeatmap)
fake_mat <- matrix(runif(nrow(seq_matrix) * ncol(seq_matrix)), nrow = nrow(seq_matrix))
top_anno <- HeatmapAnnotation(
  seq1 = anno_empty(height = unit(0.5, "cm"), border = F), 
  seq2 = anno_empty(height = unit(0.5, "cm"), border = F),
  pct = anno_barplot(height = unit(1, "cm"), col_stat_data$Pct, border = T, 
                                                 gp = gpar(fill = "#CA2D26"),
                                                 ylim = c(0, 50)), 
                              show_annotation_name = F)
right_anno <- rowAnnotation(pct2 = anno_barplot(ins_result_mean$Pct, border = T, 
                                                gp = gpar(fill = "#CA2D26"), 
                                                axis_param = c("side"="top"), bar_width = 0.5,
                                                ylim = c(0,65), width = unit(2, "cm")), 
                            show_annotation_name = F)
tmp <- ins_result_mean[,c("Pct", "se")]
openxlsx::write.xlsx(tmp, 
                     file= "~/data/project/ear_project/gene_therapy_qsw/mOTOF/Result/OTOF_ins1_row_stat_v2.xlsx", rowNames=F)

ht <- Heatmap(fake_mat, rect_gp = gpar(col = "gray"), col = c("white", "white"), 
        cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, 
        top_annotation = top_anno,
        right_annotation = right_anno,
        cell_fun = function(j, i, x, y, w, h, fill) {
          nt = pindex(seq_matrix, i, j)
          ins <- pindex(insert_matrix, i, j)
          if(ins != 0){
            col <- "#CA2D26"
            grid.rect(x , y, w, h, gp = gpar(col = col, fill = "white"))
          }
          if(j %in% c(25 : 27)){
            grid.text(nt, x, y,  gp = gpar(fontsize = 10, col = "#CA2D26"))
          }
          else{
            grid.text(nt, x, y,  gp = gpar(fontsize = 10, col = "black"))
          }

        }
        )
pdf("~/Nutstore Files/Tobin/Merged1NT/OTOF_insert1_stat_v3.pdf", width = 8, height = 1.8)
draw(ht)
decorate_annotation("pct",slice = 1, {
  od = 1 : 40
  grid.segments(seq_along(od), col_stat_data$Pct[od],
                seq_along(od),  col_stat_data$Pct[od] + col_stat_data$se[od], 
                default.units = "native", gp = gpar(col = "#CA2D26"))
  grid.segments(seq_along(od) - 0.1,col_stat_data$Pct[od] + col_stat_data$se[od],  
                seq_along(od) + 0.1, col_stat_data$Pct[od] + col_stat_data$se[od], 
                default.units = "native", gp = gpar(col = "#CA2D26")
                )
})
decorate_annotation("pct2", slice = 1, {
  od = 4:1

  grid.segments(ins_result_mean$Pct[od],seq_along(od), 
                ins_result_mean$Pct[od] + ins_result_mean$se[od],seq_along(od), default.units = "native", 
                gp = gpar(col = "#CA2D26"))
  grid.segments(ins_result_mean$Pct[od] + ins_result_mean$se[od], seq_along(od) - 0.1, 
                ins_result_mean$Pct[od] + ins_result_mean$se[od], seq_along(od) + 0.1, default.units = "native", 
                gp = gpar(col = "#CA2D26"))

})
seq_col <- rep("black", 40)
seq_col[25 :27] <- "#CA2D26"
seq_col[21] <- "#D21FF5"
decorate_annotation("seq1", {
  tg <- richtext_grob(gt_render(unlist(strsplit("CTCCCCGAGGGCGTGCCCCCCGAACGGCAGTGGGCACGGT", "*"))), 
                      rot = 0, 
                      x = unit(c(1 : 40) / 40 - 0.5 / 40, "npc"),
                      y=unit(0.5, "npc"), hjust = 0.5, 
                      gp = gpar(fontsize = 14, col = seq_col))
  grid.draw(tg)
  invisible(tg)
})

decorate_annotation("seq2", {

  tg <- textGrob(gt_render(unlist(strsplit("CTCCCCGAGGGCGTGCCCCC-GAACGGCAGTGGGCACGGT", "*"))), 
                      rot = 0, 
                      x = unit(c(1 : 40) / 40 - 0.5 / 40, "npc"),
                      y=unit(0.5, "npc"), hjust = 0.5, 
                      gp = gpar(fontsize = 14, col = seq_col))
  grid.draw(tg)
  invisible(tg)
})
dev.off()
 


#delete result
del_result <- lapply(otof_result_indel, function(x){
  tmp <- x[x$n_deleted == 2,]
  reframe <- x[(x$n_inserted - x$n_deleted) %% 3== 1, ]
  tmp$Pct <- tmp[,8] / sum(reframe[,8]) * 100
  #tmp$Pct <- tmp[,8] / sum(tmp[,8]) * 100
  tmp
})
del_result <- data.frame(do.call(rbind, del_result))
del_result_mean <- lapply(split(del_result, del_result$Aligned_Sequence), function(x){
  x$se <- plotrix::std.error(x$Pct)
  x$Pct <- mean(x$Pct)
  x[1,]
})
del_result_mean <- data.frame(do.call(rbind, del_result_mean))
del_result_mean <- del_result_mean[order(del_result_mean$Pct, decreasing = T),]
del_result_mean <- del_result_mean[1:2,]
seq_matrix <- matrix("", nrow = 2, ncol = 40)
for(i in 1 : nrow(del_result_mean)){
  seq <- unlist(strsplit(del_result_mean$Aligned_Sequence[i], "*"))
  seq_matrix[i,] <- seq
}
delete_matrix <- matrix(0, nrow = 2, ncol = 40)

for(i in 1 : nrow(del_result_mean)){
  seq <- unlist(strsplit(del_result_mean$Aligned_Sequence[i], "*"))
  seq <- which(seq == '-')
  delete_matrix[i,seq] <- 1
}

col_stat <-  lapply(otof_result_indel, function(x){
  tmp <- x[x[,1] %in% del_result_mean$Aligned_Sequence,]
  # tmp$Pct <- tmp[,8] / sum(tmp[,8]) * 100
  res <- rep(0, 40)
  for(i in 1 : nrow(tmp)){
    x <- tmp$Aligned_Sequence[i]
    x <- unlist(strsplit(x, "*"))
    x <- which(x == "-")
    for(pos in x){
      res[pos] <- res[pos] + tmp$Pct[i]
    }
  }
  res <- data.frame(Pct = res, pos = 1 : 40)
})


col_stat <- data.frame(do.call(rbind, col_stat))
col_stat <- lapply(split(col_stat,col_stat$pos),function(x){
  x$se <- plotrix::std.error(x$Pct)
  x$Pct <- mean(x$Pct)
  x[1,]
})
col_stat <- data.frame(do.call(rbind, col_stat))
col_stat <- col_stat[order(col_stat$pos, decreasing = F),]
col_stat$Pct[col_stat$se == 0] <- NA
col_stat$se[col_stat$se == 0] <- NA

openxlsx::write.xlsx(na.omit(col_stat)[,c(1,3,2)], file= "~/data/project/ear_project/gene_therapy_qsw/mOTOF/Result/OTOF_del2_col_stat_v2.xlsx", rowNames=F)
tmp <- del_result_mean[,c("Pct", "se")]
openxlsx::write.xlsx(tmp, 
                     file= "~/data/project/ear_project/gene_therapy_qsw/mOTOF/Result/OTOF_del2_row_stat_v2.xlsx", rowNames=F)


library(ComplexHeatmap)
top_anno <- HeatmapAnnotation(
  pct = anno_barplot(height = unit(1, "cm"), col_stat$Pct, border = T, 
                     gp = gpar(fill = "#5A962F"),
                     ylim = c(0, 7)), 
  show_annotation_name = F)
right_anno <- rowAnnotation(pct2 = anno_barplot(del_result_mean$Pct, border = T, 
                                                gp = gpar(fill = "#5A962F"), 
                                                axis_param = c("side"="top"), bar_width = 0.5,
                                                ylim = c(0,10), width = unit(2, "cm")), 
                            show_annotation_name = F)
fake_mat <- matrix(runif(nrow(seq_matrix) * ncol(seq_matrix)), nrow = nrow(seq_matrix))
ht <- Heatmap(fake_mat, rect_gp = gpar(col = "gray"), col = c("white", "white"), 
              cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, 
              top_annotation = top_anno,
              right_annotation = right_anno,
              cell_fun = function(j, i, x, y, w, h, fill) {
                nt = pindex(seq_matrix, i, j)
                del <- pindex(delete_matrix, i, j)
                if(del != 0){
                  col <- "#5A962F"
                  grid.rect(x , y, w, h, gp = gpar(col = col, fill = "white"))
                }
                if(j %in% c(24 : 26)){
                  grid.text(nt, x, y,  gp = gpar(fontsize = 10, col = "#CA2D26"))
                }
                else{
                  grid.text(nt, x, y,  gp = gpar(fontsize = 10, col = "black"))
                }
                
              }
)
pdf("~/Nutstore Files/Tobin/Merged1NT/OTOF_delete2_stat_v3.pdf", width = 8, height = 1)
draw(ht)
decorate_annotation("pct",slice = 1, {
  od = 1 : 40
  grid.segments(seq_along(od), col_stat$Pct[od],
                seq_along(od),  col_stat$Pct[od] + col_stat$se[od], 
                default.units = "native", gp = gpar(col = "#5A962F"))
  grid.segments(seq_along(od) - 0.1,col_stat$Pct[od] + col_stat$se[od],  
                seq_along(od) + 0.1, col_stat$Pct[od] + col_stat$se[od], 
                default.units = "native", gp = gpar(col = "#5A962F")
  )
})
decorate_annotation("pct2", slice = 1, {
  od = 2:1
  
  grid.segments(del_result_mean$Pct[od],seq_along(od), 
                del_result_mean$Pct[od] + del_result_mean$se[od],seq_along(od), default.units = "native", 
                gp = gpar(col = "#5A962F"))
  grid.segments(del_result_mean$Pct[od] + del_result_mean$se[od], seq_along(od) - 0.1, 
                del_result_mean$Pct[od] + del_result_mean$se[od], seq_along(od) + 0.1, default.units = "native", 
                gp = gpar(col = "#5A962F"))
  
})
dev.off()



####indel circle

otof_indel_circle <- lapply(otof_result, function(x){
  indel <- x$n_inserted + x$n_deleted
  x <- x[indel != 0 & indel < 15,]
  id <- paste(x[,1],x[,2], sep = "-")
  x$Pct <- x[,8]/ sum(x[,8])* 100
  indel <- x$n_inserted - x$n_deleted
  x_3 <- x[indel %% 3 == 0,]
  x_1 <- x[indel %% 3 == 1,]
  x_2 <- x[indel %% 3 == 2,]
  data.frame(type = c("3N", "3N + 1", "3N + 2"), 
             Pct = c(sum(x_3$Pct),
                     sum(x_1$Pct),
                     sum(x_2$Pct)))
})
otof_indel_circle <- data.frame(do.call(rbind, otof_indel_circle))
otof_indel_circle <- lapply(split(otof_indel_circle, otof_indel_circle$type), function(x){
  x$se <- plotrix::std.error(x$Pct)
  x$Pct <- mean(x$Pct)
  x[1,]
})
otof_indel_circle <- data.frame(do.call(rbind, otof_indel_circle))
otof_indel_circle$Pct_norm <- otof_indel_circle$Pct / sum(otof_indel_circle$Pct)

otof_indel_circle2 <- lapply(otof_result, function(x){
  indel <- x$n_inserted + x$n_deleted
  x <- x[indel != 0 & indel < 15,]
  id <- paste(x[,1],x[,2], sep = "-")
  x$Pct <- x[,8]/ sum(x[,8])* 100
  x$indel <- x$n_inserted - x$n_deleted
  x_1 <- x[x$indel %% 3 == 1,]
  res <- lapply(split(x_1$Pct, x_1$indel), sum)
  data.frame(type= names(res), Pct = unlist(res))
})
otof_indel_circle2 <- data.frame(do.call(rbind, otof_indel_circle2))
otof_indel_circle2 <- lapply(split(otof_indel_circle2, otof_indel_circle2$type), function(x){
  x$se <- plotrix::std.error(x$Pct)
  x$Pct <- mean(x$Pct)
  x[1,]
})
otof_indel_circle2 <- data.frame(do.call(rbind, otof_indel_circle2))
otof_indel_circle2$Pct_norm <- otof_indel_circle2$Pct / sum(otof_indel_circle2$Pct)
otof_indel_circle2_sel <- otof_indel_circle2[otof_indel_circle2$type %in% c("-2", "1"),]
otof_indel_circle2 <- data.frame(rbind(otof_indel_circle2_sel, 
                                       data.frame(type= "Other", Pct = 0, se = 0, Pct_norm = 1 - sum(otof_indel_circle2_sel$Pct_norm))))
otof_indel_circle2$Pct_norm <- otof_indel_circle2$Pct_norm / sum(otof_indel_circle2$Pct_norm) * otof_indel_circle$Pct_norm[2]
otof_indel_circle2 <- otof_indel_circle2[order(otof_indel_circle2$Pct_norm, decreasing = T),]
otof_indel_circle2_pct <- c(otof_indel_circle$Pct_norm[1], otof_indel_circle2$Pct_norm, otof_indel_circle$Pct_norm[3])
library(circlize)
pdf("~/Nutstore Files/Tobin/Merged1NT/OTOF_circlize.pdf", width = 6, height = 4)
circos.par("start.degree" = 100, "canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2))
circos.initialize(c(1,2,3), sector.width = otof_indel_circle$Pct_norm, xlim = c(0,1)) # 'a` just means there is one sector

circos.track(ylim = c(0,2), track.height = 0.4, bg.col = c("#535353","black", "#868686"))


par(new = TRUE) # <- magic
circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3), "start.degree" = 100)
circos.initialize(c(1,2,3,4,5), xlim = c(0, 1), sector.width =otof_indel_circle2_pct)
circos.track(ylim = c(0,2), track.height = 0.2, bg.col = c("white","#CA2D26", "#5A962F","black", "white"))
circos.clear()
dev.off()

stat_data <- data.frame(type = c(otof_indel_circle$type, otof_indel_circle2$type), 
                        Pct = c(otof_indel_circle$Pct_norm, otof_indel_circle2$Pct_norm))
openxlsx::write.xlsx(stat_data, file="~/Nutstore Files/Tobin/Merged1NT/OTOF_circlize_stat_v2.xlsx")
