library(Biostrings)
library(stringr)
source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first_second.Rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt", sep="\t")
#上面部分每次跑都很慢，不如直接改成保存下来读取
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_ins_table_data.rda")
###取出wt序列，用于后面的统计
wt_seqs <- lapply(unique(seqs_table$id), function(id){
  edit_table <- total_edit_table_rev[[id]][[1]]
  edit_table <- edit_table[edit_table$Unedited == "True",]
  getWtSeq(edit_table)
})
names(wt_seqs) <- unique(seqs_table$id)


ins1_pct_table <- lapply(unique(seqs_table$id), function(id){
  ins1_table <- seqs_table[seqs_table$id == id,]
  wt_seq <- wt_seqs[[id]]
  tmp <- getStat2(ins1_table, T)
  tmp$id <- id
  tmp
})
names(ins1_pct_table) <- unique(seqs_table$id)
ins1_pct_table_pos <- lapply(ins1_pct_table, function(tmp){
  tmp <- lapply(split(tmp, tmp$pos), function(x){
    x$Pct <- sum(x$Pct)
    x[1,]
  })
  tmp <- do.call(rbind, tmp)
}) 

###取出还成对的样本
remain_samples <- lapply(split(sgCmp$id2, sgCmp$id), function(x){
  if(sum(unlist(x) %in% names(ins1_pct_table_pos)) == 2){
    return(x)
  }
  return("")
})
remain_samples <- unlist(remain_samples)
remain_samples <- remain_samples[remain_samples != ""]


bulge_pos_sel <- bulge_pos[bulge_pos$V4 %in% sgCmp$id[sgCmp$id2 %in% remain_samples],]
bulge_pos_sel <- bulge_pos_sel[order(bulge_pos_sel$V5, decreasing = T),]
bulge_pos_order <- unlist(lapply(bulge_pos_sel$V4, function(x){
  tmp <- sgCmp[sgCmp$id == x,]
  #把Ins1放在前面
  tmp <- tmp[order(tmp$isInsert, decreasing = T),]
  tmp$id2
}))

remain_ins1_pct_table <- ins1_pct_table[bulge_pos_order]
remain_ins1_pct_table_pos <- ins1_pct_table_pos[bulge_pos_order]

plot_mat <- matrix(runif(17 * length(remain_ins1_pct_table)), ncol = 17, nrow = length(remain_ins1_pct_table))
rownames(plot_mat) <- bulge_pos_order

for(id in names(remain_ins1_pct_table_pos)){
  tmp <- remain_ins1_pct_table_pos[[id]]
  tmp <- tmp[tmp$pos %in% c(15 : 23),]
  #首先把15-20的算一下，15等于是-6～-5之间插入
  #20则是-1～1之间插入
  #之后算21-23的  21也是-1～1之间插入  所以得分开算一下
  for(i in 1 : 6){
    edit_pos <- 14 + i
    if(edit_pos %in% tmp$pos){
      plot_mat[id, i * 2] <- plot_mat[id, i * 2]  + tmp$Pct[tmp$pos == edit_pos]
    }
  }
  for(i in 7 : 8){
    edit_pos <- 14 + i
    if(edit_pos %in% tmp$pos){
      plot_mat[id, i * 2 - 2] <- plot_mat[id, i * 2 - 2]  + tmp$Pct[tmp$pos == edit_pos]
    }
  }
}

nt_mats <- lapply(remain_ins1_pct_table, function(tmp){
  plot_mat <- matrix(0, ncol = 17, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(15 : 23),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
  tmp_left <- tmp[tmp$pos <= 20,]
  tmp_right <- tmp[tmp$pos > 20,]
  
  #首先把15-20的算一下，15等于是-6～-5之间插入
  #20则是-1～1之间插入
  #之后算21-23的  21也是-1～1之间插入  所以得分开算一下
  if(nrow(tmp_left) > 0){
    for(i in 1 : nrow(tmp_left)){
      mat_pos <- (tmp_left$pos[i] - 14) * 2
      plot_mat[tmp_left$NT[i], mat_pos] <- plot_mat[tmp_left$NT[i], mat_pos] + tmp_left$Pct[i]
    }
  }
  
  if(nrow(tmp_right) > 0){
    for(i in 1 : nrow(tmp_right)){
      mat_pos <- (tmp_right$pos[i] - 14) * 2 - 2
      plot_mat[tmp_right$NT[i], mat_pos] <- plot_mat[tmp_right$NT[i], mat_pos] + tmp_right$Pct[i]
    }
  }
  
  plot_mat
  
})




seq_mat <- matrix("", ncol = 17, nrow = length(remain_ins1_pct_table))
rownames(seq_mat) <- bulge_pos_order
for(id in bulge_pos_order){
  wt_seq <- wt_seqs[[id]][15 : 23]
  for(i in 1 : length(wt_seq)){
    seq_mat[id, (i - 1) * 2 + 1] <- wt_seq[i]
  }
}




sg_info <- data.frame(id = bulge_pos_order, pair = unlist(lapply(bulge_pos_order, function(x){
  sgCmp$id[sgCmp$id2== x]
})))
sg_info$bulge <- unlist(lapply(sg_info$pair, function(x){
  bulge_pos_sel$V5[bulge_pos_sel$V4== x]
}))
sg_info$pair <- factor(sg_info$pair, levels = unique(sg_info$pair))
sg_info$bulge_in_mat <- unlist(lapply(sg_info$bulge, function(pos){
  17 - (pos - 1) * 2
}))
library(ComplexHeatmap)
library(circlize)
top_anno <- HeatmapAnnotation(position = anno_empty(border = F, width = unit(20, "mm")))
left_anno <- rowAnnotation(sg = anno_empty(border = F, width = unit(20, "mm")))
insert_color <- colorRamp2(c(0, 100), c("white", "#FF00FF"))
####最终要去掉一些低质量的结果
###首先判断Ins的比例是不是都小于1% indel reads

sg_info$lowPct <- unlist(lapply(sg_info$id, function(id){
  tmp <- remain_ins1_pct_table[[id]]
  tmp <- tmp[tmp$pos %in% c(20:21),]
  all(tmp$Pct < 1)
}))
sg_info_sel <- sg_info[!sg_info$lowPct,]
#去掉低质量的序列
sg_info_sel <- sg_info_sel[!sg_info_sel$id %in% c("Sg_21_144","Sg_6_37","Sg_6_38","Sg_7_45","Sg_17_115"),]
sg_info_sel <- sg_info_sel[sg_info_sel$pair %in% unlist(lapply(as.character(sg_info_sel$pair), function(x){
  if(sum(sg_info_sel$pair == x) == 2){
    return(x)
  }
  return("")
})),]

seq_mat_sel <- seq_mat[sg_info_sel$id,]
plot_mat_sel <- plot_mat[sg_info_sel$id,]
left_anno <- rowAnnotation(line = anno_empty(border = F, width = unit(2, "mm")), 
                           rowname = anno_text(x = rownames(plot_mat_sel)))
right_anno <- rowAnnotation(line2 = anno_empty(border = F, width = unit(2, "mm")))






###碱基Ins1与左右碱基关系图

remain_ins1_pct_table_sel <- remain_ins1_pct_table[sg_info_sel$id]
nt_counts_mats <- lapply(remain_ins1_pct_table_sel, function(tmp){
  plot_mat <- matrix(0, ncol = 17, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(20, 21),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  tmp_left <- tmp[tmp$pos <= 20,]
  tmp_right <- tmp[tmp$pos > 20,]
  
  #首先把15-20的算一下，15等于是-6～-5之间插入
  #20则是-1～1之间插入
  #之后算21-23的  21也是-1～1之间插入  所以得分开算一下
  if(nrow(tmp_left) > 0){
    for(i in 1 : nrow(tmp_left)){
      mat_pos <- (tmp_left$pos[i] - 14) * 2
      plot_mat[tmp_left$NT[i], mat_pos] <- plot_mat[tmp_left$NT[i], mat_pos] + tmp_left$Pct[i]
    }
  }
  
  if(nrow(tmp_right) > 0){
    for(i in 1 : nrow(tmp_right)){
      mat_pos <- (tmp_right$pos[i] - 14) * 2 - 2
      plot_mat[tmp_right$NT[i], mat_pos] <- plot_mat[tmp_right$NT[i], mat_pos] + tmp_right$Pct[i]
    }
  }
  
  plot_mat
  
})
#把有bulge和没有bulge的分开计算
besides_nt_stat_with_bulge <- list()

for(id in names(nt_counts_mats)[seq(1, length(nt_counts_mats), 2)]){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat[id,])
  for(i in 12){
    #还是过滤1一下的结果
    if(sum(tmp_pct_mat[,i]) > 1){
      leftnt <- tmp_seq_mat[i - 1]
      rightnt <- tmp_seq_mat[i + 1]
      if(!leftnt %in% names(besides_nt_stat_with_bulge)){
        besides_nt_stat_with_bulge[[leftnt]] <- list()
      }
      tmp_besides_nt_stat <- besides_nt_stat_with_bulge[[leftnt]]
      if(!rightnt %in% names(tmp_besides_nt_stat)){
        tmp_nt_stat <- c(0,0,0,0)
        names(tmp_nt_stat) <- c("A", "T", "C", "G")
        tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      }
      tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
      
      for(nt in rownames(tmp_counts_mat)){
        tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_counts_mat[nt, i]
      }
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      besides_nt_stat_with_bulge[[leftnt]] <- tmp_besides_nt_stat
    }
  }
}

besides_nt_stat_without_bulge <- list()

for(id in names(nt_counts_mats)[seq(2, length(nt_counts_mats), 2)]){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat[id,])
  for(i in 12){
    #还是过滤1一下的结果
    if(sum(tmp_pct_mat[,i]) > 1){
      leftnt <- tmp_seq_mat[i - 1]
      rightnt <- tmp_seq_mat[i + 1]
      if(!leftnt %in% names(besides_nt_stat_without_bulge)){
        besides_nt_stat_without_bulge[[leftnt]] <- list()
      }
      tmp_besides_nt_stat <- besides_nt_stat_without_bulge[[leftnt]]
      if(!rightnt %in% names(tmp_besides_nt_stat)){
        tmp_nt_stat <- c(0,0,0,0)
        names(tmp_nt_stat) <- c("A", "T", "C", "G")
        tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      }
      tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
      
      for(nt in rownames(tmp_counts_mat)){
        tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_counts_mat[nt, i]
      }
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      besides_nt_stat_without_bulge[[leftnt]] <- tmp_besides_nt_stat
    }
  }
}




fake_plot_mat <- matrix(runif(8*4), nrow = 8, ncol = 4)
rownames(fake_plot_mat) <- c("A","A2", "T","T2", "C","C2", "G","G2")
colnames(fake_plot_mat) <- c("A", "T", "C", "G")

ht <- Heatmap(fake_plot_mat, 
              name = "Name",
              column_split = NULL,
              column_title=NULL,
              row_title = NULL,
              row_split = rep(c("", " ", "  ", "   "), c(2,2,2,2)),
              show_column_names = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              top_annotation = HeatmapAnnotation(colname = anno_empty(height = unit(2, "cm"), border = F)),
              left_annotation = rowAnnotation(rowname = anno_empty(height = unit(2, "cm"), border = F)),
              right_annotation = NULL,
              show_row_names = F,
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                leftnt <- rownames(fake_plot_mat)[i]
                rightnt <- colnames(fake_plot_mat)[j]

                if(str_ends(leftnt, "2")){
                  leftnt <- str_sub(leftnt, 1, 1)
                  pct <- besides_nt_stat_without_bulge[[leftnt]][[rightnt]]
                } else {
                  pct <- besides_nt_stat_with_bulge[[leftnt]][[rightnt]]
                }
                labels <- names(pct)
                pct <- unlist(pct)
                col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                names(col) <- c("A", 'T','C','G')
                pct <- c(0, cumsum(pct)/sum(pct))
                dx <- diff(pct)
                nx <- length(dx)
                twopi <- 2 * pi
                init.angle <- 0
                t2xy <- function(t) {
                  t2p <- twopi * t + init.angle * pi/180
                  list(x = 0.8 * cos(t2p), y = 0.8 * sin(t2p))
                }
                maxx <- max(unlist(lapply(1 : nx, function(i){
                  n <- max(10, floor(200 * dx[i]))
                  P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                  max(abs(P$x))
                })))
                maxy <- max(unlist(lapply(1 : nx, function(i){
                  n <- max(10, floor(200 * dx[i]))
                  P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                  max(abs(P$y))
                })))
                facy <- as.numeric(h) / maxy * 0.3
                facx <- as.numeric(w) / maxx * 0.3
                for (i in 1L:nx) {
                  n <- max(10, floor(200 * dx[i]))
                  P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                  
                  grid.polygon(unit(as.numeric(x) + c(P$x * facx, 0), units = "npc"), 
                               unit(as.numeric(y) + c(P$y * facy, 0), units = "npc"), 
                               gp = gpar(fill = col[labels[i]]))
                  # print(P)
                  # polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = 45,
                  #         border = border[i], col = col[i], lty = lty[i])
                }
                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16, 
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)

pdf("~/Nutstore Files/Tobin/Merged1NT/merged_left_right_nt_pie_chart_v3_count.pdf", width = 6, height = 10.3)

draw(ht, annotation_legend_list = lgd_list)
decorate_annotation("colname", slice = 1, {
  tg <- richtext_grob(gt_render(c("A", "T", "C", "G")), 
                      rot = 0, 
                      x = unit(1 : 4 / 4 - 0.5 / 4, "npc"),
                      y=unit(0, "npc"), hjust = 0.5, 
                      gp = gpar(fontsize = 50, 
                                col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
  grid.draw(tg)
  invisible(tg)
})
for(i in 1 : 4){
  decorate_annotation("rowname", slice = i, {
    tg <- richtext_grob(gt_render(c("A", "T", "C", "G")[i]), 
                        rot = 0, 
                        x = unit(0.8, "npc"),
                        y=unit(0.5, "npc"), hjust = 0.5, 
                        gp = gpar(fontsize = 50, 
                                  col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")[i]))
    grid.draw(tg)
    invisible(tg)
    grid.segments(1.5, 0.1, 1.5, 0.9, gp = gpar(lwd = 4))
  })
}

dev.off()



besides_nt_stat <- list()

for(id in names(nt_counts_mats)){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat[id,])
  for(i in 12){
    #还是过滤1一下的结果
    if(sum(tmp_pct_mat[,i]) > 1){
      leftnt <- tmp_seq_mat[i - 1]
      rightnt <- tmp_seq_mat[i + 1]
      if(!leftnt %in% names(besides_nt_stat)){
        besides_nt_stat[[leftnt]] <- list()
      }
      tmp_besides_nt_stat <- besides_nt_stat[[leftnt]]
      if(!rightnt %in% names(tmp_besides_nt_stat)){
        tmp_nt_stat <- c(0,0,0,0)
        names(tmp_nt_stat) <- c("A", "T", "C", "G")
        tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      }
      tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
      
      for(nt in rownames(tmp_counts_mat)){
        tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_counts_mat[nt, i]
      }
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      besides_nt_stat[[leftnt]] <- tmp_besides_nt_stat
    }
  }
}


besides_2nt_stat <- list()

for(id in names(nt_counts_mats)){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat[id,])
  for(i in 12){
    #还是过滤1一下的结果
    if(sum(tmp_pct_mat[,i]) > 1){
      leftnt <- paste0(tmp_seq_mat[c(i - 3, i - 1)], collapse = "")
      rightnt <- paste0(tmp_seq_mat[c(i + 1, i + 3)], collapse = "")
      if(!leftnt %in% names(besides_2nt_stat)){
        besides_2nt_stat[[leftnt]] <- list()
      }
      tmp_besides_nt_stat <- besides_2nt_stat[[leftnt]]
      if(!rightnt %in% names(tmp_besides_nt_stat)){
        tmp_nt_stat <- c(0,0,0,0)
        names(tmp_nt_stat) <- c("A", "T", "C", "G")
        tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      }
      tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
      
      for(nt in rownames(tmp_counts_mat)){
        tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_counts_mat[nt, i]
      }
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      besides_2nt_stat[[leftnt]] <- tmp_besides_nt_stat
    }
  }
}

fake_plot_mat <- matrix(runif(4*4), nrow = 4, ncol = 4)
rownames(fake_plot_mat) <- c("A", "T", "C", "G")
colnames(fake_plot_mat) <- c("A", "T", "C", "G")
pdf("~/Nutstore Files/Tobin/Merged1NT/merged_left_right_nt_pie_chart_v4_groupPct.pdf", width = 6, height = 5)
col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(col) <-  c("A", 'T','C','G')
ht <- Heatmap(fake_plot_mat, 
              name = "Name",
              column_split = NULL,
              column_title=NULL,
              row_title = NULL,
              row_split = NULL,
              show_column_names = T, 
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_gp = gpar(fontsize = 24, col = col),
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 24, col = col),
              column_names_side = "top", top_annotation = NULL,
              left_annotation = NULL,
              right_annotation = NULL,
              show_row_names = T,
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                leftnt <- rownames(fake_plot_mat)[i]
                rightnt <- colnames(fake_plot_mat)[j]
                pct <- besides_nt_stat[[leftnt]][[rightnt]]
                labels <- names(pct)
                pct <- unlist(pct)
                col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                names(col) <- c("A", 'T','C','G')
                pct <- c(0, cumsum(pct)/sum(pct))
                dx <- diff(pct)
                nx <- length(dx)
                twopi <- 2 * pi
                init.angle <- 0
                t2xy <- function(t) {
                  t2p <- twopi * t + init.angle * pi/180
                  list(x = 0.8 * cos(t2p), y = 0.8 * sin(t2p))
                }
                maxx <- max(unlist(lapply(1 : nx, function(i){
                  n <- max(10, floor(200 * dx[i]))
                  P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                  max(abs(P$x))
                })))
                maxy <- max(unlist(lapply(1 : nx, function(i){
                  n <- max(10, floor(200 * dx[i]))
                  P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                  max(abs(P$y))
                })))
                facy <- as.numeric(h) / maxy * 0.3
                facx <- as.numeric(w) / maxx * 0.3
                for (i in 1L:nx) {
                  n <- max(10, floor(200 * dx[i]))
                  P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                  
                  grid.polygon(unit(as.numeric(x) + c(P$x * facx, 0), units = "npc"), 
                               unit(as.numeric(y) + c(P$y * facy, 0), units = "npc"), 
                               gp = gpar(fill = col[labels[i]]))
                  # print(P)
                  # polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = 45,
                  #         border = border[i], col = col[i], lty = lty[i])
                }
                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16, 
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)
draw(ht, annotation_legend_list = lgd_list)
dev.off()



besides_nt_sample <- lapply(names(nt_counts_mats), function(id){
  tmp_seq_mat <- unlist(seq_mat[id,])
  leftnt <- tmp_seq_mat[11]
  rightnt <- tmp_seq_mat[13]
  return(data.frame(id = id, left = leftnt, right = rightnt))
})
besides_nt_sample <- data.frame(do.call(rbind, besides_nt_sample))
besides_nt_sample$lr <- paste0(besides_nt_sample$left, besides_nt_sample$right)
besides_nt_sample$bulge <- T
besides_nt_sample$bulge[seq(2, nrow(besides_nt_sample), 2)] <- F
table(besides_nt_sample$lr, besides_nt_sample$bulge)


output_table <- lapply(names(besides_nt_stat), function(nt){
  x <- besides_nt_stat[[nt]]
  tmp <- lapply(x, function(y){
    y <- unlist(y)
    y / sum(y) * 100
  })
  tmp <- data.frame(do.call(rbind, tmp))
  tmp$rightNT <- rownames(tmp)
  tmp$leftNT <- nt
  tmp
})
output_table <- data.frame(do.call(rbind, output_table))
output_table <- output_table[,c(6, 5, 1,2,3,4)]
openxlsx::write.xlsx(output_table, 
                     file="~/Nutstore Files/Tobin/Merged1NT/merged_ins1_left_right_nt_pie_data.xlsx", 
                     rowNames=F, colNames=T)

save(besides_2nt_stat, file="~/Nutstore Files/Tobin/Merged1NT/ins1_beside_2nt.rda")


besides_nt_stat <- list()
besides_nt_sample_counts <- list()
for(id in names(nt_counts_mats)){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat[id,])
  for(i in 12){
    #还是过滤1一下的结果
    if(sum(tmp_pct_mat[,i]) > 1){
      leftnt <- paste0(tmp_seq_mat[i - c(3,1)], collapse = "")
      rightnt <- paste0(tmp_seq_mat[i + c(1, 3)], collapse = "")
      if(!leftnt %in% names(besides_nt_stat)){
        besides_nt_stat[[leftnt]] <- list()
        besides_nt_sample_counts[[leftnt]] <- list()
      }
      tmp_besides_nt_stat <- besides_nt_stat[[leftnt]]
      tmp_besides_nt_sample_counts <- besides_nt_sample_counts[[leftnt]]
      if(!rightnt %in% names(tmp_besides_nt_stat)){
        tmp_nt_stat <- c(0,0,0,0)
        names(tmp_nt_stat) <- c("A", "T", "C", "G")
        tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
        tmp_besides_nt_sample_counts[[rightnt]] <- 0
      }
      tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
      
      for(nt in rownames(tmp_counts_mat)){
        tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_counts_mat[nt, i]
      }
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      tmp_besides_nt_sample_counts[[rightnt]] <- tmp_besides_nt_sample_counts[[rightnt]] + 1
      besides_nt_stat[[leftnt]] <- tmp_besides_nt_stat
      besides_nt_sample_counts[[leftnt]] <- tmp_besides_nt_sample_counts
    }
  }
}


output_table <- lapply(names(besides_nt_stat), function(nt){
  x <- besides_nt_stat[[nt]]
  tmp <- lapply(x, function(y){
    y <- unlist(y)
    y / sum(y) * 100
  })
  tmp <- data.frame(do.call(rbind, tmp))
  tmp$rightNT <- rownames(tmp)
  tmp$leftNT <- nt
  tmp
})
output_table <- data.frame(do.call(rbind, output_table))
output_table <- output_table[,c(6, 5, 1,2,3,4)]
openxlsx::write.xlsx(output_table, 
                     file="~/Nutstore Files/Tobin/Merged1NT/merged_ins1_left_right_2nt_pie_data.xlsx", 
                     rowNames=F, colNames=T)
output_table <- output_table[grep("[ATCG]G",output_table$leftNT),]
output_table <- output_table[grep("[AT][ATCG]",output_table$rightNT),]
output_table$id <- paste(output_table$leftNT, output_table$rightNT)
output_table2 <- lapply(names(besides_nt_sample_counts), function(nt){
  x <- besides_nt_sample_counts[[nt]]
  tmp <- lapply(x, function(y){
    y
  })
  tmp <- data.frame(do.call(rbind, tmp))
  tmp$rightNT <- rownames(tmp)
  tmp$leftNT <- nt
  tmp
})
output_table2 <- data.frame(do.call(rbind, output_table2))
output_table2 <- output_table2[,c(3,2,1)]
output_table2 <- output_table2[grep("[ATCG]G",output_table2$leftNT),]
output_table2 <- output_table2[grep("[AT][ATCG]",output_table2$rightNT),]
output_table2$id <- paste(output_table2$leftNT, output_table2$rightNT)


output_table <- merge(output_table, output_table2, by="id")
output_table <- output_table[,-1]
output_table <- output_table[,c(-7,-8)]
colnames(output_table) <- c("leftNT", 'rightNT', "A", "T", 'C', 'G', "counts")
