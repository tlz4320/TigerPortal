library(Biostrings)
library(stringr)
load("~/data/project/ear_project/gene_therapy_ll/Result/first_sgCmp.rda")
source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_sgRNA_cmp_reDiff.txt", sep="\t")
load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first.rda")
###首先处理Del1的结果  把所有Del的编辑结果拿出来并且去掉存在Ins的结果
select_edit_pattern <- list()

for(id in names(total_edit_table_rev)){
  edit_tables <- total_edit_table_rev[[id]]
  edit_table <- lapply(edit_tables, function(edit_table){
    isindel <- edit_table$n_inserted + edit_table$n_deleted != 0
    edit_table <- edit_table[isindel,]
    edit_table$Pct <- edit_table[,7] / sum(edit_table[,7]) * 100
    edit_table <- lapply(split(edit_table, edit_table$Aligned_Sequence), 
                         function(x){
                           if(nrow(x) == 1){
                             return(x)
                           }
                           x[,9] <- sum(x[,9])
                           return(x[1,])
                         })
    edit_table <- data.frame(do.call(rbind, edit_table))
    edit_table <- edit_table[order(edit_table$Pct, decreasing = T),]
    ##选出只有Ins1的结果
    edit_table <- edit_table[edit_table$n_deleted == 1 & 
                               edit_table$n_inserted == 0& 
                               edit_table$n_mutated <= 1,]
    edit_table
  })
  select_edit_pattern[[id]] <- edit_table
}
seqs_table <- list()
##对取出来的结果算样本之间的均值并且删掉一些Ins发生的频率低于1%的
for(id in names(select_edit_pattern)){
  edit_table <- select_edit_pattern[[id]]
  
  edit_table <- data.frame(do.call(rbind, edit_table))
  if(nrow(edit_table) == 0)
    next
  if(sum(edit_table[,8]) < 0.1)
    next
  edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, edit_table$Aligned_Sequence), function(x){
    x[,7] <- mean(x[,7])
    x[,8] <- mean(x[,8])
    x[,9] <- mean(x[,9])
    x[1,]
  })))
  edit_table$id <- id
  edit_table <- edit_table[order(edit_table$Pct, decreasing = T),]
  seqs_table[[id]] <- edit_table[,c(1,2,7, 8, 9, 10)]
}

seqs_table <- data.frame(do.call(rbind, seqs_table))
###到这里只剩下148个样本了

###取出wt序列，用于后面的统计
wt_seqs <- lapply(unique(seqs_table$id), function(id){
  edit_table <- total_edit_table_rev[[id]][[1]]
  edit_table <- edit_table[edit_table$Unedited == "True",]
  getWtSeq(edit_table)
})
names(wt_seqs) <- unique(seqs_table$id)


del1_pct_table <- lapply(unique(seqs_table$id), function(id){
  del1_table <- seqs_table[seqs_table$id == id,]
  wt_seq <- wt_seqs[[id]]
  tmp <- getStat2(del1_table, F)
  tmp$id <- id
  tmp
})
names(del1_pct_table) <- unique(seqs_table$id)
del1_pct_table_pos <- lapply(del1_pct_table, function(tmp){
  tmp <- lapply(split(tmp, tmp$pos), function(x){
    x$Pct <- sum(x$Pct)
    x[1,]
  })
  tmp <- do.call(rbind, tmp)
}) 

###取出还成对的样本
remain_samples <- lapply(split(sgCmp$id2, sgCmp$id), function(x){
  if(sum(unlist(x) %in% names(del1_pct_table_pos)) == 2){
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
  #把del1放在前面
  tmp <- tmp[order(tmp$isInsert, decreasing = T),]
  tmp$id2
}))

remain_del1_pct_table <- del1_pct_table[bulge_pos_order]
remain_del1_pct_table_pos <- del1_pct_table_pos[bulge_pos_order]

plot_mat <- matrix(0, ncol = 9, nrow = length(remain_del1_pct_table))
rownames(plot_mat) <- bulge_pos_order

for(id in names(remain_del1_pct_table_pos)){
  tmp <- remain_del1_pct_table_pos[[id]]
  tmp <- tmp[tmp$pos %in% c(20:21),]
  tmp$Pct2 <- tmp$Pct / sum(tmp$Pct) * 100
  for(i in 1 : nrow(tmp)){
    pos <- tmp$pos[i] - 14
    plot_mat[id, pos] <- plot_mat[id, pos]  + tmp$Pct2[i]
    
  }

}



seq_mat <- matrix("", ncol = 17, nrow = length(remain_del1_pct_table))
rownames(seq_mat) <- bulge_pos_order
for(id in bulge_pos_order){
  wt_seq <- wt_seqs[[id]][15 : 23]
  for(i in 1 : length(wt_seq)){
    seq_mat[id, (i - 1) * 2 + 1] <- wt_seq[i]
  }
}
seq_mat <- seq_mat[,seq(1,ncol(seq_mat), 2)]



sg_info <- data.frame(id = bulge_pos_order, pair = unlist(lapply(bulge_pos_order, function(x){
  sgCmp$id[sgCmp$id2== x]
})))
sg_info$bulge <- unlist(lapply(sg_info$pair, function(x){
  bulge_pos_sel$V5[bulge_pos_sel$V4== x]
}))
sg_info$pair <- factor(sg_info$pair, levels = unique(sg_info$pair))
sg_info$bulge_in_mat <- unlist(lapply(sg_info$bulge, function(pos){
  10 - pos
}))

top_anno <- HeatmapAnnotation(position = anno_empty(border = F, width = unit(20, "mm")))
left_anno <- rowAnnotation(sg = anno_empty(border = F, width = unit(20, "mm")))
maxratio <- max(unlist(plot_mat))
library(circlize)
delet_color <- colorRamp2(c(0, maxratio), c("white", "#FF00FF"))


####最终要去掉一些低质量的结果
###首先判断Delete Reads是不是小于200条

sg_info$lowCount <- unlist(lapply(sg_info$id, function(id){
  tmp <- remain_del1_pct_table_pos[[id]]
  tmp <- tmp[tmp$pos %in% c(20:21),]
  all(tmp$counts < 200)
}))
sg_info_sel <- sg_info[!sg_info$lowCount,]
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

pdf("~/data/project/ear_project/gene_therapy_ll/Result/Del1_seq_logo_final_v3.pdf", width = 6, height = 40)
Heatmap(plot_mat_sel, 
        name = "del1 Pct",
        column_split = NULL,
        column_title=NULL,
        row_title = NULL,
        row_split = sg_info_sel$pair,
        show_column_names = F, 
        cluster_rows = F, 
        cluster_columns = F, 
        column_names_rot = 0,
        column_names_centered = T,
        column_names_gp = gpar(fontsize = 24),
        column_names_side = "top", top_annotation = top_anno,
        left_annotation = left_anno,
        right_annotation = right_anno,
        show_row_names = F,
        row_names_side = "left",
        cell_fun = function(j, i, x, y, w, h, fill) {
          
          nt = pindex(seq_mat_sel, i, j)
          col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
          names(col) <- c("A", 'T','C','G')
          sg <- rownames(seq_mat_sel)[i]
          pct <- pindex(plot_mat_sel,i, j)
          if(i %% 2== 1 & j == sg_info_sel$bulge_in_mat[sg_info_sel$id == sg]){
            grid.circle(x, y, r = 0.25, gp = gpar(fill = delet_color(pct), col = "#82328C", lwd = 3))
          } else if(pct != 0){
            grid.circle(x, y, r = 0.25, gp = gpar(fill = delet_color(pct), col = "#FF6633", lwd = 3))
          }
          if(nt != ""){
            
            grid.text(nt, x, y, 
                      gp = gpar(fontsize = 20, col = col[nt]))
          }
          
          
  
          
        },layer_fun = function(j, i, x, y, w, h, fill) {
          grid.lines(y = unit(c(0.5, 0.5), units = "npc"), x = unit(c(-0.2, 1), units = "npc"))
        }, col=c("white", "white"),show_heatmap_legend = F,
        
        border = F)


poss <- c(-6 : 3)
poss <- poss[poss != 0]

decorate_annotation("position", slice = 1, {
  tg <- richtext_grob(gt_render(poss), 
                      rot = 0, 
                      x = unit(seq(1, 17, 2) / 17 - 0.5 / 17, "npc"),
                      y=unit(0.5, "npc"), hjust = 0.5, gp = gpar(fontsize = 20))
  grid.draw(tg)
  invisible(tg)
})
for(i in 1 : length(unique(sg_info_sel$pair))){
  decorate_annotation("line", slice = i, {
    grid.lines(y = unit(c(1 / 4, 1 / 4), "npc"), 
               x=unit(c(0, 1), "npc"))
    grid.lines(y = unit(c(1 / 4, 1 - 1 / 4), "npc"), 
               x=unit(c(0, 0), "npc"))
    grid.lines(y = unit(c(1 - 1 / 4, 1 - 1 / 4), "npc"), 
               x=unit(c(0, 1), "npc"))
  })
  decorate_annotation("line2", slice = i, {
    grid.lines(y = unit(c(1 / 4, 1 / 4), "npc"), 
               x=unit(c(0, 1), "npc"))
    grid.lines(y = unit(c(1 / 4, 1 - 1 / 4), "npc"), 
               x=unit(c(1, 1), "npc"))
    grid.lines(y = unit(c(1 - 1 / 4, 1 - 1 / 4), "npc"), 
               x=unit(c(0, 1), "npc"))
  })
}



dev.off()
# ht <- draw(ht, annotation_legend_list = lgd_list)



###碱基del1与左右碱基关系图


nt_counts_mats <- lapply(remain_del1_pct_table, function(tmp){
  plot_mat <- matrix(0, ncol = 9, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(15 : 23),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]

    for(i in 1 : nrow(tmp)){
      mat_pos <- tmp$pos[i] - 14
      plot_mat[tmp$NT[i], mat_pos] <- plot_mat[tmp$NT[i], mat_pos] + tmp$counts[i]
    }
  

  
  plot_mat
  
})


seq_mat2 <- matrix("", ncol = 14, nrow = length(remain_del1_pct_table))
rownames(seq_mat2) <- bulge_pos_order
for(id in bulge_pos_order){
  wt_seq <- wt_seqs[[id]][13 : 26]
  for(i in 1 : length(wt_seq)){
    seq_mat2[id, i] <- wt_seq[i]
  }
}
nt_pct_mats <- lapply(remain_del1_pct_table, function(tmp){
  plot_mat <- matrix(0, ncol = 9, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(15 : 23),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
  
  for(i in 1 : nrow(tmp)){
    mat_pos <- tmp$pos[i] - 14
    plot_mat[tmp$NT[i], mat_pos] <- plot_mat[tmp$NT[i], mat_pos] + tmp$Pct[i]
  }
  
  
  
  plot_mat
  
})


besides_nt_stat <- list()

for(id in sg_info_sel$id){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat2[id,])
  tmp_pct_mat <- nt_pct_mats[[id]]
  leftnt <- paste(tmp_seq_mat[c(7,8)], collapse = "")
  rightnt <- paste(tmp_seq_mat[c(9, 10)], collapse = "")
  for(i in c(6,7)){
    #还是过滤1一下的结果
    if(sum(tmp_pct_mat[,i]) > 1){
      tmp_pct_mat[,i] <- tmp_pct_mat[,i] / sum(unlist(tmp_pct_mat[,c(6,7)]))
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

fake_plot_mat <- matrix(runif(16*16), nrow = 16, ncol = 16)

rownames(fake_plot_mat) <- unlist(lapply(1 : 4, function(i){
  lapply(1 : 4, function(j){
    paste(c("A", "T", "C", "G")[c(i,j)], collapse = "")
  })
  
}))
colnames(fake_plot_mat) <- rownames(fake_plot_mat) 
pdf("~/data/project/ear_project/gene_therapy_ll/Result/left_right_nt_pie_chart_del_v4.pdf", width = 20, height = 15)
top_anno <- HeatmapAnnotation("rightNT2" = anno_empty(border = F, width = unit(8, "mm")), 
                              "line1" = anno_empty(border = F, width = unit(2, "mm")),
                              "rightNT1" = anno_empty(border = F, width = unit(15, "mm")))
left_anno <- rowAnnotation("leftNT2" = anno_empty(border = F, width = unit(8, "mm")), 
                           "line2" = anno_empty(border = F, width = unit(15, "mm")),
                           "leftNT1" = anno_empty(border = F, width = unit(15, "mm")))

library(gridtext)
ht <- Heatmap(fake_plot_mat, 
              name = "Name",
              column_split = rep(c(" ", "  ", "   ", "    "), c(4,4,4,4)),
              column_title=NULL,
              row_title = NULL,
              row_split = rep(c(" ", "  ", "   ", "    "), c(4,4,4,4)),
              show_column_names = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_gp = gpar(fontsize = 24),
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 24),
              column_names_side = "top", 
              top_annotation = top_anno,
              left_annotation = left_anno,
              right_annotation = NULL,
              show_row_names = F,
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                leftnt <- rownames(fake_plot_mat)[i]
                leftnt <- paste(rev(unlist(strsplit(leftnt, "*"))), collapse = "")
                rightnt <- colnames(fake_plot_mat)[j]
                # rightnt <- paste(rev(unlist(strsplit(rightnt, "*"))), collapse = "")
                pct <- besides_nt_stat[[leftnt]][[rightnt]]
                #& ceiling(i / 4) != ceiling(j / 4)
                if(!is.null(pct) ){
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
                }
                }
                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F, rect_gp = gpar(border = T, col = "#E0E0E0"))
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "delertNT", type = "grid", pch = 16, 
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)
draw(ht, annotation_legend_list = lgd_list)
for(i in 1 : 4){
  nt <- c("A", "T", "C", "G")[i]

  decorate_annotation("rightNT1", slice = i, {
    tg <- richtext_grob(gt_render(nt), 
                        rot = 0, 
                        x = unit(0.5, "npc"),
                        y=unit(0.5, "npc"), vjust = 0.5, gp = gpar(fontsize = 40))
    grid.draw(tg)
    invisible(tg)
  })
  decorate_annotation("leftNT1", slice = i, {
    tg <- richtext_grob(gt_render(nt), 
                        rot = 0, 
                        x = unit(0.5, "npc"),
                        y=unit(0.5, "npc"), hjust = 0.9, gp = gpar(fontsize = 40))
    grid.draw(tg)
    invisible(tg)
  })
  
  decorate_annotation("rightNT2", slice = i, {
    tg <- richtext_grob(gt_render(c("A", "T", "C", "G")), 
                        rot = 0, 
                        x = unit(c(0.125, 0.375, 0.625, 0.875), "npc"),
                        y=unit(0.5, "npc"), vjust = 0.9, gp = gpar(fontsize = 20))
    grid.draw(tg)
    invisible(tg)
  })
  decorate_annotation("leftNT2", slice = i, {
    tg <- richtext_grob(gt_render(rev(c("A", "T", "C", "G"))), 
                        rot = 0, 
                        x = unit(0.5, "npc"),
                        y=unit(c(0.125, 0.375, 0.625, 0.875), "npc"), hjust = 0.1, gp = gpar(fontsize = 20))
    grid.draw(tg)
    invisible(tg)
  })
  decorate_annotation("line1", slice = i, {
    grid.lines(x = unit(c(0.1, 0.9), "npc"), 
              y=unit(c(0.5, 0.5), "npc"))
  })
  decorate_annotation("line2", slice = i, {
    grid.lines(x = unit(c(0.5, 0.5), "npc"), 
              y=unit(c(0.1, 0.9), "npc"))
  })
}

dev.off()









