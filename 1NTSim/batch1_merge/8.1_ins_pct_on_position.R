library(Biostrings)
library(stringr)
load("~/data/project/ear_project/gene_therapy_ll/Result/first_sgCmp.rda")


# bulge_pos <- lapply(split(sgCmp, sgCmp$id), function(tmp){
#   ins_one <- tmp$Cmp[tmp$isInsert == 'Insert']
#   del_one <- tmp$Cmp[tmp$isInsert != 'Insert']
#   ins_one_seq <- unlist(strsplit(ins_one, "*"))
#   del_one_seq <- unlist(strsplit(del_one, "*"))
#   ins_nt <- ins_one_seq[del_one_seq == "-"]
#   c(ins_one, del_one, ins_nt)
# 
# })
# bulge_pos <- data.frame(do.call(rbind, bulge_pos))
# bulge_pos$id <- rownames(bulge_pos)
# ToNX::write_tb(bulge_pos, "~/data/project/ear_project/gene_therapy_ll/First_sgRNA_cmp.txt")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_sgRNA_cmp_reDiff.txt", sep="\t")
load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first.rda")
###首先处理Ins1的结果  把所有Ins的编辑结果拿出来并且去掉存在Del的结果
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
    edit_table <- edit_table[edit_table$n_deleted == 0 & 
                               edit_table$n_inserted == 1& 
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
###到这里只剩下151个样本了

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

top_anno <- HeatmapAnnotation(position = anno_empty(border = F, width = unit(20, "mm")))
left_anno <- rowAnnotation(sg = anno_empty(border = F, width = unit(20, "mm")))
insert_color <- colorRamp2(c(0, 100), c("white", "#FF00FF"))
pdf("~/data/project/ear_project/gene_therapy_ll/test2.pdf", width = 6, height = 40)
Heatmap(plot_mat, 
        name = "Ins1 Pct",
        column_split = NULL,
        column_title=NULL,
        row_title = NULL,
        row_split = sg_info$pair,
        show_column_names = F, 
        cluster_rows = F, 
        cluster_columns = F, 
        column_names_rot = 0,
        column_names_centered = T,
        column_names_gp = gpar(fontsize = 24),
        column_names_side = "top", top_annotation = top_anno,
        left_annotation = NULL,
        right_annotation = NULL,
        show_row_names = T,
        row_names_side = "left",
        cell_fun = function(j, i, x, y, w, h, fill) {
          
          nt = pindex(seq_mat, i, j)
          col = c(2 : 5)
          names(col) <- c("A", 'T','C','G')
          sg <- rownames(seq_mat)[i]
          if(i %% 2== 1 & j == sg_info$bulge_in_mat[sg_info$id == sg]){
            grid.rect(x, y, w, h, gp = gpar(fill = "white", col = col[nt]))
          }
          if(nt != ""){
 
            grid.text(nt, x, y, 
                      gp = gpar(fontsize = 20, col = col[nt]))
          }


          if(nt == ""){
            nt_mat <- nt_mats[[sg]]
            if(sum(nt_mat[,j]) > 1){
              labels <- c("A", "T", "C", "G")
              pct <- unlist(nt_mat[,j])
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
            }
          }
          # if(nt == ""){
          #   col2 = "white"
          #   pct = pindex(plot_mat, i, j)
          #   col = insert_color(pct)
          #   grid.rect(x, y, w, h, gp = gpar(fill = col, col = col2))
          # }
          
        },layer_fun = function(j, i, x, y, w, h, fill) {
          grid.lines(y = unit(c(0.5, 0.5),units = "npc"))
        }, col=insert_color,show_heatmap_legend = T,
        
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
# sg_info_rev <- sg_info[nrow(sg_info) : 1,]
# for(i in 1 : length(unique(sg_info_rev$pair))){
#   decorate_annotation("sg", slice = i, {
#     tg <- richtext_grob(gt_render(sg_info_rev$id[c((i * 2) - 1,(i * 2))]), 
#                         rot = 0, 
#                         x = unit(0.5, "npc"),
#                         y=unit(c(0.25,0.75), "npc"), hjust = 0.5, gp = gpar(fontsize = 20))
#     grid.draw(tg)
#     invisible(tg)
#   })
# }


dev.off()
# ht <- draw(ht, annotation_legend_list = lgd_list)



###碱基Ins1与左右碱基关系图

nt_counts_mats <- lapply(remain_ins1_pct_table, function(tmp){
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
      plot_mat[tmp_left$NT[i], mat_pos] <- plot_mat[tmp_left$NT[i], mat_pos] + tmp_left$counts[i]
    }
  }
  
  if(nrow(tmp_right) > 0){
    for(i in 1 : nrow(tmp_right)){
      mat_pos <- (tmp_right$pos[i] - 14) * 2 - 2
      plot_mat[tmp_right$NT[i], mat_pos] <- plot_mat[tmp_right$NT[i], mat_pos] + tmp_right$counts[i]
    }
  }
  
  plot_mat
  
})
besides_nt_stat <- list()

for(id in names(nt_counts_mats)){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat[id,])
  for(i in seq(2, 16,2)){
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

fake_plot_mat <- matrix(runif(4*4), nrow = 4, ncol = 4)
rownames(fake_plot_mat) <- c("A", "T", "C", "G")
colnames(fake_plot_mat) <- c("A", "T", "C", "G")
pdf("~/data/project/ear_project/gene_therapy_ll/Result/left_right_nt_pie_chart.pdf", width = 6, height = 5)
ht <- Heatmap(fake_plot_mat, 
        name = "Name",
        column_split = NULL,
        column_title=NULL,
        row_title = NULL,
        row_split = NULL,
        show_column_names = T, 
        cluster_rows = F, 
        cluster_columns = F, 
        row_names_gp = gpar(fontsize = 24),
        column_names_rot = 0,
        column_names_centered = T,
        column_names_gp = gpar(fontsize = 24),
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
            col = c(2 : 5)
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
         legend_gp = gpar(fill = 2:5))
)
draw(ht, annotation_legend_list = lgd_list)
dev.off()







