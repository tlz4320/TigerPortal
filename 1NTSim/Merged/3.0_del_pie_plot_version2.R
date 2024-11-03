
library(Biostrings)
library(stringr)
#因为有相当一部分是forecast的数据 所以先读取这些数据
forecast_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/del1_indelphi_predict.xlsx")
sel_by_leftright <- unlist(lapply(split(forecast_info$sg, forecast_info$id), function(x){
  x <- unlist(x)
  if(str_sub(x[1], 17, 18) == str_sub(x[2], 17, 18)){
    return(x)
  }
  return("")
}))
forecast_info <- forecast_info[forecast_info$sg %in% sel_by_leftright,]
forecast_besides_nt_stat_withbulge <- list()
forecast_besides_nt_stat_withoutbulge <- list()
test <- list()
setwd("~/del1_forecast_res/")
for(id in forecast_info$id2){
  tmp <- read.table(paste0(id, "_predictedreads.txt"), header = F)
  tmp2 <- read.table(paste0(id, "_predictedindelsummary.txt"), header = F)
  tmp$Counts <- tmp2$V3
  test[[id]] <- tmp$V3[grep("D1_", tmp$V3)]
  
  tmp_del1 <- tmp[grep("D1_", tmp$V3),]
  tmp_del1$V3[grep("R0", tmp_del1$V3)] <- "Left"
  tmp_del1$V3[grep("L-1", tmp_del1$V3)] <- "Right"
  if(nrow(tmp_del1) > 2){
    print(paste0("not 2", id))
    next
  }
  tmp_del1$Pct <- tmp_del1$Counts / sum(tmp_del1$Counts) * 100

  leftnt <- str_sub(forecast_info$sg[forecast_info$id2 == id], 17,17)
  rightnt <- str_sub(forecast_info$sg[forecast_info$id2 == id], 18,18)
  if(leftnt == rightnt){
    tmp_del1 <- tmp_del1[1,]
    tmp_del1$V3 <- "Left"
  }
  if(any(!tmp_del1$V3 %in% c("Left", "Right"))){
    print(tmp_del1$V3[!tmp_del1$V3 %in% c("Left", "Right")])
  }
  if(forecast_info$bulge[forecast_info$id2 == id] == TRUE){
    if(!leftnt %in% names(forecast_besides_nt_stat_withbulge)){
      forecast_besides_nt_stat_withbulge[[leftnt]] <- list()
    }
    tmp_besides_nt_stat <- forecast_besides_nt_stat_withbulge[[leftnt]]
    if(!rightnt %in% names(tmp_besides_nt_stat)){
      tmp_nt_stat <- c(0,0,0,0)
      names(tmp_nt_stat) <- c("A", "T", "C", "G")
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
    }
    tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
    if("Left" %in% tmp_del1$V3){
      tmp_nt_stat[leftnt] <- tmp_nt_stat[leftnt] + tmp_del1$Pct[tmp_del1$V3 == "Left"]
    }
    if("Right" %in% tmp_del1$V3){
      tmp_nt_stat[rightnt] <- tmp_nt_stat[rightnt] + tmp_del1$Pct[tmp_del1$V3 == "Right"]
    }

    tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
    forecast_besides_nt_stat_withbulge[[leftnt]] <- tmp_besides_nt_stat
  } else {
    if(!leftnt %in% names(forecast_besides_nt_stat_withoutbulge)){
      forecast_besides_nt_stat_withoutbulge[[leftnt]] <- list()
    }
    tmp_besides_nt_stat <- forecast_besides_nt_stat_withoutbulge[[leftnt]]
    if(!rightnt %in% names(tmp_besides_nt_stat)){
      tmp_nt_stat <- c(0,0,0,0)
      names(tmp_nt_stat) <- c("A", "T", "C", "G")
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
    }
    tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
    if("Left" %in% tmp_del1$V3){
      tmp_nt_stat[leftnt] <- tmp_nt_stat[leftnt] + tmp_del1$Pct[tmp_del1$V3 == "Left"]
    }
    if("Right" %in% tmp_del1$V3){
      tmp_nt_stat[rightnt] <- tmp_nt_stat[rightnt] + tmp_del1$Pct[tmp_del1$V3 == "Right"]
    }
    
    tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
    forecast_besides_nt_stat_withoutbulge[[leftnt]] <- tmp_besides_nt_stat
  }
}
rm(test)

forecast_besides_nt_stat <- list()
for(i in c("A", "T", "C", "G")){
  forecast_besides_nt_stat[[i]] <- list()
  for(j in c("A", "T", "C", "G")){
    if(i %in% names(forecast_besides_nt_stat_withbulge)){
      if(j %in% names(forecast_besides_nt_stat_withbulge[[i]]))
      forecast_besides_nt_stat[[i]][[j]] <- forecast_besides_nt_stat_withbulge[[i]][[j]]
    }
    if(i %in% names(forecast_besides_nt_stat_withoutbulge)){
      if(j %in% names(forecast_besides_nt_stat_withoutbulge[[i]])){
        if(j %in% names(forecast_besides_nt_stat[[i]])){
          forecast_besides_nt_stat[[i]][[j]] <- forecast_besides_nt_stat[[i]][[j]]+ forecast_besides_nt_stat_withoutbulge[[i]][[j]]
        } else{
          forecast_besides_nt_stat[[i]][[j]] <- forecast_besides_nt_stat_withoutbulge[[i]][[j]]
          
        }
      }
    }
  }

}

source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first_second.Rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt", sep="\t")
###首先处理Del1的结果  把所有Del的编辑结果拿出来并且去掉存在Ins的结果#上面部分每次跑都很慢，不如直接改成保存下来读取
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_del_table_data.rda")
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
###首先判断Delete比例是不是都小于1% indel reads

sg_info$lowCount <- unlist(lapply(sg_info$id, function(id){
  tmp <- remain_del1_pct_table_pos[[id]]
  tmp <- tmp[tmp$pos %in% c(20:21),]
  all(tmp$Pct < 1)
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
sel_by_leftright <- unlist(lapply(as.character(sg_info_sel$pair), function(x){
  ids <- sg_info_sel$id[sg_info_sel$pair == x]
  if(sum(wt_seqs[[ids[1]]][c(20,21)] == wt_seqs[[ids[2]]][c(20,21)]) == 2){
    return(x)
  }
  return("")
  
}))
sg_info_sel <- sg_info_sel[sg_info_sel$pair %in% sel_by_leftright,]
seq_mat_sel <- seq_mat[sg_info_sel$id,]
plot_mat_sel <- plot_mat[sg_info_sel$id,]
left_anno <- rowAnnotation(line = anno_empty(border = F, width = unit(2, "mm")), 
                           rowname = anno_text(x = rownames(plot_mat_sel)))
right_anno <- rowAnnotation(line2 = anno_empty(border = F, width = unit(2, "mm")))




remain_del1_pct_table_sel <- remain_del1_pct_table[sg_info_sel$id]

#重新统计一个切割位点左右1NT的情况



nt_counts_mats <- lapply(remain_del1_pct_table_sel, function(tmp){
  plot_mat <- matrix(0, ncol = 9, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(20 : 21),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  for(i in 1 : nrow(tmp)){
    mat_pos <- tmp$pos[i] - 14
    plot_mat[tmp$NT[i], mat_pos] <- plot_mat[tmp$NT[i], mat_pos] + tmp$Pct[i]
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
nt_pct_mats <- lapply(remain_del1_pct_table_sel, function(tmp){
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
sg_info_sel <- sg_info[!sg_info$lowCount,]
#去掉低质量的序列
for(id in sg_info_sel$id){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat2[id,])
  tmp_pct_mat <- nt_pct_mats[[id]]
  leftnt <- tmp_seq_mat[8]
  rightnt <- tmp_seq_mat[9]
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
                     file="~/Nutstore Files/Tobin/Merged1NT/merged_del1_left_right_nt_pie_data_add_forecast.xlsx", 
                     rowNames=F, colNames=T)


fake_plot_mat <- matrix(runif(4*4), nrow = 4, ncol = 4)

rownames(fake_plot_mat) <- c("A", "T", "C", "G")
colnames(fake_plot_mat) <- rownames(fake_plot_mat) 
top_anno <- HeatmapAnnotation("rightNT1" = anno_empty(border = F, width = unit(15, "mm")))
left_anno <- rowAnnotation("leftNT1" = anno_empty(border = F, width = unit(15, "mm")))

library(gridtext)
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
pdf("~/Nutstore Files/Tobin/Merged1NT/merged_del1_left_right_nt_pie_add_forecast.pdf", width = 6, height = 5)
draw(ht, annotation_legend_list = lgd_list)
dev.off()

nt_counts_mats <- lapply(remain_del1_pct_table_sel, function(tmp){
  plot_mat <- matrix(0, ncol = 9, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(20 : 21),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  for(i in 1 : nrow(tmp)){
    mat_pos <- tmp$pos[i] - 14
    plot_mat[tmp$NT[i], mat_pos] <- plot_mat[tmp$NT[i], mat_pos] + tmp$Pct[i]
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

###用于保存已经包含的左右两个碱基的sgRNA情况  用于补充forecast的数据
# left_right_nt_table <- lapply(sg_info_sel$id, function(x){
#   paste(unlist(seq_mat2[x,c(7,8,9,10)]), collapse = "")
# })
# left_right_nt_table <- data.frame(id = sg_info_sel$id, nt = unlist(left_right_nt_table))
# ToNX::write_tb(left_right_nt_table, "~/data/project/ear_project/gene_therapy_ll/left_right_nt_table_after_filtered.txt")
#把有bulge和没有bulge的分开计算
besides_nt_stat_with_bulge <- list()


for(id in names(nt_counts_mats)[seq(1, length(nt_counts_mats), 2)]){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat2[id,])
  tmp_pct_mat <- nt_counts_mats[[id]]
  leftnt <- tmp_seq_mat[8]
  rightnt <- tmp_seq_mat[9]
  if(leftnt == "A" & rightnt == "C"){
    print(id)
  }
  for(i in c(6,7)){
    #还是过滤1一下的结果
    if(sum(tmp_pct_mat[,i]) > 1){
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
  tmp_seq_mat <- unlist(seq_mat2[id,])
  tmp_pct_mat <- nt_counts_mats[[id]]
  leftnt <- tmp_seq_mat[8]
  rightnt <- tmp_seq_mat[9]
  if(leftnt == "A" & rightnt == "C"){
    print(id)
  }
  for(i in c(6,7)){
    #还是过滤1一下的结果
    if(sum(tmp_pct_mat[,i]) > 1){
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


for(name in c("A",'T', "C", "G")){
  tmp <- forecast_besides_nt_stat_withbulge[[name]]
  for(name2 in c("A",'T', "C", "G")){
    if(name == name2){
      next
    }
    if(name2 %in% names(besides_nt_stat_with_bulge[[name]])){
      besides_nt_stat_with_bulge[[name]][[name2]] <- besides_nt_stat_with_bulge[[name]][[name2]] + tmp[[name2]]
    } else {
      besides_nt_stat_with_bulge[[name]][[name2]] <- tmp[[name2]]
    }
    
  }
  tmp <- forecast_besides_nt_stat_withoutbulge[[name]]
  for(name2 in c("A",'T', "C", "G")){
    if(name == name2){
      next
    }
    if(name2 %in% names(besides_nt_stat_without_bulge[[name]])){
      besides_nt_stat_without_bulge[[name]][[name2]] <- besides_nt_stat_without_bulge[[name]][[name2]] + tmp[[name2]]
    } else {
      besides_nt_stat_without_bulge[[name]][[name2]] <- tmp[[name2]]
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
                if(!is.null(pct)){
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
                
                }
                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16, 
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)

pdf("~/Nutstore Files/Tobin/Merged1NT/merged_left_right_nt_pie_chart_v4_add_forecast.pdf", width = 6, height = 10.3)

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
