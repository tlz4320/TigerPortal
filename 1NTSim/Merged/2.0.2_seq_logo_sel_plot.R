library(Biostrings)
library(stringr)
source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first_second.Rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt", sep="\t")
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
top_anno <- HeatmapAnnotation(position = anno_empty(border = F, width = unit(1, "inches")))
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

###输出ins1的比例
ins1_as_left_pct <- lapply(rownames(plot_mat_sel), function(id){
  mat <- nt_mats[[id]]
  left_nt <- seq_mat_sel[id, 11]
  pct <- mat[left_nt, 12]
  group_pct <- pct / sum(mat[,12]) * 100
  data.frame(sg = id, nt = left_nt, pct_in_indel = pct, pct_in_cut_site_ins1 = group_pct)
})
ins1_as_left_pct <- data.frame(do.call(rbind, ins1_as_left_pct))
openxlsx::write.xlsx(ins1_as_left_pct, file="~/Nutstore Files/Tobin/Merged1NT/ins1_is_leftNT_Percentage.xlsx", rowNames=F, colNames=T)


for(i in 1 : 15){
  tmp <- get(paste0("sel", i))
  print(all(tmp %in% rownames(plot_mat_sel)))
  library(ToNX)
  sg_info_sel2 <- sg_info_sel[match(tmp, sg_info_sel$id),]
  #为了删除空隙，所以去掉中间的空白
  plot_mat_sel2 <- plot_mat_sel[tmp, c(1,3,5,7,9,11,12,13,15,17)]
  seq_mat_sel2 <- seq_mat_sel[tmp, ]
  sg_info_sel2$pair <- as.character(sg_info_sel2$pair)
  sg_info_sel2$pair <- factor(sg_info_sel2$pair, levels = unique(sg_info_sel2$pair))
  left_anno2 <- rowAnnotation(line = anno_empty(border = F, width = unit(2, "mm")), 
                              rowname = anno_text(x = rownames(plot_mat_sel2)))
  right_anno2 <- rowAnnotation(line2 = anno_empty(border = F, width = unit(2, "mm")))
  
  
  ht <- Heatmap(plot_mat_sel2, 
                name = "Ins1 Pct",
                column_split = NULL,
                column_title=NULL,
                row_title = NULL,
                row_split = sg_info_sel2$pair,
                show_column_names = F, 
                cluster_rows = F, 
                cluster_columns = F, 
                column_names_rot = 0,
                column_names_centered = T,
                column_names_gp = gpar(fontsize = 24),
                column_names_side = "top", top_annotation = top_anno,
                left_annotation = left_anno2,
                right_annotation = right_anno2,
                show_row_names = F,
                row_names_side = "left",
                cell_fun = function(j, i, x, y, w, h, fill) {
                  #因为去掉了空白，和之前的统计坐标可能对不上 需要重新对一下坐标
                  if(j <= 6){
                    j <- j * 2 - 1
                  } else if(j == 7){
                    j <- 12
                  } else if(j > 7){
                    j <- j * 2 - 3
                  }
                  nt = pindex(seq_mat_sel2, i, j)
                  col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                  names(col) <- c("A", 'T','C','G')
                  sg <- rownames(seq_mat_sel2)[i]
                  if(i %% 2 == 1 & j == sg_info_sel2$bulge_in_mat[sg_info_sel2$id == sg]){
                    grid.rect(x, y,width = as.numeric(w) * 0.8, height = as.numeric(h) * 0.8, 
                              gp = gpar(fill = "white", col = "#82328C", lwd = 1))
                  }
                  if(nt != ""){
                    
                    grid.text(nt, x, y, 
                              gp = gpar(fontsize = 26, col = col[nt]))
                  }
                  
                  
                  if(nt == ""){
                    left_nt <- pindex(seq_mat_sel2, i, j - 1)
                    nt_mat <- nt_mats[[sg]]
                    if(sum(nt_mat[,j]) > 1 & j  == 12){
                      pct <- unlist(nt_mat[,j])
                      which_nt <- which(left_nt == names(pct))
                      pct <- c(sort(pct[-which_nt], decreasing = F), pct[which_nt])
                      ig <- 270
                      if(i %% 2 == 0){
                        pct <- rev(pct)
                        ig <- 90
                      }
                      grid.pie(pct,x = x,y = y, w= w, h=h, 
                               colors = col,size.scales = 0.45,
                               init.angle = ig)
                    }
                  }
                },layer_fun = function(j, i, x, y, w, h, fill) {
                  grid.lines(y = unit(c(0.5, 0.5), units = "npc"), x = unit(c(-0.15, 1), units = "npc"))
                  c}, col=c("white", "white"),show_heatmap_legend = F,
                
                border = F)
  col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
  
  pdf(paste0("~/Nutstore Files/Tobin/Merged1NT/Final_seq_logo_ins1_sFig4b_", i, "_fix.pdf"), 
      width =5, height = 6.14 * length(tmp) / 16 + 1)
  ht <- draw(ht, annotation_legend_list = list(  
    Legend(labels = c("A", 'T','C','G'), title = "NT", type = "grid", pch = 16, 
           legend_gp = gpar(fill = col), nrow = 1)), 
    annotation_legend_side = "top")
  poss <- c(-6 : 3)
  poss <- poss[poss != 0]
  library(gridtext)
  # decorate_annotation("position", slice = 1, {
  #   tg <- richtext_grob(gt_render(poss), 
  #                       rot = 0, 
  #                       x = unit(c(1 : 10 / 10 - 0.5 / 10)[-7], "npc"),
  #                       y=unit(0.5, "npc"), hjust = 0.5, gp = gpar(fontsize = 20))
  #   grid.draw(tg)
  #   invisible(tg)
  # })
  for(i in 1 : length(unique(sg_info_sel2$pair))){
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
}
