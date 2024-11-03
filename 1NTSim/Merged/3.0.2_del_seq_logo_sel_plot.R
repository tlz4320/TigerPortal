library(Biostrings)
library(stringr)
source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first_second.Rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt", sep="\t")
###首先处理Del1的结果  把所有Del的编辑结果拿出来并且去掉存在Ins的结果
#上面部分每次跑都很慢，不如直接改成保存下来读取
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

top_anno <- HeatmapAnnotation(position = anno_empty(border = F, width = unit(1, "inches")))
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

seq_mat_sel <- seq_mat[sg_info_sel$id,]
plot_mat_sel <- plot_mat[sg_info_sel$id,]
left_anno <- rowAnnotation(line = anno_empty(border = F, width = unit(2, "mm")), 
                           rowname = anno_text(x = rownames(plot_mat_sel)))
right_anno <- rowAnnotation(line2 = anno_empty(border = F, width = unit(2, "mm")))

###输出del1的比例
del1_as_left_pct <- lapply(rownames(plot_mat_sel), function(id){
  left_nt <- seq_mat_sel[id, 6]
  right_nt <- seq_mat_sel[id, 7]
  tmp <- remain_del1_pct_table_pos[[id]]
  tmp <- tmp[tmp$pos %in% c(20:21),]
  pct_left <- sum(tmp$Pct[tmp$pos == 20])
  pct_right <- sum(tmp$Pct[tmp$pos == 21])
  if(left_nt == right_nt){
    pct_left <- pct_right <- (pct_left + pct_right) / 2
  }
  data.frame(sg = id, nt = left_nt, left_pct_in_indel = pct_left, 
             right_pct_in_indel = pct_right, 
             left_pct_in_cut_site_del1 = pct_left / (pct_left + pct_right) * 100)
})
del1_as_left_pct <- data.frame(do.call(rbind, del1_as_left_pct))
del1_as_left_pct$right_pct_in_cut_site_del1 <- 100 - del1_as_left_pct$left_pct_in_cut_site_del1
openxlsx::write.xlsx(del1_as_left_pct, file="~/Nutstore Files/Tobin/Merged1NT/del1_is_leftNT_Percentage.xlsx", rowNames=F, colNames=T)



for(i in 1 : 11){
  sel <- get(paste0("sel", i))
  print(all(sel %in% rownames(plot_mat_sel)))
  
  sg_info_sel2 <- sg_info_sel[match(sel, sg_info_sel$id),]
  tmp <- plot_mat_sel[sel, ]
  plot_mat_sel2 <- matrix(0, nrow = nrow(tmp), ncol = 10)
  plot_mat_sel2[,1:6] <- tmp[,1:6]
  plot_mat_sel2[,8 : 10] <- tmp[,7:9]
  rownames(plot_mat_sel2) <- sel
  seq_mat_sel2 <- seq_mat_sel[sel, ]
  
  sg_info_sel2$pair <- as.character(sg_info_sel2$pair)
  sg_info_sel2$pair <- factor(sg_info_sel2$pair, levels = unique(sg_info_sel2$pair))
  left_anno2 <- rowAnnotation(line = anno_empty(border = F, width = unit(2, "mm")), 
                              rowname = anno_text(x = rownames(plot_mat_sel2)))
  right_anno2 <- rowAnnotation(line2 = anno_empty(border = F, width = unit(2, "mm")))
  
  filename <- rename$name[rename$sel == paste0("sel", i)]
  pdf(paste0("~/Nutstore Files/Tobin/Merged1NT/Final_seq_logo_del1_",filename,".pdf"), 
      width =5, height = 5.87 * length(sel) / length(sel1) + 1)
  ht <- Heatmap(plot_mat_sel2, 
                name = "del1 Pct",
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
                  #因为插入了一个空列
                  if(j == 7){
                    left_nt = pindex(seq_mat_sel2, i, 6)
                    right_nt = pindex(seq_mat_sel2, i, 7)
                    left_pct <- pindex(plot_mat_sel2,i, 6)
                    right_pct <- pindex(plot_mat_sel2,i, 8)
                    point_pct <- left_pct / (left_pct + right_pct)
                    if(left_nt == right_nt){
                      point_pct <- 0.5
                    }
                    top_pos <- as.numeric(y) + as.numeric(h) * 0.4
                    bot_pos <- as.numeric(y) - as.numeric(h) * 0.4
                    left_pos <- as.numeric(x) - as.numeric(w) * 0.5
                    right_pos <- as.numeric(x) + as.numeric(w) * 0.5
                  
                    center_pos <- left_pos + point_pct * as.numeric(w)
      
                    col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                    names(col) <- c("A", 'T','C','G')
                    grid.polygon(x = c(left_pos, left_pos, center_pos), 
                                 y = c(top_pos, bot_pos, y), gp = gpar(fill = col[left_nt]))
                    grid.polygon(x = c(right_pos, right_pos, center_pos), 
                                 y = c(top_pos, bot_pos, y), gp = gpar(fill = col[right_nt]))
                    
                  }else{
                    if(j > 6){
                      j <- j - 1
                    }
                    nt = pindex(seq_mat_sel2, i, j)
                    col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                    names(col) <- c("A", 'T','C','G')
                    sg <- rownames(seq_mat_sel2)[i]
                    if(i %% 2== 1 & j == sg_info_sel2$bulge_in_mat[sg_info_sel2$id == sg]){
                      grid.rect(x, y, as.numeric(w)*0.8, as.numeric(h)*0.8, gp = gpar(fill = "white", col = "#82328C", lwd = 1))
                    } 
                    if(nt != ""){
                      
                      grid.text(nt, x, y, 
                                gp = gpar(fontsize = 26, col = col[nt]))
                    }
                    
                  }
                  
                  
                  
                  
                },layer_fun = function(j, i, x, y, w, h, fill) {
                  grid.lines(y = unit(c(0.5, 0.5), units = "npc"), x = unit(c(-0.15, 1), units = "npc"))
                }, col=c("white", "white"),show_heatmap_legend = F,
                
                border = F)
  
  ht <- draw(ht)
  poss <- c(-6 : 3)
  poss <- poss[poss != 0]
  
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
