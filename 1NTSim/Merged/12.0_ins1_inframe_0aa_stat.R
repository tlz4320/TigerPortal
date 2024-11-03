library(Biostrings)
library(stringr)
source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first_second.Rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt", sep="\t")
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

seq_mat_sel <- seq_mat[sg_info_sel$id,]
plot_mat_sel <- plot_mat[sg_info_sel$id,]
sel_ids <- bulge_pos_order[seq(2,length(bulge_pos_order), 2)]
noAA_change_seq <- list()
seqs_table <- list()
sel_ids <- sel_ids[sel_ids%in% sg_info_sel$id]
for(id in sel_ids){
  
  edit_tables <- total_edit_table_rev[[id]]
  wt_seq <- sgCmp$id[sgCmp$id2 == id]
  wt_seq <- sgCmp$Cmp[sgCmp$id == wt_seq]
  sg <- sgCmp$Cmp[sgCmp$id2 == id][1]
  wt_seq <- wt_seq[wt_seq != sg]
  edit_table <- lapply(edit_tables, function(edit_table){
    edit_table <- edit_table[edit_table$Unedited == "False",]
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
    edit_table$sgCmp <- sg
    edit_table$wtCmp <- wt_seq
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}



for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  remain_table <- list()
  for(rep in names(edit_table)){
    if(nrow(edit_table[[rep]]) == 0){
      next
    }
    if(sum(edit_table[[rep]]$Pct) < 1){
      next
    }
    remain_table[[rep]] <- edit_table[[rep]]
  }
  if(length(remain_table) == 0){
    print(id)
  }
  noAA_change_seq[[id]] <- remain_table
}



####接下来要计算来自组织的Left+Right预测Insert的结果
wt_seqs <- lapply(unique(seqs_table$id), function(id){
  edit_table <- total_old_edit_table_rev[[id]][[1]]
  edit_table <- edit_table[edit_table$Unedited == "True",]
  getWtSeq(edit_table)
})
names(wt_seqs) <- unique(seqs_table$id)


sample_info <- data.frame(sg = bulge_pos_order)
sample_info$id <- unlist(lapply(sample_info$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
$pair



seqs_table_stat <- lapply(unique(sg_info_sel$pair), function(id){
  sgs <- sample_info$sg[sample_info$id == id]
  ins_one_seq <- wt_seqs[[sgs[1]]]
  del_one_seq <- wt_seqs[[sgs[2]]]
  bulgepos <- bulge_pos$V5[bulge_pos$V4 == id]
  
  seq_tables <- noAA_change_seq[[sgs[2]]]
  seq_tables <- lapply(seq_tables, function(seq_table){
    del_one_stat <- getStat2(seq_table, T)
    del_one_stat$id <- sgs[2]
    ref <- ins_one_seq
    ins_target <- del_one_stat
    mut <- ref[- (20 + 4 - bulgepos)]
    
    ins_target$diff <- 0
    for(i in 1 : nrow(ins_target)){
      pos <- ins_target$pos[i]
      if(pos > 20){
        pos <- pos - 1
      }
      tmp <- c(mut[1 : (pos - 1)], ins_target$NT[i], mut[pos : length(mut)])
      ins_target$diff[i] <- sum(tmp != ref)
    }
    ins_target
  })

})
###改成了分开样本处理，接下来是要先根据位置+NT进行一个合并用于后续计算Se

names(seqs_table_stat) <- names(noAA_change_seq)

###789三列进行加和
seqs_table_merge <- lapply(seqs_table_stat, function(seq_tables){
  lapply(seq_tables, function(seq_table){
    pos_nt <- paste(seq_table$pos, seq_table$NT, by="-")
    seq_table_split <- lapply(split(seq_table, pos_nt), function(sp){
      sp[,3] <- sum(sp[,3])
      sp[,5] <- sum(sp[,5])
      sp[1,]
    })
    data.frame(do.call(rbind, seq_table_split))
  })
})

seqs_table_mean <- lapply(names(seqs_table_merge), function(id){
  seq_tables <- seqs_table_merge[[id]]
  counts <- length(seq_tables)
  seq_table <- data.frame(do.call(rbind, seq_tables))
  pos_nt <- paste(seq_table$pos, seq_table$NT, by="-")
  seq_table_split <- lapply(split(seq_table, pos_nt), function(sp){
    sp[,3] <- sum(sp[,3]) / counts
    sp[,5] <- sum(sp[,5])/ counts
    sp[1,]
  })
  seq_table <- data.frame(do.call(rbind, seq_table_split))
  seq_table$id <- id
  seq_table
})
seqs_table_mean <- data.frame(do.call(rbind, seqs_table_mean))
openxlsx::write.xlsx(seqs_table_mean, file="~/Nutstore Files/Tobin/Merged1NT/total_ins1_0AA_stat.xlsx",rowNames=F, colNames=T)
###然后基于diff进行合并后 计算se
seqs_table_se <- lapply(names(seqs_table_merge), function(id){
  seq_tables <- seqs_table_merge[[id]]
  counts <- length(seq_tables)
  seq_tables <- lapply(seq_tables, function(seq_table){
    seq_table_split <- lapply(split(seq_table, seq_table$diff), function(sp){
      sp[,3] <- sum(sp[,3])
      sp[,5] <- sum(sp[,5])
      sp[1,]
    })
    data.frame(do.call(rbind, seq_table_split))
  })
  
  seq_table <- data.frame(do.call(rbind, seq_tables))
  
  seq_table$se <- NULL
  seq_table$groupSe <- NULL
  seq_table_split <- lapply(split(seq_table, seq_table$diff), function(sp){
    se <- sp$Pct
    groupse <- sp$GroupPct
    if(length(se) < counts){
      se <- unlist(c(se, rep(0, counts - length(se))))
      groupse <- unlist(c(groupse, rep(0, counts - length(se))))
    }
    sp$se <- plotrix::std.error(unlist(se))
    sp$groupSe <- plotrix::std.error(unlist(groupse))
    sp[,3] <- sum(sp[,3]) / counts
    sp[,5] <- sum(sp[,5])/ counts
    sp[1,]
  })
  seq_table <- data.frame(do.call(rbind, seq_table_split))
  seq_table$id <- id
  seq_table
})
names(seqs_table_se) <- names(seqs_table_merge)
seqs_table_se <- data.frame(do.call(rbind, seqs_table_se))
openxlsx::write.xlsx(seqs_table_se, file="~/Nutstore Files/Tobin/Merged1NT/total_ins1_0AA_stat_Se.xlsx",rowNames=F, colNames=T)

seqs_table_mean <- lapply(names(seqs_table_merge), function(id){
  seq_tables <- seqs_table_merge[[id]]
  counts <- length(seq_tables)
  seq_tables <- lapply(seq_tables, function(seq_table){
    if(any(seq_table$diff > 6)){
      seq_table_6 <- seq_table[seq_table$diff > 6,]
      
      seq_table_split <- lapply(split(seq_table_6, seq_table_6$NT), function(sp){
        sp[,3] <- sum(sp[,3]) 
        sp[,5] <- sum(sp[,5])
        sp[1,]
      })
      seq_table_6 <- data.frame(do.call(rbind, seq_table_split))
      #单独把大于6的拿出来
      #原因是因为本身我们关心6范围内的数据
      #另外如果把全范围画出来会非常乱
      seq_table_6$pos <- 9999
      seq_table_6$diff <- 7
      seq_table <- data.frame(rbind(seq_table_6, seq_table[seq_table$diff <= 6,]))
    }
    seq_table
  })
  
  seq_table <- data.frame(do.call(rbind, seq_tables))
  pos_nt <- paste(seq_table$pos, seq_table$NT, by="-")
  seq_table_split <- lapply(split(seq_table, pos_nt), function(sp){
    sp[,3] <- sum(sp[,3]) / counts
    sp[,5] <- sum(sp[,5])/ counts
    sp[1,]
  })
  seq_table <- data.frame(do.call(rbind, seq_table_split))
  seq_table$id <- id
  seq_table
})
seqs_table_mean <- data.frame(do.call(rbind, seqs_table_mean))
tmp <- seqs_table_mean
tmp$diff[tmp$diff == 7] <- ">6"
openxlsx::write.xlsx(tmp, file="~/Nutstore Files/Tobin/Merged1NT/total_ins1_0AA_stat_merge_large6.xlsx",rowNames=F, colNames=T)



seqs_table_se <- lapply(names(seqs_table_merge), function(id){
  seq_tables <- seqs_table_merge[[id]]
  counts <- length(seq_tables)
  seq_tables <- lapply(seq_tables, function(seq_table){
    seq_table$diff[seq_table$diff > 6] <- 7
    seq_table_split <- lapply(split(seq_table, seq_table$diff), function(sp){
      sp[,3] <- sum(sp[,3]) 
      sp[,5] <- sum(sp[,5])
      sp[1,]
    })
    data.frame(do.call(rbind, seq_table_split))
  })
  
  seq_table <- data.frame(do.call(rbind, seq_tables))
  
  seq_table$se <- NULL
  seq_table$groupSe <- NULL
  seq_table_split <- lapply(split(seq_table, seq_table$diff), function(sp){
    se <- sp$Pct
    groupse <- sp$GroupPct
    if(length(se) < counts){
      se <- unlist(c(se, rep(0, counts - length(se))))
      groupse <- unlist(c(groupse, rep(0, counts - length(se))))
    }
    sp$se <- plotrix::std.error(unlist(se))
    sp$groupSe <- plotrix::std.error(unlist(groupse))
    sp[,3] <- sum(sp[,3]) / counts
    sp[,5] <- sum(sp[,5])/ counts
    sp[1,]
  })
  seq_table <- data.frame(do.call(rbind, seq_table_split))
  seq_table$id <- id
  seq_table
})
names(seqs_table_se) <- names(seqs_table_merge)
seqs_table_se <- data.frame(do.call(rbind, seqs_table_se))
tmp <- seqs_table_se
tmp$diff[tmp$diff == 7] <- ">6"
openxlsx::write.xlsx(tmp, file="~/Nutstore Files/Tobin/Merged1NT/total_ins1_0AA_stat_Se_merge_large6.xlsx",rowNames=F, colNames=T)





cols = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(cols) <- c("A", 'T','C','G')
cols2 <- c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")
names(cols2) <- c("-2", "-1", "0", "1", "2", "Others")

plot_list <- list()

for(id in unique(seqs_table_mean$id)){
  plot_data <- inframe_result_processed[inframe_result_processed$id == id,]
  if(nrow(plot_data) != 6){
    append <- data.frame("indel" = 0,
                         "fq" = 0,
                         "pct" = 0,
                         "aa" = c("-2", "-1", "0", "1", "2", "Others"), 
                         "se" = 0,
                         "pct2" = 0,
                         "id" = id)
    append <- append[!append$aa %in% plot_data$aa,]
    plot_data <- data.frame(rbind(append, plot_data))
  }
  plot_data$aa  <- factor(plot_data$aa, 
                          levels = c("-2", "-1", "0", "1", "2", "Others"))
  inframe_plot <- ggplot(plot_data, aes(x = aa, y = pct2, fill = aa)) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymin=pct2-se, ymax=pct2+se), width=.2, color = "black") + 
    ylab("Percent in Indel") + xlab("number of AA change") + 
    scale_x_discrete(breaks = c("-2", "-1", "0", "1", "2", "Others")) + 
    scale_fill_manual(values = cols2) + 
    theme_bw() + theme(panel.grid = element_blank())
  plot_data_0aa <- seqs_table_mean[seqs_table_mean$id == id,]
  plot_data_se <- seqs_table_se[seqs_table_se$id == id,]
  if(!all(0 : 7 %in% plot_data_0aa$diff)){
    needappend <- setdiff(c(0 : 7), plot_data_0aa$diff)
    append <- lapply(needappend, function(d){
      tmp <- plot_data_0aa[1,]
      tmp$Pct <- tmp$GroupPct <- 0
      tmp$diff <- d
      tmp
    })
    append <- data.frame(do.call(rbind, append))
    plot_data_0aa <- data.frame(rbind(append, plot_data_0aa))
    append <- lapply(needappend, function(d){
      tmp <- plot_data_se[1,]
      tmp$Pct <- tmp$GroupPct <- 0
      tmp$diff <- d
      tmp$se <- 0
      tmp
    })
    append <- data.frame(do.call(rbind, append))
    plot_data_se <- data.frame(rbind(append, plot_data_se))
    
    
  }
  plot_data_0aa <- plot_data_0aa[plot_data_0aa$diff %in% 0 : 7,]
  plot_data_se <- plot_data_se[plot_data_se$diff %in% 0 : 7,]
  
  plot_data_0aa$diff <- factor(plot_data_0aa$diff, levels = c(0 : 7))
  plot_data_se$diff <- factor(plot_data_se$diff, levels = c(0 : 7))
  plot_data_0aa$NT <- stringr::str_to_upper(plot_data_0aa$NT)
  noaa_plot <- ggplot(plot_data_0aa, aes(x = diff, y = GroupPct, fill = NT)) + 
    geom_bar(stat = "identity", position = "stack") + 
    geom_errorbar(data = plot_data_se, 
                  aes(x = diff, y = GroupPct, 
                      ymin=GroupPct-se, ymax=GroupPct+se), width=.2, color = "black") + 
    ylab("Percent in 0AA") + xlab("number of NT identity change") + 
    scale_x_discrete(breaks = c(0:7), labels = c(0:6, ">6")) + 
    scale_fill_manual(values = cols) + 
    theme_bw() + theme(panel.grid = element_blank())
  dis <- unique(bulge_pos$V5[bulge_pos$V4 == sg_info_sel$pair[sg_info_sel$id == id]])
  if(dis >= 4){
    dis <- 3 - dis
  } else {
    dis <- 4 - dis
  }
  gene_name <- id
  plot_list[[id]] <- (inframe_plot + ggtitle(gene_name)) * 
    (noaa_plot + ggtitle(paste0("Bulge Position:", dis)))
}



pdf("~/Nutstore Files/Tobin/Merged1NT/total_ins1_inframe_0aa_stat.pdf", width = 8, height = 4)
for(plot in plot_list){
  print(plot)
}
dev.off()





