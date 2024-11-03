###太惨了 找到的规律都被人发表了 接着2.0
high_cor_pair <- total_pair_sg_cor[total_pair_sg_cor$cor > 0.5,]
high_cor_pair_sg <- unlist(lapply(high_cor_pair$pair, function(x){
  unlist(strsplit(x, "[-]"))
}))


setwd("~/data/project/ear_project/gene_therapy/find1NTSim/")
tmp <- read.table("human_sg_final_anno2_addid.txt", sep="\t")
tmp2 <- read.table("threshold20000/human_sg_final_anno2_addid.txt", sep="\t")
tmp3 <- read.table("human_sg_lastNTSim_final_addid2.txt", sep = "\t")
tmp3$V14 <- lapply(1 : nrow(tmp3), function(x){
  if(str_sub(tmp3[x,2], 2, 20) == str_sub(tmp3[x, 8], 1, 19)){
    return(paste0(paste0(tmp3[x,2], "-"), " ", paste0("-",tmp3[x,8])))
  }
  return(paste0(paste0(tmp3[x,8], "-"), " ", paste0("-",tmp3[x,2])))
})
tmp4 <-unlist(c(tmp[,14], tmp2[,14], tmp3[,14]))
tmp <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[1]
}))
tmp2 <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[2]
}))
tmp4 <- data.frame(sgRNA = str_remove(c(tmp, tmp2), "[-]"), 
                   Cmp = c(tmp, tmp2))
tmp4 <- tmp4[!duplicated(tmp4$sgRNA),]
sgCmp <- tmp4
sgCmp <- sgCmp[sgCmp$sgRNA %in% sgRNA$sgRNA,]

high_cor_pair$cmp1 <- unlist(lapply(high_cor_pair$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[1]
  x <- sgRNA$sgRNA[sgRNA$id2 == x]
  sgCmp$Cmp[sgCmp$sgRNA == x]
}))

high_cor_pair$cmp2 <- unlist(lapply(high_cor_pair$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[2]
  x <- sgRNA$sgRNA[sgRNA$id2 == x]
  sgCmp$Cmp[sgCmp$sgRNA == x]
}))
nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_diffNT_Ins1_all.txt", sep = "\t")
colnames(nt_diff) <- c("Align", "Ref", "CmpRef", "CmpWT", "ID","Percent", "Pos","Pos2", "Diff", "editNT", "changeNTs", "editPos")
nt_diff <- nt_diff[nt_diff$Pos2 %in% c(-1, 1),]
nt_diff$side <- unlist(lapply(1 : nrow(nt_diff), function(i){
  cmp <- nt_diff$CmpRef[i]
  cmppair <- nt_diff$CmpWT[i]
  cmp <- unlist(strsplit(cmp, "*"))
  cmppair <- unlist(strsplit(cmppair, "*"))
  if(which(cmp == '-') == 1 | which(cmppair == '-') == 1){
    return("First")
  }
  return("Last")
}))
nt_diff <- nt_diff[order(nt_diff$Percent, decreasing = T),]
nt_diff <- nt_diff[nt_diff$ID %in% high_cor_pair_sg,]
nt_diff <- nt_diff[nt_diff$editNT !='N',]


pair_nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_Ins1_paired_output.txt", sep = "\t")
colnames(pair_nt_diff) <- c("Align", "Ref","CmpRef", "CmpWT","ID","Percent",  "editNT")

pair_nt_diff <- pair_nt_diff[pair_nt_diff$ID %in% high_cor_pair_sg,]
pair_nt_diff <- pair_nt_diff[pair_nt_diff$editNT != 'N',]

####在1号位置 并且最后补位的结果
nt_diff_pos1_last <-  nt_diff[nt_diff$Pos2 == 1 & nt_diff$side == 'Last',]
high_cor_pair_pos1_last <- high_cor_pair[unlist(lapply(high_cor_pair$pair, function(x){
  sum(unlist(strsplit(x, "[-]")) %in% nt_diff_pos1_last$ID) != 0
})) ,]
total_plots <- list()
for(pairids in high_cor_pair_pos1_last$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  #put insert one at front
  sg_seqs <- unlist(high_cor_pair[high_cor_pair$pair == pairids, 
                                     c("cmp1", "cmp2")])
  cmp <- sg_seqs[1]
  cmp_pair <- sg_seqs[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  #在第一位有一个插入情况下，说明后面可能是last或者非last的情况
  if(this_pos == 1 | pair_pos == 1){
    
    #说明是last位有插入 且在这个sgRNA插入
    if(pair_pos == length(cmp)){
      ids <- ids
    }
    #说明是last位有插入 且不在这个sgRNA插入
    if(this_pos == length(cmp)){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
    #说明不是last位有插入，那就看谁在第一位有-那就是后面有插入
    if(this_pos == 1){
      ids <- ids
    }else{
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
  } else {
    #说明是最后一个碱基出现了shift，那就看最早的-是出现在谁那边谁就有delete
    if(this_pos < pair_pos){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    } else{
      ids <- ids
    }
  }
  edit_one <- ids[ids %in% nt_diff_pos1_last$ID]
  paired_one <- ids[!ids %in% nt_diff_pos1_last$ID]
  
  
  edit_one_res <- nt_diff[nt_diff$ID == edit_one,]
  paired_one_res <- pair_nt_diff[pair_nt_diff$ID == paired_one,]
  
  edit_one_res_plot <- lapply(split(edit_one_res$Percent,edit_one_res$editNT), sum)
  edit_one_res_plot <- data.frame(editNT = names(edit_one_res_plot), Percent = unlist(edit_one_res_plot))
  
  paired_one_res_plot <- lapply(split(paired_one_res$Percent,paired_one_res$editNT), sum)
  paired_one_res_plot <- data.frame(editNT = names(paired_one_res_plot), Percent = unlist(paired_one_res_plot))

  
  if(edit_one == ids[1]){
    leftone <- edit_one_res_plot
  } else{
    leftone <-paired_one_res_plot
  }
  if(edit_one != ids[1]){
    rightone <- edit_one_res_plot
  } else{
    rightone <-paired_one_res_plot
  }
  max_ratio <- max(c(leftone$Percent, rightone$Percent))
  for_cor <- merge(leftone, rightone, by="editNT", all=T)
  for_cor$Percent.x[is.na(for_cor$Percent.x)] <- 0
  for_cor$Percent.y[is.na(for_cor$Percent.y)] <- 0
  
  p1 <- ggplot(leftone, aes(y = editNT, x = Percent)) + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    xlab(paste0("% in indel reads of ", ids[1])) + ylab("") + 
    scale_x_reverse(expand = c(0,0), limits = c(max_ratio, 0)) + 
    scale_y_discrete(position = "right") + 
    theme_classic2() + 
    theme(axis.text.y = element_blank(),
          plot.margin =  unit(c(0.2,0,0,0.2), "cm"))
  
   p2 <- ggplot(rightone, aes(y = editNT, x = Percent)) + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    scale_x_continuous(expand = c(0,0),  limits = c( 0, max_ratio)) + 
    scale_y_discrete() + 
    xlab(paste0( " % in indel reads of ", ids[2])) + ylab("") + 
    theme_classic2() + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5) , 
          plot.margin =  unit(c(0.2,0.75,0,-1), "cm"))
   merge_plot2 <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T)
   sg_seqs <- paste(c(paste0(c(ids[1],sg_seqs[1]), collapse = ":"), 
                      paste0(c(ids[2],sg_seqs[2]), collapse = ":")), 
                    collapse = "\n")
   total_plots[[pairids]] <- annotate_figure(merge_plot2, top = text_grob(sg_seqs, 
                                                                          color = "black", face = "bold", size = 14))
}
pdf("~/data/project/ear_project/gene_therapy_ll/NT_stat/Pos1_last_append.pdf", width = 12, height = 6)
ggpubr::ggarrange(plotlist = total_plots, nrow = 2, ncol = 3)
dev.off()


####在-1号位置 并且最开始补位的结果
nt_diff_pos_1_first <-  nt_diff[nt_diff$Pos2 == -1 & nt_diff$side == 'First',]
high_cor_pair_pos_1_first <- high_cor_pair[unlist(lapply(high_cor_pair$pair, function(x){
  sum(unlist(strsplit(x, "[-]")) %in% nt_diff_pos_1_first$ID) != 0
})) ,]
total_plots <- list()
for(pairids in high_cor_pair_pos_1_first$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  #put insert one at front
  sg_seqs <- unlist(high_cor_pair[high_cor_pair$pair == pairids, 
                                  c("cmp1", "cmp2")])
  cmp <- sg_seqs[1]
  cmp_pair <- sg_seqs[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  #在第一位有一个插入情况下，说明后面可能是last或者非last的情况
  if(this_pos == 1 | pair_pos == 1){
    
    #说明是last位有插入 且在这个sgRNA插入
    if(pair_pos == length(cmp)){
      ids <- ids
    }
    #说明是last位有插入 且不在这个sgRNA插入
    if(this_pos == length(cmp)){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
    #说明不是last位有插入，那就看谁在第一位有-那就是后面有插入
    if(this_pos == 1){
      ids <- ids
    }else{
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
  } else {
    #说明是最后一个碱基出现了shift，那就看最早的-是出现在谁那边谁就有delete
    if(this_pos < pair_pos){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    } else{
      ids <- ids
    }
  }
  edit_one <- ids[ids %in% nt_diff_pos_1_first$ID]
  paired_one <- ids[!ids %in% nt_diff_pos_1_first$ID]
  
  
  edit_one_res <- nt_diff[nt_diff$ID == edit_one,]
  paired_one_res <- pair_nt_diff[pair_nt_diff$ID == paired_one,]
  
  edit_one_res_plot <- lapply(split(edit_one_res$Percent,edit_one_res$editNT), sum)
  edit_one_res_plot <- data.frame(editNT = names(edit_one_res_plot), Percent = unlist(edit_one_res_plot))
  
  paired_one_res_plot <- lapply(split(paired_one_res$Percent,paired_one_res$editNT), sum)
  paired_one_res_plot <- data.frame(editNT = names(paired_one_res_plot), Percent = unlist(paired_one_res_plot))
  
  
  if(edit_one == ids[1]){
    leftone <- edit_one_res_plot
  } else{
    leftone <-paired_one_res_plot
  }
  leftone <- leftone[leftone$editNT %in% c("A", "T", "C", "G"),]
  if(nrow(leftone) != 4){
    nts <- setdiff(c("A", "T", "C", "G"), leftone$editNT)
    leftone <- data.frame(rbind(leftone, 
                                data.frame(editNT = nts, Percent = 0)))
  }
  if(edit_one != ids[1]){
    rightone <- edit_one_res_plot
  } else{
    rightone <-paired_one_res_plot
  }
  rightone <- rightone[rightone$editNT %in% c("A", "T", "C", "G"),]
  
  if(nrow(rightone) != 4){
    nts <- setdiff(c("A", "T", "C", "G"), rightone$editNT)
    rightone <- data.frame(rbind(rightone, 
                                data.frame(editNT = nts, Percent = 0)))
  }
  max_ratio <- max(c(leftone$Percent, rightone$Percent))
  for_cor <- merge(leftone, rightone, by="editNT", all=T)
  for_cor$Percent.x[is.na(for_cor$Percent.x)] <- 0
  for_cor$Percent.y[is.na(for_cor$Percent.y)] <- 0
  
  p1 <- ggplot(leftone, aes(y = editNT, x = Percent)) + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    xlab(paste0("% in indel reads of ", ids[1])) + ylab("") + 
    scale_x_reverse(expand = c(0,0), limits = c(max_ratio, 0)) + 
    scale_y_discrete(position = "right") + 
    theme_classic2() + 
    theme(axis.text.y = element_blank(),
          plot.margin =  unit(c(0.2,0,0,0.2), "cm"))
  
  p2 <- ggplot(rightone, aes(y = editNT, x = Percent)) + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    scale_x_continuous(expand = c(0,0),  limits = c( 0, max_ratio)) + 
    scale_y_discrete() + 
    xlab(paste0( " % in indel reads of ", ids[2])) + ylab("") + 
    theme_classic2() + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5) , 
          plot.margin =  unit(c(0.2,0.75,0,-1), "cm"))
  merge_plot2 <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T)
  sg_seqs <- paste(c(paste0(c(ids[1],sg_seqs[1]), collapse = ":"), 
                     paste0(c(ids[2],sg_seqs[2]), collapse = ":")), 
                   collapse = "\n")
  total_plots[[pairids]] <- annotate_figure(merge_plot2, top = text_grob(sg_seqs, 
                                                                         color = "black", face = "bold", size = 14))
}
pdf("~/data/project/ear_project/gene_therapy_ll/NT_stat/Pos-1_first_append.pdf", width = 12, height = 6)

ggpubr::ggarrange(plotlist = total_plots, nrow = 4, ncol = 2)
dev.off()

##第三种情况
####在1号位置 并且最开始补位的结果
nt_diff_pos1_first <-  nt_diff[nt_diff$Pos2 == 1 & nt_diff$side == 'First',]
high_cor_pair_pos1_first <- high_cor_pair[unlist(lapply(high_cor_pair$pair, function(x){
  sum(unlist(strsplit(x, "[-]")) %in% nt_diff_pos1_first$ID) != 0
})) ,]
total_plots <- list()
for(pairids in high_cor_pair_pos1_first$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  #put insert one at front
  sg_seqs <- unlist(high_cor_pair[high_cor_pair$pair == pairids, 
                                  c("cmp1", "cmp2")])
  cmp <- sg_seqs[1]
  cmp_pair <- sg_seqs[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  #在第一位有一个插入情况下，说明后面可能是last或者非last的情况
  if(this_pos == 1 | pair_pos == 1){
    
    #说明是last位有插入 且在这个sgRNA插入
    if(pair_pos == length(cmp)){
      ids <- ids
    }
    #说明是last位有插入 且不在这个sgRNA插入
    if(this_pos == length(cmp)){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
    #说明不是last位有插入，那就看谁在第一位有-那就是后面有插入
    if(this_pos == 1){
      ids <- ids
    }else{
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
  } else {
    #说明是最后一个碱基出现了shift，那就看最早的-是出现在谁那边谁就有delete
    if(this_pos < pair_pos){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    } else{
      ids <- ids
    }
  }
  edit_one <- ids[ids %in% nt_diff_pos1_first$ID]
  paired_one <- ids[!ids %in% nt_diff_pos1_first$ID]
  
  
  edit_one_res <- nt_diff[nt_diff$ID == edit_one,]
  paired_one_res <- pair_nt_diff[pair_nt_diff$ID == paired_one,]
  
  edit_one_res_plot <- lapply(split(edit_one_res$Percent,edit_one_res$editNT), sum)
  edit_one_res_plot <- data.frame(editNT = names(edit_one_res_plot), Percent = unlist(edit_one_res_plot))
  
  paired_one_res_plot <- lapply(split(paired_one_res$Percent,paired_one_res$editNT), sum)
  paired_one_res_plot <- data.frame(editNT = names(paired_one_res_plot), Percent = unlist(paired_one_res_plot))
  
  
  if(edit_one == ids[1]){
    leftone <- edit_one_res_plot
  } else{
    leftone <-paired_one_res_plot
  }
  leftone <- leftone[leftone$editNT %in% c("A", "T", "C", "G"),]
  if(nrow(leftone) != 4){
    nts <- setdiff(c("A", "T", "C", "G"), leftone$editNT)
    leftone <- data.frame(rbind(leftone, 
                                data.frame(editNT = nts, Percent = 0)))
  }
  if(edit_one != ids[1]){
    rightone <- edit_one_res_plot
  } else{
    rightone <-paired_one_res_plot
  }
  rightone <- rightone[rightone$editNT %in% c("A", "T", "C", "G"),]
  
  if(nrow(rightone) != 4){
    nts <- setdiff(c("A", "T", "C", "G"), rightone$editNT)
    rightone <- data.frame(rbind(rightone, 
                                 data.frame(editNT = nts, Percent = 0)))
  }
  max_ratio <- max(c(leftone$Percent, rightone$Percent))
  for_cor <- merge(leftone, rightone, by="editNT", all=T)
  for_cor$Percent.x[is.na(for_cor$Percent.x)] <- 0
  for_cor$Percent.y[is.na(for_cor$Percent.y)] <- 0
  
  p1 <- ggplot(leftone, aes(y = editNT, x = Percent)) + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    xlab(paste0("% in indel reads of ", ids[1])) + ylab("") + 
    scale_x_reverse(expand = c(0,0), limits = c(max_ratio, 0)) + 
    scale_y_discrete(position = "right") + 
    theme_classic2() + 
    theme(axis.text.y = element_blank(),
          plot.margin =  unit(c(0.2,0,0,0.2), "cm"))
  
  p2 <- ggplot(rightone, aes(y = editNT, x = Percent)) + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    scale_x_continuous(expand = c(0,0),  limits = c( 0, max_ratio)) + 
    scale_y_discrete() + 
    xlab(paste0( " % in indel reads of ", ids[2])) + ylab("") + 
    theme_classic2() + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5) , 
          plot.margin =  unit(c(0.2,0.75,0,-1), "cm"))
  merge_plot2 <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T)
  sg_seqs <- paste(c(paste0(c(ids[1],sg_seqs[1]), collapse = ":"), 
                     paste0(c(ids[2],sg_seqs[2]), collapse = ":")), 
                   collapse = "\n")
  total_plots[[pairids]] <- annotate_figure(merge_plot2, top = text_grob(sg_seqs, 
                                                                         color = "black", face = "bold", size = 14))
}
pdf("~/data/project/ear_project/gene_therapy_ll/NT_stat/Pos1_first_append.pdf", width = 12, height = 6)

ggpubr::ggarrange(plotlist = total_plots, nrow = 2, ncol = 2)
dev.off()



##第四种情况
####在-1号位置 并且最后补位的结果
nt_diff_pos_1_last <-  nt_diff[nt_diff$Pos2 == -1 & nt_diff$side == 'Last',]
high_cor_pair_pos_1_last <- high_cor_pair[unlist(lapply(high_cor_pair$pair, function(x){
  sum(unlist(strsplit(x, "[-]")) %in% nt_diff_pos_1_last$ID) != 0
})) ,]
total_plots <- list()
for(pairids in high_cor_pair_pos_1_last$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  #put insert one at front
  sg_seqs <- unlist(high_cor_pair[high_cor_pair$pair == pairids, 
                                  c("cmp1", "cmp2")])
  cmp <- sg_seqs[1]
  cmp_pair <- sg_seqs[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  #在第一位有一个插入情况下，说明后面可能是last或者非last的情况
  if(this_pos == 1 | pair_pos == 1){
    
    #说明是last位有插入 且在这个sgRNA插入
    if(pair_pos == length(cmp)){
      ids <- ids
    }
    #说明是last位有插入 且不在这个sgRNA插入
    if(this_pos == length(cmp)){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
    #说明不是last位有插入，那就看谁在第一位有-那就是后面有插入
    if(this_pos == 1){
      ids <- ids
    }else{
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
  } else {
    #说明是最后一个碱基出现了shift，那就看最早的-是出现在谁那边谁就有delete
    if(this_pos < pair_pos){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    } else{
      ids <- ids
    }
  }
  edit_one <- ids[ids %in% nt_diff_pos_1_last$ID]
  paired_one <- ids[!ids %in% nt_diff_pos_1_last$ID]
  
  
  edit_one_res <- nt_diff[nt_diff$ID == edit_one,]
  paired_one_res <- pair_nt_diff[pair_nt_diff$ID == paired_one,]
  
  edit_one_res_plot <- lapply(split(edit_one_res$Percent,edit_one_res$editNT), sum)
  edit_one_res_plot <- data.frame(editNT = names(edit_one_res_plot), Percent = unlist(edit_one_res_plot))
  
  paired_one_res_plot <- lapply(split(paired_one_res$Percent,paired_one_res$editNT), sum)
  paired_one_res_plot <- data.frame(editNT = names(paired_one_res_plot), Percent = unlist(paired_one_res_plot))
  
  
  if(edit_one == ids[1]){
    leftone <- edit_one_res_plot
  } else{
    leftone <-paired_one_res_plot
  }
  leftone <- leftone[leftone$editNT %in% c("A", "T", "C", "G"),]
  if(nrow(leftone) != 4){
    nts <- setdiff(c("A", "T", "C", "G"), leftone$editNT)
    leftone <- data.frame(rbind(leftone, 
                                data.frame(editNT = nts, Percent = 0)))
  }
  if(edit_one != ids[1]){
    rightone <- edit_one_res_plot
  } else{
    rightone <-paired_one_res_plot
  }
  rightone <- rightone[rightone$editNT %in% c("A", "T", "C", "G"),]
  
  if(nrow(rightone) != 4){
    nts <- setdiff(c("A", "T", "C", "G"), rightone$editNT)
    rightone <- data.frame(rbind(rightone, 
                                 data.frame(editNT = nts, Percent = 0)))
  }
  max_ratio <- max(c(leftone$Percent, rightone$Percent))
  for_cor <- merge(leftone, rightone, by="editNT", all=T)
  for_cor$Percent.x[is.na(for_cor$Percent.x)] <- 0
  for_cor$Percent.y[is.na(for_cor$Percent.y)] <- 0
  
  p1 <- ggplot(leftone, aes(y = editNT, x = Percent)) + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    xlab(paste0("% in indel reads of ", ids[1])) + ylab("") + 
    scale_x_reverse(expand = c(0,0), limits = c(max_ratio, 0)) + 
    scale_y_discrete(position = "right") + 
    theme_classic2() + 
    theme(axis.text.y = element_blank(),
          plot.margin =  unit(c(0.2,0,0,0.2), "cm"))
  
  p2 <- ggplot(rightone, aes(y = editNT, x = Percent)) + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    scale_x_continuous(expand = c(0,0),  limits = c( 0, max_ratio)) + 
    scale_y_discrete() + 
    xlab(paste0( " % in indel reads of ", ids[2])) + ylab("") + 
    theme_classic2() + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5) , 
          plot.margin =  unit(c(0.2,0.75,0,-1), "cm"))
  merge_plot2 <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T)
  sg_seqs <- paste(c(paste0(c(ids[1],sg_seqs[1]), collapse = ":"), 
                     paste0(c(ids[2],sg_seqs[2]), collapse = ":")), 
                   collapse = "\n")
  total_plots[[pairids]] <- annotate_figure(merge_plot2, top = text_grob(sg_seqs, 
                                                                         color = "black", face = "bold", size = 14))
}
pdf("~/data/project/ear_project/gene_therapy_ll/NT_stat/Pos-1_last_append.pdf", width = 12, height = 6)

ggpubr::ggarrange(plotlist = total_plots, nrow = 1, ncol = 2)
dev.off()

