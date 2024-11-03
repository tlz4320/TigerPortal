load("~/data/project/ear_project/gene_therapy_ll/Result/first_total_edit_table.rda")
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
nt_diff <- nt_diff[str_length(nt_diff$editNT) == 1,]

nt_diff$type <- paste(nt_diff$Pos2, nt_diff$side, sep = " ")

pair_nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_Ins1_paired_output.txt", sep = "\t")
colnames(pair_nt_diff) <- c("Align", "Ref","CmpRef", "CmpWT","ID","Percent",  "editNT", "editPos")

pair_nt_diff <- pair_nt_diff[pair_nt_diff$ID %in% high_cor_pair_sg,]
pair_nt_diff <- pair_nt_diff[pair_nt_diff$editNT != 'N',]
pair_nt_diff <- pair_nt_diff[str_length(pair_nt_diff$editNT) == 1,]

###上面处理好原始数据后 下面开始统计了
high_cor_pair <- total_pair_sg_cor[total_pair_sg_cor$cor > 0.5,]
high_cor_pair_sg <- unlist(lapply(high_cor_pair$pair, function(x){
  unlist(strsplit(x, "[-]"))
}))
high_cor_pair_sel <- high_cor_pair[unlist(lapply(high_cor_pair$pair, function(x){
  sum(unlist(strsplit(x, "[-]")) %in% nt_diff$ID) != 0
})) ,]
total_stat <- list()
for(pairids in high_cor_pair_sel$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  #put insert one at front
  

  edit_one <- ids[ids %in% nt_diff$ID]
  paired_one <- ids[!ids %in% nt_diff$ID]
  edit_one_edit_table <- total_edit_table_rev[[edit_one]][[1]]
  pair_one_edit_table <- total_edit_table_rev[[paired_one]][[1]]
  edit_one_wt <- getWtSeq(edit_one_edit_table[edit_one_edit_table$Unedited == "True",])
  pair_one_wt <- getWtSeq(pair_one_edit_table[pair_one_edit_table$Unedited == "True",])
  
  edit_left <- edit_one_wt[20]
  edit_right <- edit_one_wt[21]
  pair_left <- pair_one_wt[20]
  pair_right <- pair_one_wt[21]
  edit_one_res <- nt_diff[nt_diff$ID == edit_one,]
  type <- unique(edit_one_res$type)[1]
  paired_one_res <- pair_nt_diff[pair_nt_diff$ID == paired_one,]
  
  edit_one_res_plot <- lapply(split(edit_one_res$Percent,edit_one_res$editNT), sum)
  edit_one_res_plot <- data.frame(editNT = names(edit_one_res_plot), Percent = unlist(edit_one_res_plot))
  
  paired_one_res_plot <- lapply(split(paired_one_res$Percent,paired_one_res$editNT), sum)
  paired_one_res_plot <- data.frame(editNT = names(paired_one_res_plot), Percent = unlist(paired_one_res_plot))
  
  

  leftone <- edit_one_res_plot

  leftone <- leftone[leftone$editNT %in% c("A", "T", "C", "G"),]
  if(nrow(leftone) != 4){
    nts <- setdiff(c("A", "T", "C", "G"), leftone$editNT)
    leftone <- data.frame(rbind(leftone, 
                                data.frame(editNT = nts, Percent = 0)))
  }

  rightone <- paired_one_res_plot

  rightone <- rightone[rightone$editNT %in% c("A", "T", "C", "G"),]
  
  if(nrow(rightone) != 4){
    nts <- setdiff(c("A", "T", "C", "G"), rightone$editNT)
    rightone <- data.frame(rbind(rightone, 
                                 data.frame(editNT = nts, Percent = 0)))
  }
  
  leftone$norm <- leftone$Percent / sum(leftone$Percent) * 100
  rightone$norm <- rightone$Percent / sum(rightone$Percent) * 100
  
  leftone$leftNT <- edit_left
  leftone$rightNT <- edit_right
  rightone$leftNT <- edit_left
  rightone$rightNT <- edit_right
  leftone$group <- "edit"
  rightone$group <- "ref"
  leftone$pair <- pairids
  rightone$pair <- pairids
  leftone$type <- type
  rightone$type <- type
  total_stat[[length(total_stat) + 1]] <- data.frame(rbind(leftone, rightone))
}

total_stat <- data.frame(do.call(rbind, total_stat))

plot_data_left <- total_stat

pdf("~/data/project/ear_project/gene_therapy_ll/NT_stat/insert_point_plot_1_last.pdf", width = 6, height = 4)
plot_data_left <- plot_data_left[plot_data_left$type %in% c("1 Last", "-1 First"),]
plot_data_left$editNT2 <- paste0(plot_data_left$editNT, plot_data_left$group)
plot_data_tmp <- plot_data_left[plot_data_left$type == "1 Last",]
plot_data_tmp$editNT2 <- factor(plot_data_tmp$editNT2, levels = c("Aedit", "Aref", "Cedit","Cref", "Gedit","Gref", "Tedit","Tref"))
plot_data_tmp$editNT2 <- as.numeric(plot_data_tmp$editNT2)
ggplot(plot_data_tmp, aes(x = editNT2, y = norm, color = pair, shape = group, group = paste0(pair, editNT))) + 
  geom_line(aes(x = editNT2, y = norm)) + geom_point(size = 2) + 
  scale_x_continuous(breaks = c(1.5,3.5, 5.5, 7.5), labels = c('A','C','G','T')) + 
  facet_wrap(~leftNT, nrow = 1) + theme_classic2()
  
dev.off()

pdf("~/data/project/ear_project/gene_therapy_ll/NT_stat/insert_point_plot_-1_first.pdf", width = 6, height = 4)
plot_data_left <- plot_data_left[plot_data_left$type %in% c("1 Last", "-1 First"),]
plot_data_left$editNT2 <- paste0(plot_data_left$editNT, plot_data_left$group)
plot_data_tmp <- plot_data_left[plot_data_left$type == "-1 First",]
plot_data_tmp$editNT2 <- factor(plot_data_tmp$editNT2, levels = c("Aedit", "Aref", "Cedit","Cref", "Gedit","Gref", "Tedit","Tref"))
plot_data_tmp$editNT2 <- as.numeric(plot_data_tmp$editNT2)
ggplot(plot_data_tmp, aes(x = editNT2, y = norm, color = pair, shape = group, group = paste0(pair, editNT))) + 
  geom_line(aes(x = editNT2, y = norm)) + geom_point(size = 2) + 
  scale_x_continuous(breaks = c(1.5,3.5, 5.5, 7.5), labels = c('A','C','G','T')) + 
  facet_wrap(~leftNT, nrow = 1) + theme_classic2()

dev.off()
