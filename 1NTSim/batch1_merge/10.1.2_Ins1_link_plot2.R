#这个脚本没用 只是用来备份的
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

total_plot_data <- list()
for(id in unique(sample_info$id)){
  sgs <- sample_info$sg[sample_info$id == id]
  ins_one_seq <- wt_seqs[[sgs[1]]]
  del_one_seq <- wt_seqs[[sgs[2]]]
  ins_one_left <- ins_one_seq[20]
  del_one_left <- del_one_seq[20]
  ins_one_right <- ins_one_seq[21]
  del_one_right <- del_one_seq[21]
  
  ins_one_stat <- remain_ins1_pct_table[[sgs[1]]]
  del_one_stat <- remain_ins1_pct_table[[sgs[2]]]
  
  ifremain1 <- sum(ins_one_stat$Pct[ins_one_stat$pos %in% c(20, 21)]) / sum(ins_one_stat$Pct) * 100
  ifremain2 <- sum(del_one_stat$Pct[del_one_stat$pos %in% c(20, 21)]) / sum(del_one_stat$Pct) * 100
  if(ifremain1 < 1 | ifremain2 < 1){
    next
  }
  
  if(ins_one_left != del_one_left){
    ins_one_stat <- ins_one_stat[ins_one_stat$pos %in% c(20, 21),]
    del_one_stat <- del_one_stat[del_one_stat$pos %in% c(20, 21),]
    
    
    ins_one_pct <- sum(ins_one_stat$Pct[ins_one_stat$NT == del_one_left]) / 
      sum(ins_one_stat$Pct) * 100
    del_one_pct <- sum(del_one_stat$Pct[del_one_stat$NT == del_one_left]) / 
      sum(del_one_stat$Pct) * 100
  }else{
    ins_one_stat <- ins_one_stat[ins_one_stat$pos %in% c(20, 21),]
    del_one_stat <- del_one_stat[del_one_stat$pos %in% c(20, 21),]
    
    
    ins_one_pct <- sum(ins_one_stat$Pct[ins_one_stat$NT == ins_one_left]) / 
      sum(ins_one_stat$Pct) * 100
    del_one_pct <- sum(del_one_stat$Pct[del_one_stat$NT == del_one_left]) / 
      sum(del_one_stat$Pct) * 100
  }
  
  total_plot_data[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), type=c("+1 strand", "-1 strand"),
                                      sg = sgs, NT = ins_one_left, id = id, sameRight = (ins_one_right == del_one_right), 
                                      sameLeft = (ins_one_left == del_one_left))
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data$type <- factor(total_plot_data$type, levels = c("-1 strand", "+1 strand"))
total_plot_data$NT <- paste0(total_plot_data$NT, "/", total_plot_data$NT)


total_plot_data$ltp <- ifelse(total_plot_data$sameRight, "SameRight", "DiffRight")
total_plot_data$ltp <- factor(total_plot_data$ltp, levels = c("SameRight", "DiffRight"))

total_plot_data_sel <- total_plot_data[total_plot_data$sg %in% sg_info_sel$id,]
total_plot_data <- total_plot_data[total_plot_data$sameLeft,]
# pdf("~/data/project/ear_project/gene_therapy_ll/Result/Left_Ins1_point_link_plot_v2.pdf", width = 6, height = 4)
# ggplot(total_plot_data, aes(x = type, y = Pct, color = NT, shape = type, group = id)) + 
#   geom_line(aes(x = type, y = Pct, linetype = ltp)) + geom_point(size = 2) + xlab("") + ylab("Normalized(%)") + 
#   facet_wrap(~NT, nrow = 1) + 
#   scale_color_manual(values = c("A/A" = "#2F89FC",
#                                    "T/T" = "#30E3CA", 
#                                    "C/C" = "#66CD00", 
#                                    "G/G" = "#98ABEF")) + 
#   theme_classic2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dev.off()

total_plot_data_sel$NT[!total_plot_data_sel$sameLeft] <- "N/N"
total_plot_data_sel$NT <- factor(total_plot_data_sel$NT, levels = c("A/A", "C/C", "G/G", "T/T", "N/N"))
total_plot_data_sel$shape <- ifelse(total_plot_data_sel$type == "-1 strand", 24, 25)


pdf("~/data/project/ear_project/gene_therapy_ll/Result/Left_Ins1_point_link_plot_v4_bulge_cause.pdf", width = 6, height = 4)
ggplot(total_plot_data_sel2, aes(x = type, y = Pct, color = NT, fill = NT, group = id)) +
  geom_line(aes(x = type, y = Pct, linetype = ltp)) + geom_point(size = 2, 
                                                                 shape = total_plot_data_sel2$shape, 
                                                                 color ="black") + xlab("") + ylab("Normalized(%)") +
  facet_wrap(~NT, nrow = 1) +
  scale_fill_manual(values = c("A/A" = "#2F89FC",
                               "T/T" = "#30E3CA",
                               "C/C" = "#66CD00",
                               "G/G" = "#98ABEF", 
                               "N/N" = "#666666")) + 
  scale_color_manual(values = c("A/A" = "#2F89FC",
                                "T/T" = "#30E3CA",
                                "C/C" = "#66CD00",
                                "G/G" = "#98ABEF", 
                                "N/N" = "#666666")) +
  theme_classic2() + theme(axis.text.x = element_blank())
dev.off()