###在2.0脚本之后运行

sample_info <- data.frame(sg = bulge_pos_order)
sample_info$id <- unlist(lapply(sample_info$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
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
  if(ins_one_right == del_one_right & ins_one_left != del_one_left){
    ins_one_stat <- ins_one_stat[ins_one_stat$pos %in% c(20, 21),]
    del_one_stat <- del_one_stat[del_one_stat$pos %in% c(20, 21),]
    
    
    ins_one_pct <- sum(ins_one_stat$Pct[ins_one_stat$NT == ins_one_right]) / 
      sum(ins_one_stat$Pct) * 100
    del_one_pct <- sum(del_one_stat$Pct[del_one_stat$NT == ins_one_right]) / 
      sum(del_one_stat$Pct) * 100
  }
  total_plot_data[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), type=c("+1 strand", "-1 strand"),
                                      sg = sgs, NT = ins_one_left, id = id, sameRight = (ins_one_right == del_one_right), 
                                      sameLeft = (ins_one_left == del_one_left), 
                                      sameLR = c(ins_one_left == ins_one_right, del_one_left == del_one_right))
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
total_plot_data_sel$shape <- ifelse(total_plot_data_sel$type == "-1 strand", 25, 24)
total_plot_data_sel$type <- factor(total_plot_data_sel$type, levels = c("+1 strand", "-1 strand"))
total_plot_data_sel <- total_plot_data_sel[total_plot_data_sel$NT != "N/N" |
                                             (!total_plot_data_sel$sameLR & total_plot_data_sel$sameRight) 
                                             ,]
tmp_remain <- lapply(split(total_plot_data_sel$id, total_plot_data_sel$id), length)
tmp_remain <- data.frame(id = names(tmp_remain), remain = unlist(tmp_remain))
total_plot_data_sel <- total_plot_data_sel[total_plot_data_sel$id %in% tmp_remain$id[tmp_remain$remain == 2],]
pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_Left_Ins1_point_link_plot_v4.pdf", width = 6, height = 4)
ggplot(total_plot_data_sel, aes(x = type, y = Pct, color = NT, fill = NT, group = id)) +
  geom_line(aes(x = type, y = Pct, linetype = ltp)) + geom_point(size = 2, 
                                                                 shape = total_plot_data_sel$shape, 
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

total_plot_data_sel <- merge(total_plot_data_sel, sg_info_sel, by.x="sg", by.y="id")
total_plot_data_sel2 <- total_plot_data_sel[order(total_plot_data_sel$id),]
total_plot_data_sel2 <- total_plot_data_sel[total_plot_data_sel$bulge %in% c(3,4),]
pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_Left_Ins1_point_link_plot_v6_bulge_cause.pdf", width = 6, height = 4)
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
total_plot_data_sel3 <- total_plot_data_sel2[total_plot_data_sel2$NT != "N/N",]
pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_Left_Ins1_point_link_plot_rmNN_bulge_cause.pdf", width = 6, height = 4)
ggplot(total_plot_data_sel3, aes(x = type, y = Pct, color = NT, fill = NT, group = id)) +
  geom_line(aes(x = type, y = Pct, linetype = ltp)) + geom_point(size = 2, 
                                                                 shape = total_plot_data_sel3$shape, 
                                                                 color ="black") + xlab("") + ylab("Normalized(%)") +
  facet_wrap(~NT, nrow = 1) +
  scale_fill_manual(values = c("A/A" = "#2F89FC",
                               "T/T" = "#30E3CA",
                               "C/C" = "#66CD00",
                               "G/G" = "#98ABEF")) + 
  scale_color_manual(values = c("A/A" = "#2F89FC",
                                "T/T" = "#30E3CA",
                                "C/C" = "#66CD00",
                                "G/G" = "#98ABEF")) +
  theme_classic2() + theme(axis.text.x = element_blank())
dev.off()

##写出输入给indelphi用
for_indelphi <- sgCmp
for_indelphi$cutpos <- unlist(lapply(1 : nrow(for_indelphi), function(i){
  sg <- for_indelphi$sgRNA[i]
  genome <- for_indelphi$genomeSeq[i]
  unlist(str_locate(genome, sg)[1]) + 16
}))
for_indelphi <- for_indelphi[,c("id2", "genomeSeq", "cutpos")]
for_indelphi$genomeSeq <- unlist(lapply(1 : nrow(for_indelphi), function(i){
  str_sub(for_indelphi$genomeSeq[i], for_indelphi$cutpos[i] - 49, for_indelphi$cutpos[i] + 50)
}))
for_indelphi$cutpos <- 50
write.table(for_indelphi, "~/for_indelphi.csv", col.names = T, row.names = F,quote = T,sep = ",")



sample_info <- data.frame(sg = bulge_pos_order)
sample_info$id <- unlist(lapply(sample_info$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))

sample_info$leftSame <- unlist(lapply(sample_info$sg, function(x){
  x <- unlist(x)
  sg <- wt_seqs[[x]]
  res <- 0
  while(sg[20 - res] == sg[20]){
    res <- res + 1
    if(res == 3)
      break
  }
  return(res)
}))
sample_info$rightSame <- unlist(lapply(sample_info$sg, function(x){
  x <- unlist(x)
  sg <- wt_seqs[[x]]
  res <- 0
  while(sg[21 + res] == sg[21]){
    res <- res + 1
    if(res == 3)
      break
  }
  return(res)
}))
sample_info$same34 <- unlist(lapply(sample_info$sg, function(x){
  x <- unlist(x)
  sg <- wt_seqs[[x]]
  return(sg[20] == sg[21])
}))
sample_info$leftNT <- unlist(lapply(sample_info$sg, function(x){
  x <- unlist(x)
  sg <- wt_seqs[[x]]
  return(sg[20])
}))
sameContext <- unlist(lapply(split(sample_info$sg, sample_info$id), function(x){
  x <- unlist(x)
  sg1 <- wt_seqs[[x[1]]]
  sg2 <- wt_seqs[[x[2]]]
  for(i in 20 : 21){
    if(sg1[i] != sg2[i]){
      return("")
    }
  }
  return(x)
}))
sel_sample <- sample_info[!sample_info$same34,]
sel_sample <- sel_sample[sel_sample$leftSame != 1 | sel_sample$rightSame != 1,]
sel_sample <- sel_sample[sel_sample$id %in% unlist(lapply(split(sel_sample$id, sel_sample$id), function(x){
  if(length(x) == 2){
    return(x[1])
  }
  return("")
})),]
sel_sample <- sel_sample[sel_sample$sg %in% sameContext,]
sel_sample <- sel_sample[sel_sample$sg %in% sg_info_sel$id,]
sel_sample$Pct <- 0
for(index in 1 : nrow(sel_sample)){
  sg <- sel_sample$sg[index]
  wt_seq <- wt_seqs[[sg]]
  tmp_stat <- remain_ins1_pct_table[[sg]]
  leftnt <- sel_sample$leftNT[index]
  ifremain <- sum(tmp_stat$Pct[tmp_stat$pos %in% c(20, 21)]) 
  tmp_stat <- tmp_stat[tmp_stat$pos %in% c(20, 21),]
  tmp_pct <- sum(tmp_stat$Pct[tmp_stat$NT == leftnt]) / sum(tmp_stat$Pct) * 100
  sel_sample$Pct[index] <- tmp_pct
}
sel_sample$type <- rep(c("Edit","Ref"), nrow(sel_sample) / 2)
# sel_sample$leftSame <- as.numeric(sel_sample$leftSame)
# sel_sample$rightSame <- as.numeric(sel_sample$rightSame)


sel_sample$right[seq(1, nrow(sel_sample), 2)] <- sel_sample$rightSame[seq(1, nrow(sel_sample), 2)] - 0.3
sel_sample$right[seq(2, nrow(sel_sample), 2)] <- sel_sample$rightSame[seq(2, nrow(sel_sample), 2)] + 0.3
sel_sample <- sel_sample[order(sel_sample$leftSame, sel_sample$rightSame),]

sel_sample$left <- c(1.3, 1.3, 1.1, 1.1, 0.9, 0.9, 0.7, 0.7, 1, 
                     2.3, 2.3, 2.1, 2.1, 1.9, 1.9, 1.7, 1.7, 
                     1.8, 1.8, 2.2, 2,
                     3.2, 3.2, 2.8, 2.8, 3)
sel_sample$shape <- ifelse(sel_sample$type == "Ref", 24, 25)


sel_sample$seq <- unlist(lapply(sel_sample$sg, function(x){
  sgCmp$sgRNA[sgCmp$id2 == x]
}))
pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_Tandem_Ins1_point_link_plot_v1.pdf", width = 5.5, height = 4)

ggplot(sel_sample, aes(x = right, y = left)) + 
  geom_line(aes(group = id)) + 
  geom_point(aes(fill = Pct), size = 3, shape = sel_sample$shape) + 
  geom_vline(xintercept = c(1.5, 2.5)) + geom_hline(yintercept = c(1.5, 2.5)) + 
  xlab("Right NT tandem number") + ylab("Left NT tandem number") + 
  scale_x_continuous(expand = c(0,0),limits = c(0.5, 3.5), breaks = c(1, 2, 3), sec.axis = ~.) + 
  scale_y_continuous(expand = c(0,0), limits = c(0.5,3.5), breaks = c(1,2,3), sec.axis = ~.) + 
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(0,100), name = "Insert same as Left(%)") +
  theme_classic2() + theme(
    axis.text.x.top = element_blank(), 
    axis.text.y.right = element_blank(), 
    axis.ticks.x.top = element_blank(), 
    axis.ticks.y.right = element_blank()
  )
dev.off()
