sample_info <- data.frame(sg = bulge_pos_order)
rownames(bulge_pos) <- bulge_pos$V4

sample_info$id <- unlist(lapply(sample_info$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
sample_info$bulge <- bulge_pos[sample_info$id, "V5"]
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
selCond <- unlist(lapply(split(sample_info, sample_info$id), function(xx){
  x <- unlist(xx$sg)
  sg1 <- wt_seqs[[x[1]]]
  sg2 <- wt_seqs[[x[2]]]
  id <- sample_info$id[sample_info$sg == x[1]]
  pos <- bulge_pos$V5[bulge_pos$V4 == id]
  if(pos == 4 & xx$leftSame[1] != 1){
    return(x)
  }
  return("")
}))
sameLeft <- unlist(lapply(split(sample_info, sample_info$id), function(x){
  left <- unlist(x$leftNT)
  if(left[1] == left[2]){
    return(x$sg)
  }
  return("")
}))
sel_sample <- sample_info[sample_info$leftSame != 1 | sample_info$sg %in% selCond,]
sel_sample <- sel_sample[sel_sample$sg %in% sameLeft,]
sel_sample <- sel_sample[sel_sample$id %in% unlist(lapply(split(sel_sample$id, sel_sample$id), function(x){
  if(length(x) == 2){
    return(x[1])
  }
  return("")
})),]
sel_sample <- sel_sample[sel_sample$sg %in% sg_info_sel$id,]

sel_sample$tandem <- unlist(lapply(sel_sample$sg, function(x){
  id <- sel_sample$id[sel_sample$sg == x]
  x <- sel_sample[sel_sample$id == id,]
  x$leftSame[1]
}))
sel_sample$bulgeInTandem <- unlist(lapply(sel_sample$sg, function(x){
  id <- sel_sample$id[sel_sample$sg == x]
  x <- sel_sample[sel_sample$id == id,]
  x$leftSame[1] != x$leftSame[2]
}))

sgCmp$pos[sgCmp$pos == 20] <- 1
sample_info <- data.frame(sg = bulge_pos_order)
sample_info$id <- unlist(lapply(sample_info$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
total_plot_data <- list()
for(id in unique(sample_info$id)){
  sgs <- sample_info$sg[sample_info$id == id]
  ins_one_seq <- wt_seqs[[sgs[1]]]
  del_one_seq <- wt_seqs[[sgs[2]]]
  bulgepos <- bulge_pos$V5[bulge_pos$V4 == id]
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
  
  
  ins_one_stat <- ins_one_stat[ins_one_stat$pos %in% c(20, 21),]
  del_one_stat <- del_one_stat[del_one_stat$pos %in% c(20, 21),]
  tmp_res <- calDiff(ins_one_seq, as.numeric(bulgepos), ins_one_stat,del_one_stat)
  ins_one_diff <- tmp_res$ref
  ins_one_diff <- sum(ins_one_diff$Pct / 
                        sum(ins_one_diff$Pct) * ins_one_diff$diff)
  del_one_diff <- tmp_res$target
  del_one_diff <- sum(del_one_diff$Pct / 
                        sum(del_one_diff$Pct) * del_one_diff$diff)
  total_plot_data[[id]] <- data.frame(diff = c(ins_one_diff, del_one_diff), type=c("Ref", "Target"),
                                      sg = sgs, NT = ins_one_left, id = id, sameRight = (ins_one_right == del_one_right), 
                                      sameLeft = (ins_one_left == del_one_left), 
                                      sameLR = c(ins_one_left == ins_one_right, del_one_left == del_one_right))
}


total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data$type <- factor(total_plot_data$type, levels = c("Ref", "Target"))
total_plot_data$ltp <- ifelse(total_plot_data$bulgeInTandem, "bulgeInTandem", "bulgeOutTandem")
total_plot_data$ltp <- factor(total_plot_data$ltp, levels = c("bulgeOutTandem", "bulgeInTandem"))
total_plot_data <- merge(total_plot_data, sel_sample, by="sg")
colnames(total_plot_data)[5] <- "id"
total_plot_data$shape <- ifelse(total_plot_data$type == "Target", 25, 24)
pdf("~/Nutstore Files/Tobin/Merged1NT/left_tandem_link_plot.pdf", width = 6, height = 4)
ggplot(total_plot_data, aes(x = type, y = diff, color = NT, fill = NT, group = id)) +
  geom_line(aes(x = type, y = diff, linetype = ltp)) + 
  geom_point(size = 2, 
             shape = total_plot_data$shape, 
             color ="black") + xlab("") + ylab("Normalized(%)") +
  facet_wrap(~tandem, nrow = 1) +
  scale_fill_manual(values = c("A" = "#2F89FC",
                               "T" = "#30E3CA",
                               "C" = "#66CD00",
                               "G" = "#98ABEF")) + 
  scale_color_manual(values = c("A" = "#2F89FC",
                                "T" = "#30E3CA",
                                "C" = "#66CD00",
                                "G" = "#98ABEF")) + 
  scale_y_reverse() + 
  theme_classic2()  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
