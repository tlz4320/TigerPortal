calDiff <- function(ref, pos, ins_ref, ins_target){
  mut <- ref[- (20 + 4 - pos)]
  pos <- 20
  ins_ref$diff <- 0
  for(i in 1 : nrow(ins_ref)){
    tmp <- c(mut[1 : (pos - 1)], ins_ref$NT[i], mut[pos : length(mut)])
    ins_ref$diff[i] <- sum(tmp != ref)
  }
  ins_target$diff <- 0
  for(i in 1 : nrow(ins_target)){
    tmp <- c(mut[1 : (pos - 1)], ins_target$NT[i], mut[pos : length(mut)])
    ins_target$diff[i] <- sum(tmp != ref)
  }
  return(list(ref = ins_ref, target = ins_target))
}



###在2.0脚本之后运行
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
rownames(bulge_pos) <- bulge_pos$V4
total_plot_data$bulge <- bulge_pos[total_plot_data$id, "V5"]
total_plot_data <- total_plot_data[total_plot_data$sg %in% sg_info_sel$id,]

total_plot_data$leftSame <- unlist(lapply(total_plot_data$sg, function(x){
  x <- unlist(x)
  sg <- wt_seqs[[x]]
  res <- 0
  while(sg[20 - res] == sg[20]){
    res <- res + 1
    if(res == 4)
      break
  }
  return(res)
}))
total_plot_data$tandem <- unlist(lapply(total_plot_data$sg, function(x){
  id <- total_plot_data$id[total_plot_data$sg == x]
  x <- total_plot_data[total_plot_data$id == id,]
  x$leftSame[1]
}))
tandem <- unlist(lapply(split(total_plot_data, total_plot_data$id), function(x){
  left <- unlist(x$leftSame)
  if(left[1] != 1){
    return(x$sg)
  }
  return("")
}))
total_plot_data$isTandem <- total_plot_data$sg %in% tandem

total_plot_data$bulge <- 3 - total_plot_data$bulge
# total_plot_data_tandem <- total_plot_data[total_plot_data$isTandem & 
#                                             total_plot_data$type == "Target",]
# total_plot_data_tandem$tandem <- as.character(total_plot_data_tandem$tandem)
# total_plot_data_not_tandem <- total_plot_data[!total_plot_data$isTandem & 
#                                                 total_plot_data$type == "Target",]

total_plot_data$shape <- 21
total_plot_data_target <- total_plot_data[total_plot_data$type == "Target",]

total_plot_data_target$shape[total_plot_data_target$isTandem & total_plot_data_target$tandem == 2] <- 22
total_plot_data_target$shape[total_plot_data_target$isTandem & total_plot_data_target$tandem == 3] <- 23
total_plot_data_target$shape[total_plot_data_target$isTandem & total_plot_data_target$tandem == 4] <- 24
total_plot_data_target$shape <- factor(total_plot_data_target$shape, levels = c(21:24))

pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_new_volcano_plot_v3.pdf", width = 9, height = 6)
ggplot(total_plot_data_target) + geom_point(aes(x = bulge, y =diff, 
                                                shape = shape), 
                                                position = position_jitter(seed = 1), 
                                                size = 3, fill =  '#DDDDDD') + 
  ylim(c(0,5)) + 

  scale_x_continuous(breaks = c(-5 : 2), labels = c("-5", "-4", "-3", "-2", "-1", "1", "2", "3")) +
  scale_y_continuous(breaks = c(0 : 5), limits = c(0, 5)) + 
  scale_color_manual(values = c("black", "red")) + theme_bw() + xlab("") + ylab("")+
  scale_shape_manual(values = c(21,22,23,24), labels = c("No Tandem", "2 Tandem", "3 Tandem", "4 Tandem")) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_line(colour = "#AAAAAA"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 15), 
        axis.ticks.x = element_blank())
dev.off()

# pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_new_volcano_plot.pdf", width = 8, height = 6)
# ggplot(total_plot_data_not_tandem) + geom_point(aes(x = bulge, y =diff), 
#                                position = position_jitter(seed = 1), 
#                                shape = 21, size = 3, fill =  '#DDDDDD') + 
#   ylim(c(0,5)) + 
#   scale_x_continuous(breaks = c(-5 : 2), labels = c("-5", "-4", "-3", "-2", "-1", "1", "2", "3")) +
#   scale_color_manual(values = c("black", "red")) + theme_bw() + xlab("") + ylab("")+
#   ggtitle("No Tandem") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor.x = element_line(colour = "#AAAAAA"),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         axis.text = element_text(size = 15), 
#         axis.ticks.x = element_blank())
# 
# 
# ggplot(total_plot_data_tandem) + geom_point(aes(x = bulge, y =diff, shape = tandem), 
#                                                 position = position_jitter(seed = 1), 
#                                                 size = 3, fill =  '#DDDDDD') + 
#   ylim(c(0,5)) + 
#   ggtitle("With Tandem") +
#   scale_x_continuous(breaks = c(-5 : 2), labels = c("-5", "-4", "-3", "-2", "-1", "1", "2", "3")) +
#   scale_color_manual(values = c("black", "red")) + theme_bw() + xlab("") + ylab("")+
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor.x = element_line(colour = "#AAAAAA"),
#         panel.grid.major.y = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         axis.text = element_text(size = 15), 
#         axis.ticks.x = element_blank())
# dev.off()
