#####更加细致的分析在3.0的分析之后

###最后看一下连续碱基的情况
####画连线图的版本
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
sel_sample <- data.frame(sg = bulge_pos_order)
sel_sample$id <- unlist(lapply(sel_sample$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))

sel_sample <- sel_sample[sel_sample$sg %in% sg_info_sel$id,]
sel_sample <- sample_info[!sample_info$same34,]
sel_sample <- sel_sample[sel_sample$leftSame != 1 | sel_sample$rightSame != 1,]
sel_sample <- sel_sample[sel_sample$id %in% unlist(lapply(split(sel_sample$id, sel_sample$id), function(x){
  if(length(x) == 2){
    return(x[1])
  }
  return("")
})),]
sel_sample <- sel_sample[sel_sample$sg %in% sameContext,]

sel_sample$Pct <- 0
for(index in 1 : nrow(sel_sample)){
  sg <- sel_sample$sg[index]
  wt_seq <- wt_seqs[[sg]]
  tmp_stat <- remain_del1_pct_table[[sg]]
  
  ifremain <- sum(tmp_stat$Pct[tmp_stat$pos %in% c(20, 21)]) 
  tmp_stat <- tmp_stat[tmp_stat$pos %in% c(20, 21),]
  tmp_pct <- sum(tmp_stat$Pct[tmp_stat$pos == 20]) / sum(tmp_stat$Pct) * 100
  sel_sample$Pct[index] <- tmp_pct
}
sel_sample$type <- rep(c("Edit","Ref"), nrow(sel_sample) / 2)
# sel_sample$leftSame <- as.numeric(sel_sample$leftSame)
# sel_sample$rightSame <- as.numeric(sel_sample$rightSame)


sel_sample$right[seq(1, nrow(sel_sample), 2)] <- sel_sample$rightSame[seq(1, nrow(sel_sample), 2)] - 0.3
sel_sample$right[seq(2, nrow(sel_sample), 2)] <- sel_sample$rightSame[seq(2, nrow(sel_sample), 2)] + 0.3
sel_sample <- sel_sample[order(sel_sample$leftSame, sel_sample$rightSame),]
sel_sample$left <- c(1.3,1.3,1.1,1.1, 0.9,0.9, 0.7, 0.7, 1, 2.4, 2.4, 2.2, 2.2, 2, 2, 1.8, 1.8, 1.6, 1.6,
                     2.3, 1.8, 1.8, 2.3, 2.3, 2, 2, 3.2, 3.2, 2.8, 2.8, 3, 3)
sel_sample$right[24] <- 1.7
sel_sample$shape <- ifelse(sel_sample$type == "Ref", 24, 25)

sel_sample$seq <- unlist(lapply(sel_sample$sg, function(x){
  sgCmp$sgRNA[sgCmp$id2 == x]
}))
pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_Tandem_Del1_point_link_plot_v1.pdf", width = 5.5, height = 4)

ggplot(sel_sample, aes(x = right, y = left)) + 
  geom_line(aes(group = id)) + 
  geom_point(aes(fill = Pct), size = 3, shape = sel_sample$shape) + 
  geom_vline(xintercept = c(1.5, 2.5)) + geom_hline(yintercept = c(1.5, 2.5)) + 
  xlab("Right NT tandem number") + ylab("Left NT tandem number") + 
  scale_x_continuous(expand = c(0,0),limits = c(0.5, 3.5), breaks = c(1, 2, 3), sec.axis = ~.) + 
  scale_y_continuous(expand = c(0,0), limits = c(0.5,3.5), breaks = c(1,2,3), sec.axis = ~.) + 
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(0,100), name = "Delete in Left(%)") +
  theme_classic2() + theme(
    axis.text.x.top = element_blank(), 
    axis.text.y.right = element_blank(), 
    axis.ticks.x.top = element_blank(), 
    axis.ticks.y.right = element_blank()
  )
dev.off()



###周围2个碱基都相同的情况

sel_34same_sample <- data.frame(sg = bulge_pos_order)
sel_34same_sample$id <- unlist(lapply(sel_34same_sample$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
sample_info <- data.frame(sg = bulge_pos_order)
sample_info$id <- unlist(lapply(sample_info$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
sameContext <- unlist(lapply(unique(sample_info$id), function(x){
  x <- sample_info$sg[sample_info$id == x]
  x <- unlist(x)
  sg1 <- wt_seqs[[x[1]]]
  sg2 <- wt_seqs[[x[2]]]
  if(sg1[20] == sg2[20]){
    if(sg1[21] != sg2[21]){
      return("Right Diff")
    }
    return("Same")
  } 
  if(sg1[21] == sg2[21]){
    return("Left Diff")
  }
  
  return("Total Diff")
}))
sameContext <- data.frame(id = unique(sample_info$id), sameContext = sameContext)
sample_info <- merge(sample_info, sameContext, by="id")
sample_info <- sample_info[sample_info$sameContext != "Total Diff",]



total_plot_data <- list()
for(id in unique(sel_34same_sample$id)){
  sgs <- sel_34same_sample$sg[sel_34same_sample$id == id]
  ins_one_seq <- wt_seqs[[sgs[1]]]
  del_one_seq <- wt_seqs[[sgs[2]]]
  ins_one_left <- ins_one_seq[20]
  del_one_left <- del_one_seq[20]
  ins_one_stat <- remain_del1_pct_table[[sgs[1]]]
  del_one_stat <- remain_del1_pct_table[[sgs[2]]]
  
  ifremain1 <- sum(ins_one_stat$Pct[ins_one_stat$pos %in% c(20, 21)]) 
  ifremain2 <- sum(del_one_stat$Pct[del_one_stat$pos %in% c(20, 21)])
  if(ifremain1 < 1 | ifremain2 < 1){
    next
  }
  
  
  ins_one_stat <- ins_one_stat[ins_one_stat$pos %in% c(20, 21),]
  del_one_stat <- del_one_stat[del_one_stat$pos %in% c(20, 21),]
  
  
  ins_one_pct <- sum(ins_one_stat$Pct[ins_one_stat$pos == 20]) / 
    sum(ins_one_stat$Pct[ins_one_stat$pos %in% c(20, 21)]) * 100
  del_one_pct <- sum(del_one_stat$Pct[del_one_stat$pos == 20]) / 
    sum(del_one_stat$Pct[del_one_stat$pos %in% c(20, 21)]) * 100
  total_plot_data[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), Pct2 = c(100 - ins_one_pct, 
                                                                                  100 - del_one_pct),
                                      type=c("+1 strand", "-1 strand"),
                                      sg = sgs, NT = c(ins_one_left, del_one_left))
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data <- merge(total_plot_data, sample_info, by="sg")
total_plot_data$type <- factor(total_plot_data$type, levels = c("-1 strand", "+1 strand"))
total_plot_data$shape <- ifelse(total_plot_data$type == "-1 strand", 24, 25)
total_plot_data$sameContext <- factor(total_plot_data$sameContext, levels = c("Same", "Right Diff", "Left Diff"))
total_plot_data$cond <- ""
total_plot_data$cond[total_plot_data$sameContext == "Same"] <- 
  "<span style = 'color:#8EC22F;'>■ </span>|<span style = 'color:#34A4DD;'>■</span>  <span style = 'color:#8EC22F;'>■ </span>|<span style = 'color:#34A4DD;'>■</span>"
total_plot_data$cond[total_plot_data$sameContext == "Right Diff"] <- 
  "<span style = 'color:#8EC22F;'>■ </span>|<span style = 'color:#34A4DD;'>■</span>  <span style = 'color:#8EC22F;'>■ </span>|<span style = 'color:#910981;'>■</span>"
total_plot_data$cond[total_plot_data$sameContext == "Left Diff"] <- 
  "<span style = 'color:#910981;'>■ </span>|<span style = 'color:#34A4DD;'>■</span>  <span style = 'color:#8EC22F;'>■ </span>|<span style = 'color:#34A4DD;'>■</span>"


total_plot_data$Pct2 <- total_plot_data$Pct
total_plot_data$deleteLeft <- "Left"
total_plot_data$deleteLeft[total_plot_data$Pct < 50] <- "Right"

total_plot_data$deleteLeft[total_plot_data$Pct < 50 & 
                             total_plot_data$type == "+1 strand" & 
                             total_plot_data$sameContext == "Right Diff"] <- "Third"

total_plot_data$deleteLeft[total_plot_data$Pct > 50 & 
                             total_plot_data$type == "-1 strand" & 
                             total_plot_data$sameContext == "Left Diff"] <- "Third"

total_plot_data$Pct2[total_plot_data$Pct2 < 50] <- 100 - total_plot_data$Pct2[total_plot_data$Pct2 < 50]

# total_plot_data$NT <- paste0(total_plot_data$NT, "/", total_plot_data$NT)
library(ggtext)

total_plot_data_sel <- merge(total_plot_data, sg_info_sel, by.x="sg", by.y='id')

grDevices::cairo_pdf("~/data/project/ear_project/gene_therapy_ll/Result/Context_oneNT_Del1_point_link_plot_changeColor_v2.pdf",
                     width = 6, height = 4)
ggplot(total_plot_data_sel, aes(x = type, y = Pct, color = sameContext, group = id)) +
  geom_line() +
  geom_point(aes(fill = deleteLeft),
             color = "black",
             size = 2,
             shape = total_plot_data_sel$shape) +
  xlab("") +
  scale_fill_manual(values = c(Left = "#8EC22F", Right = "#34A4DD", Third = "#910981")) +
  scale_y_continuous(
    name = "<span style = 'color:#8EC22F;'>■ </span>Normalized Delete Left(%)",
    sec.axis = sec_axis( trans=~abs(100 - .),
                         name="<span style = 'color:#34A4DD;'>■ </span>Normalized Delete Right(%)")
  ) +
  facet_wrap(~cond, nrow = 1) + theme_classic2() +
  theme(axis.text.x = element_blank(), axis.title.y = element_markdown(),
        axis.title.y.right = element_markdown(),
        strip.text = element_markdown())
dev.off()


total_plot_data_sel2 <- total_plot_data_sel
total_plot_data_sel2$isAppendLast <- unlist(lapply(total_plot_data_sel2$pair, function(x){
  if(!is.na(unlist(str_match(x, "last")))){
    return(F)
  }
  tmp <- sgCmp[sgCmp$id == x,]
  
  a <- unlist(strsplit(tmp$Cmp[1], "*"))
  b <- unlist(strsplit(tmp$Cmp[2], "*"))
  if(a[length(a)] == '-' | b[length(b)] == '-'){
    return(T)
  }
  return(F)
}))
total_plot_data_sel2 <- total_plot_data_sel2[(total_plot_data_sel2$bulge >= 6 & 
                                                !total_plot_data_sel2$isAppendLast) |
                                               (total_plot_data_sel2$bulge == 3 & 
                                                  total_plot_data_sel2$isAppendLast & 
                                                  total_plot_data_sel2$sameContext != "Same") | 
                                               (total_plot_data_sel2$bulge == 4 & 
                                                  !total_plot_data_sel2$isAppendLast & 
                                                  total_plot_data_sel2$sameContext != "Same"),]

total_plot_data_sel2 <- total_plot_data_sel2[order(total_plot_data_sel2$bulge),]

grDevices::cairo_pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_Context_oneNT_Del1_point_link_plot_changeColor_v2_bulge_cause.pdf",
                     width = 6, height = 4)
ggplot(total_plot_data_sel2, aes(x = type, y = Pct, color = sameContext, group = id)) +
  geom_line() +
  geom_point(aes(fill = deleteLeft),
             color = "black",
             size = 2,
             shape = total_plot_data_sel2$shape) +
  xlab("") +
  scale_fill_manual(values = c(Left = "#8EC22F", Right = "#34A4DD", Third = "#910981")) +
  scale_y_continuous(
    name = "<span style = 'color:#8EC22F;'>■ </span>Normalized Delete Left(%)",
    sec.axis = sec_axis( trans=~abs(100 - .),
                         name="<span style = 'color:#34A4DD;'>■ </span>Normalized Delete Right(%)")
  ) +
  facet_wrap(~cond, nrow = 1) + theme_classic2() +
  theme(axis.text.x = element_blank(), axis.title.y = element_markdown(),
        axis.title.y.right = element_markdown(),
        strip.text = element_markdown())
dev.off()