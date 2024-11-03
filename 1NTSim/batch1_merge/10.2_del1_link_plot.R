#####更加细致的分析在8.2的分析之后

###最后看一下连续碱基的情况
###根据连续性来找结果
####这个画热图的版本先不用了
sample_info <- data.frame(sg = bulge_pos_order)
sample_info$id <- unlist(lapply(sample_info$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))

sel_result <- list()
for(left in 3 : 1){
  for(right in 3 : 1){
    sel <- unlist(lapply(split(sample_info$sg, sample_info$id), function(x, left, right){
      x <- unlist(x)
      sg1 <- wt_seqs[[x[1]]]
      sg2 <- wt_seqs[[x[2]]]
      for(i in (21 - left) : (20 + right)){
        if(sg1[i] != sg2[i]){
          return("")
        }
      }
      for(i in (21 - left) : 20){
        if(sg1[i] != sg1[20]){
          return("")
        }
      }
      for(i in 21 : (20 + right)){
        if(sg1[i] != sg1[21]){
          return("")
        }
      }
      if(sg1[20] == sg1[21]){
        return("")
      }
      return(x)
    }, left, right))
    sel <- sel[sel != ""]
    sel <- sel[!sel %in% unlist(sel_result)]
    if(length(sel) != 0){
      sel_result[[paste0(left, "_", right)]] <- sel
    }
  }
}
sel_sample <- data.frame(sg = bulge_pos_order)
sel_sample$id <- unlist(lapply(sel_sample$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))

plot_mat <- matrix(0, nrow = 3, ncol = 6)
rownames(plot_mat) <- as.character(3 : 1)
colnames(plot_mat) <- paste0(rep(as.character(1 : 3),c(2,2,2)), rep(c("-1", "+1"), c(3)))
for(left in 1 : 3){
  for(right in 1 : 3){
    if(left == 1 & right == 1){
      next
    }
    if(!paste0(left, "_", right) %in% names(sel_result)){
      next
    }
    sel <- sel_result[[paste0(left, "_", right)]]
    left <- as.character(left)
    right <- as.character(right)
    sel <- sel_sample[sel_sample$sg %in% sel,]
    for(id in unique(sel$id)){
      sgs <- sel_sample$sg[sel_sample$id == id]
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
      plot_mat[left, paste0(right, "-1")] <- plot_mat[left, paste0(right, "-1")] + del_one_pct
      plot_mat[left, paste0(right, "+1")] <- plot_mat[left, paste0(right, "+1")] + ins_one_pct
    }
    plot_mat[left, paste0(right, "-1")] <- plot_mat[left, paste0(right, "-1")] / length(unique(sel$id))
    plot_mat[left, paste0(right, "+1")] <- plot_mat[left, paste0(right, "+1")] / length(unique(sel$id))
  }
}
library(ComplexHeatmap)
Heatmap(plot_mat, 
        cluster_rows = F, 
        cluster_columns = F)
####这个是画连线图的版本
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
sel_sample$left <- c(1,1,1,1.8, 1.8,2.2, 1.8, 1.8, 2.2, 2.2, 2, 2, 2.8, 2.8, 3.2, 3.2, 3, 3)
sel_sample$right[10] <- 1.7
sel_sample$shape <- ifelse(sel_sample$type == "Ref", 24, 25)

sel_sample$seq <- unlist(lapply(sel_sample$sg, function(x){
  sgCmp$sgRNA[sgCmp$id2 == x]
}))
pdf("~/data/project/ear_project/gene_therapy_ll/Result/Tandem_Del1_point_link_plot_v3.pdf", width = 5.5, height = 4)

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



###周围四个碱基都相同的情况
sel_34same_sample <- data.frame(sg = bulge_pos_order)
sel_34same_sample$id <- unlist(lapply(sel_34same_sample$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))

sample_info <- data.frame(sg = bulge_pos_order)
sample_info$id <- unlist(lapply(sample_info$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
sameContext <- unlist(lapply(split(sample_info$sg, sample_info$id), function(x){
  x <- unlist(x)
  sg1 <- wt_seqs[[x[1]]]
  sg2 <- wt_seqs[[x[2]]]
  for(i in 19 : 22){
    if(sg1[i] != sg2[i]){
      return("")
    }
  }
  return(x)
}))
sample_info$sameContext <- sample_info$sg %in% sameContext



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
  total_plot_data[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), Pct2 = c(1 - ins_one_pct, 
                                                                                  1 - del_one_pct),
                                      type=c("+1 strand", "-1 strand"),
                                      sg = sgs, NT = c(ins_one_left, del_one_left))
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data <- merge(total_plot_data, sample_info, by="sg")
total_plot_data$type <- factor(total_plot_data$type, levels = c("-1 strand", "+1 strand"))
# total_plot_data$NT <- paste0(total_plot_data$NT, "/", total_plot_data$NT)
pdf("~/data/project/ear_project/gene_therapy_ll/Result/Context_Del1_point_link_plot.pdf", width = 6, height = 4)
ggplot(total_plot_data, aes(x = type, y = Pct, color = sameContext, shape = type, group = id)) + 
  geom_line(aes(x = type, y = Pct)) + geom_point(size = 2) + xlab("") + ylab("Normalized Delete Left(%)") + 
  facet_wrap(~sameContext, nrow = 1) + theme_classic2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
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

grDevices::cairo_pdf("~/data/project/ear_project/gene_therapy_ll/Result/Context_oneNT_Del1_point_link_plot_changeColor_v2_bulge_cause.pdf",
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



# grDevices::cairo_pdf("~/data/project/ear_project/gene_therapy_ll/Result/Context_oneNT_Del1_point_link_plot_changeColor.pdf", 
#     width = 6, height = 4)
# ggplot(total_plot_data, aes(x = type, y = Pct, color = sameContext, group = id)) + 
#   geom_line() + 
#   geom_point(aes(fill = deleteLeft),
#              color = "black", 
#              size = 2, 
#              shape = total_plot_data$shape) + 
#   xlab("") + 
#   scale_fill_manual(values = c(Left = "#8EC22F", Right = "#34A4DD", Third = "#910981")) + 
#   scale_y_continuous(
#     name = "<span style = 'color:#8EC22F;'>■ </span>Normalized Delete Left(%)",
#     sec.axis = sec_axis( trans=~abs(100 - .), 
#                          name="<span style = 'color:#34A4DD;'>■ </span>Normalized Delete Right(%)")
#   ) + 
#   facet_wrap(~cond, nrow = 1) + theme_classic2() + 
#   theme(axis.text.x = element_blank(), axis.title.y = element_markdown(), 
#         axis.title.y.right = element_markdown(), 
#         strip.text = element_markdown())
# dev.off()




####34不同但是 33 44 相同情况
sel_34same_sample <- data.frame(sg = bulge_pos_order)
sel_34same_sample$id <- unlist(lapply(sel_34same_sample$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
sel <- lapply(split(sel_34same_sample$sg, sel_34same_sample$id), function(x){
  x <- unlist(x)
  sg1 <- wt_seqs[[x[1]]]
  sg2 <- wt_seqs[[x[2]]]
  return(sg1[20] == sg2[20] & sg1[21] == sg2[21] & sg1[20] != sg1[21])
})
sel <- data.frame(id = names(sel), sel = unlist(sel))
sel <- sel[sel$sel,]
sel_34same_sample <- sel_34same_sample[sel_34same_sample$id %in% sel$id,]



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
  total_plot_data[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), type=c("+1 strand", "-1 strand"),
                                      sg = sgs, NT = ins_one_left, id = id)
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data$type <- factor(total_plot_data$type, levels = c("-1 strand", "+1 strand"))
total_plot_data$NT <- paste0(total_plot_data$NT, "/", total_plot_data$NT)
pdf("~/data/project/ear_project/gene_therapy_ll/Result/Left_Del1_point_link_plot.pdf", width = 6, height = 4)
ggplot(total_plot_data, aes(x = type, y = Pct, color = NT, shape = type, group = id)) + 
  geom_line(aes(x = type, y = Pct)) + geom_point(size = 2) + xlab("") + ylab("Normalized(%)") + 
  facet_wrap(~NT, nrow = 1) + theme_classic2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




sel_34same_sample <- data.frame(sg = bulge_pos_order)
sel_34same_sample$id <- unlist(lapply(sel_34same_sample$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
sel <- lapply(split(sel_34same_sample$sg, sel_34same_sample$id), function(x){
  x <- unlist(x)
  sg1 <- wt_seqs[[x[1]]]
  sg2 <- wt_seqs[[x[2]]]
  return(sg1[20] == sg2[20] & sg1[21] == sg2[21] & sg1[20] != sg1[21])
})
sel <- data.frame(id = names(sel), sel = unlist(sel))
sel <- sel[sel$sel,]
sel_34same_sample <- sel_34same_sample[sel_34same_sample$id %in% sel$id,]



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
  total_plot_data[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), type=c("+1 strand", "-1 strand"),
                                      sg = sgs, NT = ins_one_left, id = id)
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data$type <- factor(total_plot_data$type, levels = c("-1 strand", "+1 strand"))
total_plot_data$NT <- paste0(total_plot_data$NT, "/", total_plot_data$NT)
pdf("~/data/project/ear_project/gene_therapy_ll/Result/Left_Del1_point_link_plot.pdf", width = 6, height = 4)
ggplot(total_plot_data, aes(x = type, y = Pct, color = NT, shape = type, group = id)) + 
  geom_line(aes(x = type, y = Pct)) + geom_point(size = 2) + xlab("") + ylab("Normalized(%)") + 
  facet_wrap(~NT, nrow = 1) + theme_classic2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




##首先统计左边为A/T/C/G时各自删除左侧的情况
###要先选出pair sgRNA的34碱基都相同的
###并且只统计左边是A[TCG]右边不是A[TCG]的情况

sel_34same_sample <- data.frame(sg = bulge_pos_order)
sel_34same_sample$id <- unlist(lapply(sel_34same_sample$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
sel <- lapply(split(sel_34same_sample$sg, sel_34same_sample$id), function(x){
  x <- unlist(x)
  sg1 <- wt_seqs[[x[1]]]
  sg2 <- wt_seqs[[x[2]]]
  return(sg1[20] == sg2[20] & sg1[21] == sg2[21] & sg1[20] != sg1[21])
})
sel <- data.frame(id = names(sel), sel = unlist(sel))
sel <- sel[sel$sel,]
sel_34same_sample <- sel_34same_sample[sel_34same_sample$id %in% sel$id,]



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
  total_plot_data[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), type=c("+1 strand", "-1 strand"),
                                      sg = sgs, NT = ins_one_left, id = id)
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data$type <- factor(total_plot_data$type, levels = c("-1 strand", "+1 strand"))
total_plot_data$NT <- paste0(total_plot_data$NT, "/", total_plot_data$NT)
pdf("~/data/project/ear_project/gene_therapy_ll/Result/Left_Del1_point_link_plot.pdf", width = 6, height = 4)
ggplot(total_plot_data, aes(x = type, y = Pct, color = NT, shape = type, group = id)) + 
  geom_line(aes(x = type, y = Pct)) + geom_point(size = 2) + xlab("") + ylab("Normalized(%)") + 
  facet_wrap(~NT, nrow = 1) + theme_classic2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

###然后是右侧的情况

total_plot_data2 <- list()
for(id in unique(sel_34same_sample$id)){
  sgs <- sel_34same_sample$sg[sel_34same_sample$id == id]
  ins_one_seq <- wt_seqs[[sgs[1]]]
  del_one_seq <- wt_seqs[[sgs[2]]]
  ins_one_left <- ins_one_seq[21]
  del_one_left <- del_one_seq[21]
  ins_one_stat <- remain_del1_pct_table[[sgs[1]]]
  del_one_stat <- remain_del1_pct_table[[sgs[2]]]
  
  ifremain1 <- sum(ins_one_stat$Pct[ins_one_stat$pos %in% c(20, 21)]) 
  ifremain2 <- sum(del_one_stat$Pct[del_one_stat$pos %in% c(20, 21)])
  if(ifremain1 < 1 | ifremain2 < 1){
    next
  }
  
  
  ins_one_stat <- ins_one_stat[ins_one_stat$pos %in% c(20, 21),]
  del_one_stat <- del_one_stat[del_one_stat$pos %in% c(20, 21),]
  
  
  ins_one_pct <- sum(ins_one_stat$Pct[ins_one_stat$pos == 21]) / 
    sum(ins_one_stat$Pct[ins_one_stat$pos %in% c(20, 21)]) * 100
  del_one_pct <- sum(del_one_stat$Pct[del_one_stat$pos == 21]) / 
    sum(del_one_stat$Pct[del_one_stat$pos %in% c(20, 21)]) * 100
  total_plot_data2[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), type=c("+1 strand", "-1 strand"),
                                       sg = sgs, NT = ins_one_left, id = id)
}
total_plot_data2 <- data.frame(do.call(rbind, total_plot_data2))
total_plot_data2$type <- factor(total_plot_data2$type, levels = c("-1 strand", "+1 strand"))
total_plot_data2$NT <- paste0(total_plot_data2$NT, "/", total_plot_data2$NT)
pdf("~/data/project/ear_project/gene_therapy_ll/Result/Right_Del1_point_link_plot.pdf", width = 6, height = 4)
ggplot(total_plot_data2, aes(x = type, y = Pct, color = NT, shape = type, group = id)) + 
  geom_line(aes(x = type, y = Pct)) + geom_point(size = 2) + xlab("") + ylab("Normalized(%)") + 
  facet_wrap(~NT, nrow = 1) + theme_classic2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



##在之后统计左边为AA/TT/CC/GG时各自删除左侧的情况
###要先选出pair sgRNA的34碱基都相同的
###并且只统计左边是A[TCG]右边不是A[TCG]的情况


sel <- lapply(split(sel_34same_sample$sg, sel_34same_sample$id), function(x){
  x <- unlist(x)
  sg1 <- wt_seqs[[x[1]]]
  sg2 <- wt_seqs[[x[2]]]
  return(sg1[19] == sg2[20] & sg2[19] == sg2[20] & sg1[21] != sg2[22] & sg2[21] != sg2[22])
})
sel <- data.frame(id = names(sel), sel = unlist(sel))
sel <- sel[sel$sel,]
sel_34same_sample2 <- sel_34same_sample[sel_34same_sample$id %in% sel$id,]



total_plot_data <- list()
for(id in unique(sel_34same_sample2$id)){
  sgs <- sel_34same_sample2$sg[sel_34same_sample2$id == id]
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
  total_plot_data[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), type=c("+1 strand", "-1 strand"),
                                      sg = sgs, NT = ins_one_left, id = id)
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data$type <- factor(total_plot_data$type, levels = c("-1 strand", "+1 strand"))
total_plot_data$NT <- paste0(total_plot_data$NT, total_plot_data$NT, "/", total_plot_data$NT)
pdf("~/data/project/ear_project/gene_therapy_ll/Result/Left2_Del1_point_link_plot.pdf", width = 6, height = 4)
ggplot(total_plot_data, aes(x = type, y = Pct, color = NT, shape = type, group = id)) + 
  geom_line(aes(x = type, y = Pct)) + geom_point(size = 2) + xlab("") + ylab("Normalized(%)") + 
  facet_wrap(~NT, nrow = 1) + theme_classic2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


sel <- lapply(split(sel_34same_sample$sg, sel_34same_sample$id), function(x){
  x <- unlist(x)
  sg1 <- wt_seqs[[x[1]]]
  sg2 <- wt_seqs[[x[2]]]
  return(sg1[19] != sg2[20] & sg2[19] != sg2[20] & sg1[21] == sg2[22] & sg2[21] == sg2[22])
})

sel <- data.frame(id = names(sel), sel = unlist(sel))
sel <- sel[sel$sel,]
sel_34same_sample2 <- sel_34same_sample[sel_34same_sample$id %in% sel$id,]

total_plot_data2 <- list()
for(id in unique(sel_34same_sample2$id)){
  sgs <- sel_34same_sample2$sg[sel_34same_sample2$id == id]
  ins_one_seq <- wt_seqs[[sgs[1]]]
  del_one_seq <- wt_seqs[[sgs[2]]]
  ins_one_left <- ins_one_seq[21]
  del_one_left <- del_one_seq[21]
  ins_one_stat <- remain_del1_pct_table[[sgs[1]]]
  del_one_stat <- remain_del1_pct_table[[sgs[2]]]
  
  ifremain1 <- sum(ins_one_stat$Pct[ins_one_stat$pos %in% c(20, 21)])
  ifremain2 <- sum(del_one_stat$Pct[del_one_stat$pos %in% c(20, 21)])
  if(ifremain1 < 1 | ifremain2 < 1){
    next
  }
  
  
  ins_one_stat <- ins_one_stat[ins_one_stat$pos %in% c(20, 21),]
  del_one_stat <- del_one_stat[del_one_stat$pos %in% c(20, 21),]
  
  
  ins_one_pct <- sum(ins_one_stat$Pct[ins_one_stat$pos == 21]) / 
    sum(ins_one_stat$Pct[ins_one_stat$pos %in% c(20, 21)]) * 100
  del_one_pct <- sum(del_one_stat$Pct[del_one_stat$pos == 21]) / 
    sum(del_one_stat$Pct[del_one_stat$pos %in% c(20, 21)]) * 100
  total_plot_data2[[id]] <- data.frame(Pct = c(ins_one_pct, del_one_pct), type=c("+1 strand", "-1 strand"),
                                       sg = sgs, NT = ins_one_left, id = id)
}
total_plot_data2 <- data.frame(do.call(rbind, total_plot_data2))
total_plot_data2$type <- factor(total_plot_data2$type, levels = c("-1 strand", "+1 strand"))
total_plot_data2$NT <- paste0(total_plot_data2$NT,total_plot_data2$NT, "/", total_plot_data2$NT)
pdf("~/data/project/ear_project/gene_therapy_ll/Result/Right2_Del1_point_link_plot.pdf", width = 6, height = 4)
ggplot(total_plot_data2, aes(x = type, y = Pct, color = NT, shape = type, group = id)) + 
  geom_line(aes(x = type, y = Pct)) + geom_point(size = 2) + xlab("") + ylab("Normalized(%)") + 
  facet_wrap(~NT, nrow = 1) + theme_classic2() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


####最后一种，计算34位置不同 然后2345反向的情况删除左右是不是也会反过来

sel_34diff_sample <- data.frame(sg = bulge_pos_order)
sel_34diff_sample$id <- unlist(lapply(sel_34diff_sample$sg, function(x){
  sgCmp$id[sgCmp$id2 == x]
}))
sel <- lapply(sel_34diff_sample$sg, function(x){
  x <- unlist(x)
  sg1 <- wt_seqs[[x]]
  return(sg1[21] != sg1[20]) # & sg1[19] != sg1[20] & sg1[21] != sg1[22] 
})
sel <- data.frame(sg = sel_34diff_sample$sg, sel = unlist(sel))
sel <- sel[sel$sel,]
sel_34diff_sample <- sel_34diff_sample[sel_34diff_sample$sg %in% sel$sg,]



total_plot_data <- list()
for(sg in sel_34diff_sample$sg){
  wt_seq <- wt_seqs[[sg]]
  nts <- paste(wt_seq[c(19,20,21,22)], collapse = "")

  if(!nts %in% names(total_plot_data)){
    tmp_stat <- c(A = 0, 'T' = 0, 'C' = 0, 'G' = 0)
  }
  else{
    tmp_stat <- total_plot_data[[nts]]
  }
  tmp_edit_table <- remain_del1_pct_table[[sg]]
  tmp_edit_table <- tmp_edit_table[tmp_edit_table$pos %in% c(20, 21),]
  left_pct <- sum(tmp_edit_table$counts[tmp_edit_table$pos == 20])
  right_pct <- sum(tmp_edit_table$counts[tmp_edit_table$pos == 21])
  if(sum(tmp_edit_table$Pct[tmp_edit_table$pos == 20]) >= 1){
    tmp_stat[wt_seq[20]] <- tmp_stat[wt_seq[20]] + left_pct
    total_plot_data[[nts]] <- tmp_stat
  }
  if(sum(tmp_edit_table$Pct[tmp_edit_table$pos == 21]) >= 1){
    tmp_stat[wt_seq[21]] <- tmp_stat[wt_seq[21]] + right_pct
    total_plot_data[[nts]] <- tmp_stat
  }
 
}

reversed_remain <- lapply(names(total_plot_data), function(x){
  y <- paste(rev(unlist(strsplit(x, "*"))), collapse = "")
  if(y %in% names(total_plot_data)){
    return(c(x, y))
  }
  return(c("",""))
})
reversed_remain <- data.frame(do.call(rbind, reversed_remain))
reversed_remain <- reversed_remain[reversed_remain$X1 != "",]
reversed_remain$id <- unlist(lapply(1 : nrow(reversed_remain), function(i){
  x <- unlist(reversed_remain[i,])
  paste(sort(x), collapse = "-")
}))
reversed_remain <- reversed_remain[!duplicated(reversed_remain$id),]

fake_plot_matrix <- matrix(runif(2 * nrow(reversed_remain)), 
                           nrow = nrow(reversed_remain), ncol = 2)
rownames(fake_plot_matrix) <- unlist(lapply(1 : nrow(reversed_remain), function(i){
  x <- unlist(reversed_remain[i,c(1,2)])
  paste(x, collapse = "-")
}))



left_name <- rowAnnotation(left = anno_empty(border = F, width = unit(10, "mm")))
right_name <- rowAnnotation(right = anno_empty(border = F, width = unit(10, "mm")))
pdf("~/data/project/ear_project/gene_therapy_ll/Result/compare_reverse_del1_pie_chart_v2.pdf", width = 3.5, height = 10)
ht <- Heatmap(fake_plot_matrix, 
              name = "Name",
              cluster_rows = F, 
              cluster_columns = F,
              show_row_names = F,
              left_annotation = left_name,
              right_annotation = right_name,
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                nts <- reversed_remain[i, j]
                if(nts %in% names(total_plot_data)){
                  pct <- total_plot_data[[nts]]
                  pct <- unlist(pct)
                  col = c(2 : 5)
                  names(col) <- c("A", 'T', "C", "G")
                  grid.pie(pct,x = x,y = y, w= w, h=h, colors = col)
                }
                
                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16, 
         legend_gp = gpar(fill = 2:5))
)
draw(ht, annotation_legend_list = lgd_list)
decorate_annotation("left", slice = 1, {
  tg <- richtext_grob(gt_render(rev(reversed_remain$X1)), 
                      rot = 0, 
                      x = unit(0.5, "npc"),
                      y=unit(1 / nrow(reversed_remain) * (1 : nrow(reversed_remain)) - 1 / nrow(reversed_remain) / 2, "npc"), hjust = 0.5)
  grid.draw(tg)
  invisible(tg)
})

decorate_annotation("right", slice = 1, {
  tg <- richtext_grob(gt_render(rev(reversed_remain$X2)), 
                      rot = 0, 
                      x = unit(0.5, "npc"),
                      y=unit(1 / nrow(reversed_remain) * (1 : nrow(reversed_remain)) - 1 / nrow(reversed_remain) / 2, "npc"), hjust = 0.5)
  grid.draw(tg)
  invisible(tg)
})
dev.off()
