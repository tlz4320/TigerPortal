#3.0之后

sel_ids <- bulge_pos_order[seq(1,length(bulge_pos_order), 2)]
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
    edit_table <- edit_table[edit_table$n_deleted == 1 & 
                               edit_table$n_inserted == 0& 
                               edit_table$n_mutated <= 1,]
    edit_table$sgCmp <- sg
    edit_table$wtCmp <- wt_seq
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}
seqs_table <- list()
for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  
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
  seqs_table[[id]] <- edit_table[,c(1,2, 10, 11, 12, 9)]
}

seqs_table <- data.frame(do.call(rbind, seqs_table))
seqs_table <- seqs_table[seqs_table$id %in% sg_info_sel$id,]

# ToNX::write_tb(seqs_table,
#                file="~/data/project/ear_project/gene_therapy_ll/merged_seqs_table_Del1_for_heatmap.txt")



nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/merged_seqs_table_Del1_for_heatmap_output.txt", sep = "\t")
colnames(nt_diff) <- c("Align", "Ref", "CmpRef", "CmpWT", "ID","Percent", "Pos","Pos2", "Diff", "editNT", "changeNTs","editBegin")
plot_data <- lapply(split(nt_diff, nt_diff$ID), function(x){
  bulge <- unique(x$Pos)
  pos <- unlist(lapply(1 : nrow(x), function(index){
    tmp <- unlist(strsplit(x$Align[index], "*"))
    poss <- which(tmp == "-")
    poss <- poss[which.min(abs(poss - 20))]
    if(poss <= 20){
      poss <- poss - 1
    }
    poss - 20
  }))
  x <- x[abs(pos) == 1,]
  diff <- sum(x$Diff * x$Percent / sum(x$Percent))
  pos <- pos[abs(pos) == 1]
  diffpos <- unlist(lapply(1 : length(pos), function(index){
    pos2 <- x$Pos2[index]
    pos1 <- pos[index]
    ##如果3  4碱基相同时要注意可能delete和bulge其实可以算作同一个碱基
    if(pos1 != pos2){
      wt_seq <- wt_seqs[[x$ID[1]]]
      if(wt_seq[20] == wt_seq[21] & abs(pos2) == 1)
        return(0)
      
    }
    if(pos1 * pos2 < 0){
      pos <- pos2 - pos1
      if(pos < 0){
        return(pos + 1)
      }
      return(pos - 1)
    }else{
      return(pos2 - pos1)
    }
  }))
  pos <- sum(diffpos * x$Percent / sum(x$Percent))
  
  data.frame(pos = pos, diff = diff, bulge = bulge)
})


plot_data <- data.frame(do.call(rbind, plot_data))
# plot_data$pos[plot_data$pos > 0] <- plot_data$pos[plot_data$pos > 0]  - 1
table(plot_data$Pos2)
plot_data$id <- rownames(plot_data)
# plot_data$pos <- as.integer(plot_data$pos)
# plot_data <- merge(plot_data, sgCmp, by.x = "id", by.y="id2")
# colnames(plot_data)[2] <- "pos"
# pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_volcano_plot_of_bulge_distance_ntdiff_del_has0_v1.pdf", width = 8, height = 6)
# ggplot(plot_data) + geom_point(aes(x = pos, y =diff),  
#                                shape = 21, size = 3, fill =  '#DDDDDD') + 
#   scale_x_continuous(breaks = (-5 : 3)) +
#   scale_color_manual(values = c("black", "red")) + theme_bw() + xlab("") + ylab("")+
#   theme(panel.grid = element_blank(), axis.text = element_text(size = 15))
# dev.off()
plot_data$leftSame <- unlist(lapply(plot_data$id, function(x){
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
plot_data$rightSame <- unlist(lapply(plot_data$id, function(x){
  x <- unlist(x)
  sg <- wt_seqs[[x]]
  res <- 0
  while(sg[21 + res] == sg[21]){
    res <- res + 1
    if(res == 4)
      break
  }
  return(res)
}))
plot_data$shape <- 21

plot_data$shape[(plot_data$bulge > 3 & plot_data$leftSame == 2) | (plot_data$bulge < 4 & plot_data$rightSame == 2)] <- 22
plot_data$shape[(plot_data$bulge > 3 & plot_data$leftSame == 3) | (plot_data$bulge < 4 & plot_data$rightSame == 3)] <- 23
plot_data$shape[(plot_data$bulge > 3 & plot_data$leftSame == 4) | (plot_data$bulge < 4 & plot_data$rightSame == 4)] <- 24
plot_data$shape <- factor(plot_data$shape, levels = c(21:24))


pdf("~/Nutstore Files/Tobin/Merged1NT/del1_new_volcano_plot_v3.pdf", width = 8, height = 6)
ggplot(plot_data) + geom_point(aes(x = -bulge, y =diff, shape = shape),  
                               size = 3, fill =  '#DDDDDD',
                               position = position_jitter(seed = 1)) + 
  scale_x_continuous(breaks = -8 : -1, labels = c("-5", "-4", "-3", "-2", "-1", "1", "2", "3")) +
  scale_y_continuous(breaks = c(0 : 5), limits = c(0, 5)) + 
  scale_shape_manual(values = c(21,22,23,24), labels = c("No Tandem", "2 Tandem", "3 Tandem", "4 Tandem")) + 
  scale_color_manual(values = c("black", "red")) + theme_bw() + xlab("") + ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_line(colour = "#AAAAAA"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 15), 
        axis.ticks.x = element_blank())

dev.off()
