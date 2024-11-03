#8.1之后

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

#这个是为了8.1数据的结果
# ToNX::write_tb(seqs_table,
               # file="~/data/project/ear_project/gene_therapy_ll/seqs_table_Ins1_for_heatmap.txt")



nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_diffNT_Ins1_heatmap2.txt", sep = "\t")
colnames(nt_diff) <- c("Align", "Ref", "CmpRef", "CmpWT", "ID","Percent", "Pos","Pos2", "Diff", "editNT", "changeNTs","editBegin")
plot_data <- lapply(split(nt_diff, nt_diff$ID), function(x){
  pos <- unlist(lapply(1 : nrow(x), function(index){
    tmp <- unlist(strsplit(x$Ref[index], "*"))
    poss <- which(tmp == "-")
    poss <- poss[which.min(abs(poss - 20))]
    if(poss > 20){
      poss <- poss - 1
    }
    poss - 20
  }))
  x <- x[pos == 0,]
  diff <- sum(x$Diff * x$Percent / sum(x$Percent))
  pos <- x$Pos2[1]
  
  data.frame(pos = pos, diff = diff)
})

plot_data <- data.frame(do.call(rbind, plot_data))
# plot_data$pos[plot_data$pos > 0] <- plot_data$pos[plot_data$pos > 0]  - 1
table(plot_data$Pos2)
plot_data$id <- rownames(plot_data)
plot_data$pos <- as.integer(plot_data$pos)
# plot_data <- merge(plot_data, sgCmp, by.x = "id", by.y="id2")
# colnames(plot_data)[2] <- "pos"
pdf("~/data/project/ear_project/gene_therapy_ll/Result/volcano_plot_of_bulge_distance_ntdiff_ins_has0_v3.pdf", width = 8, height = 6)
ggplot(plot_data) + geom_point(aes(x = pos, y =diff), 
                               position = position_jitter(seed = 1), 
                               shape = 21, size = 3, fill =  '#DDDDDD') + 
  scale_x_continuous(breaks = (-4 : 3)) +
  scale_color_manual(values = c("black", "red")) + theme_bw() + xlab("") + ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_line(colour = "#AAAAAA"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 15), 
        axis.ticks.x = element_blank())
dev.off()

###下面都是被淘汰的  没用的代码

nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_diffNT_Ins1_heatmap2.txt", sep = "\t")
colnames(nt_diff) <- c("Align", "Ref", "CmpRef", "CmpWT", "ID","Percent", "Pos","Pos2", "Diff", "editNT", "changeNTs","editBegin")
plot_data <- lapply(split(nt_diff, nt_diff$ID), function(x){
  diff <- sum(x$Diff * x$Percent / sum(x$Percent))
  data.frame(pos = x$Pos2[1], diff = diff)
})

plot_data <- data.frame(do.call(rbind, plot_data))
plot_data$pos[plot_data$pos > 0] <- plot_data$pos[plot_data$pos > 0]  - 1
table(plot_data$Pos2)
plot_data$id <- rownames(plot_data)
# plot_data <- merge(plot_data, sgCmp, by.x = "id", by.y="id2")
# colnames(plot_data)[2] <- "pos"
pdf("~/data/project/ear_project/gene_therapy_ll/Result/volcano_plot_of_bulge_distance_ntdiff_ins_mv2.pdf", width = 8, height = 6)
ggplot(plot_data) + geom_point(aes(x = pos, y =diff), 
                               position = position_jitter(seed = 1), shape = 1, size = 3) + 
  scale_x_continuous(breaks = (-5 : 2), labels = c(-5 : -1, 1 : 3)) +
  scale_color_manual(values = c("black", "red")) + theme_bw() + theme(panel.grid.major = element_blank())
dev.off()

plot_data <- lapply(split(nt_diff, nt_diff$ID), function(x){
  x <- x[order(x$Percent, decreasing = T),]
  x <- x[1,]
  x$Pos3 <- unlist(lapply(1 : nrow(x), function(i){
    align_seq <- unlist(strsplit(x$Ref[i], "*"))
    insertpos <- which(align_seq == "-")
    insertpos <- insertpos[which.min(abs(insertpos - 20))]
    insertpos <- insertpos - 20
    if(insertpos > 0){
      insertpos <- insertpos - 1
    }
    bulgepos <- x$Pos2[i]
    res <- bulgepos - insertpos
    if(res == 0){
      insertpos <- which(align_seq == "-")
      insertpos <- insertpos[which.min(abs(insertpos - 20))]
      if(insertpos > 0){
        return(1)
      }
      return(-1)
    }
    if(res < 0){
      return(res - 1)
    }
    return(res + 1)
  }))
 x
})
plot_data <- data.frame(do.call(rbind, plot_data))
plot_data$Pos4 <- unlist(lapply(plot_data$Ref, function(x){
  align_seq <- unlist(strsplit(x, "*"))
  insertpos <- which(align_seq == "-")
  insertpos <- insertpos[which.min(abs(insertpos - 20))]
  insertpos <- insertpos - 20
}))
plot_data$Pos3[plot_data$Pos3 > 0] <- plot_data$Pos3[plot_data$Pos3 > 0]  - 1
table(plot_data$Pos3)
plot_data$id <- rownames(plot_data)
# plot_data <- merge(plot_data, sgCmp, by.x = "id", by.y="id2")
# colnames(plot_data)[2] <- "pos"
pdf("~/data/project/ear_project/gene_therapy_ll/Result/volcano_plot_of_bulge_distance_ntdiff_ins_test.pdf", width = 8, height = 6)
ggplot(plot_data) + geom_point(aes(x = Pos3, y =Diff), 
                               position = position_jitter(seed = 1), shape = 1, size = 3) + 
  scale_x_continuous(breaks = (-6 : 3), labels = c(-6 : -1, 1 : 4)) +
  scale_color_manual(values = c("black", "red")) + theme_bw() + theme(panel.grid.major = element_blank())
dev.off()



plot_data$Pos2[plot_data$Pos2 > 0] <- plot_data$Pos2[plot_data$Pos2 > 0]  - 1
table(plot_data$Pos2)
plot_data$id <- rownames(plot_data)
# plot_data <- merge(plot_data, sgCmp, by.x = "id", by.y="id2")
# colnames(plot_data)[2] <- "pos"
pdf("~/data/project/ear_project/gene_therapy_ll/Result/volcano_plot_of_bulge_distance_ntdiff_ins_test.pdf", width = 8, height = 6)
ggplot(plot_data) + geom_point(aes(x = Pos2, y =Diff), 
                               position = position_jitter(seed = 1), shape = 1, size = 3) + 
  scale_x_continuous(breaks = (-6 : 3), labels = c(-6 : -1, 1 : 4)) +
  scale_color_manual(values = c("black", "red")) + theme_bw() + theme(panel.grid.major = element_blank())
dev.off()



#####最后选择了这个部分的结果
nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_diffNT_Ins1_heatmap2.txt", sep = "\t")
colnames(nt_diff) <- c("Align", "Ref", "CmpRef", "CmpWT", "ID","Percent", "Pos","Pos2", "Diff", "editNT", "changeNTs","editBegin")
plot_data <- lapply(split(nt_diff, nt_diff$ID), function(x){
  sel_indel_in_cutsite <- unlist(lapply(1 : nrow(x), function(i){
    align_seq <- unlist(strsplit(x$Ref[i], "*"))
    insertpos <- which(align_seq == "-")
    insertpos <- insertpos[which.min(abs(insertpos - 20))]
    return(insertpos %in% c(20, 21))
  }))
  x <- x[sel_indel_in_cutsite,]
  diff <- sum(x$Diff * x$Percent / sum(x$Percent))
  data.frame(pos = x$Pos2[1], diff = diff)
})

plot_data <- data.frame(do.call(rbind, plot_data))
# plot_data$pos[plot_data$pos > 0] <- plot_data$pos[plot_data$pos > 0]  - 1
table(plot_data$Pos2)
plot_data$id <- rownames(plot_data)
# plot_data <- merge(plot_data, sgCmp, by.x = "id", by.y="id2")
# colnames(plot_data)[2] <- "pos"
pdf("~/data/project/ear_project/gene_therapy_ll/Result/volcano_plot_of_bulge_distance_ntdiff_ins_has0.pdf", width = 8, height = 6)
ggplot(plot_data) + geom_point(aes(x = pos, y =diff), 
                               position = position_jitter(seed = 1), shape = 1, size = 3) + 
  scale_x_continuous(breaks = (-4 : 3)) +
  scale_color_manual(values = c("black", "red")) + theme_bw() + theme(panel.grid.major = element_blank())
dev.off()


plot_data <- lapply(split(nt_diff, nt_diff$ID), function(x){
  x <- x[order(x$Percent, decreasing = T),]
  x <- x[1,]
  x$Pos3 <- unlist(lapply(1 : nrow(x), function(i){
    align_seq <- unlist(strsplit(x$Ref[i], "*"))
    insertpos <- which(align_seq == "-")
    insertpos <- insertpos[which.min(abs(insertpos - 20))]
    insertpos <- insertpos - 20
    if(insertpos > 0){
      insertpos <- insertpos - 1
    }
    bulgepos <- x$Pos2[i]
    res <- bulgepos - insertpos
    if(res == 0){
      insertpos <- which(align_seq == "-")
      insertpos <- insertpos[which.min(abs(insertpos - 20))]
      if(insertpos > 0){
        return(1)
      }
      return(-1)
    }
    if(res < 0){
      return(res - 1)
    }
    return(res + 1)
  }))
  x
})
plot_data <- data.frame(do.call(rbind, plot_data))
plot_data$Pos4 <- unlist(lapply(plot_data$Ref, function(x){
  align_seq <- unlist(strsplit(x, "*"))
  insertpos <- which(align_seq == "-")
  insertpos <- insertpos[which.min(abs(insertpos - 20))]
  insertpos <- insertpos - 20
}))
plot_data$Pos3[plot_data$Pos3 > 0] <- plot_data$Pos3[plot_data$Pos3 > 0]  - 1
table(plot_data$Pos3)
plot_data$id <- rownames(plot_data)
# plot_data <- merge(plot_data, sgCmp, by.x = "id", by.y="id2")
# colnames(plot_data)[2] <- "pos"
pdf("~/data/project/ear_project/gene_therapy_ll/Result/volcano_plot_of_bulge_distance_ntdiff_ins_test.pdf", width = 8, height = 6)
ggplot(plot_data) + geom_point(aes(x = Pos3, y =Diff), 
                               position = position_jitter(seed = 1), shape = 1, size = 3) + 
  scale_x_continuous(breaks = (-6 : 3), labels = c(-6 : -1, 1 : 4)) +
  scale_color_manual(values = c("black", "red")) + theme_bw() + theme(panel.grid.major = element_blank())
dev.off()



plot_data$Pos2[plot_data$Pos2 > 0] <- plot_data$Pos2[plot_data$Pos2 > 0]  - 1
table(plot_data$Pos2)
plot_data$id <- rownames(plot_data)
# plot_data <- merge(plot_data, sgCmp, by.x = "id", by.y="id2")
# colnames(plot_data)[2] <- "pos"
pdf("~/data/project/ear_project/gene_therapy_ll/Result/volcano_plot_of_bulge_distance_ntdiff_ins_test.pdf", width = 8, height = 6)
ggplot(plot_data) + geom_point(aes(x = Pos2, y =Diff), 
                               position = position_jitter(seed = 1), shape = 1, size = 3) + 
  scale_x_continuous(breaks = (-6 : 3), labels = c(-6 : -1, 1 : 4)) +
  scale_color_manual(values = c("black", "red")) + theme_bw() + theme(panel.grid.major = element_blank())
dev.off()





