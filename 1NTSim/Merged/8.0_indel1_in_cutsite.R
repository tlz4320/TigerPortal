#2.0.2之后  需要用到过滤后的sg_info_sel
#统计Ins1 Pct and Se
ins1_stat <- lapply(sg_info_sel$id, function(id){
  tmp <- total_edit_table_rev[[id]]
  res <- unlist(lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted != 0,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    ins_pct <- x[x$n_inserted == 1 & x$n_deleted == 0, ]
    return(sum(ins_pct$Pct))
  }))
  
  data.frame(Pct = mean(res), se = plotrix::std.error(res), id = id)
})
ins1_stat <- data.frame(do.call(rbind, ins1_stat))
ins1_stat <- ins1_stat[order(ins1_stat$Pct, decreasing = T),]
ins1_stat$id <- factor(ins1_stat$id, levels = ins1_stat$id)
pdf("~/Nutstore Files/Tobin/Merged1NT/Ins1_used_sgRNA_ins1_pct.pdf", width = 18, height = 4)
ggplot(ins1_stat,aes(x = id,y = Pct)) + geom_bar(stat = "identity", fill = "#06a0c9") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), width = .2, color = "black") + 
  xlab("") + ylab("%Ins1 in Indel Reads") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
dev.off()


#3.0.2之后  需要用到过滤后的sg_info_sel
#统计del1 Pct and Se
del1_stat <- lapply(sg_info_sel$id, function(id){
  tmp <- total_edit_table_rev[[id]]
  res <- unlist(lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted != 0,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    ins_pct <- x[x$n_inserted == 0 & x$n_deleted == 1, ]
    return(sum(ins_pct$Pct))
  }))
  
  data.frame(Pct = mean(res), se = plotrix::std.error(res), id = id)
})
del1_stat <- data.frame(do.call(rbind, del1_stat))
del1_stat <- del1_stat[order(del1_stat$Pct, decreasing = T),]
del1_stat$id <- factor(del1_stat$id, levels = del1_stat$id)
pdf("~/Nutstore Files/Tobin/Merged1NT/del1_used_sgRNA_del1_pct.pdf", width = 18, height = 4)
ggplot(del1_stat,aes(x = id,y = Pct)) + geom_bar(stat = "identity", fill = "#06a0c9") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), width = .2, color = "black") + 
  xlab("") + ylab("%del1 in Indel Reads") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
dev.off()

###之后统计在InDel1发生在切割位点的比例
total_sample <- unique(c(as.character(del1_stat$id), as.character(ins1_stat$id)))
indel1_stat <- lapply(total_sample, function(id){
  tmp <- total_edit_table_rev[[id]]
  res <- lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted == 1,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    x_ins <- x[x$n_inserted == 1,]
    x_del <- x[x$n_deleted == 1,]
    x_ins$pos <- unlist(lapply(x_ins$Reference_Sequence, function(seq){
      seq <- unlist(strsplit(seq, "*"))
      pos <- which(seq == "-")
      pos[which.min(abs(pos - 20))]
    }))
    x_del$pos <- unlist(lapply(x_del$Aligned_Sequence, function(seq){
      seq <- unlist(strsplit(seq, "*"))
      pos <- which(seq == "-")
      pos[which.min(abs(pos - 20))]
    }))
    x_ins_cutsite <- x_ins[x_ins$pos %in% c(20, 21),]
    x_del_cutsite <- x_del[x_del$pos %in% c(20, 21),]
    return(data.frame(ins = sum(x_ins_cutsite$Pct), del = sum(x_del_cutsite$Pct)))
  })
  res <- data.frame(do.call(rbind, res))
  indel_sum <- rowSums(res)
  data.frame(insPct = mean(res$ins), delPct = mean(res$del), indelPct = mean(indel_sum), 
             insSe = plotrix::std.error(res$ins), delSe = plotrix::std.error(res$del), 
             indelSe = plotrix::std.error(indel_sum), id = id)
})
indel1_stat <- data.frame(do.call(rbind, indel1_stat))
indel1_stat <- indel1_stat[order(indel1_stat$indelPct, decreasing = T),]
indel1_stat$id <- factor(indel1_stat$id, levels = indel1_stat$id)
pdf("~/Nutstore Files/Tobin/Merged1NT/total_used_sgRNA_indel1_in_cut.pdf", width = 20, height = 4)
ggplot(indel1_stat,aes(x = id,y = indelPct)) + geom_bar(stat = "identity", fill = "#06a0c9") + 
  geom_errorbar(aes(ymin = indelPct - indelSe, ymax = indelPct + indelSe), width = .2, color = "black") + 
  xlab("") + ylab("%InDel1 in Cut Site") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
dev.off()


###之后统计在InDel1发生在切割位点的比例 但是分母是总Indel reads
total_sample <- unique(c(as.character(del1_stat$id), as.character(ins1_stat$id)))
indel1_stat <- lapply(total_sample, function(id){
  tmp <- total_edit_table_rev[[id]]
  res <- lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted != 0,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    x_ins <- x[x$n_inserted == 1 & x$n_deleted == 0,]
    x_del <- x[x$n_deleted == 1 & x$n_inserted == 0,]
    x_ins$pos <- unlist(lapply(x_ins$Reference_Sequence, function(seq){
      seq <- unlist(strsplit(seq, "*"))
      pos <- which(seq == "-")
      pos[which.min(abs(pos - 20))]
    }))
    x_del$pos <- unlist(lapply(x_del$Aligned_Sequence, function(seq){
      seq <- unlist(strsplit(seq, "*"))
      pos <- which(seq == "-")
      pos[which.min(abs(pos - 20))]
    }))
    x_ins_cutsite <- x_ins[x_ins$pos %in% c(20, 21),]
    x_del_cutsite <- x_del[x_del$pos %in% c(20, 21),]
    return(data.frame(ins = sum(x_ins_cutsite$Pct), del = sum(x_del_cutsite$Pct)))
  })
  res <- data.frame(do.call(rbind, res))
  indel_sum <- rowSums(res)
  data.frame(insPct = mean(res$ins), delPct = mean(res$del), indelPct = mean(indel_sum), 
             insSe = plotrix::std.error(res$ins), delSe = plotrix::std.error(res$del), 
             indelSe = plotrix::std.error(indel_sum), id = id)
})
indel1_stat <- data.frame(do.call(rbind, indel1_stat))
indel1_stat <- indel1_stat[order(indel1_stat$indelPct, decreasing = T),]
indel1_stat$id <- factor(indel1_stat$id, levels = indel1_stat$id)
pdf("~/Nutstore Files/Tobin/Merged1NT/total_used_sgRNA_indel1_in_cut.pdf", width = 20, height = 4)
ggplot(indel1_stat,aes(x = id,y = indelPct)) + geom_bar(stat = "identity", fill = "#06a0c9") + 
  geom_errorbar(aes(ymin = indelPct - indelSe, ymax = indelPct + indelSe), width = .2, color = "black") + 
  xlab("") + ylab("%InDel1 in Cut Site") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
dev.off()


ins1_stat2 <- indel1_stat[,c("insPct", "insSe", "id")]
colnames(ins1_stat2) <- colnames(ins1_stat)
ins1_stat2$type <- "indel"
ins1_stat2$id <- as.character(ins1_stat2$id)
ins1_stat2 <- ins1_stat2[ins1_stat2$id %in% ins1_stat$id,]
tmp <- ins1_stat
tmp$type <- "ins"
tmp$id <- as.character(tmp$id)
ins1_stat2 <- data.frame(rbind(tmp, ins1_stat2))
ins1_stat2$id <- factor(as.character(ins1_stat2$id), levels = as.character(ins1_stat$id))
ins1_stat2$type <- factor(ins1_stat2$type, levels = c("ins", "indel"), labels = c("Total Ins1", "Ins1 in Cut Site"))


pdf("~/Nutstore Files/Tobin/Merged1NT/Ins1_used_sgRNA_ins1_pct_v2.pdf", width = 20, height = 4)
ggplot(ins1_stat2,aes(x = id,y = Pct, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), 
                color = "black", position = "dodge", width = 0.9) + 
  xlab("") + ylab("%del1 in Indel Reads") + 
  scale_fill_manual(values = c("#06a0c9" ,"#E1FFFF")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
dev.off()

del1_stat2 <- indel1_stat[,c("delPct", "delSe", "id")]
colnames(del1_stat2) <- colnames(del1_stat)
del1_stat2$type <- "indel"
del1_stat2$id <- as.character(del1_stat2$id)
del1_stat2 <- del1_stat2[del1_stat2$id %in% del1_stat$id,]
tmp <- del1_stat
tmp$type <- "del"
tmp$id <- as.character(tmp$id)
del1_stat2 <- data.frame(rbind(tmp, del1_stat2))
del1_stat2$id <- factor(as.character(del1_stat2$id), levels = as.character(del1_stat$id))
del1_stat2$type <- factor(del1_stat2$type, levels = c("del", "indel"), labels = c("Total Del1", "Del1 in Cut Site"))
pdf("~/Nutstore Files/Tobin/Merged1NT/del1_used_sgRNA_del1_pct_v2.pdf", width = 20, height = 4)
ggplot(del1_stat2,aes(x = id,y = Pct, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), 
                color = "black", position = "dodge", width = 0.9) + 
  xlab("") + ylab("%del1 in Indel Reads") + 
  scale_fill_manual(values = c("#06a0c9" ,"#E1FFFF")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
dev.off()


ins1_stat2 <- merge(ins1_stat2, sgCmp, by.x="id", by.y="id2")
ins1_order <- ins1_stat2[ins1_stat2$type == "Total Ins1",]
ins1_order <- lapply(split(ins1_order$Pct, ins1_order$id.y), function(x){
  mean(x)
})
ins1_order <- data.frame(id = names(ins1_order), pct = unlist(ins1_order))
ins1_order <- ins1_order[order(ins1_order$pct,decreasing = T),]
ins1_stat2 <- ins1_stat2[order(ins1_stat2$id),]
ins1_stat2_B <- ins1_stat2[ins1_stat2$isInsert == "Insert",]
ins1_stat2_C <- ins1_stat2[ins1_stat2$isInsert != "Insert",]
id_order <- unlist(lapply(as.character(ins1_stat2_B$id), function(x){
  ids <- sgCmp$id2[sgCmp$id == sgCmp$id[sgCmp$id2 == x]]
  ids[ids != x]
}))
ins1_stat2_C$id <- factor(as.character(ins1_stat2_C$id), levels = unique(id_order))
p1 <- ggplot(ins1_stat2_B,aes(x = id,y = Pct, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), 
                color = "black", position = "dodge", width = 0.9) + 
  xlab("") + ylab("%del1 in Indel Reads") + scale_x_discrete(position = "top")  + 
  scale_fill_manual(values = c("#06a0c9" ,"#E1FFFF")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
p2 <- ggplot(ins1_stat2_C,aes(x = id,y = Pct, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), 
                color = "black", position = "dodge", width = 0.9) + 
  xlab("") + ylab("%del1 in Indel Reads") + 
  scale_fill_manual(values = c("#06a0c9" ,"#E1FFFF")) + 
  scale_y_reverse() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
pdf("~/Nutstore Files/Tobin/Merged1NT/Ins1_used_sgRNA_ins1_pct_split_BC.pdf", width = 20, height = 10)
ggpubr::ggarrange(p1, p2, nrow = 2)
dev.off()
ins1_stat2_B$id.y <- factor(ins1_stat2_B$id.y, levels = ins1_order$id)
ins1_stat2_B <- ins1_stat2_B[order(ins1_stat2_B$id.y),]
ins1_stat2_B$id <- factor(as.character(ins1_stat2_B$id), levels = unique(as.character(ins1_stat2_B$id)))

ins1_stat2_C$id.y <- factor(ins1_stat2_C$id.y, levels = ins1_order$id)
ins1_stat2_C <- ins1_stat2_C[order(ins1_stat2_C$id.y),]
ins1_stat2_C$id <- factor(as.character(ins1_stat2_C$id), levels = unique(as.character(ins1_stat2_C$id)))


p1 <- ggplot(ins1_stat2_B,aes(x = id,y = Pct, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), 
                color = "black", position = "dodge", width = 0.9) + 
  xlab("") + ylab("%del1 in Indel Reads") + scale_x_discrete(position = "top")  + 
  scale_fill_manual(values = c("#06a0c9" ,"#E1FFFF")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
p2 <- ggplot(ins1_stat2_C,aes(x = id,y = Pct, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), 
                color = "black", position = "dodge", width = 0.9) + 
  xlab("") + ylab("%del1 in Indel Reads") + 
  scale_fill_manual(values = c("#06a0c9" ,"#E1FFFF")) + 
  scale_y_reverse() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
pdf("~/Nutstore Files/Tobin/Merged1NT/Ins1_used_sgRNA_ins1_pct_split_BC_sort_by_mean.pdf", width = 20, height = 10)
ggpubr::ggarrange(p1, p2, nrow = 2)
dev.off()
