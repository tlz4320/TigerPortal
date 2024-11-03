print(load("~/data/project/ear_project/gene_therapy_ll/Previews/Result/Rep123_indel_stat.rda"))
sample_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)

names(old_tissue) <- str_remove(names(old_tissue), "CRISPResso_on_")
names(old_cell) <- str_remove(names(old_cell), "CRISPResso_on_")
old_cell_order <- old_cell[names(old_cell) %in% sample_info$`Cell-w-W`]
split_name <- data.frame(cellid = sample_info$`Cell-w-W`, 
                         tissueid = sample_info$`Tissue-w-W(Mix)`,
                         mutation = sample_info$GRNAa)
split_name$mutation <- unlist(lapply(split_name$mutation, function(x){
  x <- unlist(strsplit(x, "[-]"))
  paste(x[-length(x)], collapse = "-")
}))

split_name_cell <- split_name[split_name$cellid %in% names(old_cell),]
split_name_cell <- na.omit(split_name_cell)

split_name_tissue <- split_name[split_name$tissueid %in% names(old_tissue),]
split_name_tissue <- na.omit(split_name_tissue)



split_name_cell$mutation <- unlist(lapply(split_name_cell$mutation, function(x){
  x <- unlist(str_split(x, "[-]"))
  x[2] <- stringr::str_to_title(x[2])
  paste(x, collapse = "-")
}))
split_name_cell$type <- unlist(lapply(split_name_cell$mutation, function(x){
  str_sub(x, str_length(x) - 1,  str_length(x) - 1)
}))


split_name_tissue$mutation <- unlist(lapply(split_name_tissue$mutation, function(x){
  x <- unlist(str_split(x, "[-]"))
  x[2] <- stringr::str_to_title(x[2])
  paste(x, collapse = "-")
}))
split_name_tissue$type <- unlist(lapply(split_name_tissue$mutation, function(x){
  str_sub(x, str_length(x) - 1,  str_length(x) - 1)
}))

###读取新测的4个sgRNA结果 然后合并到之前的结果里面去
load(file="~/data/project/ear_project/gene_therapy_ll/Previews/Result/old_rep123_edit.rda")
print(load("~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data_newer_newer.rda"))

rep1_new_otof <- otof_cell_sg12_15_edit[seq(1,8, 2)]
rep2_new_otof <- otof_cell_sg12_15_edit[seq(2,8, 2)]
for(name in names(rep1_new_otof)){
  name2 <- str_remove(name, "[-].*")
  old_rep1_edit[[name2]] <- rep1_new_otof[[name]]
}
for(name in names(rep2_new_otof)){
  name2 <- str_remove(name, "[-].*")
  old_rep2_edit[[name2]] <- rep2_new_otof[[name]]
}


otof_tissue_sg12_15_rep1 <- otof_tissue_sg12_15_resote[grep("_1",names(otof_tissue_sg12_15_resote))]
otof_tissue_sg12_15_rep2 <- otof_tissue_sg12_15_resote[grep("_2",names(otof_tissue_sg12_15_resote))]
otof_tissue_sg12_15_rep3 <- otof_tissue_sg12_15_resote[grep("_3",names(otof_tissue_sg12_15_resote))]
otof_tissue_sg12_15_rep4 <- otof_tissue_sg12_15_resote[grep("_4",names(otof_tissue_sg12_15_resote))]
names(otof_tissue_sg12_15_rep1) <- paste0(str_to_lower(str_remove(names(otof_tissue_sg12_15_rep1), "[_].*")), '-tissue')
names(otof_tissue_sg12_15_rep2) <- paste0(str_to_lower(str_remove(names(otof_tissue_sg12_15_rep2), "[_].*")), '-tissue')
names(otof_tissue_sg12_15_rep3) <- paste0(str_to_lower(str_remove(names(otof_tissue_sg12_15_rep3), "[_].*")), '-tissue')
names(otof_tissue_sg12_15_rep4) <- paste0(str_to_lower(str_remove(names(otof_tissue_sg12_15_rep4), "[_].*")), '-tissue')



for(name in names(otof_tissue_sg12_15_rep1)){
  old_rep1_edit[[name]] <- otof_tissue_sg12_15_rep1[[name]]
}
for(name in names(otof_tissue_sg12_15_rep2)){
  old_rep2_edit[[name]] <- otof_tissue_sg12_15_rep2[[name]]
}
for(name in names(otof_tissue_sg12_15_rep3)){
  old_rep3_edit[[name]] <- otof_tissue_sg12_15_rep3[[name]]
}
old_rep4 <- otof_tissue_sg12_15_rep4
total_old_edit_table <- list(Rep1 = old_rep1_edit, Rep2 = old_rep2_edit, Rep3 = old_rep3_edit, 
                             Rep4 = old_rep4)
for(name in names(total_old_edit_table)){
  names(total_old_edit_table[[name]]) <- str_remove(names(total_old_edit_table[[name]]), "CRISPResso_on_")
}



###接下来根据是Insert突变还是Delete突变来计算Del1与Ins1的编辑比例
split_name_tissue <- data.frame(rbind(split_name_tissue, 
                                      data.frame("cellid" = paste0("m", 12 : 15, "-tissue"),
                                                 "tissueid" = paste0("m", 12 : 15, "-tissue"),
                                                 "mutation" = 'w-Otof-1236dC', 
                                                 "type" = 'd')))
split_name_cell <- data.frame(rbind(split_name_cell, 
                                    data.frame("cellid" = paste0("m", 12 : 15),
                                               "tissueid" = paste0("m", 12 : 15),
                                               "mutation" = 'w-Otof-1236dC', 
                                               "type" = 'd')))
split_name_tissue$id <- str_remove(split_name_tissue$tissueid, "[-].*")
split_name_cell$id <- str_remove(split_name_cell$cellid, "[-].*")
split_name <- merge(split_name_tissue[,-1], split_name_cell[,-2], by="id")
split_name <- split_name[,c(-1, -3,-4)]
colnames(split_name) <- str_remove(colnames(split_name), "[.].*")


tissue_ins1_stat <- lapply(split_name$tissueid, function(id){
  tmp <- list()
  for(name in names(total_old_edit_table)){
    if(id %in% names(total_old_edit_table[[name]])){
      tmp[[name]] <- total_old_edit_table[[name]][[id]]
    }
  }
  res <- unlist(lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted != 0,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    ins_pct <- x[x$n_inserted == 1 & x$n_deleted == 0, ]
    return(sum(ins_pct$Pct))
  }))
  
  data.frame(Pct = mean(res), se = plotrix::std.error(res), id = id)
})
tissue_ins1_stat <- data.frame(do.call(rbind, tissue_ins1_stat))
tissue_ins1_stat <- tissue_ins1_stat[order(tissue_ins1_stat$Pct, decreasing = T),]
tissue_ins1_stat$id <- factor(tissue_ins1_stat$id, levels = tissue_ins1_stat$id)
ggplot(tissue_ins1_stat,aes(x = id,y = Pct)) + geom_bar(stat = "identity", fill = "#06a0c9") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), width = .2, color = "black") + 
  xlab("") + ylab("%Ins1 in Indel Reads") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())


#统计del1 Pct and Se
tissue_del1_stat <- lapply(split_name$tissueid, function(id){
  tmp <- list()
  for(name in names(total_old_edit_table)){
    if(id %in% names(total_old_edit_table[[name]])){
      tmp[[name]] <- total_old_edit_table[[name]][[id]]
    }
  }
  res <- unlist(lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted != 0,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    ins_pct <- x[x$n_inserted == 0 & x$n_deleted == 1, ]
    return(sum(ins_pct$Pct))
  }))
  
  data.frame(Pct = mean(res), se = plotrix::std.error(res), id = id)
})
tissue_del1_stat <- data.frame(do.call(rbind, tissue_del1_stat))
tissue_del1_stat <- tissue_del1_stat[order(tissue_del1_stat$Pct, decreasing = T),]
tissue_del1_stat$id <- factor(tissue_del1_stat$id, levels = tissue_del1_stat$id)
ggplot(tissue_del1_stat,aes(x = id,y = Pct)) + geom_bar(stat = "identity", fill = "#06a0c9") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), width = .2, color = "black") + 
  xlab("") + ylab("%del1 in Indel Reads") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())


tissue_indel1_stat <- lapply(split_name$tissueid, function(id){
  tmp <- list()
  for(name in names(total_old_edit_table)){
    if(id %in% names(total_old_edit_table[[name]])){
      tmp[[name]] <- total_old_edit_table[[name]][[id]]
    }
  }
  res <- unlist(lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted != 0,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    ins_pct <- x[xor(x$n_inserted == 1,  x$n_deleted == 1), ]
    return(sum(ins_pct$Pct))
  }))
  
  data.frame(Pct = mean(res), se = plotrix::std.error(res), id = id)
})
tissue_indel1_stat <- data.frame(do.call(rbind, tissue_indel1_stat))
tissue_indel1_stat <- tissue_indel1_stat[order(tissue_indel1_stat$Pct, decreasing = T),]
tissue_indel1_stat$id <- factor(tissue_indel1_stat$id, levels = tissue_indel1_stat$id)

cell_ins1_stat <- lapply(split_name$cellid, function(id){
  tmp <- list()
  for(name in names(total_old_edit_table)){
    if(id %in% names(total_old_edit_table[[name]])){
      tmp[[name]] <- total_old_edit_table[[name]][[id]]
    }
  }
  if(length(tmp) == 0){
    print(id)
  }
  res <- unlist(lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted != 0,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    ins_pct <- x[x$n_inserted == 1 & x$n_deleted == 0, ]
    return(sum(ins_pct$Pct))
  }))
  
  data.frame(Pct = mean(res), se = plotrix::std.error(res), id = id)
})
cell_ins1_stat <- data.frame(do.call(rbind, cell_ins1_stat))
cell_ins1_stat <- cell_ins1_stat[order(cell_ins1_stat$Pct, decreasing = T),]
cell_ins1_stat$id <- factor(cell_ins1_stat$id, levels = cell_ins1_stat$id)
ggplot(cell_ins1_stat,aes(x = id,y = Pct)) + geom_bar(stat = "identity", fill = "#06a0c9") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), width = .2, color = "black") + 
  xlab("") + ylab("%Ins1 in Indel Reads") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())


#统计del1 Pct and Se
cell_del1_stat <- lapply(split_name$cellid, function(id){
  tmp <- list()
  for(name in names(total_old_edit_table)){
    if(id %in% names(total_old_edit_table[[name]])){
      tmp[[name]] <- total_old_edit_table[[name]][[id]]
    }
  }
  res <- unlist(lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted != 0,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    ins_pct <- x[x$n_inserted == 0 & x$n_deleted == 1, ]
    return(sum(ins_pct$Pct))
  }))
  
  data.frame(Pct = mean(res), se = plotrix::std.error(res), id = id)
})
cell_del1_stat <- data.frame(do.call(rbind, cell_del1_stat))
cell_del1_stat <- cell_del1_stat[order(cell_del1_stat$Pct, decreasing = T),]
cell_del1_stat$id <- factor(cell_del1_stat$id, levels = cell_del1_stat$id)
ggplot(cell_del1_stat,aes(x = id,y = Pct)) + geom_bar(stat = "identity", fill = "#06a0c9") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), width = .2, color = "black") + 
  xlab("") + ylab("%del1 in Indel Reads") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())

cell_indel1_stat <- lapply(split_name$cellid, function(id){
  tmp <- list()
  for(name in names(total_old_edit_table)){
    if(id %in% names(total_old_edit_table[[name]])){
      tmp[[name]] <- total_old_edit_table[[name]][[id]]
    }
  }
  res <- unlist(lapply(tmp, function(x){
    x <- x[x$n_deleted + x$n_inserted != 0,]
    x$Pct <- x[,7] / sum(x[,7]) * 100
    ins_pct <- x[xor(x$n_inserted == 1 , x$n_deleted == 1), ]
    return(sum(ins_pct$Pct))
  }))
  
  data.frame(Pct = mean(res), se = plotrix::std.error(res), id = id)
})
cell_indel1_stat <- data.frame(do.call(rbind, cell_indel1_stat))
cell_indel1_stat <- cell_indel1_stat[order(cell_indel1_stat$Pct, decreasing = T),]
cell_indel1_stat$id <- factor(cell_indel1_stat$id, levels = cell_indel1_stat$id)

###之后统计在InDel1发生在切割位点的比例 但是分母是总Indel reads
tissue_indel1_stat <- lapply(split_name$tissueid, function(id){
  tmp <- list()
  for(name in names(total_old_edit_table)){
    if(id %in% names(total_old_edit_table[[name]])){
      tmp[[name]] <- total_old_edit_table[[name]][[id]]
    }
  }
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
tissue_indel1_stat <- data.frame(do.call(rbind, tissue_indel1_stat))
tissue_indel1_stat <- tissue_indel1_stat[order(tissue_indel1_stat$indelPct, decreasing = T),]
tissue_indel1_stat$id <- factor(tissue_indel1_stat$id, levels = tissue_indel1_stat$id)



cell_indel1_stat <- lapply(split_name$cellid, function(id){
  tmp <- list()
  for(name in names(total_old_edit_table)){
    if(id %in% names(total_old_edit_table[[name]])){
      tmp[[name]] <- total_old_edit_table[[name]][[id]]
    }
  }
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
cell_indel1_stat <- data.frame(do.call(rbind, cell_indel1_stat))
cell_indel1_stat <- cell_indel1_stat[order(cell_indel1_stat$indelPct, decreasing = T),]
cell_indel1_stat$id <- factor(cell_indel1_stat$id, levels = cell_indel1_stat$id)


total_indel_stat <- split_name
rownames(tissue_del1_stat) <- tissue_del1_stat$id
colnames(tissue_del1_stat) <- paste("tissue_del", colnames(tissue_del1_stat), sep = "_")
rownames(tissue_ins1_stat) <- tissue_ins1_stat$id
colnames(tissue_ins1_stat) <- paste("tissue_ins", colnames(tissue_ins1_stat), sep = "_")
rownames(tissue_indel1_stat) <- tissue_indel1_stat$id
colnames(tissue_indel1_stat) <- paste("tissue", colnames(tissue_indel1_stat), "in_cut", sep = "_")
tissue_del1_stat <- tissue_del1_stat[total_indel_stat$tissueid,]
tissue_ins1_stat <- tissue_ins1_stat[total_indel_stat$tissueid,]
tissue_indel1_stat <- tissue_indel1_stat[total_indel_stat$tissueid,]

rownames(cell_del1_stat) <- cell_del1_stat$id
colnames(cell_del1_stat) <- paste("cell_del", colnames(cell_del1_stat), sep = "_")
rownames(cell_ins1_stat) <- cell_ins1_stat$id
colnames(cell_ins1_stat) <- paste("cell_ins", colnames(cell_ins1_stat), sep = "_")
rownames(cell_indel1_stat) <- cell_indel1_stat$id
colnames(cell_indel1_stat) <- paste("cell", colnames(cell_indel1_stat), "in_cut", sep = "_")
cell_del1_stat <- cell_del1_stat[total_indel_stat$cellid,]
cell_ins1_stat <- cell_ins1_stat[total_indel_stat$cellid,]
cell_indel1_stat <- cell_indel1_stat[total_indel_stat$cellid,]

total_indel_stat <- cbind(total_indel_stat, tissue_ins1_stat[,c(1,2)], 
                          tissue_del1_stat[,c(1,2)],
                          tissue_indel1_stat[,-ncol(tissue_indel1_stat)],
                          cell_ins1_stat[,c(1,2)], 
                          cell_del1_stat[,c(1,2)],
                          cell_indel1_stat[,-ncol(cell_indel1_stat)])
total_indel_stat <- total_indel_stat[gtools::mixedorder(total_indel_stat$tissueid),]
openxlsx::write.xlsx(total_indel_stat, "~/Nutstore Files/Tobin/Previous/total_indel1_stat_result.xlsx", rowNames=F, colNames=T)



total_indel_stat  <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/total_indel1_stat_result.xlsx")

tmp <- total_indel_stat[,c("tissue_indelPct_in_cut", "tissue_indelSe_in_cut", "tissueid")]
tmp$type <- "cut"
tissue_indel1_plot <- tissue_indel1_stat
tissue_indel1_plot$type <- "indel1"
colnames(tmp) <- colnames(tissue_indel1_plot)
tissue_indel1_plot <- data.frame(rbind(tissue_indel1_plot, tmp))
tissue_indel1_plot$type <- factor(tissue_indel1_plot$type, levels = c("indel1", "cut"), labels = c("Total Indel1", 
                                                                                                   "Indel1 in cutsite"))
pdf("~/Nutstore Files/Tobin/Previous/Indel1_in_118_tissue.pdf", width = 20, height = 4)
ggplot(tissue_indel1_plot,aes(x = id,y = Pct, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), 
                color = "black", position = "dodge", width = 0.9) + 
  xlab("") + ylab("%Indel1 in Indel Reads") + 
  scale_fill_manual(values = c("#06a0c9" ,"#E1FFFF")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
dev.off()



tmp <- total_indel_stat[,c("cell_indelPct_in_cut", "cell_indelSe_in_cut", "cellid")]
tmp$type <- "cut"
cell_indel1_plot <- cell_indel1_stat
cell_indel1_plot$type <- "indel1"
colnames(tmp) <- colnames(cell_indel1_plot)
cell_indel1_plot <- data.frame(rbind(cell_indel1_plot, tmp))
cell_indel1_plot$type <- factor(cell_indel1_plot$type, levels = c("indel1", "cut"), labels = c("Total Indel1", 
                                                                                                   "Indel1 in cutsite"))
pdf("~/Nutstore Files/Tobin/Previous/Indel1_in_118_cell.pdf", width = 20, height = 4)
ggplot(cell_indel1_plot,aes(x = id,y = Pct, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = Pct - se, ymax = Pct + se), 
                color = "black", position = "dodge", width = 0.9) + 
  xlab("") + ylab("%Indel1 in Indel Reads") + 
  scale_fill_manual(values = c("#06a0c9" ,"#E1FFFF")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
dev.off()
