#total_pair_sg_cor 来自1.0
high_cor_sample <- total_pair_sg_cor[total_pair_sg_cor$cor > 0.8,]

sgRNA_sel <- sgRNA[sgRNA$id %in% high_cor_sample$id,]
setwd("~/data/project/ear_project/gene_therapy/find1NTSim/")
tmp <- read.table("human_sg_final_anno2_addid.txt", sep="\t")
tmp2 <- read.table("threshold20000/human_sg_final_anno2_addid.txt", sep="\t")
tmp3 <- read.table("human_sg_lastNTSim_final_addid2.txt", sep = "\t")
tmp3$V14 <- lapply(1 : nrow(tmp3), function(x){
  if(str_sub(tmp3[x,2], 2, 20) == str_sub(tmp3[x, 8], 1, 19)){
    return(paste0(paste0(tmp3[x,2], "-"), " ", paste0("-",tmp3[x,8])))
  }
  return(paste0(paste0(tmp3[x,8], "-"), " ", paste0("-",tmp3[x,2])))
})
tmp4 <-unlist(c(tmp[,14], tmp2[,14], tmp3[,14]))
tmp <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[1]
}))
tmp2 <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[2]
}))
tmp4 <- data.frame(sgRNA = str_remove(c(tmp, tmp2), "[-]"), 
                   Cmp = c(tmp, tmp2))
tmp4 <- tmp4[!duplicated(tmp4$sgRNA),]
sgCmp <- tmp4
sgCmp <- sgCmp[sgCmp$sgRNA %in% sgRNA_sel$sgRNA,]

high_cor_sample$cmp1 <- unlist(lapply(high_cor_sample$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[1]
  x <- sgRNA$sgRNA[sgRNA$id2 == x]
  sgCmp$Cmp[sgCmp$sgRNA == x]
}))

high_cor_sample$cmp2 <- unlist(lapply(high_cor_sample$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[2]
  x <- sgRNA$sgRNA[sgRNA$id2 == x]
  sgCmp$Cmp[sgCmp$sgRNA == x]
}))
ToNX::write_tb(high_cor_sample, file="~/data/project/ear_project/gene_therapy_ll/high_cor_sample_reDiff.txt")

total_edit_table <- list()
setwd("~/data/project/ear_project/gene_therapy_ll/high_cor/")
for(ids in high_cor_sample$pair){
  ids <- unlist(strsplit(ids, "[-]"))
  for(sample in c("Rep1", "Rep2", "Rep3")){
    setwd(sample)
    if(!sample %in% names(total_edit_table)){
      total_edit_table[[sample]] <- list()
    }
    for(id in ids){
      if(file.exists(id)){
        setwd(id)
        filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
        filename <- filename[which.max(str_length(filename))]
        edit_table <- read.table(filename, 
                                 sep="\t", header = T, comment.char = "")
        total_edit_table[[sample]][[id]] <- edit_table
        setwd("..")
      }
    }
    setwd("..")
  }
}


total_indel_table_aln <- lapply(total_edit_table, function(x){
  lapply(x, function(y){
    y$indel <- y[,5] - y[,4]
    forwt <- y[y$indel == 0,]
    wtSeq <- getWtSeq(forwt)
    y <- y[y[,4] != 0 | y[,5] != 0,]
    y <- y[order(y$indel),]
    align_edit_table(y, wtSeq)
  })
})

# save(total_indel_table_aln, file="total_indel_table_aln_new.rda")
setwd("~/data/project/ear_project/gene_therapy_ll/high_cor/")
load("~/data/project/ear_project/gene_therapy_ll/high_cor/total_indel_table_aln_new.rda")

high_cor_fixDiff <- read.table("~/data/project/ear_project/gene_therapy_ll/high_cor_sample_reDiff_output.txt", 
                               sep="\t")
high_cor_fixDiff <- high_cor_fixDiff[high_cor_fixDiff$V8 %in% c(3, 4),]
colnames(high_cor_fixDiff) <- c("pair", "cor",  "pval", "id",   "pos",  "cmp1", "cmp2", "diffPos")
pdf("indel_in_indel_output_version2.pdf", width = 24, height = 12)
for(pairids in high_cor_fixDiff$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  
  tmp_res <- lapply(ids, function(id){
    tmp <- list()
    for(sample in names(total_indel_table_aln)){
      if(id %in% names(total_indel_table_aln[[sample]])){
        sel <- total_indel_table_aln[[sample]][[id]]
        sel$sample <- sample
        sel$Pct <- sel[,7] / sum(sel[,7])
        tmp[[sample]] <- sel
      }
    }
    tmp <- data.frame(do.call(rbind, tmp))
    
    tmp_mean <- lapply(split(tmp, paste(tmp$Aligned_Sequence, tmp$Reference_Sequence)), function(x){
      x$Pct <- mean(x$Pct)
      x[1,]
    })
    tmp_mean <- data.frame(do.call(rbind, tmp_mean))
    tmp_stat <- lapply(split(tmp_mean$Pct, tmp_mean$indel), sum)
    tmp_stat <- data.frame(indel = names(tmp_stat), Pct = unlist(tmp_stat))
    return(list(mean_table = tmp_mean, stat = tmp_stat, orig_table = tmp))
  })
  names(tmp_res) <- ids
  pair_stat <- data.frame(rbind(tmp_res[[1]]$stat, tmp_res[[2]]$stat))
  pair_stat <- lapply(split(pair_stat$Pct, pair_stat$indel), mean)
  pair_stat <- data.frame(indel = names(pair_stat), Pct = unlist(pair_stat))
  pair_stat <- pair_stat[order(pair_stat$Pct, decreasing = T),]
  pair_stat <- pair_stat[as.numeric(pair_stat$indel) %in% c(-3 : 3),]
  pair_stat <- pair_stat[pair_stat$indel != "0",]
  ####统计折线图需要的百分比信息  需要算三组之间的均值
  mean_stat <- lapply(as.numeric(pair_stat[,1]), function(x){

    sel_pos1 <- tmp_res[[1]]$orig_table
    sel_pos1 <- sel_pos1[sel_pos1$indel == x,]
    sel_pos2 <- tmp_res[[2]]$orig_table
    sel_pos2 <- sel_pos2[sel_pos2$indel == x,]
    
    sel_stat1 <- lapply(split(sel_pos1, sel_pos1$sample), function(y){
      res <- getStat(y, x > 0)
      res$sample <- unique(y$sample)
      res
    })
    sel_stat2 <- lapply(split(sel_pos2, sel_pos2$sample), function(y){
      res <- getStat(y, x > 0)
      res$sample <- unique(y$sample)
      res
    })
    sel_stat1 <- data.frame(do.call(rbind, sel_stat1))
    sel_stat2 <- data.frame(do.call(rbind, sel_stat2))
    sel_stat1_mean <- lapply(split(sel_stat1, sel_stat1$pos), function(y){
      y$se <- plotrix::std.error(y[,2])
      y$se[is.na(y$se)] <- 0
      y$type <- paste0(ifelse(x > 0, "Ins", "Del"), x)
      y[,2] <- mean(y[,2])
      y[,3] <- mean(y[,3])
      y[1,]
    })
    sel_stat1_mean <- data.frame(do.call(rbind, sel_stat1_mean))
    
    sel_stat2_mean <- lapply(split(sel_stat2, sel_stat2$pos), function(y){
      y$se <- plotrix::std.error(y[,2])
      y$se[is.na(y$se)] <- 0
      y$type <- paste0(ifelse(x > 0, "Ins", "Del"), x)
      y[,2] <- mean(y[,2])
      y[,3] <- mean(y[,3])
      y[1,]
    })
    sel_stat2_mean <- data.frame(do.call(rbind, sel_stat2_mean))
    result <- list(sg1 = sel_stat1_mean, sg2 = sel_stat2_mean)
    names(result) <- names(tmp_res)
    result
  })
  names(mean_stat) <- unlist(lapply(as.numeric(pair_stat[,1]), function(y){
    paste0(ifelse(y > 0, "Ins", "Del"), y)
  }))
  mean_stat1 <- data.frame(do.call(rbind, lapply(mean_stat, function(y){
    y[[ids[1]]]
  })))
  mean_stat2 <- data.frame(do.call(rbind, lapply(mean_stat, function(y){
    y[[ids[2]]]
  })))
  
  colors <- hue_pal()(6)
  names(colors) <- c("Ins1", "Ins2", "Ins3", "Del-1", "Del-2", "Del-3")
  
  wt_seq1 <- total_edit_table$Rep3[[ids[1]]]
  wt_seq1 <- getWtSeq(wt_seq1[wt_seq1[,3] == "True",])
  maxPct <- max(c(mean_stat1$Pct + mean_stat1$se, mean_stat2$Pct+ mean_stat2$se))
  plot_data <- mean_stat1
  plot_data <- plot_data[plot_data$pos %in% as.character(c(7 : 26)),]
  plot_data$pos <- factor(plot_data$pos, 
                          levels = gtools::mixedsort(unique(plot_data$pos)))
  plot_data$type <- factor(plot_data$type, levels =names(mean_stat))
  se_data <- plot_data[plot_data$Pct != 0 & plot_data$se > 0.01,]
  p1 <- ggplot(plot_data, aes(x=pos, y=Pct, color=type, group = type)) + 
    geom_line() +  
    geom_vline(xintercept = 14.5, linetype = "dashed") + 
    geom_errorbar(data = se_data,
                  aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black")  + 
    ggplot2::facet_wrap(~type, nrow = nrow(pair_stat), ncol = 1, ) + 
    scale_x_discrete(labels = wt_seq1[7:26]) + 
    scale_y_continuous(limits = c(0, ceiling(maxPct))) + 
    theme_minimal() + theme(strip.background = element_blank(),
                            strip.text.x = element_blank())
  
  wt_seq2 <- total_edit_table$Rep3[[ids[2]]]
  wt_seq2 <- getWtSeq(wt_seq2[wt_seq2[,3] == "True",])
  
  plot_data <- mean_stat2
  plot_data <- plot_data[plot_data$pos %in% as.character(c(7 : 26)),]
  plot_data$pos <- factor(plot_data$pos, 
                          levels = gtools::mixedsort(unique(plot_data$pos)))
  plot_data$type <- factor(plot_data$type, levels =names(mean_stat))
  se_data <- plot_data[plot_data$Pct != 0 & plot_data$se > 0.01,]
  p2 <- ggplot(plot_data, aes(x=pos, y=Pct, color=type, group = type)) + 
    geom_line() +  
    geom_vline(xintercept = 14.5, linetype = "dashed") + 
    geom_errorbar(data = se_data,
                  aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black")  + 
    ggplot2::facet_wrap(~type, nrow = nrow(pair_stat), ncol = 1, ) + 
    scale_x_discrete(labels = wt_seq2[7:26]) + 
    scale_y_continuous(limits = c(0, ceiling(maxPct))) + 
    theme_minimal() + theme(  strip.background = element_blank(),
                              strip.text.x = element_blank())
  
  sg_seqs <- unlist(high_cor_fixDiff[high_cor_fixDiff$pair == pairids, 
                              c("cmp1", "cmp2")])
  
  sg_seqs <- paste(c(paste0(c(ids[1],sg_seqs[1]), collapse = ":"), 
                     paste0(c(ids[2],sg_seqs[2]), collapse = ":")), 
                     collapse = "\n")
  merge_plot <- ggarrange(p1, p2, nrow = 1, legend = "top", common.legend = T)
  merge_plot <- annotate_figure(merge_plot, top = text_grob(sg_seqs, 
                                        color = "black", face = "bold", size = 14))
  ####画柱状图
  
  p_l1 <- lapply(names(total_data[[ids[1]]]), function(n){
    region <- c(-19: 19)
    sel_sg <- total_data[[ids[1]]][[n]]
    long_del <- sum(sel_sg[sel_sg[,1] < min(region), 2])
    long_ins <- sum(sel_sg[sel_sg[,1] > max(region), 2])
    long_res <- data.frame(indel_size = c("<=-20", ">=20"), fq = c(long_del, long_ins))
    tmp_region <- region[!region %in% sel_sg[,1]]
    sel_sg <- sel_sg[sel_sg[,1] %in% region,]
    if(length(tmp_region) != 0){
      sel_sg <- data.frame(rbind(sel_sg, data.frame(indel_size = tmp_region, fq = 0)))
      sel_sg <- sel_sg[order(sel_sg[,1]),]
    }
    sel_sg <- data.frame(rbind(sel_sg, long_res))
    sel_sg[sel_sg[,1] == 0, 2] <- 0 
    sel_sg$sample <- n
    sel_sg$pct <- sel_sg[,2] / sum(sel_sg[,2]) * 100
    sel_sg
  })
  p_l1 <- data.frame(do.call(rbind, p_l1))
  se1 <- lapply(split(p_l1$pct, p_l1$indel_size), function(x){
    plotrix::std.error(unlist(x))
  })
  se1 <- data.frame(indel_size = names(se1), se = unlist(se1))
  se1$se[is.na(se1$se)] <- 0
  p_l1_mean <- lapply(split(p_l1$pct, p_l1$indel_size), mean)
  p_l1_mean <- data.frame(indel_size = names(p_l1_mean), Pct = unlist(p_l1_mean))
  p_l1_mean <- merge(p_l1_mean, se1, by="indel_size")
  p_l1_mean$indel_size <- factor(p_l1_mean$indel_size, 
                                 levels = c("<=-20", -19 : 19, ">=20"))
 
  p_l2 <- lapply(names(total_data[[ids[2]]]), function(n){
    region <- c(-19: 19)
    sel_sg <- total_data[[ids[2]]][[n]]
    long_del <- sum(sel_sg[sel_sg[,1] < min(region), 2])
    long_ins <- sum(sel_sg[sel_sg[,1] > max(region), 2])
    long_res <- data.frame(indel_size = c("<=-20", ">=20"), fq = c(long_del, long_ins))
    tmp_region <- region[!region %in% sel_sg[,1]]
    sel_sg <- sel_sg[sel_sg[,1] %in% region,]
    if(length(tmp_region) != 0){
      sel_sg <- data.frame(rbind(sel_sg, data.frame(indel_size = tmp_region, fq = 0)))
      sel_sg <- sel_sg[order(sel_sg[,1]),]
    }
    sel_sg <- data.frame(rbind(sel_sg, long_res))
    sel_sg[sel_sg[,1] == 0, 2] <- 0 
    sel_sg$sample <- n
    sel_sg$pct <- sel_sg[,2] / sum(sel_sg[,2]) * 100
    sel_sg
  })
  p_l2 <- data.frame(do.call(rbind, p_l2))
  se2 <- lapply(split(p_l2$pct, p_l2$indel_size), function(x){
    plotrix::std.error(unlist(x))
  })
  se2 <- data.frame(indel_size = names(se2), se = unlist(se2))
  se2$se[is.na(se2$se)] <- 0
  p_l2_mean <- lapply(split(p_l2$pct, p_l2$indel_size), mean)
  p_l2_mean <- data.frame(indel_size = names(p_l2_mean), Pct = unlist(p_l2_mean))
  p_l2_mean <- merge(p_l2_mean, se2, by="indel_size")
  p_l2_mean$indel_size <- factor(p_l2_mean$indel_size, 
                                 levels = c("<=-20", -19 : 19, ">=20"))

  
  max_ratio <- max(c(p_l2_mean$Pct + p_l2_mean$se, p_l1_mean$Pct + p_l1_mean$se))
  x_label <- c("<=-20", -19 : 19, ">=20")
  x_label[seq(2, 40, 2)] <- ""
  p1 <- ggplot(p_l1_mean, aes(x = indel_size, y = Pct)) + 
    geom_bar(stat="identity", fill = "#06a0c9") + 
    geom_errorbar(aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black")  + 
    ylab(paste0("% in indel reads of ", ids[1])) + xlab("") + 
    scale_x_discrete(labels = x_label) + 
    scale_y_continuous(expand = c(0,0), limits = c(0, max_ratio)) + 
    theme_classic2() + 
    theme(plot.margin =  unit(c(0,0,0,0.2), "cm"))
  
  p2 <- ggplot(p_l2_mean, aes(x = indel_size, y = Pct)) + 
    geom_bar(stat="identity", fill = "#06a0c9") + 
    geom_errorbar(aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black")  + 
    scale_y_reverse(expand = c(0,0),  limits = c(max_ratio, 0)) + 
    scale_x_discrete(position = "top") + 
    ylab(paste0( " % in indel reads of ", ids[2])) + xlab("") + 
    theme_classic2() + 
    theme(axis.text.x = element_blank(), 
          plot.margin =  unit(c(-1,0,0,0.2), "cm"))
  merge_plot2 <- ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T)
  merge_plot2 <- annotate_figure(merge_plot2, top = text_grob(paste0("Cor: ", round(high_cor_fixDiff$cor[high_cor_fixDiff$pair == pairids], 3)), 
                                                            color = "black", face = "bold", size = 14))
  rm(p1, p2, p_l1, p_l2)
  print(aplot::insert_right(merge_plot2, merge_plot, width = 1))
  
  
  # bar_stat <- lapply(as.numeric(pair_stat[,1]), function(x){
  #   sel_pos1 <- tmp_res[[1]]$orig_table
  #   sel_pos1 <- sel_pos1[sel_pos1$indel == x,]
  #   sel_pos1 <- lapply(split(sel_pos1$Pct, sel_pos1$sample), sum)
  #   sel_pos2 <- tmp_res[[2]]$orig_table
  #   sel_pos2 <- sel_pos2[sel_pos2$indel == x,]
  #   sel_pos2 <- lapply(split(sel_pos2$Pct, sel_pos2$sample), sum)
  #   se1 <- plotrix::std.error(unlist(sel_pos1) * 100)
  #   se2 <- plotrix::std.error(unlist(sel_pos2) * 100)
  #   
  #   sel_pos1 <- data.frame(id = ids[1], 
  #                          Pct = mean(unlist(sel_pos1)) * 100, 
  #                          se = se1,
  #                          type = paste0(ifelse(x > 0, "Ins", "Del"), x))
  #   sel_pos2 <- data.frame(id = ids[2], 
  #                          Pct = mean(unlist(sel_pos2)) * 100, 
  #                          se = se2,
  #                          type = paste0(ifelse(x > 0, "Ins", "Del"), x))
  #   data.frame(rbind(sel_pos1, sel_pos2))
  # })
  # bar_stat <- data.frame(do.call(rbind, bar_stat))
  # bar_stat$type <- factor(bar_stat$type, levels = names(mean_stat))
  # barplot1 <- ggplot(bar_stat[bar_stat$id == ids[1],], aes(x = type, y = Pct)) + 
  #   geom_bar(aes(fill = type), stat = "identity") + 
  #   geom_errorbar(aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black") + 
  #   theme_classic2() + xlab("")+
  #   theme(plot.margin =  unit(c(0,0,0,0.2), "cm"))
  # barplot2 <- ggplot(bar_stat[bar_stat$id == ids[2],], aes(x = type, y = Pct)) + 
  #   geom_bar(aes(fill = type), stat = "identity") + 
  #   geom_errorbar(aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black") + 
  #   xlab("") +
  #   scale_y_reverse() + scale_x_discrete(position = "top") + theme_classic2() + 
  #   theme(axis.text.x = element_blank(), 
  #         plot.margin =  unit(c(-1,0,0,0.2), "cm"))
  # merge_plot2 <- ggpubr::ggarrange(barplot1, barplot2, ncol = 1, legend = "none", heights = c(0.55, 0.45))
  

}
dev.off()


ins1_table <- lapply(total_indel_table_aln, function(x){
  lapply(x, function(y){
    y$Pct <- y[,7] / sum(y[,7])
    y[y$indel == 1,]
  })
})
ins1_stat <- lapply(ins1_table, function(x){
  lapply(x, function(y){
    getStat(y, T)
  })
})
del1_stat <- lapply(total_indel_table_aln, function(x){
  lapply(x, function(y){
    y$Pct <- y[,7] / sum(y[,7])
    y <- y[y$indel == -1,]
    if(nrow(y) != 0){
      getStat(y, F)
    }
    else{
      data.frame(pos = -1, Pct = 0, GroupPct = 0)
    }
  })
})
del2_stat <- lapply(total_indel_table_aln, function(x){
  lapply(x, function(y){
    y$Pct <- y[,7] / sum(y[,7])
    y <- y[y$indel == -2,]
    if(nrow(y) != 0){
      getStat(y, F)
    }
    else{
      data.frame(pos = -1, Pct = 0, GroupPct = 0)
    }
  })
})

ins1_stat_mean <- list()
del1_stat_mean <- list()
del2_stat_mean <- list()
for(rep in names(ins1_stat)){
  tmp_data <- ins1_stat[[rep]]
  for(sg in names(tmp_data)){
    tmp <- tmp_data[[sg]]
    tmp$Rep <- rep
    tmp$sg <- sg
    tmp_data[[sg]] <- tmp
    ins1_stat_mean[[length(ins1_stat_mean) + 1]] <- tmp
  }
  ins1_stat[[rep]] <- tmp_data
}

for(rep in names(del1_stat)){
  tmp_data <- del1_stat[[rep]]
  for(sg in names(tmp_data)){
    tmp <- tmp_data[[sg]]
    tmp$Rep <- rep
    tmp$sg <- sg
    tmp_data[[sg]] <- tmp
    del1_stat_mean[[length(del1_stat_mean) + 1]] <- tmp
  }
  del1_stat[[rep]] <- tmp_data
}

for(rep in names(del2_stat)){
  tmp_data <- del2_stat[[rep]]
  for(sg in names(tmp_data)){
    tmp <- tmp_data[[sg]]
    tmp$Rep <- rep
    tmp$sg <- sg
    tmp_data[[sg]] <- tmp
    del2_stat_mean[[length(del2_stat_mean) + 1]] <- tmp
  }
  del2_stat[[rep]] <- tmp_data
}

ins1_stat_mean <- data.frame(do.call(rbind, ins1_stat_mean))
ins1_stat_mean <- ins1_stat_mean[ins1_stat_mean$pos != -1,]
del1_stat_mean <- data.frame(do.call(rbind, del1_stat_mean))
del1_stat_mean <- del1_stat_mean[del1_stat_mean$pos != -1,]

del2_stat_mean <- data.frame(do.call(rbind, del2_stat_mean))
del2_stat_mean <- del2_stat_mean[del2_stat_mean$pos != -1,]

ins1_stat_mean <- lapply(split(ins1_stat_mean, ins1_stat_mean$sg),function(x){
  x <- lapply(split(x, x$pos), function(y){
    y$se <- plotrix::std.error(y[,2])
    y$type <- "Ins1"
    y[,2] <- mean(y[,2])
    y[,3] <- mean(y[,3])
    y[1,]
  })
  data.frame(do.call(rbind, x))
})

del1_stat_mean <- lapply(split(del1_stat_mean, del1_stat_mean$sg),function(x){
  x <- lapply(split(x, x$pos), function(y){
    y$se <- plotrix::std.error(y[,2])
    y$type <- "Del1"
    y[,2] <- mean(y[,2])
    y[,3] <- mean(y[,3])
    y[1,]
  })
  data.frame(do.call(rbind, x))
})

del2_stat_mean <- lapply(split(del2_stat_mean, del2_stat_mean$sg),function(x){
  x <- lapply(split(x, x$pos), function(y){
    y$se <- plotrix::std.error(y[,2])
    y$type <- "Del2"
    y[,2] <- mean(y[,2])
    y[,3] <- mean(y[,3])
    y[1,]
  })
  data.frame(do.call(rbind, x))
})


high_cor_sample


id <- "Sg_1_2"
plot_data <- data.frame(rbind(ins1_stat_mean[[id]], del1_stat_mean[[id]], 
                              del2_stat_mean[[id]]))
plot_data <- plot_data[plot_data$pos %in% as.character(c(7 : 26)),]
plot_data$pos <- factor(plot_data$pos, 
                        levels = gtools::mixedsort(unique(plot_data$pos)))
plot_data$type <- factor(plot_data$type, levels = c("Ins1", "Del1", "Del2"))
se_data <- plot_data[plot_data$Pct != 0 & plot_data$se > 0.01,]
ggplot(plot_data, aes(x=pos, y=Pct, color=type, group = type)) + 
  geom_line() +  
  geom_errorbar(data = se_data,
                aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black")  + 
  ggplot2::facet_wrap(~type, nrow = 3, ncol = 1) + 
  theme_minimal()









###output result to file
setwd("~/data/project/ear_project/gene_therapy_ll/edit_seqlist/")
lapply(names(total_indel_table_aln), function(xn){
  x <- total_indel_table_aln[[xn]]
  lapply(names(x), function(yn){
    print(paste0(xn, "_", yn))
    y <- x[[yn]]
    y$Pct <- y[,7] / sum(y[,7])
    result <- list()
    for(indel in unique(y$indel)){
      tmp <- y[y$indel == indel,]
      stat <- getStat(tmp, indel > 0)
      seq <- unlist(lapply(tmp[,10], function(x){
        paste(unlist(strsplit(x, "*")), collapse = "  ")
      }))
      tmp_res <- data.frame(indel = indel, name = "seq", seq = seq, Pct = tmp$Pct * 100)
      tmp_res[nrow(tmp_res) + 1, ] <- c(indel, "score_table", 
                                        paste(stat[,2], collapse = "  "), "None")
      tmp_res[nrow(tmp_res) + 1, ] <- c(indel, "score_group", 
                                        paste(stat[,3], collapse = "  "), "None")
      result[[as.character(indel)]] <- tmp_res
    }
    indel_order <- gtools::mixedsort(as.character(unique(y$indel)))
    result <- result[indel_order]
    result <- data.frame(do.call(rbind, result))
    ToNX::write_tb(result, file=paste0(xn, "/", yn, "_Pure_Editing_product_seqlist.tsv"))
  })
})


for(rep in c("Rep1", "Rep2", "Rep3")){
  setwd(rep)
  for(pair in high_cor_sample$pair){
    files <- unlist(strsplit(pair, "[-]"))
    exist_files <- unlist(lapply(files, function(file){
      unlist(list.files(pattern = file))
    }))
    if(length(exist_files) == 0)
      next
    tmp_res <- list()
    for(file in exist_files){
      tmp <- read.table(file)
      tmp_res[[str_remove(file, "_Pure_Editing_product_seqlist.tsv")]] <- tmp

    }
    openxlsx::write.xlsx(tmp_res, file=paste0(pair, ".xlsx"), 
                         rowNames=F, colNames=F)
  }
  setwd("..")
}



















