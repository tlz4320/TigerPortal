setwd(
  "~/data/project/ear_project/gene_therapy_ll/Second/Batch1/"
)
print(load("indel_stat.rda"))
print(load("processed_cor_result.rda"))
print(load("../second_sg.rda"))

high_cor_sample <- second_batch1_pair_cor[second_batch1_pair_cor$cor > 0.8,]


high_cor_sample$cmp1 <- unlist(lapply(high_cor_sample$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[1]
  second_sg$Cmp[second_sg$ID2 == x]
}))

high_cor_sample$cmp2 <- unlist(lapply(high_cor_sample$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[2]
  second_sg$Cmp[second_sg$ID2 == x]
}))

second_batch1_edit_rev <- lapply(names(second_batch1_indel[[1]]), function(x){
  lapply(second_batch1_indel, function(y){
    y[[x]]
  })
})
names(second_batch1_edit_rev) <- names(second_batch1_indel[[1]])

second_batch1_edit_aln <- lapply(second_batch1_edit_table, function(x){
  lapply(x, function(y){
    y$indel <- y[,5] - y[,4]
    forwt <- y[y$indel == 0,]
    wtSeq <- getWtSeq(forwt)
    y <- y[y[,4] != 0 | y[,5] != 0,]
    y <- y[order(y$indel),]
    align_edit_table(y, wtSeq)
  })
})


pdf("Result/indel_in_indel_second_batch1_result.pdf", width = 20, height = 12)
for(pairids in high_cor_sample$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  #put insert one at front
  sg_seqs <- unlist(high_cor_sample[high_cor_sample$pair == pairids, 
                                     c("cmp1", "cmp2")])
  cmp <- sg_seqs[1]
  cmp_pair <- sg_seqs[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  #在第一位有一个插入情况下，说明后面可能是last或者非last的情况
  if(this_pos == 1 | pair_pos == 1){
    
    #说明是last位有插入 且在这个second_sg插入
    if(pair_pos == length(cmp)){
      ids <- ids
    }
    #说明是last位有插入 且不在这个second_sg插入
    if(this_pos == length(cmp)){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
    #说明不是last位有插入，那就看谁在第一位有-那就是后面有插入
    if(this_pos == 1){
      ids <- ids
    }else{
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
  } else {
    #说明是最后一个碱基出现了shift，那就看最早的-是出现在谁那边谁就有delete
    if(this_pos < pair_pos){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    } else{
      ids <- ids
    }
  }
  
  
  
  
  tmp_res <- lapply(ids, function(id){
    tmp <- list()
    for(sample in names(second_batch1_edit_aln)){
      if(id %in% names(second_batch1_edit_aln[[sample]])){
        sel <- second_batch1_edit_aln[[sample]][[id]]
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
  rownames(pair_stat) <- pair_stat$indel
  pair_stat <- pair_stat[c("-1", "-2", "-3", "1", "2", "3"),]
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
      y$type <- paste0(ifelse(x > 0, "Ins", "Del"), abs(x))
      y[,2] <- mean(y[,2])
      y[,3] <- mean(y[,3])
      y[1,]
    })
    sel_stat1_mean <- data.frame(do.call(rbind, sel_stat1_mean))
    
    sel_stat2_mean <- lapply(split(sel_stat2, sel_stat2$pos), function(y){
      y$se <- plotrix::std.error(y[,2])
      y$se[is.na(y$se)] <- 0
      y$type <- paste0(ifelse(x > 0, "Ins", "Del"), abs(x))
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
    paste0(ifelse(y > 0, "Ins", "Del"), abs(y))
  }))
  mean_stat1 <- data.frame(do.call(rbind, lapply(mean_stat, function(y){
    y[[ids[1]]]
  })))
  mean_stat2 <- data.frame(do.call(rbind, lapply(mean_stat, function(y){
    y[[ids[2]]]
  })))
  mean_stat1$sg <- "sg1"
  mean_stat2$sg <- "sg2"
  merged_mean_stat <- data.frame(rbind(mean_stat1, mean_stat2))
  colors <- hue_pal()(6)
  names(colors) <- c("Ins1", "Ins2", "Ins3", "Del1", "Del2", "Del3")
  wt_seq1 <- second_batch1_edit_table[[1]][[ids[1]]]
  wt_seq1 <- getWtSeq(wt_seq1[wt_seq1[,3] == "True",])
  wt_seq2 <- second_batch1_edit_table[[1]][[ids[2]]]
  wt_seq2 <- getWtSeq(wt_seq2[wt_seq2[,3] == "True",])
  plot_list <- lapply(names(mean_stat), function(x){
    plot_data <- merged_mean_stat[merged_mean_stat$type == x,]
    maxPct <- ceiling(max(c(plot_data$Pct + plot_data$se)))
    plot_data <- plot_data[plot_data$pos %in% as.character(c(7 : 26)),]
    plot_data$pos <- factor(plot_data$pos, 
                            levels = gtools::mixedsort(unique(plot_data$pos)))
    se_data <- plot_data[plot_data$Pct != 0 & plot_data$se > 0.01,]
    
    color_seq <- mkColorfulSeq(wt_seq1[7:26], wt_seq2[7:26], sg_seqs[1], sg_seqs[2])
    
    ggplot(plot_data, aes(x=pos, y=Pct,group = type)) + 
      geom_vline(xintercept = 14.5, linetype = "dashed") + 
      geom_line(color = colors[x]) +  
      geom_errorbar(data = se_data,
                    aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black")  + 
      ggplot2::facet_wrap(~sg, nrow = 1, ncol = 2, scales = 'free_x') + 
      ggh4x::facetted_pos_scales(
        x = list(
          sg == "sg1" ~ scale_x_discrete(labels = color_seq[[1]]),
          sg == "sg2" ~ scale_x_discrete(labels = color_seq[[2]])
        )
      ) + 
      scale_y_continuous(limits = c(0, ceiling(maxPct))) + 
      ylab(x) + 
      theme_minimal() + theme(strip.background = element_blank(),
                              strip.text.x = element_blank(), axis.title.x = element_blank(), 
                              axis.text.x = ggtext::element_markdown())
    
  })
  
  
  sg_seqs <- paste(c(paste0(c(ids[1],sg_seqs[1]), collapse = ":"), 
                     paste0(c(ids[2],sg_seqs[2]), collapse = ":")), 
                   collapse = "\n")
  legend_plot <- get_legend(ggplot(data.frame(name = factor(names(colors), levels = names(colors)))) + 
                              geom_point(aes(x=name, y = 1, color = name)) + 
                              scale_color_manual(values = colors) + guides( colour = guide_legend("", nrow = 1, ncol = 6)))
  merge_plot <- ggarrange(plotlist = plot_list, ncol = 1, legend = "top", common.legend = T, legend.grob = legend_plot)
  merge_plot <- annotate_figure(merge_plot, top = text_grob(sg_seqs, 
                                                            color = "black", face = "bold", size = 14))
  
  ####画柱状图
  p_l1 <- lapply(names(second_batch1_edit_rev[[ids[1]]]), function(n){
    region <- c(-19: 19)
    sel_sg <- second_batch1_edit_rev[[ids[1]]][[n]]
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
  
  p_l2 <- lapply(names(second_batch1_edit_rev[[ids[2]]]), function(n){
    region <- c(-19: 19)
    sel_sg <- second_batch1_edit_rev[[ids[2]]][[n]]
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
  merge_plot2 <- annotate_figure(merge_plot2, top = text_grob(paste0("Cor: ", round(high_cor_sample$cor[high_cor_sample$pair == pairids], 3)), 
                                                              color = "black", face = "bold", size = 14))
  rm(p1, p2, p_l1, p_l2)
  print(aplot::insert_right(merge_plot2, merge_plot, width = 1))
  
}
dev.off()


####输出所有补位在前面且Bulge在34位置的结果

bulge_append_first <- second_batch1_pair_cor


bulge_append_first$cmp1 <- unlist(lapply(bulge_append_first$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[1]
  second_sg$Cmp[second_sg$ID2 == x]
}))

bulge_append_first$cmp2 <- unlist(lapply(bulge_append_first$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[2]
  second_sg$Cmp[second_sg$ID2 == x]
}))

bulge_append_first <- bulge_append_first[unlist(lapply(1 : nrow(bulge_append_first), function(i){
  cmps <- unlist(bulge_append_first[i, c("cmp1", "cmp2")])
  cmp <- cmps[1]
  cmp_pair <- cmps[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  return(which(cmp == '-') == 1 | which(cmp_pair=='-') == 1)
})),]
pdf("~/data/project/ear_project/gene_therapy_ll/Second/Batch1/Result/indel_in_indel_second_batch1_append_first_result.pdf", width = 20, height = 12)
for(pairids in bulge_append_first$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  #put insert one at front
  sg_seqs <- unlist(bulge_append_first[bulge_append_first$pair == pairids, 
                                    c("cmp1", "cmp2")])
  cmp <- sg_seqs[1]
  cmp_pair <- sg_seqs[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  #在第一位有一个插入情况下，说明后面可能是last或者非last的情况
  if(this_pos == 1 | pair_pos == 1){
    
    #说明是last位有插入 且在这个second_sg插入
    if(pair_pos == length(cmp)){
      ids <- ids
    }
    #说明是last位有插入 且不在这个second_sg插入
    if(this_pos == length(cmp)){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
    #说明不是last位有插入，那就看谁在第一位有-那就是后面有插入
    if(this_pos == 1){
      ids <- ids
    }else{
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    }
    
  } else {
    #说明是最后一个碱基出现了shift，那就看最早的-是出现在谁那边谁就有delete
    if(this_pos < pair_pos){
      ids <- rev(ids)
      sg_seqs <- rev(sg_seqs)
    } else{
      ids <- ids
    }
  }
  
  
  
  
  tmp_res <- lapply(ids, function(id){
    tmp <- list()
    for(sample in names(second_batch1_edit_aln)){
      if(id %in% names(second_batch1_edit_aln[[sample]])){
        sel <- second_batch1_edit_aln[[sample]][[id]]
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
  rownames(pair_stat) <- pair_stat$indel
  pair_stat <- pair_stat[c("-1", "-2", "-3", "1", "2", "3"),]
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
      y$type <- paste0(ifelse(x > 0, "Ins", "Del"), abs(x))
      y[,2] <- mean(y[,2])
      y[,3] <- mean(y[,3])
      y[1,]
    })
    sel_stat1_mean <- data.frame(do.call(rbind, sel_stat1_mean))
    
    sel_stat2_mean <- lapply(split(sel_stat2, sel_stat2$pos), function(y){
      y$se <- plotrix::std.error(y[,2])
      y$se[is.na(y$se)] <- 0
      y$type <- paste0(ifelse(x > 0, "Ins", "Del"), abs(x))
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
    paste0(ifelse(y > 0, "Ins", "Del"), abs(y))
  }))
  mean_stat1 <- data.frame(do.call(rbind, lapply(mean_stat, function(y){
    y[[ids[1]]]
  })))
  mean_stat2 <- data.frame(do.call(rbind, lapply(mean_stat, function(y){
    y[[ids[2]]]
  })))
  mean_stat1$sg <- "sg1"
  mean_stat2$sg <- "sg2"
  merged_mean_stat <- data.frame(rbind(mean_stat1, mean_stat2))
  colors <- hue_pal()(6)
  names(colors) <- c("Ins1", "Ins2", "Ins3", "Del1", "Del2", "Del3")
  wt_seq1 <- second_batch1_edit_table[[1]][[ids[1]]]
  wt_seq1 <- getWtSeq(wt_seq1[wt_seq1[,3] == "True",])
  wt_seq2 <- second_batch1_edit_table[[1]][[ids[2]]]
  wt_seq2 <- getWtSeq(wt_seq2[wt_seq2[,3] == "True",])
  plot_list <- lapply(names(mean_stat), function(x){
    plot_data <- merged_mean_stat[merged_mean_stat$type == x,]
    maxPct <- ceiling(max(c(plot_data$Pct + plot_data$se)))
    plot_data <- plot_data[plot_data$pos %in% as.character(c(7 : 26)),]
    plot_data$pos <- factor(plot_data$pos, 
                            levels = gtools::mixedsort(unique(plot_data$pos)))
    se_data <- plot_data[plot_data$Pct != 0 & plot_data$se > 0.01,]
    
    color_seq <- mkColorfulSeq(wt_seq1[7:26], wt_seq2[7:26], sg_seqs[1], sg_seqs[2])
    
    ggplot(plot_data, aes(x=pos, y=Pct,group = type)) + 
      geom_vline(xintercept = 14.5, linetype = "dashed") + 
      geom_line(color = colors[x]) +  
      geom_errorbar(data = se_data,
                    aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black")  + 
      ggplot2::facet_wrap(~sg, nrow = 1, ncol = 2, scales = 'free_x') + 
      ggh4x::facetted_pos_scales(
        x = list(
          sg == "sg1" ~ scale_x_discrete(labels = color_seq[[1]]),
          sg == "sg2" ~ scale_x_discrete(labels = color_seq[[2]])
        )
      ) + 
      scale_y_continuous(limits = c(0, ceiling(maxPct))) + 
      ylab(x) + 
      theme_minimal() + theme(strip.background = element_blank(),
                              strip.text.x = element_blank(), axis.title.x = element_blank(), 
                              axis.text.x = ggtext::element_markdown())
    
  })
  
  
  sg_seqs <- paste(c(paste0(c(ids[1],sg_seqs[1]), collapse = ":"), 
                     paste0(c(ids[2],sg_seqs[2]), collapse = ":")), 
                   collapse = "\n")
  legend_plot <- get_legend(ggplot(data.frame(name = factor(names(colors), levels = names(colors)))) + 
                              geom_point(aes(x=name, y = 1, color = name)) + 
                              scale_color_manual(values = colors) + guides( colour = guide_legend("", nrow = 1, ncol = 6)))
  merge_plot <- ggarrange(plotlist = plot_list, ncol = 1, legend = "top", common.legend = T, legend.grob = legend_plot)
  merge_plot <- annotate_figure(merge_plot, top = text_grob(sg_seqs, 
                                                            color = "black", face = "bold", size = 14))
  
  ####画柱状图
  p_l1 <- lapply(names(second_batch1_edit_rev[[ids[1]]]), function(n){
    region <- c(-19: 19)
    sel_sg <- second_batch1_edit_rev[[ids[1]]][[n]]
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
  
  p_l2 <- lapply(names(second_batch1_edit_rev[[ids[2]]]), function(n){
    region <- c(-19: 19)
    sel_sg <- second_batch1_edit_rev[[ids[2]]][[n]]
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
  merge_plot2 <- annotate_figure(merge_plot2, top = text_grob(paste0("Cor: ", round(bulge_append_first$cor[bulge_append_first$pair == pairids], 3)), 
                                                              color = "black", face = "bold", size = 14))
  rm(p1, p2, p_l1, p_l2)
  print(aplot::insert_right(merge_plot2, merge_plot, width = 1))
  
}
dev.off()

names(second_batch1_edit_aln) <- c("Rep1", "Rep2", "Rep3")
lapply(names(second_batch1_edit_aln), function(xn){
  x <- second_batch1_edit_aln[[xn]]
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





