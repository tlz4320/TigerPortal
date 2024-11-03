#total_pair_sg_cor 来自1.0 total_data也是来自1.0
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
# ToNX::write_tb(high_cor_sample, file="~/data/project/ear_project/gene_therapy_ll/high_cor_sample_reDiff.txt")

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
pdf("~/data/project/ear_project/gene_therapy_ll/high_cor/indel_in_indel_output_version4.pdf", width = 20, height = 12)
for(pairids in high_cor_fixDiff$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  #put insert one at front
  sg_seqs <- unlist(high_cor_fixDiff[high_cor_fixDiff$pair == pairids, 
                                     c("cmp1", "cmp2")])
  cmp <- sg_seqs[1]
  cmp_pair <- sg_seqs[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  #在第一位有一个插入情况下，说明后面可能是last或者非last的情况
  if(this_pos == 1 | pair_pos == 1){
    
    #说明是last位有插入 且在这个sgRNA插入
    if(pair_pos == length(cmp)){
      ids <- ids
    }
    #说明是last位有插入 且不在这个sgRNA插入
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
  mean_stat1$sg <- 'sg1'
  mean_stat2$sg <- 'sg2'
  merged_mean_stat <- data.frame(rbind(mean_stat1, mean_stat2))
  colors <- hue_pal()(6)
  names(colors) <- c("Ins1", "Ins2", "Ins3", "Del1", "Del2", "Del3")
  wt_seq1 <- total_edit_table$Rep3[[ids[1]]]
  wt_seq1 <- getWtSeq(wt_seq1[wt_seq1[,3] == "True",])
  wt_seq2 <- total_edit_table$Rep3[[ids[2]]]
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
  
}
dev.off()




####输出所有补位在前面且Bulge在34位置的结果

bulge_append_first_res <- total_pair_sg_cor

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
sgCmp <- sgCmp[sgCmp$sgRNA %in% sgRNA$sgRNA,]

bulge_append_first_res$cmp1 <- unlist(lapply(bulge_append_first_res$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[1]
  x <- sgRNA$sgRNA[sgRNA$id2 == x]
  sgCmp$Cmp[sgCmp$sgRNA == x]
}))

bulge_append_first_res$cmp2 <- unlist(lapply(bulge_append_first_res$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))[2]
  x <- sgRNA$sgRNA[sgRNA$id2 == x]
  sgCmp$Cmp[sgCmp$sgRNA == x]
}))
ToNX::write_tb(bulge_append_first_res, file="~/data/project/ear_project/gene_therapy_ll/cor_sample_reDiff_all.txt")
bulge_append_first_res <- read.table("~/data/project/ear_project/gene_therapy_ll/cor_sample_reDiff_all_output.txt", 
                               sep="\t")
bulge_append_first_res <- bulge_append_first_res[bulge_append_first_res$V8 %in% c(3, 4),]
colnames(bulge_append_first_res) <- c("pair", "cor",  "pval", "id",   "pos",  "cmp1", "cmp2", "diffPos")

bulge_append_first_res <- bulge_append_first_res[unlist(lapply(1 : nrow(bulge_append_first_res), function(i){
  cmps <- unlist(bulge_append_first_res[i, c("cmp1", "cmp2")])
  cmp <- cmps[1]
  cmp_pair <- cmps[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  return(which(cmp == '-') == 1 | which(cmp_pair=='-') == 1)
})),]

####Read sg28 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/batch1/")
samples <- list.files(pattern = "^S")
total_edit_table <- list()
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  sgs <- list.files(pattern = "^sg")
  for(sg in sgs){
    setwd(sg)
    setwd("CRISPResso_on_nhej")
    filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
    filename <- filename[which.max(str_length(filename))]
    edit_table <- read.table(filename, 
                             sep="\t", header = T, comment.char = "")
    
    sgnum <- str_remove(str_remove(sg, "sg"), "_split")
    sg <- sgRNA$id2[as.numeric(sgnum)]
    total_edit_table[[sample]][[sg]] <- edit_table
    setwd("../..")
  }
  setwd("..")
  rm(edit_table, sg, sgs)
}

####Read batch1_1 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
sgRNA_tmp <- read.table("sg_info.txt")
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    if(sample == "CRISPRessoPooled_on_B29" & sg %in% c("Sg_12_79", 
                                                       "Sg_12_80", 
                                                       "Sg_12_81")){
      print("rm sample")
      next
    }
    else{
      setwd(paste0("CRISPResso_on_",sg))
      filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
      filename <- filename[which.max(str_length(filename))]
      edit_table <- read.table(filename, 
                               sep="\t", header = T, comment.char = "")
      total_edit_table[[sample]][[sg]] <- edit_table
      setwd("..")
    }
    
  }
  setwd("..")
}
###Read batch1_2 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
sgRNA_tmp <- read.table("sg_info.txt")
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    if(sg %in% c("Sg_23_159", "Sg_19_133")){
      print("rm sample")
      next
    }
    if(sample == "CRISPRessoPooled_on_A99" & sg %in% c("Sg_23_157")){
      print("rm sample")
      next
    }
    if(sample == "CRISPRessoPooled_on_B99" & sg %in% c("Sg_20_137", 
                                                       "Sg_20_138", 
                                                       "Sg_20_139")){
      print("rm sample")
      next
    }
    if(sample == "CRISPRessoPooled_on_C99" & sg %in% c("Sg_19_132")){
      print("rm sample")
      next
    }
    setwd(paste0("CRISPResso_on_",sg))
    filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
    filename <- filename[which.max(str_length(filename))]
    edit_table <- read.table(filename, 
                             sep="\t", header = T, comment.char = "")
    total_edit_table[[sample]][[sg]] <- edit_table
    setwd("..")
  }
  setwd("..")
}
rm(sample, samples, sg, sgRNA_tmp, edit_table)

total_edit_table_tmp <- total_edit_table
total_edit_table <- list()
for(i in 1 : 9){
  rep <- paste0("Rep", rep(1:3,3)[i])
  if(!rep %in% names(total_edit_table)){
    total_edit_table[[rep]] <- list()
  }
  for(sg in names(total_edit_table_tmp[[i]])){
    total_edit_table[[rep]][[sg]] <- total_edit_table_tmp[[i]][[sg]]
  }
}
sg_sel <- unlist(lapply(bulge_append_first_res$pair, function(x){
  unlist(strsplit(x, "[-]"))
}))
total_edit_table_sel <- lapply(total_edit_table, function(x){
  x[names(x) %in% sg_sel]
})
total_edit_table_tmp <- total_edit_table
total_edit_table <- total_edit_table_sel
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


pdf("~/data/project/ear_project/gene_therapy_ll/high_cor/indel_in_indel_append_first_result.pdf", width = 20, height = 12)
for(pairids in bulge_append_first_res$pair){
  ids <- unlist(strsplit(pairids, "[-]"))
  #put insert one at front
  sg_seqs <- unlist(bulge_append_first_res[bulge_append_first_res$pair == pairids, 
                                     c("cmp1", "cmp2")])
  cmp <- sg_seqs[1]
  cmp_pair <- sg_seqs[2]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  #在第一位有一个插入情况下，说明后面可能是last或者非last的情况
  if(this_pos == 1 | pair_pos == 1){
    
    #说明是last位有插入 且在这个sgRNA插入
    if(pair_pos == length(cmp)){
      ids <- ids
    }
    #说明是last位有插入 且不在这个sgRNA插入
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
  mean_stat1$sg <- 'sg1'
  mean_stat2$sg <- 'sg2'
  merged_mean_stat <- data.frame(rbind(mean_stat1, mean_stat2))
  colors <- hue_pal()(6)
  names(colors) <- c("Ins1", "Ins2", "Ins3", "Del1", "Del2", "Del3")
  wt_seq1 <- total_edit_table$Rep3[[ids[1]]]
  wt_seq1 <- getWtSeq(wt_seq1[wt_seq1[,3] == "True",])
  wt_seq2 <- total_edit_table$Rep3[[ids[2]]]
  wt_seq2 <- getWtSeq(wt_seq2[wt_seq2[,3] == "True",])
  plot_list <- lapply(names(mean_stat), function(x){
    plot_data <- merged_mean_stat[merged_mean_stat$type == x,]
    maxPct <- ceiling(max(c(plot_data$Pct + plot_data$se)))
    plot_data <- plot_data[plot_data$pos %in% as.character(c(7 : 26)),]
    plot_data$pos <- factor(plot_data$pos, 
                            levels = gtools::mixedsort(unique(plot_data$pos)))
    se_data <- plot_data[plot_data$Pct != 0 & plot_data$se > 0.01,]
    
    color_seq <- mkColorfulSeq(wt_seq1[7:26], wt_seq2[7:26], sg_seqs[1], sg_seqs[2])
    
    p <- ggplot(plot_data, aes(x=pos, y=Pct,group = type)) + 
      geom_vline(xintercept = 14.5, linetype = "dashed") + 
      geom_line(color = colors[x])
    if(nrow(se_data) > 0){
      p <- p +       geom_errorbar(data = se_data,
                                   aes(ymin=Pct-se, ymax=Pct+se), width=.2, color = "black")
    }
 
      p + ggplot2::facet_wrap(~sg, nrow = 1, ncol = 2, scales = 'free_x') + 
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
  merge_plot2 <- annotate_figure(merge_plot2, top = text_grob(paste0("Cor: ", round(bulge_append_first_res$cor[bulge_append_first_res$pair == pairids], 3)), 
                                                              color = "black", face = "bold", size = 14))
  rm(p1, p2, p_l1, p_l2)
  print(aplot::insert_right(merge_plot2, merge_plot, width = 1))
  
}
dev.off()



