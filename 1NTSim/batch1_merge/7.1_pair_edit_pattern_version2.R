#pdf("~/data/project/ear_project/gene_therapy_ll/high_cor/edit_pattern_version1.pdf", width = 6, height = 4)
sel_pair <- high_cor_fixDiff[high_cor_fixDiff$cor > 0.95,]
for(pairids in sel_pair$pair){
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

  
  ####画柱状图
  p_l1 <- lapply(names(total_data[[ids[1]]]), function(n){
    region <- c(-18: 12)
    sel_sg <- total_data[[ids[1]]][[n]]
    long_del <- sum(sel_sg[sel_sg[,1] < min(region), 2])
    long_ins <- sum(sel_sg[sel_sg[,1] > max(region), 2])
    long_res <- data.frame(indel_size = c("<=-19", ">=13"), fq = c(long_del, long_ins))
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
  p_l1_mean <- p_l1_mean[p_l1_mean$indel_size %in% as.character(-18:12),]
  p_l1_mean$indel_size <- factor(p_l1_mean$indel_size, 
                                 levels = as.character(-18:12))
  
  p_l2 <- lapply(names(total_data[[ids[2]]]), function(n){
    region <- c(-18: 12)
    sel_sg <- total_data[[ids[2]]][[n]]
    long_del <- sum(sel_sg[sel_sg[,1] < min(region), 2])
    long_ins <- sum(sel_sg[sel_sg[,1] > max(region), 2])
    long_res <- data.frame(indel_size = c("<=-19", ">=13"), fq = c(long_del, long_ins))
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
  p_l2_mean <- p_l2_mean[p_l2_mean$indel_size %in% as.character(-18:12),]
  p_l2_mean$indel_size <- factor(p_l2_mean$indel_size, 
                                 levels = as.character(-18:12))
  
  
  max_ratio <- max(c(p_l2_mean$Pct + p_l2_mean$se, p_l1_mean$Pct + p_l1_mean$se))
  x_label <- c(-18 : 12)
  x_label[seq(3, length(x_label), 3)] <- ""
  x_label[seq(2, length(x_label), 3)] <- ""
  p1 <- ggplot(p_l1_mean, aes(y = indel_size, x = Pct)) + 
    geom_errorbar(aes(xmin=Pct-se, xmax=Pct+se), width=.2, color = "black")  + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    xlab(paste0("% in indel reads of ", ids[1])) + ylab("") + 
    scale_x_reverse(expand = c(0,0), limits = c(max_ratio, 0)) + 
    scale_y_discrete(position = "right") + 
    theme_classic2() + 
    theme(axis.text.y = element_blank(),
          plot.margin =  unit(c(0.2,0,0,0.2), "cm"))
  
  p2 <- ggplot(p_l2_mean, aes(y = indel_size, x = Pct)) + 
    geom_errorbar(aes(xmin=Pct-se, xmax=Pct+se), width=.2, color = "black")  + 
    geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
    scale_x_continuous(expand = c(0,0),  limits = c( 0, max_ratio)) + 
    scale_y_discrete(labels = x_label) + 
    xlab(paste0( " % in indel reads of ", ids[2])) + ylab("") + 
    theme_classic2() + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5) , 
          plot.margin =  unit(c(0.2,0.75,0,-1), "cm"))

  merge_plot2 <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T)
  merge_plot2 <- annotate_figure(merge_plot2, top = text_grob(paste0(pairids, "  Cor: ", round(high_cor_fixDiff$cor[high_cor_fixDiff$pair == pairids], 3)), 
                                                              color = "black", face = "bold", size = 14))
  pdf(paste0("~/data/project/ear_project/gene_therapy_ll/high_cor/edit_pattern_version2_", pairids,".pdf"), width = 6, height = 4)
  print(merge_plot2)
  dev.off()
}
dev.off()

