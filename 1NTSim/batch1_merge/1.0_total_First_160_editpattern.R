print(load("~/data/project/ear_project/gene_therapy_ll/batch1/Result/indel_stat_total.rda"))
names(sg28) <- str_remove(names(sg28), "_.*")
sgs <- paste0("sg", c(1 : 7, 10, 11, 12, 14 : 17, 19:21,24,25, 27))

sg28_single_sel <- sg28[sgs]
sg28_single_sel <- lapply(sg28_single_sel, function(x){
    samples <- str_remove(names(x), "./")
    names(x) <- samples
    multi_samples <- samples[grep("^B", samples)]
    single_samples <- x[!samples %in% multi_samples]
    names(single_samples) <- c("A", "B", "C")
    single_samples[unlist(lapply(single_samples, function(y){
      sum(y[,2]) != 0
    }))]
})
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/First1NT/2024_1_12_integrated_design_result.xlsx")
sgRNA$id2 <- paste("Sg", sgRNA$new_id, sep = "-")
sgRNA$id2 <- str_replace_all(sgRNA$id2, "-", "_")
sg28 <- read.xlsx("~/data/project/ear_project/gene_therapy/find1NTSim/sg28_redesign_fixbad_fix_sgprimer.xlsx")
sg28$order <- 1 : 28
table(sg28$sgRNA == sgRNA$sgRNA[1:28])
for(i in 1 : length(sg28_single_sel)){
  n <- as.numeric(str_remove(names(sg28_single_sel)[i], "sg"))
  names(sg28_single_sel)[i] <- sgRNA$id2[n]
}



setwd("~/data/project/ear_project/gene_therapy_ll/batch2/")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/indel_stat2.rda")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599//indel_stat.rda")
max_region <- c(-200, 30)
sg28_mean <- lapply(sg28_single_sel, function(x){
  if(length(x) == 1){
    return(x[[1]])
  }
  fq <- unlist(rowMeans(data.frame(do.call(cbind, lapply(x, function(y){y[,2]})))))
  x[[1]][,2] <- fq
  x[[1]]
})
batch2_total_mean <- batch2_mean
for(name in names(batch2_mean2)){
  batch2_total_mean[[name]] <- batch2_mean2[[name]]
}
total_mean <- batch2_total_mean
for(name in names(sg28_mean)){
  total_mean[[name]] <- sg28_mean[[name]]
}
###去掉那些不匹配的样本
sgRNA_pair <- split(sgRNA$id2, sgRNA$id)
sgRNA_pair_rm <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(total_mean)) != 2)
}))]
sgRNA_pair_remain <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(total_mean)) == 2)
}))]
###计算均值的相关性，因为本次数据每次样本是分开测得 不好AA BB CC一一对比了
total_pair_sg_cor <- lapply(sgRNA_pair_remain, function(y){
  ids <- unlist(y)
  sel_sg1 <- total_mean[[ids[1]]]
  sel_sg2 <- total_mean[[ids[2]]]
  max_del <- min(c(sel_sg1[,1], sel_sg2[,1]))
  max_del <- max(c(max_del, max_region[1]))
  max_in <- max(c(sel_sg1[,1], sel_sg2[,1]))
  max_in <- min(c(max_in, max_region[2]))
  region <- max_del : max_in
  tmp_region <- region[!region %in% sel_sg1[,1]]
  if(length(tmp_region) != 0){
    sel_sg1 <- data.frame(rbind(sel_sg1, data.frame(indel_size = tmp_region, fq = 0)))
  }
  tmp_region <- region[!region %in% sel_sg2[,1]]
  if(length(tmp_region) != 0){
    sel_sg2 <- data.frame(rbind(sel_sg2, data.frame(indel_size = tmp_region, fq = 0)))
  }
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  common_indel <- intersect(sel_sg1[,1], sel_sg2[,1])
  sel_sg1 <- sel_sg1[sel_sg1[,1] %in% common_indel,]
  sel_sg2 <- sel_sg2[sel_sg2[,1] %in% common_indel,]
  sel_sg1 <- sel_sg1[order(sel_sg1[,1]),]
  sel_sg2 <- sel_sg2[order(sel_sg2[,1]),]
  # rm0 <- sel_sg1[,2] == 0 & sel_sg2[,2] == 0
  #cor.test(sel_sg1[!rm0,2], sel_sg2[!rm0,2])
  cor.test(sel_sg1[,2], sel_sg2[,2])
})
names(total_pair_sg_cor) <- unlist(lapply(sgRNA_pair_remain, function(x){
  paste(x, collapse = "-")
}))
total_pair_sg_cor <- data.frame(pair = names(total_pair_sg_cor), 
                                cor = unlist(lapply(total_pair_sg_cor, function(x){x$estimate})),
                                pval = unlist(lapply(total_pair_sg_cor, function(x){x$p.value})))

sg_name <- data.frame(id = names(sgRNA_pair_remain), 
                      pair = unlist(lapply(sgRNA_pair_remain, function(x){
                        paste(x, collapse = "-")
                      })))
sg_name$pos <- unlist(lapply(sg_name$id, function(x){
  if(!is.na(str_match(x, "last"))){
    return(1)
  }
  unlist(strsplit(x, "[_-]"))[2]
}))
total_pair_sg_cor <- merge(total_pair_sg_cor, sg_name, by = "pair")

total_pair_sg_cor$pos <- as.integer(total_pair_sg_cor$pos)
plot_data <- total_pair_sg_cor[order(total_pair_sg_cor$pos, total_pair_sg_cor$cor),]
plot_data$pair <- factor(plot_data$pair, levels = plot_data$pair)
plot_data$pos <- factor(as.character(plot_data$pos), levels = as.character(unique(plot_data$pos)))
sample_counts <- data.frame(table(plot_data$pos))

pdf("Result/cor_mean_point_plot_add_sg28.pdf", width = 10, height = 4)
ggplot(plot_data, aes(x = pos, y = cor)) +
  geom_boxplot(outlier.size = 0, aes(color = pos)) + 
  geom_point(aes(color = pos), position = position_dodge2(width = 0.7)) + 
  geom_text(data = sample_counts, 
            aes(x = Var1, y = 1, label = paste0("Counts:", Freq))) + 
  theme_bw() + 
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        text = element_text(size=15))
dev.off()


setwd("~/data/project/ear_project/gene_therapy_ll/batch2/")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/indel_stat2.rda")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599//indel_stat.rda")

max_region <- c(-200, 30)
batch2 <- batch2_1
for(name in names(batch2_2)){
  batch2[[name]] <- batch2_2[[name]]
}
total_sg <- unique(unlist(lapply(batch2, names)))
batch2_rev <- lapply(total_sg, function(sg){
  res <- list()
  for(name in names(batch2)){
    n <- str_remove(name, "CRISPRessoPooled_on_")
    n <- str_remove(n, "[0-9]+")
    if(sg %in% names(batch2[[name]])){
      res[[n]] <- batch2[[name]][[sg]]
    }
  }
  res
})
names(batch2_rev) <- total_sg
total_data <- batch2_rev
for(name in names(sg28_single_sel)){
  total_data[[name]] <- sg28_single_sel[[name]]
}


sg_name <- data.frame(id = names(sgRNA_pair_remain), 
                      pair = unlist(lapply(sgRNA_pair_remain, function(x){
                        paste(x, collapse = "-")
                      })))
sg_name$pos <- unlist(lapply(sg_name$id, function(x){
  if(!is.na(str_match(x, "last"))){
    return(1)
  }
  unlist(strsplit(x, "[_-]"))[2]
}))
sg_name <- sg_name[order(sg_name$pos, sg_name$pair),]
sg_name <- sg_name[sg_name$id %in% names(sgRNA_pair_remain),]

total_pair_sg_cor <- total_pair_sg_cor[order(total_pair_sg_cor$pos, total_pair_sg_cor$cor),]

####绘制所有结果
pdf("Result/all_pair_sample_result_add_sg28.pdf", width = 10, height = 8)
for(name in total_pair_sg_cor$id){
  ids <- sgRNA_pair_remain[[name]]
  pos <- sg_name$pos[sg_name$id == name]
  id <- paste(ids, collapse = "-")
  cols <- scales::hue_pal()(3)
  names(cols) <- c("A", "B", "C")
  
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
    sel_sg$pct <- sel_sg[,2] / sum(sel_sg[,2])
    sel_sg
  })
  p_l1 <- data.frame(do.call(rbind, p_l1))
  p_l1$indel_size <- factor(p_l1$indel_size, 
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
    sel_sg$pct <- sel_sg[,2] / sum(sel_sg[,2])
    sel_sg
  })
  p_l2 <- data.frame(do.call(rbind, p_l2))
  p_l2$indel_size <- factor(p_l2$indel_size, 
                            levels = c("<=-20", -19 : 19, ">=20"))
  
  max_ratio <- max(c(p_l1$pct, p_l2$pct))
  p1 <- 
    ggplot(p_l1, aes(x = indel_size, y = pct * 100, fill = sample)) + 
    geom_bar(stat="identity", position = "dodge") + 
    ylab(paste0("% in indel reads of ", ids[1])) + xlab("") + 
    scale_y_continuous(expand = c(0,0), limits = c(0, max_ratio * 100)) + 
    scale_fill_manual(values = cols) + 
    ggtitle(paste0(ids[1], "-", ids[2], " Bulge Pos: ", pos, 
                   " Cor: ", 
                   round(total_pair_sg_cor$cor[total_pair_sg_cor$pair == id], 2))) + 
    theme_classic2() + 
    theme(plot.margin =  unit(c(0,0,0,0.2), "cm"))
  
  p2 <- ggplot(p_l2, aes(x = indel_size, y = pct * 100, fill = sample)) + 
    geom_bar(stat="identity", position = "dodge") + 
    scale_y_reverse(expand = c(0,0),  limits = c(max_ratio * 100, 0)) + 
    scale_x_discrete(position = "top") + 
    scale_fill_manual(values = cols) + 
    ylab(paste0( " % in indel reads of ", ids[2])) + xlab("") + 
    theme_classic2() + 
    theme(axis.text.x = element_blank(), 
          plot.margin =  unit(c(-1,0,0,0.2), "cm"))
  print(ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T))
  rm(p1, p2, p_l1, p_l2)
}
dev.off()







