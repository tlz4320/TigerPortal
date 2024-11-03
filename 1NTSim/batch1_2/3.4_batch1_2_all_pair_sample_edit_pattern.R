#需要来自3.1的total_pair_sg_cor对象

setwd("~/data/project/ear_project/gene_therapy_ll/batch2/")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/indel_stat2.rda")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599//indel_stat.rda")
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/First1NT/2024_1_12_integrated_design_result.xlsx")
sgRNA <- sgRNA[c(29 : nrow(sgRNA)),]
sgRNA$id2 <- paste("Sg", sgRNA$new_id, sep = "-")
sgRNA$id2 <- str_replace_all(sgRNA$id2, "-", "_")
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

###去掉那些不匹配的样本
sgRNA_pair <- split(sgRNA$id2, sgRNA$id)
sgRNA_pair_rm <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(batch2_rev)) != 2)
}))]
sgRNA_pair_remain <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(batch2_rev)) == 2)
}))]

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
pdf("Result/all_pair_sample_result_fix.pdf", width = 10, height = 8)
for(name in total_pair_sg_cor$id){
  ids <- sgRNA_pair_remain[[name]]
  pos <- sg_name$pos[sg_name$id == name]
  id <- paste(ids, collapse = "-")
  cols <- scales::hue_pal()(3)
  names(cols) <- c("A", "B", "C")
  
  p_l1 <- lapply(names(batch2_rev[[ids[1]]]), function(n){
    region <- c(-19: 19)
    sel_sg <- batch2_rev[[ids[1]]][[n]]
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
  
  p_l2 <- lapply(names(batch2_rev[[ids[2]]]), function(n){
    region <- c(-19: 19)
    sel_sg <- batch2_rev[[ids[2]]][[n]]
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
















