setwd("~/data/project/ear_project/gene_therapy_ll/batch2/")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/indel_stat2.rda")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599//indel_stat.rda")
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/First1NT/2024_1_12_integrated_design_result.xlsx")
sgRNA <- sgRNA[c(29 : nrow(sgRNA)),]
sgRNA$id2 <- paste("Sg", sgRNA$new_id, sep = "-")
sgRNA$id2 <- str_replace_all(sgRNA$id2, "-", "_")
max_region <- c(-200, 30)
batch2_total_mean <- batch2_mean
for(name in names(batch2_mean2)){
  batch2_total_mean[[name]] <- batch2_mean2[[name]]
}
###去掉那些不匹配的样本
sgRNA_pair <- split(sgRNA$id2, sgRNA$id)
sgRNA_pair_rm <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(batch2_total_mean)) != 2)
}))]
sgRNA_pair_remain <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(batch2_total_mean)) == 2)
}))]
###计算均值的相关性，因为本次数据每次样本是分开测得 不好AA BB CC一一对比了
total_pair_sg_cor <- lapply(sgRNA_pair_remain, function(y){
    ids <- unlist(y)
    sel_sg1 <- batch2_total_mean[[ids[1]]]
    sel_sg2 <- batch2_total_mean[[ids[2]]]
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
order_name <- unlist(lapply(sgRNA_pair_remain, function(x){
  paste(x, collapse = "-")
}))
names(sgRNA_pair_remain) <- order_name
order_name <- gtools::mixedsort(order_name)
pdf("Result/sgRNA_pair_edit_pattern_mean_cor.pdf", width = 10, height = 6)
for(name in order_name){
  ids <- unlist(sgRNA_pair_remain[[name]])
  id <- paste(ids, collapse = "-")
  plot_data <- lapply(ids, function(x){
    region <- c(-19: 19)
    sel_sg <- batch2_total_mean[[x]]
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
    sel_sg$sg <- x
    sel_sg$pct <- sel_sg[,2] / sum(sel_sg[,2])
    sel_sg
  })

  plot_data <- data.frame(do.call(rbind, plot_data))
  plot_data$indel_size <- factor(plot_data$indel_size, levels = c("<=-20", -19 : 19, ">=20"))
  print(ggplot(plot_data, aes(x = indel_size, y = pct * 100, fill = sg)) + 
          geom_bar(stat="identity", position = "dodge") + 
          ylab("% in indel reads") + xlab("") + 
          ggtitle(paste0(id, " Cor: ", 
                         round(total_pair_sg_cor$cor[total_pair_sg_cor$pair == id], 2))))
  
}
dev.off()

total_pair_sg_cor$pos <- as.integer(total_pair_sg_cor$pos)
plot_data <- total_pair_sg_cor[order(total_pair_sg_cor$pos, total_pair_sg_cor$cor),]
plot_data$pair <- factor(plot_data$pair, levels = plot_data$pair)
plot_data$pos <- factor(as.character(plot_data$pos), levels = as.character(unique(plot_data$pos)))
sample_counts <- data.frame(table(plot_data$pos))

pdf("Result/cor_mean_point_plot_fix.pdf", width = 10, height = 4)
ggplot(plot_data, aes(x = pos, y = cor)) +
  geom_boxplot(outlier.size = 0, aes(color = pos)) + 
  geom_point(aes(color = pos), position = position_dodge2(width = 0.7)) + 
  geom_text(data = sample_counts, 
            aes(x = Var1, y = 1, label = paste0("Counts:", Freq))) + 
  theme_bw() + 
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        text = element_text(size=15))
dev.off()
