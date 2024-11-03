load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/merged_indel_table_first_second.Rda")
###去掉那些不匹配的样本
sgCmp_sel <- sgCmp[!sgCmp$id2 %in%c("Sg_21_144","Sg_6_37","Sg_6_38","Sg_7_45","Sg_17_115"), ]
sgRNA_pair <- split(sgCmp_sel$id2, sgCmp_sel$id)

sgRNA_pair_rm <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(total_mean_first_second)) != 2)
}))]
sgRNA_pair_remain <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(total_mean_first_second)) == 2)
}))]
# sgCmp132 <- sgCmp[sgCmp$id2 %in% unlist(sgRNA_pair_remain),]
# save(sgCmp132 , file = "~/data/project/ear_project/gene_therapy_ll/Result/sgCmp132.rda")
c("Sg_21_144","Sg_6_37","Sg_6_38","Sg_7_45","Sg_17_115") %in% unlist(sgRNA_pair_remain)
###计算均值的相关性，因为本次数据每次样本是分开测得 不好AA BB CC一一对比了
max_region <- c(-200, 30)
total_pair_sg_cor <- lapply(sgRNA_pair_remain, function(y){
  ids <- unlist(y)
  sel_sg1 <- total_mean_first_second[[ids[1]]]
  sel_sg2 <- total_mean_first_second[[ids[2]]]
  max_del <- min(c(sel_sg1[,1], sel_sg2[,1]))
  max_in <- max(c(sel_sg1[,1], sel_sg2[,1]))
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
table(plot_data$cor >= 0.6)


plot_region <- c(-18 : 12)
with_bulge_sg <- unlist(lapply(sgRNA_pair_remain, function(x){
  x <- unlist(x)
  x <- sgCmp[sgCmp$id2 %in% x,]
  x$id2[x$isInsert == "Insert"]
}))
without_bulge_sg <- unlist(lapply(sgRNA_pair_remain, function(x){
  x <- unlist(x)
  x <- sgCmp[sgCmp$id2 %in% x,]
  x$id2[x$isInsert != "Insert"]
}))

mean_indel_table_sel_with_bulge <- lapply(with_bulge_sg, function(x){
  res <- total_mean_first_second[[x]]
  res[,3] <- res[,2] / sum(res[,2]) * 100
  res[res[,1] %in% plot_region,]
})
names(mean_indel_table_sel_with_bulge) <- with_bulge_sg
mean_indel_table_sel_without_bulge <- lapply(without_bulge_sg, function(x){
  res <- total_mean_first_second[[x]]
  res[,3] <- res[,2] / sum(res[,2]) * 100
  res[res[,1] %in% plot_region,]
})
names(mean_indel_table_sel_without_bulge) <- without_bulge_sg
total_pair_sg_cor$bulgeOne <- unlist(lapply(total_pair_sg_cor$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))
  x <- sgCmp[sgCmp$id2 %in% x,]
  x$id2[x$isInsert == "Insert"]
}))
total_pair_sg_cor$notBulgeOne <- unlist(lapply(total_pair_sg_cor$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))
  x <- sgCmp[sgCmp$id2 %in% x,]
  x$id2[x$isInsert != "Insert"]
}))
total_pair_sg_cor <- total_pair_sg_cor[order(total_pair_sg_cor$cor, decreasing = T),]

plot_mat <- matrix(0, ncol = length(plot_region) * 2 - 1, nrow = nrow(total_pair_sg_cor))
for(i in 1 : nrow(total_pair_sg_cor)){
  tmp1 <- mean_indel_table_sel_with_bulge[[total_pair_sg_cor$bulgeOne[i]]]
  tmp1 <- tmp1[tmp1[,1] != 0,]
  tmp2 <- mean_indel_table_sel_without_bulge[[total_pair_sg_cor$notBulgeOne[i]]]
  tmp2 <- tmp2[tmp2[,1] != 0,]
  plot_mat[i, 1 : (length(plot_region) - 1)] <- rev(tmp1[,3]) / 100
  plot_mat[i, length(plot_region)] <- total_pair_sg_cor$cor[i]
  plot_mat[i, (length(plot_region) + 1) : ncol(plot_mat)] <- tmp2[,3] / 100
}

col_acc <- colSums(plot_mat)
col_acc[length(plot_region)] <- NA
plot_data <- data.frame(Pct = c(col_acc[-length(plot_region)], 0, 0), 
                        indel_size = c(12 : 1, -1 : -18, -18 : -1, 1 : 12, 0, 0), 
                        type = c(rep(c("NoBulge", "withBulge"), c(30, 30)), "NoBulge", "withBulge"))
max_ratio <- max(plot_data$Pct)
x_label <- c(-18 : 12)
x_label[seq(3, length(x_label), 3)] <- ""
x_label[seq(2, length(x_label), 3)] <- ""
plot_data$indel_size <- factor(plot_data$indel_size, levels = -18 : 12)
p1 <- ggplot(plot_data[plot_data$type == "NoBulge",], aes(y = indel_size, x = Pct)) + 
  geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
  xlab("Sum % in indel reads of noBulge") + ylab("") + 
  scale_x_reverse(expand = c(0,0), limits = c(max_ratio, 0)) + 
  scale_y_discrete(position = "right") + 
  theme_classic2() + 
  theme(axis.text.y = element_blank(),
        plot.margin =  unit(c(0.2,0,0,0.2), "cm"))

p2 <- ggplot(plot_data[plot_data$type != "NoBulge",], aes(y = indel_size, x = Pct)) + 
  geom_bar(stat="identity", fill = "#A8A8A8", color = "black", width = 0.75) + 
  scale_x_continuous(expand = c(0,0),  limits = c( 0, max_ratio)) + 
  scale_y_discrete(labels = x_label) + 
  xlab("Sum % in indel reads of withBulge") + ylab("") + 
  theme_classic2() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5) , 
        plot.margin =  unit(c(0.2,0.75,0,-1), "cm"))

merge_plot2 <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T)
pdf("~/Nutstore Files/Tobin/Merged1NT/paired_indel_size_hist.pdf", width = 6, height = 4)
print(merge_plot2)
dev.off()
