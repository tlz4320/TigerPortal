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


plot_region <- c(-8: 5)
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

library(ComplexHeatmap)
library(circlize)
max(unlist(plot_mat[,-length(plot_region)]))
indel_color <- colorRamp2(c(0, 1), c("white", "blue"))
cor_color <- colorRamp2(c(0, 1), c("white", "red"))

indel_reads_with_bulge <- lapply(with_bulge_sg, function(x){
  tmp <- total_indel_table_first_second_rev[[x]]
  tmp <- unlist(lapply(tmp, function(y){
    sum(y[y[,1] %in% plot_region & y[,1] != 0,2])
  }))
  data.frame(count = log(mean(tmp)), se = plotrix::std.error(log(tmp)))
})
indel_reads_with_bulge <- data.frame(do.call(rbind, indel_reads_with_bulge))
rownames(indel_reads_with_bulge) <- with_bulge_sg
indel_reads_with_bulge <- indel_reads_with_bulge[total_pair_sg_cor$bulgeOne,]

indel_reads_without_bulge <- lapply(without_bulge_sg, function(x){
  tmp <- total_indel_table_first_second_rev[[x]]
  tmp <- unlist(lapply(tmp, function(y){
    sum(y[y[,1] %in% plot_region & y[,1] != 0,2])
  }))
  data.frame(count = log(mean(tmp)), se = plotrix::std.error(log(tmp)))
})
indel_reads_without_bulge <- data.frame(do.call(rbind, indel_reads_without_bulge))
rownames(indel_reads_without_bulge) <- without_bulge_sg
indel_reads_without_bulge <- indel_reads_without_bulge[total_pair_sg_cor$notBulgeOne,]

plot_mat <- plot_mat[,c(5,6,22,23)]
rownames(plot_mat) <- total_pair_sg_cor$id
total_pair_sg_cor$ins1 <- plot_mat[,4]
total_pair_sg_cor$del1 <- plot_mat[,2]
library(ggbeeswarm)
pdf("~/Nutstore Files/Tobin/Merged1NT/Indel1_pct_132_v2.pdf", width = 8,height = 6)
ggplot(total_pair_sg_cor) + geom_boxplot(aes(x = "Ins1", y = ins1), outlier.shape = NA) + 
  geom_quasirandom(aes(x = "Ins1", y = ins1), shape = 21, size =2 ,fill = "gray90") + 
  geom_boxplot(aes(x = "Del1", y = del1), outlier.shape = NA) + 
  geom_quasirandom(aes(x = "Del1", y = del1),shape = 21, size = 2, fill = "gray90") + xlab("") + ylab("Percent(%)") + 
  geom_hline(yintercept = c(0.09,0.2, 0.3, 0.4)) + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()
total_pair_sg_cor <- total_pair_sg_cor[,c(-4,-6,-7)]
openxlsx::write.xlsx(total_pair_sg_cor, file="~/Nutstore Files/Tobin/Merged1NT/Indel1_stat_data.xlsx", rowNames=F,colNames=T)
