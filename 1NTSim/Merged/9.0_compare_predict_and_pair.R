ins1_lr_nt_stat <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Merged1NT/merged_ins1_left_right_nt_pie_data.xlsx")
###在2.0脚本之后运行
remain_ins1_pct_table_sel <- remain_ins1_pct_table[sg_info_sel$id]

nt_counts_mats <- lapply(remain_ins1_pct_table_sel, function(tmp){
  plot_mat <- matrix(0, ncol = 17, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(20, 21),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  tmp_left <- tmp[tmp$pos <= 20,]
  tmp_right <- tmp[tmp$pos > 20,]
  
  #首先把15-20的算一下，15等于是-6～-5之间插入
  #20则是-1～1之间插入
  #之后算21-23的  21也是-1～1之间插入  所以得分开算一下
  if(nrow(tmp_left) > 0){
    for(i in 1 : nrow(tmp_left)){
      mat_pos <- (tmp_left$pos[i] - 14) * 2
      plot_mat[tmp_left$NT[i], mat_pos] <- plot_mat[tmp_left$NT[i], mat_pos] + tmp_left$Pct[i]
    }
  }
  
  if(nrow(tmp_right) > 0){
    for(i in 1 : nrow(tmp_right)){
      mat_pos <- (tmp_right$pos[i] - 14) * 2 - 2
      plot_mat[tmp_right$NT[i], mat_pos] <- plot_mat[tmp_right$NT[i], mat_pos] + tmp_right$Pct[i]
    }
  }
  
  plot_mat
  
})


ins1_real_vs_pie <- sg_info_sel
ins1_real_vs_pie$diff <- unlist(lapply(ins1_real_vs_pie$id, function(id){
  left_nt <- seq_mat[id, 11]
  right_nt <- seq_mat[id, 13]
  tmp_counts_mat <- nt_counts_mats[[id]]
  ins1_pct <- tmp_counts_mat[,12]
  pie_pred <- unlist(ins1_lr_nt_stat[ins1_lr_nt_stat$leftNT == left_nt & 
                                ins1_lr_nt_stat$rightNT == right_nt, c(-1,-2)])
  pie_pred <- pie_pred[names(ins1_pct)]
  ins1_pct <- ins1_pct / sum(ins1_pct) * 100
  sum(abs(ins1_pct - pie_pred))
}))

paired_predict <- sg_info_sel
paired_predict$pair <- as.character(paired_predict$pair)
paired_predict <- lapply(split(paired_predict$id, paired_predict$pair), function(ids){
  ids <- unlist(ids)
  id1 <- ids[1]
  left_nt1 <- seq_mat[id1, 11]
  right_nt1 <- seq_mat[id1, 13]
  tmp_counts_mat <- nt_counts_mats[[id1]]
  ins1_pct1 <- tmp_counts_mat[,12]
  ins1_pct1 <- ins1_pct1 / sum(ins1_pct1) * 100
  
  id2 <- ids[2]
  left_nt2 <- seq_mat[id2, 11]
  right_nt2 <- seq_mat[id2, 13]
  tmp_counts_mat <- nt_counts_mats[[id2]]
  ins1_pct2 <- tmp_counts_mat[,12]
  ins1_pct2 <- ins1_pct2 / sum(ins1_pct2) * 100
  
  c(sum(abs(ins1_pct1 - ins1_pct2)), left_nt1 == left_nt2, right_nt1 == right_nt2)
})
paired_predict <- data.frame(do.call(rbind, paired_predict))
colnames(paired_predict) <- c("pair_diff", "sameLeft", "sameRight")
paired_predict <- paired_predict[paired_predict$sameLeft == 1,]
paired_predict$id <- rownames(paired_predict)
###后面是读取预测的结果
setwd('~/indelphi_res2/')
indelphi_132_res <- list()
for(file in list.files()){
  indelphi_132_res[[file]] <- read.table(file,sep = ",", header = T)
}
names(indelphi_132_res) <- str_remove(names(indelphi_132_res), ".csv")
indelphi_res_ins1 <- lapply(indelphi_132_res, function(x){
  x <- x[x$Category == "ins" & x$Length == 1,]
  x$Pct <- x$Predicted.frequency / sum(x$Predicted.frequency) * 100
  colnames(x)[4] <- "NT"
  x
})
# save(indelphi_132_res, file="~/data/project/ear_project/gene_therapy_ll/Result/132_indelphi.rda")

setwd('~/forecast132_result/')
forcast_file <- read.table("~/132_for_forecast.txt", header = F)
forecast_132_res_ins1 <- list()
for(file in forcast_file$V1){
  tmp1 <- read.table(paste0(file, "_predictedindelsummary.txt"), 
                     sep = "\t", header = T)
  colnames(tmp1) <- c("V1", "NT", "Count")
  tmp2 <- read.table(paste0(file, "_predictedreads.txt"), 
                     sep = "\t", header = F)
  sg <- sgCmp$sgRNA[sgCmp$id2 == file]
  cutpos <- unlist(str_locate(tmp2[1,2], sg)[1]) + 17
  tmp2 <- merge(tmp2, tmp1, by.y='V1', by.x = "V3")
  type <- unlist(lapply(tmp2[,1], function(x){
    unlist(strsplit(x, "[_]"))[1]
  }))
  tmp2 <- tmp2[type == "I1",]
  tmp2$NT <- unlist(lapply(tmp2$V2, function(x){
    str_sub(x, cutpos, cutpos)
  }))
  forecast_132_res_ins1[[file]] <- tmp2
}


ins1_real_vs_pie$indelphi_diff <- unlist(lapply(ins1_real_vs_pie$id, function(id){
  left_nt <- seq_mat[id, 11]
  right_nt <- seq_mat[id, 13]

  tmp <- indelphi_res_ins1[[id]]
  ins1_pct <- rep(0, 4)
  names(ins1_pct) <- c("A", "T", "C", "G")
  for(nt in names(ins1_pct)){
    if(nt %in% tmp$NT){
      ins1_pct[nt] <- tmp$Pct[tmp$NT == nt]
    }
  }
  pie_pred <- unlist(ins1_lr_nt_stat[ins1_lr_nt_stat$leftNT == left_nt & 
                                       ins1_lr_nt_stat$rightNT == right_nt, c(-1,-2)])
  pie_pred <- pie_pred[names(ins1_pct)]
  ins1_pct <- ins1_pct / sum(ins1_pct) * 100
  sum(abs(ins1_pct - pie_pred))
}))


ins1_real_vs_pie$forecast_diff <- unlist(lapply(ins1_real_vs_pie$id, function(id){
  left_nt <- seq_mat[id, 11]
  right_nt <- seq_mat[id, 13]
  
  tmp <- forecast_132_res_ins1[[id]]
  ins1_pct <- rep(0, 4)
  names(ins1_pct) <- c("A", "T", "C", "G")
  for(nt in names(ins1_pct)){
    if(nt %in% tmp$NT){
      ins1_pct[nt] <- tmp$Count[tmp$NT == nt]
    }
  }
  pie_pred <- unlist(ins1_lr_nt_stat[ins1_lr_nt_stat$leftNT == left_nt & 
                                       ins1_lr_nt_stat$rightNT == right_nt, c(-1,-2)])
  pie_pred <- pie_pred[names(ins1_pct)]
  ins1_pct <- ins1_pct / sum(ins1_pct) * 100
  sum(abs(ins1_pct - pie_pred))
}))
plot_data <- melt(ins1_real_vs_pie, id.vars = "id", 
                  measure.vars = c("diff", "indelphi_diff", "forecast_diff"))
median(ins1_real_vs_pie$diff)
median(ins1_real_vs_pie$indelphi_diff)
median(ins1_real_vs_pie$forecast_diff)
pdf("~/Nutstore Files/Tobin/Merged1NT/compare_real_predict_ins1.pdf", width = 6, height = 4)
ggplot(plot_data, aes(x= variable, y = value)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.35),
             shape= 21,
             color = "black", fill = "gray") + 
  stat_summary(fun = "median", colour = "black", size = 4,
               geom = "text", aes(label = paste0("Median=",round(after_stat(y), 2))),
               vjust = -16) +
  stat_compare_means(comparisons = list(c("diff", "indelphi_diff"), 
                                        c("diff", "forecast_diff")), paired = T) + 
  xlab("") + ylab("Real Ins1 vs Predict") +
  theme_bw()
dev.off()

plot_data_add_pair <- data.frame(id = paired_predict$id, 
                                 variable = "pair_diff", 
                                 value = paired_predict$pair_diff)
plot_data_add_pair <- data.frame(rbind(plot_data_add_pair, plot_data))

pdf("~/Nutstore Files/Tobin/Merged1NT/compare_real_predict_ins1_addpair.pdf", width = 7, height = 4)
ggplot(plot_data_add_pair, aes(x= variable, y = value)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.35),
             shape= 21,
             color = "black", fill = "gray") + 
  stat_summary(fun = "median", colour = "black", size = 4,
               geom = "text", aes(label = paste0("Median=",round(after_stat(y), 2))),
               vjust = -16) +
  stat_compare_means(comparisons = list(c("diff", "indelphi_diff"), 
                                        c("diff", "forecast_diff")), paired = T) + 
  stat_compare_means(label.y = c(190, 210, 230), comparisons = list(c("diff", "pair_diff"), 
                                        c("indelphi_diff", "pair_diff"),
                                        c("forecast_diff", "pair_diff")
                                        )) + 
  xlab("") + ylab("Real Ins1 vs Predict") +
  theme_bw()
dev.off()


###试试看只画Cstrand效果如何
plot_data <- ins1_real_vs_pie[seq(2, nrow(ins1_real_vs_pie), 2),]
rownames(plot_data) <- as.character(plot_data$pair)
plot_data <- plot_data[paired_predict$id,]
plot_data <- melt(plot_data, id.vars = "pair", 
                  measure.vars = c("diff", "indelphi_diff", "forecast_diff"))
plot_data_add_pair <- data.frame(pair = paired_predict$id, 
                                 variable = "pair_diff", 
                                 value = paired_predict$pair_diff)
plot_data_add_pair <- data.frame(rbind(plot_data_add_pair, plot_data))
pdf("~/Nutstore Files/Tobin/Merged1NT/compare_real_predict_ins1_addpair_only_could_pred.pdf", width = 7, height = 4)
ggplot(plot_data_add_pair, aes(x= variable, y = value)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.35),
             shape= 21,
             color = "black", fill = "gray") + 
  stat_summary(fun = "median", colour = "black", size = 4,
               geom = "text", aes(label = paste0("Median=",round(after_stat(y), 2))),
               vjust = -16) +
  stat_compare_means(comparisons = list(c("diff", "indelphi_diff"), 
                                        c("diff", "forecast_diff"), 
                                        c("diff", "pair_diff"), 
                                        c("indelphi_diff", "pair_diff"),
                                        c("forecast_diff", "pair_diff")), paired = T) + 
  xlab("") + ylab("Real Ins1 vs Predict") +
  theme_bw()
dev.off()
