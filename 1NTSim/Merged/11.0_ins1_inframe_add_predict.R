#这里用于计算回码
#首先考虑是Ins1的情况 所以基于2.0的结果
library(Biostrings)
library(stringr)


tmp <- lapply(sg_info_sel$id, function(id){
  tables <- total_edit_table_rev[[id]]
  res <- lapply(tables, function(tmp){
    indel <- tmp$n_deleted + tmp$n_inserted
    tmp <- tmp[indel != 0,]
    if(nrow(tmp) == 0){
      return(data.frame(fq = 0, pct = 0))
    }
    tmp$inframe <- tmp$n_inserted - tmp$n_deleted
    total_reads <- sum(tmp[,7])
    tmp <- tmp[tmp$inframe %% 3 == 1, ]
    tmp$inframe <- tmp$inframe - 1

    if(nrow(tmp) == 0){
      return(data.frame(fq = 0, pct = 0))
    }
    return(data.frame(fq = sum(tmp[,7]), pct = sum(tmp[,7]) / total_reads) * 100)
  })
  res <- data.frame(do.call(rbind, res))
  res <- res[res$pct != 0,]
  if(nrow(res) == 0){
    return(data.frame(fq = 0, pct = 0, se = 0))
  }
  res$se <- plotrix::std.error(res$pct)
  res$fq <- sum(res$fq) / length(tables)
  res$pct <- sum(res$pct) / length(tables)
  res[1,]
})
tmp <- data.frame(do.call(rbind, tmp))
sg_info_sel$inframePct <- tmp$pct
sg_info_sel$inframeSe <- tmp$se
plot_data <- hist(sg_info_sel$inframePct[seq(2, nrow(sg_info_sel), 2)], breaks = seq(from=0, to=100, by=3))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
plot_data$color <- "#FAB0B0"
plot_data$color[plot_data$inframePct > 30] <-"#FA1E1E" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 1.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", nrow(sg_info_sel) / 2), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")
#pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_hist.pdf", width = 6, height = 4)
# ggplot(plot_data) + geom_bar(aes(x = x, y = counts), 
#                              color = "black",
#                              stat="identity", fill = plot_data$color) + 
#   geom_text(x = 80, y = 8, label = anno_info, data = data.frame(a=1)) + 
#   xlab("Inframe Frequency") + ylab("Counts") + 
#   scale_y_continuous(breaks = 1 : 10) + 
#   scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90)) + 
#   theme_bw() + 
#   theme(panel.grid = element_blank())
# dev.off()

# openxlsx::write.xlsx(sg_info_sel, file="~/Nutstore Files/Tobin/Merged1NT/Ins1_inframe_ratio_info.xlsx", rowNames=F, colNames=T)




inframe_result <- lapply(sg_info_sel$id, function(id){
  tables <- total_edit_table_rev[[id]]
  res <- lapply(tables, function(tmp){
    indel <- tmp$n_deleted + tmp$n_inserted
    tmp <- tmp[indel != 0,]
    if(nrow(tmp) == 0){
      return(data.frame(indel = 0, fq = 0, pct = 0))
    }
    tmp$inframe <- tmp$n_inserted - tmp$n_deleted
    total_reads <- sum(tmp[,7])
    tmp <- tmp[tmp$inframe %% 3 == 1, ]
    tmp$inframe <- tmp$inframe - 1

    if(nrow(tmp) == 0){
      return(data.frame(indel = 0, fq = 0, pct = 0))
    }
    tmp <- lapply(split(tmp[,7], tmp$inframe), sum)
    tmp <- data.frame(indel = as.numeric(names(tmp)), fq = unlist(tmp))
    tmp$pct <- tmp$fq / total_reads * 100
    tmp <- lapply(split(tmp, tmp$indel), function(x){
      x[,2] <- sum(x[,2])
      x[,3] <- sum(x[,3])
      x[1,]
    })
    data.frame(do.call(rbind,tmp))
  })
  res <- data.frame(do.call(rbind, res))
  res <- res[res$pct != 0,]
  res$aa <- as.numeric(res$indel) / 3
  if(nrow(res) == 0){
    res$se <- NULL
  }else{
    res$se <- 0
  }
  
  
  res <- lapply(split(res, res$aa), function(tmp){
    tmp[,2] <- sum(tmp[,2]) / length(tables)
    se <- tmp[,3]
    if(length(se) < length(tables)){
      se <- unlist(c(se, rep(0, length(tables) - length(se))))
    }
    se <- plotrix::std.error(se)
    tmp[,3] <- sum(tmp[,3]) / length(tables)
    tmp[,5] <- se
    tmp[1,]
  })
  res <- data.frame(do.call(rbind, res))
  res
})
names(inframe_result) <- sg_info_sel$id
inframe_result_processed <- lapply(inframe_result, function(x){
  if(nrow(x) == 0)
    return(x)
  x <- x[order(abs(x$aa)),]
  x2 <- x
  x2$aa[abs(x2$aa) > 5] <- "Others"
  x2$fq[abs(x$aa) > 5] <- sum(x2$fq[abs(x$aa) > 5])
  x2$pct[abs(x$aa) > 5] <- sum(x2$pct[abs(x$aa) > 5])
  if(any(abs(x$aa) > 5)){
    x <- x2[c(which(abs(x$aa) <= 5), which(abs(x$aa) > 5)[1]),] 
  }
  x
})


inframe_result_processed <- lapply(sg_info_sel$id, function(id){
    tmp <- inframe_result_processed[[id]]
    if(nrow(tmp) == 0){
      tmp <- data.frame(indel=0,
                        fq=0,
                        pct=0,
                        aa=0, 
                        se = 0,
                        pct2 = 0)
    }
    tmp$pct2 <- tmp$pct / sum(tmp$pct) * 100
    tmp$id <- id
    tmp
  })
names(inframe_result_processed) <- sg_info_sel$id
inframe_result_processed <- data.frame(do.call(rbind, inframe_result_processed))
# save(inframe_result_processed, file="~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_spec_result.rda")
ids <- gtools::mixedsort(unique(inframe_result_processed$id))
inframe_result_processed$id <- factor(inframe_result_processed$id, levels = rev(ids))
plot_data <- inframe_result_processed
colnames(plot_data)[4] <- "Number of AAChange"


# pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_aa_change.pdf", width = 6, height = 12)
# ggplot(plot_data) + 
  # geom_bar(aes(x = pct, y = id, fill = `Number of AAChange`), 
           # stat = "identity", position = "stack") + 
  # theme_bw() + theme(panel.grid = element_blank())
# dev.off()


#接下来把预测的结果的Inframe也进行统计

load("~/data/project/ear_project/gene_therapy_ll/Result/132_forecast.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/indelphi_132_result.rda")


tmp <- lapply(sg_info_sel$id, function(id){
  tmp <- indelphi_132_res[[id]]
  total_reads <- sum(tmp$fq)
  tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 1, ]
  res <- data.frame(fq = sum(tmp$fq), pct = sum(tmp$fq) / total_reads * 100)
  res <- res[res$pct != 0,]
  res
})
tmp <- data.frame(do.call(rbind, tmp))
sg_info_sel$inframe_indelphi <- tmp$pct

tmp <- lapply(sg_info_sel$id, function(id){
  tmp <- forecast_132_res[[id]]
  total_reads <- sum(tmp$fq)
  tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 1, ]
  res <- data.frame(fq = sum(tmp$fq), pct = sum(tmp$fq) / total_reads * 100)
  res <- res[res$pct != 0,]
  res
})
tmp <- data.frame(do.call(rbind, tmp))
sg_info_sel$inframe_forecast <- tmp$pct


inframe_result_indelphi <- lapply(sg_info_sel$id, function(id){
  tmp <- indelphi_132_res[[id]]
  total_reads <- sum(tmp$fq)
  tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 1, ]
  tmp$inframe <- as.numeric(tmp$indel_size) - 1

  tmp <- lapply(split(tmp$fq, tmp$inframe), sum)
  tmp <- data.frame(indel = as.numeric(names(tmp)), fq = unlist(tmp))
  tmp$pct <- tmp$fq / total_reads * 100
  res <- tmp
  res <- res[res$pct != 0,]
  res$aa <- as.numeric(res$indel) / 3
  res
})
names(inframe_result_indelphi) <- sg_info_sel$id
inframe_result_indelphi <- lapply(inframe_result_indelphi, function(x){
  if(nrow(x) == 0)
    return(x)
  x <- x[order(abs(x$aa)),]
  x2 <- x
  x2$aa[abs(x2$aa) > 5] <- "Others"
  x2$fq[abs(x$aa) > 5] <- sum(x2$fq[abs(x$aa) > 5])
  x2$pct[abs(x$aa) > 5] <- sum(x2$pct[abs(x$aa) > 5])
  if(any(abs(x$aa) > 5)){
    x <- x2[c(which(abs(x$aa) <= 5), which(abs(x$aa) > 5)[1]),] 
  }
  x
})


inframe_result_indelphi <- lapply(sg_info_sel$id, function(id){
  tmp <- inframe_result_indelphi[[id]]
  if(nrow(tmp) == 0){
    tmp <- data.frame(indel=0,
                      fq=0,
                      pct=0,
                      aa=0, 
                      pct2 = 0)
  }
  tmp$pct2 <- tmp$pct / sum(tmp$pct) * 100
  tmp$id <- id
  tmp
})
names(inframe_result_indelphi) <- sg_info_sel$id
inframe_result_indelphi <- data.frame(do.call(rbind, inframe_result_indelphi))



inframe_result_forecast <- lapply(sg_info_sel$id, function(id){
  tmp <- forecast_132_res[[id]]
  total_reads <- sum(tmp$fq)
  tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 1, ]
  tmp$inframe <- as.numeric(tmp$indel_size) - 1
  
  tmp <- lapply(split(tmp$fq, tmp$inframe), sum)
  tmp <- data.frame(indel = as.numeric(names(tmp)), fq = unlist(tmp))
  tmp$pct <- tmp$fq / total_reads * 100
  res <- tmp
  res <- res[res$pct != 0,]
  res$aa <- as.numeric(res$indel) / 3
  res
})
names(inframe_result_forecast) <- sg_info_sel$id
inframe_result_forecast <- lapply(inframe_result_forecast, function(x){
  if(nrow(x) == 0)
    return(x)
  x <- x[order(abs(x$aa)),]
  x2 <- x
  x2$aa[abs(x2$aa) > 5] <- "Others"
  x2$fq[abs(x$aa) > 5] <- sum(x2$fq[abs(x$aa) > 5])
  x2$pct[abs(x$aa) > 5] <- sum(x2$pct[abs(x$aa) > 5])
  if(any(abs(x$aa) > 5)){
    x <- x2[c(which(abs(x$aa) <= 5), which(abs(x$aa) > 5)[1]),] 
  }
  x
})


inframe_result_forecast <- lapply(sg_info_sel$id, function(id){
  tmp <- inframe_result_forecast[[id]]
  if(nrow(tmp) == 0){
    tmp <- data.frame(indel=0,
                      fq=0,
                      pct=0,
                      aa=0, 
                      pct2 = 0)
  }
  tmp$pct2 <- tmp$pct / sum(tmp$pct) * 100
  tmp$id <- id
  tmp
})
names(inframe_result_forecast) <- sg_info_sel$id
inframe_result_forecast <- data.frame(do.call(rbind, inframe_result_forecast))
# sg_info_sel <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Merged1NT/Ins1_inframe_ratio_info_add_predict.xlsx")
sg_info_sel_noBulge <- sg_info_sel[seq(2, nrow(sg_info_sel), 2),]
sg_info_sel_withBulge <- sg_info_sel[seq(1, nrow(sg_info_sel), 2),]
sg_info_sel_output <- sg_info_sel_noBulge
sg_info_sel_output$inframe_reference <- sg_info_sel_withBulge$inframePct
sg_info_sel_output$inframeSe_reference <- sg_info_sel_withBulge$inframeSe
sg_info_sel_output <- sg_info_sel_output[,-4]
colnames(sg_info_sel_output)[c(5,6)] <- c("inframe_target", "inframeSe_target")
openxlsx::write.xlsx(sg_info_sel_output, file="~/Nutstore Files/Tobin/Merged1NT/Ins1_inframe_ratio_info_add_predict_target.xlsx", rowNames=F, colNames=T)



plot_data <- data.frame(pair = sg_info_sel$pair[seq(2, nrow(sg_info_sel), 2)],
                        bulgePct = sg_info_sel$inframePct[seq(1, nrow(sg_info_sel), 2)],
                        noBulgePct = sg_info_sel$inframePct[seq(2, nrow(sg_info_sel), 2)], 
                        indelPhiPct = sg_info_sel$inframe_indelphi[seq(2, nrow(sg_info_sel), 2)],
                        ForecastPct = sg_info_sel$inframe_forecast[seq(2, nrow(sg_info_sel), 2)])
plot_data <- melt(plot_data, id.vars = "pair", measure.vars = c("bulgePct", "noBulgePct", "indelPhiPct", "ForecastPct"))
plot_data$variable <- factor(as.character(plot_data$variable), levels = c("noBulgePct", "bulgePct", "indelPhiPct", "ForecastPct"))
pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_violine.pdf", width = 6, height = 4)
ggplot(plot_data, aes(x= variable, y = value)) + 
  geom_violin(aes(fill = variable)) + 
  geom_point(shape= 21,
             color = "black", fill = "gray") + 
  geom_line(aes(group = pair), color = '#DDDDDD') + 
  stat_compare_means(comparisons = list(c("bulgePct", "noBulgePct"), 
                                        c("noBulgePct", "indelPhiPct"), 
                                        c("noBulgePct", "ForecastPct")), paired = T) + 
  xlab("") + ylab("Inframe% in Indel") +
  theme_bw()
dev.off()
pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_violine_v2.pdf", width = 6, height = 4)
ggplot(plot_data, aes(x= variable, y = value)) + 
  geom_violin(aes(fill = variable)) + 
  geom_beeswarm(shape= 21, 
                color = "black", fill = "gray") +
  stat_compare_means(comparisons = list(c("bulgePct", "noBulgePct"), 
                                        c("noBulgePct", "indelPhiPct"), 
                                        c("noBulgePct", "ForecastPct")), paired = T) + 
  xlab("") + ylab("Inframe% in Indel") +
  theme_bw()
dev.off() 

bulge_one <- sg_info_sel$id[seq(1, nrow(sg_info_sel), 2)]
nobulge_one <- sg_info_sel$id[seq(2, nrow(sg_info_sel), 2)]
plot_data_0aa <- inframe_result_processed[inframe_result_processed$aa == "0",]
rownames(plot_data_0aa) <- plot_data_0aa$id
tmp2 <-inframe_result_indelphi[inframe_result_indelphi$aa == "0",]
rownames(tmp2) <- tmp2$id
tmp2 <- tmp2[nobulge_one,]
tmp3 <- inframe_result_forecast[inframe_result_forecast$aa == "0",]
rownames(tmp3) <- tmp3$id
tmp3 <- tmp3[nobulge_one,]
tmp2$type <- "IndelPhi"
tmp3$type <- "forecast"
tmp <- plot_data_0aa[bulge_one,]
tmp1 <- plot_data_0aa[nobulge_one,]
tmp$type <- "Bulge"
tmp1$type <- "noBulge"
tmp$id <- tmp1$id
tmp <- tmp[,-5]
tmp1 <- tmp1[,-5]
plot_data_0aa <- data.frame(rbind(tmp, tmp1, tmp2, tmp3))
plot_data_0aa$type <- factor(plot_data_0aa$type, levels = c("noBulge", "Bulge", "IndelPhi", "forecast"), 
                             labels = c("target", "reference", "InDelphi", "FOREcasT"))


pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_0AA_violine_of_inframe_all.pdf", width = 6, height = 4)
ggplot(plot_data_0aa, aes(x= type, y = pct2)) + 
  geom_violin(aes(fill = type)) + 
  geom_point(shape= 21,
             color = "black", fill = "gray") + 
  geom_line(aes(group = id), color = '#DDDDDD', linewidth=0.1) + 
  stat_compare_means(comparisons = list(c("target", "reference"), 
                                        c("target", "InDelphi"), 
                                        c("target", "FOREcasT")), paired = T) + 
  xlab("") + ylab("0AA/Inframe") +
  theme_bw()
dev.off()

inframe30_noBulge <- sg_info_sel[sg_info_sel$id %in% nobulge_one,]
inframe30_noBulge <- inframe30_noBulge[inframe30_noBulge$inframePct > 30,]
inframe30_withBulge <- sg_info_sel[seq(1, nrow(sg_info_sel), 2),]
inframe30_withBulge <- inframe30_withBulge[inframe30_withBulge$inframePct > 30,]
plot_data_0aa$reference30 <- plot_data_0aa$id %in% sg_info_sel$id[which(sg_info_sel$id %in% inframe30_withBulge$id) + 1]
plot_data_0aa$target30 <- plot_data_0aa$id%in% inframe30_noBulge$id
plot_data_0aa$large30 <- "No"
plot_data_0aa$large30[plot_data_0aa$target30] <- "Target>30"
plot_data_0aa$large30[plot_data_0aa$target30 & plot_data_0aa$reference30] <- "Both>30"
plot_data_0aa$large30[!plot_data_0aa$target30 & plot_data_0aa$reference30] <- "Reference>30"
plot_data_0aa$large30 <- factor(plot_data_0aa$large30, levels = c("Both>30", "Target>30", "Reference>30", "No"))
plot_data_0aa_sel <- plot_data_0aa[plot_data_0aa$large30 != 'No',]


plot_data_0aa_sel2 <- plot_data_0aa_sel[plot_data_0aa_sel$large30 == "Both>30",]
library(rstatix)
stat <- plot_data_0aa_sel2 %>%
  wilcox_test(pct2 ~ type, paired = T)
stat <- add_xy_position(stat, x = "type")
stat <- stat[stat$group1 == "target",]
stat$y.position <- c(106, 118, 130)

tmp <- plot_data_0aa
tmp$type <- as.numeric(tmp$type)
tmp <- beeswarm::beeswarm(pct2 ~ type, tmp, corral = "none", corralWidth = 0.2, spacing = 0.4)
plot_data_0aa$beex <- tmp$x
pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_0AA_violine_of_inframe_all2.pdf", width = 8, height = 5)
ggplot(plot_data_0aa, aes(x= type, y = pct2)) + 
  geom_violin(aes(fill = type), alpha = 0.5) + 
  ggnewscale::new_scale("fill") +
  geom_point(aes(x = beex, fill=large30),shape= 21,
             color = "black") + 
  geom_line(aes(x = beex, group = id, color=large30), 
            linewidth=0.1) + 
  scale_fill_manual(values = c("#FF4500", "#FF4500", "#FF4500", "white")) + 
  scale_color_manual(values = c("#FF4500", "#F08080", "#F08080", "#DDDDDD"))+ 
  stat_compare_means(comparisons = list(c("target", "reference"), 
                                        c("target", "InDelphi"), 
                                        c("target", "FOREcasT")), paired = T) +
  stat_pvalue_manual(stat,  label = "Pval Both>30 : {p}", tip.length = 0.02) +
  xlab("") + ylab("0AA/Inframe") +
  theme_bw()
dev.off()
export_info <- plot_data_0aa[,c("id", "pct2", "large30", "type")]
export_info_tmp <- export_info[export_info$type == "target",]
colnames(export_info_tmp)[2] <- "target"
export_info_tmp$reference <- export_info$pct2[export_info$type == "reference"]
export_info_tmp$InDelphi <- export_info$pct2[export_info$type == "InDelphi"]
export_info_tmp$FOREcasT <- export_info$pct2[export_info$type == "FOREcasT"]
export_info_tmp <- export_info_tmp[,c(1,3,2,5,6,7)]
openxlsx::write.xlsx(export_info_tmp, file="~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_0aa_plot_data.xlsx", rowNames=F, colNames=T)


pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_0AA_violine_of_inframe_large_30.pdf", width = 6, height = 4)
ggplot(plot_data_0aa_sel, aes(x= type, y = pct2)) + 
  geom_violin(aes(fill = type)) + 
  ggnewscale::new_scale("fill") +
  geom_point(aes(fill=large30),shape= 21,
             color = "black") + 
  geom_line(aes(group = id, color=large30), linewidth=0.1) + 
  stat_compare_means(comparisons = list(c("target", "reference"), 
                                        c("target", "InDelphi"), 
                                        c("target", "FOREcasT")), paired = T) + 
  xlab("") + ylab("0AA/Inframe") +
  theme_bw()
dev.off()
plot_data_0aa_sel2 <- plot_data_0aa_sel[plot_data_0aa_sel$large30 == "Both>30",]
pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_0AA_violine_of_inframe_Both_30.pdf", width = 6, height = 4)
ggplot(plot_data_0aa_sel2, aes(x= type, y = pct2)) + 
  geom_violin(aes(fill = type)) + 
  geom_point(shape= 21,color = "black", fill = "gray") + 
  geom_line(aes(group = id), colour = '#DDDDDD', linewidth = 0.1) + 
  stat_compare_means(comparisons = list(c("target", "reference"), 
                                        c("target", "InDelphi"), 
                                        c("target", "FOREcasT")), paired = T) + 
  xlab("") + ylab("0AA/Inframe") +
  theme_bw()
dev.off()



pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_0AA_violine_of_indel.pdf", width = 6, height = 4)
ggplot(plot_data_0aa, aes(x= type, y = pct)) + 
  geom_violin(aes(fill = type)) + 
  geom_point(shape= 21,
             color = "black", fill = "gray") + 
  geom_line(aes(group = id), color = '#DDDDDD') + 
  stat_compare_means(comparisons = list(c("target", "reference"), 
                                        c("target", "InDelphi"), 
                                        c("target", "FOREcasT")), paired = T) + 
  xlab("") + ylab("0AA/Inframe") +
  theme_bw()
dev.off()
pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_0AA_violine_v2.pdf", width = 6, height = 4)
ggplot(plot_data_0aa, aes(x= type, y = pct)) + 
  geom_violin(aes(fill = type)) + 
  geom_beeswarm(shape= 21, 
             color = "black", fill = "gray") + 
  stat_compare_means(comparisons = list(c("noBulge", "Bulge"), 
                                        c("noBulge", "IndelPhi"), 
                                        c("noBulge", "forecast")), paired = T) + 
  xlab("") + ylab("0AA in Indel") +
  theme_bw()
dev.off()
