print(load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda"))
print(load("~/data/project/ear_project/gene_therapy_ll/Result/merged_indel_table_first_second.Rda"))
print(load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first_second.Rda"))
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
sgCmp_sel <- sgCmp_sel[sgCmp_sel$id %in% names(sgRNA_pair_remain),]

library(Biostrings)
library(stringr)


tmp <- lapply(sgCmp_sel$id2, function(id){
  bc_strand <- "Insert"
  tables <- total_edit_table_rev[[id]]
  res <- lapply(tables, function(tmp){
    indel <- tmp$n_deleted + tmp$n_inserted
    tmp <- tmp[indel != 0,]
    if(nrow(tmp) == 0){
      return(data.frame(fq = 0, pct = 0))
    }
    tmp$inframe <- tmp$n_inserted - tmp$n_deleted
    total_reads <- sum(tmp[,7])
    if(bc_strand == "Delete"){
      tmp <- tmp[tmp$inframe %% 3 == 1, ]
      tmp$inframe <- tmp$inframe - 1
    } else {
      tmp <- tmp[tmp$inframe %% 3 == 2, ]
      tmp$inframe <- tmp$inframe + 1
    }
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
sgCmp_sel$inframePct <- tmp$pct
sgCmp_sel$inframeSe <- tmp$se


inframe_result <- lapply(sgCmp_sel$id2, function(id){
  bc_strand <- "Insert"
  tables <- total_edit_table_rev[[id]]
  res <- lapply(tables, function(tmp){
    indel <- tmp$n_deleted + tmp$n_inserted
    tmp <- tmp[indel != 0,]
    if(nrow(tmp) == 0){
      return(data.frame(indel = 0, fq = 0, pct = 0))
    }
    tmp$inframe <- tmp$n_inserted - tmp$n_deleted
    total_reads <- sum(tmp[,7])
    if(bc_strand == "Delete"){
      tmp <- tmp[tmp$inframe %% 3 == 1, ]
      tmp$inframe <- tmp$inframe - 1
    } else {
      tmp <- tmp[tmp$inframe %% 3 == 2, ]
      tmp$inframe <- tmp$inframe + 1
    }
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
names(inframe_result) <- sgCmp_sel$id2
inframe_result_processed <- lapply(inframe_result, function(x){
  if(nrow(x) == 0)
    return(x)
  x <- x[order(abs(x$aa)),]
  x2 <- x
  x2$aa[abs(x2$aa) > 2] <- "Others"
  x2$fq[abs(x$aa) > 2] <- sum(x2$fq[abs(x$aa) > 2])
  x2$pct[abs(x$aa) > 2] <- sum(x2$pct[abs(x$aa) > 2])
  if(any(abs(x$aa) > 2)){
    x <- x2[c(which(abs(x$aa) <= 2), which(abs(x$aa) > 2)[1]),] 
  }
  x
})


inframe_result_processed <- lapply(sgCmp_sel$id2, function(id){
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
names(inframe_result_processed) <- sgCmp_sel$id2
inframe_result_processed <- data.frame(do.call(rbind, inframe_result_processed))
inframe_0aa <- inframe_result_processed[inframe_result_processed$aa == 0,]
rownames(inframe_0aa) <- inframe_0aa$id
inframe_0aa <- inframe_0aa[sgCmp_sel$id2,]
sgCmp_sel[,"0AA"] <- inframe_0aa$pct
sgCmp_sel[,"0AASe"] <- inframe_0aa$se

orderlist <- lapply(split(sgCmp_sel$inframePct, sgCmp_sel$id), function(x){
  x <- unlist(x)
  mean(x) - abs(x[1] - x[2]) * 2
})
orderlist <- data.frame(id = names(orderlist), pct = unlist(orderlist))
orderlist <- orderlist[order(orderlist$pct, decreasing = T),]
sgCmp_sel$idorder <- factor(sgCmp_sel$id, levels = orderlist$id)
sgCmp_sel <- sgCmp_sel[order(sgCmp_sel$idorder),]
sgCmp_bstrand <- sgCmp_sel[sgCmp_sel$isInsert == "Insert",]
sgCmp_bstrand$id2 <- factor(sgCmp_bstrand$id2, levels = sgCmp_bstrand$id2)
sgCmp_cstrand <- sgCmp_sel[sgCmp_sel$isInsert != "Insert",]
sgCmp_cstrand$id2 <- factor(sgCmp_cstrand$id2, levels = sgCmp_cstrand$id2)

sgCmp_bstrand$`0AA2` <- sgCmp_bstrand$`0AA` / sgCmp_bstrand$inframePct * 100
sgCmp_cstrand$`0AA2` <- sgCmp_cstrand$`0AA` / sgCmp_cstrand$inframePct * 100
p1 <- ggplot(sgCmp_bstrand) + 
  geom_point(stat = "identity",aes(x = inframePct,y = `0AA2`)) + 
  geom_abline(slope = 1) + 
  labs(x = "", y = "BStrand 0AA% in Inframe") + 
  scale_x_continuous(limits = c(0, 100))  +
  scale_y_continuous(limits = c(0, 100)) + 
  theme_classic2() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
p2 <- ggplot(sgCmp_cstrand) + 
  geom_point(stat = "identity",aes(x = inframePct,y = `0AA2`)) + 
  geom_abline(slope = -1) + 
  labs(x = "Inframe%", y = "CStrand 0AA% in Inframe") + 
  scale_x_continuous(position = "top", limits = c(0, 100))  +
  scale_y_reverse(limits = c(100, 0)) +
  theme_classic2() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())




pdf("~/Nutstore Files/Tobin/Merged1NT/132_inframe_0aa_pointPlot_for_BStrand_v2.pdf", width = 5, height = 6)
ggpubr::ggarrange(p1, p2, nrow = 2, common.legend = T)
dev.off()


plot_data2 <- sgCmp_sel
plot_data2$`0AA2` <- plot_data2$`0AA` / plot_data2$inframePct * 100
plot_data2$`0AA2`[plot_data2$isInsert != 'Insert'] <- -plot_data2$`0AA2`[plot_data2$isInsert != 'Insert']
plot_data2$type <- factor(plot_data2$isInsert, levels = c("Insert", "Delete"), labels = c("BStrand", "CStrand"))

plot_data2$inframePct2 <- plot_data2$inframePct
plot_data2$inframePct2[plot_data2$isInsert != 'Insert'] <- -plot_data2$inframePct2[plot_data2$isInsert != 'Insert']

slope_cal <-lapply(split(plot_data2[,c("inframePct", "0AA2")], plot_data2$id), function(x){
  suby <- abs(x[1,1]) + abs(x[2,1])
  subx <- abs(x[1,2]) - abs(x[2,2])
  if(subx == 0){
    return(99999)
  }
  return(abs(suby / subx))
})
slope_cal <- data.frame(id = names(slope_cal), slope = unlist(slope_cal))

plot_data2 <- merge(plot_data2, slope_cal, by="id")
plot_data2$large60 <- plot_data2$slope > sqrt(3)


pdf("~/Nutstore Files/Tobin/Merged1NT/132_inframe_0aa_pointPlot_for_BStrand_alpha_v5.pdf", width = 6, height = 4)
ggplot() + 
  geom_line(data = plot_data2 ,aes(x = abs(`0AA2`), y = inframePct2, 
                                   group = id, alpha = large60), linewidth = 0.1) + 
  geom_point(data = plot_data2 ,aes(x = abs(`0AA2`), y = inframePct2, color = type, 
                                    alpha = large60)) + 
  geom_abline(slope = c(1, -1), intercept = 0) + 
  geom_hline(yintercept = 0, linewidth = 1) + 
  labs(y = "Inframe%", x = "0AA% in Inframe") + 
  scale_x_continuous(limits = c(0, 100), expand = c(0,0))  +
  scale_y_continuous(limits = c(-100, 100)) + 
  scale_color_manual(values = c("#FF6A6A", "#FF6A6A")) + 
  theme_classic2() + theme(axis.text.x = element_text(angle = 90, 
                                                      hjust = 0, 
                                                      vjust = 0.5), 
                           panel.grid  = element_blank())
dev.off()


plot_data2$alpha <- stats::approxfun(c(sqrt(3), 1/ tan(15 / 180 * pi),15), c(0.1,0.7,1))(plot_data2$slope)
plot_data2$alpha[plot_data2$slope < sqrt(3)] <- 0.1
plot_data2$alpha[plot_data2$slope > 15] <- 1
pdf("~/Nutstore Files/Tobin/Merged1NT/132_inframe_0aa_pointPlot_for_BStrand_gradient_alpha_v8.pdf", width = 5, height = 4)
ggplot() + 
  geom_line(data = plot_data2 ,aes(x = abs(`0AA2`), y = inframePct2, 
                                   group = id), linewidth = 0.1, alpha = plot_data2$alpha) + 
  geom_point(data = plot_data2 ,aes(x = abs(`0AA2`), y = inframePct2, 
                                    color = type), alpha = plot_data2$alpha) +
  
  # geom_abline(slope = c(1, -1), intercept = 0) + 
  geom_hline(yintercept = 0, linewidth = 1) + 
  labs(y = "Inframe%", x = "0AA% in Inframe") + 
  scale_x_continuous(limits = c(0, 100), expand = c(0,0))  +
  scale_color_manual(values = c("#0A5028", "#ab1d25")) + 
  theme_classic2() + theme(axis.text.x = element_text(angle = 90, 
                                                      hjust = 0, 
                                                      vjust = 0.5), 
                           panel.grid  = element_blank()) + 
  coord_flip() + 
  scale_y_reverse(limits = c(100, -100))
dev.off()

pdf("~/Nutstore Files/Tobin/Merged1NT/132_inframe_0aa_BStrand_legend_v3.pdf", width = 3, height = 7)

legend_plot <- data.frame(id = 1 : 132, 
                          color = rev(alpha("#0A5028", sort(plot_data2$alpha[plot_data2$isInsert == "Insert"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data2$slope[plot_data2$isInsert == "Insert"] < tan(x / 180 * pi))
})) + 0.5
p1 <- ggplot(legend_plot) + 
  geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 132) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.3) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 132) + 0.5), 
            aes(x = x, y = 0.5, label = round(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.3) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.3) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 132.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()


legend_plot <- data.frame(id = 1 : 132, 
                          color = rev(alpha("#ab1d25", sort(plot_data2$alpha[plot_data2$isInsert != "Insert"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data2$slope[plot_data2$isInsert != "Insert"] < tan(x / 180 * pi))
})) + 0.5

p2 <- ggplot(legend_plot) + 
  geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 132) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.3) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 132) + 0.5), 
            aes(x = x, y = 0.5, label = round(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.3) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.3) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 132.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()
ggpubr::ggarrange(p1, p2, nrow = 1)
dev.off()



p1 <- ggplot(sgCmp_bstrand) + 
  geom_bar(stat = "identity",aes(x = id2,y = inframePct, fill = "Inframe")) + 
  geom_bar(stat = "identity",aes(x = id2,y = `0AA`, fill = "0AA"), width = 0.5) + 
  geom_errorbar(aes(x = id2, ymin = inframePct - inframeSe, ymax = inframePct + inframeSe), 
                color = "black") + 
  geom_errorbar(aes(x = id2, ymin = `0AA` - `0AASe`, ymax =  `0AA` + `0AASe`), 
                color = "black") + 
  labs(x = "", y = "BStrand % in Indel Reads", fill = "Type") + 
  scale_fill_manual(values = c("Inframe" = "#06a0c9", "0AA" = "#E1FFFF")) + 
  scale_x_discrete(position = "top")  +
  scale_y_continuous(limits = c(0, 100)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
p2 <- ggplot(sgCmp_cstrand) + 
  geom_bar(stat = "identity",aes(x = id2,y = inframePct, fill = "Inframe")) + 
  geom_bar(stat = "identity",aes(x = id2,y = `0AA`, fill = "0AA"), width = 0.5) + 
  geom_errorbar(aes(x = id2, ymin = inframePct - inframeSe, ymax = inframePct + inframeSe), 
                color = "black") + 
  geom_errorbar(aes(x = id2, ymin = `0AA` - `0AASe`, ymax =  `0AA` + `0AASe`), 
                color = "black") +
  labs(x = "", y = "CStrand % in Indel Reads", fill = "Type") + 
  scale_fill_manual(values = c("Inframe" = "#06a0c9", "0AA" = "#E1FFFF")) + 
  scale_y_reverse(limits = c(100, 0)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
pdf("~/Nutstore Files/Tobin/Merged1NT/132_inframe_0aa_barplot_for_BStrand_Del1.pdf", width = 20, height = 10)
ggpubr::ggarrange(p1, p2, nrow = 2, common.legend = T)
dev.off()
vcd::assocstats()
#fisher test
test_result <- lapply(1 : nrow(sgCmp_bstrand), function(i){
  mat <- matrix(c(sgCmp_bstrand$inframePct[i], sgCmp_bstrand$`0AA`[i], 
                  sgCmp_cstrand$inframePct[i], sgCmp_cstrand$`0AA`[i]
                  ), nrow = 2)
  chisq.test(mat)$p.value
})

test_result <- unlist(test_result)
sgCmp_sel <- sgCmp_sel[order(sgCmp_sel$id, sgCmp_sel$isInsert),]
sd_b <- sd(sgCmp_sel$inframePct[sgCmp_sel$isInsert == "Insert"])
sd_c <- sd(sgCmp_sel$inframePct[sgCmp_sel$isInsert != "Insert"])
sd_0aa_b <- sd(sgCmp_sel$`0AA`[sgCmp_sel$isInsert == "Insert"])
sd_0aa_c <- sd(sgCmp_sel$`0AA`[sgCmp_sel$isInsert != "Insert"])
orderlist <- lapply(split(sgCmp_sel[,c("inframePct", "0AA2")], sgCmp_sel$id), function(x){
  p1 <- pnorm(x[2,1],mean = x[1,1], sd = sd_c)
  p2 <- pnorm(x[2,2],mean = x[1,2], sd = sd_0aa_c)
  if(p1 < 0.5){
    p1 <- 2 * p1
  } else {
    p1 <- 2 * (1 - p1)
  }
  if(p2 < 0.5){
    p2 <- 2 * p2
  } else {
    p2 <- 2 * (1 - p2)
  }
  
  
  data.frame(pval = p1 * p2, 
             mean = mean(unlist(x[,1])))
})
tmp <- names(orderlist)
orderlist <- data.frame(do.call(rbind, orderlist))
orderlist$id <- tmp
orderlist <- orderlist[order(orderlist$pval, orderlist$mean, decreasing = T),]
sgCmp_sel$idorder <- factor(sgCmp_sel$id, levels = orderlist$id)
sgCmp_sel <- sgCmp_sel[order(sgCmp_sel$idorder),]
sgCmp_bstrand <- sgCmp_sel[sgCmp_sel$isInsert == "Insert",]
sgCmp_bstrand$id2 <- factor(sgCmp_bstrand$id2, levels = sgCmp_bstrand$id2)
sgCmp_cstrand <- sgCmp_sel[sgCmp_sel$isInsert != "Insert",]
sgCmp_cstrand$id2 <- factor(sgCmp_cstrand$id2, levels = sgCmp_cstrand$id2)


p1 <- ggplot(sgCmp_bstrand) + 
  geom_bar(stat = "identity",aes(x = id2,y = inframePct, fill = "Inframe")) + 
  geom_bar(stat = "identity",aes(x = id2,y = `0AA`, fill = "0AA"), width = 0.5) + 
  geom_errorbar(aes(x = id2, ymin = inframePct - inframeSe, ymax = inframePct + inframeSe), 
                color = "black") + 
  geom_errorbar(aes(x = id2, ymin = `0AA` - `0AASe`, ymax =  `0AA` + `0AASe`), 
                color = "black") + 
  labs(x = "", y = "BStrand % in Indel Reads", fill = "Type") + 
  scale_fill_manual(values = c("Inframe" = "#06a0c9", "0AA" = "#E1FFFF")) + 
  scale_x_discrete(position = "top")  +
  scale_y_continuous(limits = c(0, 100)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
p2 <- ggplot(sgCmp_cstrand) + 
  geom_bar(stat = "identity",aes(x = id2,y = inframePct, fill = "Inframe")) + 
  geom_bar(stat = "identity",aes(x = id2,y = `0AA`, fill = "0AA"), width = 0.5) + 
  geom_errorbar(aes(x = id2, ymin = inframePct - inframeSe, ymax = inframePct + inframeSe), 
                color = "black") + 
  geom_errorbar(aes(x = id2, ymin = `0AA` - `0AASe`, ymax =  `0AA` + `0AASe`), 
                color = "black") +
  labs(x = "", y = "CStrand % in Indel Reads", fill = "Type") + 
  scale_fill_manual(values = c("Inframe" = "#06a0c9", "0AA" = "#E1FFFF")) + 
  scale_y_reverse(limits = c(100, 0)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
pdf("~/Nutstore Files/Tobin/Merged1NT/132_inframe_0aa_barplot_for_BStrand_Del1_reorder.pdf", width = 20, height = 10)
ggpubr::ggarrange(p1, p2, nrow = 2, common.legend = T)
dev.off()
openxlsx::write.xlsx(orderlist, file="~/Nutstore Files/Tobin/Merged1NT/132_inframe_0aa_barplot_for_BStrand_Del1_reorder.xlsx", rowNames=F, colNames=T)


sgCmp_sel$`0AA2` <- sgCmp_sel$`0AA` / sgCmp_sel$inframePct * 100
sgCmp_sel <- sgCmp_sel[order(sgCmp_sel$id, sgCmp_sel$isInsert),]
sd_b <- sd(sgCmp_sel$inframePct[sgCmp_sel$isInsert == "Insert"])
sd_c <- sd(sgCmp_sel$inframePct[sgCmp_sel$isInsert != "Insert"])
sd_0aa_b <- sd(sgCmp_sel$`0AA`[sgCmp_sel$isInsert == "Insert"])
sd_0aa_c <- sd(sgCmp_sel$`0AA`[sgCmp_sel$isInsert != "Insert"])
orderlist <- lapply(split(sgCmp_sel[,c("inframePct", "0AA2")], sgCmp_sel$id), function(x){
  p1 <- pnorm(x[2,1],mean = x[1,1], sd = sd_c)
  p2 <- pnorm(x[2,2],mean = x[1,2], sd = sd_0aa_c)
  if(p1 < 0.5){
    p1 <- 2 * p1
  } else {
    p1 <- 2 * (1 - p1)
  }
  if(p2 < 0.5){
    p2 <- 2 * p2
  } else {
    p2 <- 2 * (1 - p2)
  }
  
  
  data.frame(pval = p1 * p2, 
             mean = mean(unlist(x[,1])))
})
tmp <- names(orderlist)
orderlist <- data.frame(do.call(rbind, orderlist))
orderlist$id <- tmp
orderlist <- orderlist[order(orderlist$pval, orderlist$mean, decreasing = T),]
sgCmp_sel$idorder <- factor(sgCmp_sel$id, levels = orderlist$id)
sgCmp_sel <- sgCmp_sel[order(sgCmp_sel$idorder),]
sgCmp_bstrand <- sgCmp_sel[sgCmp_sel$isInsert == "Insert",]
sgCmp_bstrand$id2 <- factor(sgCmp_bstrand$id2, levels = sgCmp_bstrand$id2)
sgCmp_cstrand <- sgCmp_sel[sgCmp_sel$isInsert != "Insert",]
sgCmp_cstrand$id2 <- factor(sgCmp_cstrand$id2, levels = sgCmp_cstrand$id2)


p1 <- ggplot(sgCmp_bstrand) + 
  # geom_bar(stat = "identity",aes(x = id2,y = inframePct, fill = "Inframe")) + 
  geom_bar(stat = "identity",aes(x = id2,y = `0AA2`, fill = "0AA")) + 
  # geom_errorbar(aes(x = id2, ymin = inframePct - inframeSe, ymax = inframePct + inframeSe), 
  #               color = "black") + 
  geom_errorbar(aes(x = id2, ymin = `0AA2` - `0AASe`, ymax =  `0AA2` + `0AASe`), 
                color = "black") + 
  labs(x = "", y = "BStrand 0AA/Inframe", fill = "Type") + 
  scale_fill_manual(values = c("Inframe" = "#06a0c9", "0AA" = "#06a0c9")) + 
  scale_x_discrete(position = "top")  +
  scale_y_continuous(limits = c(0, 100)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
p2 <- ggplot(sgCmp_cstrand) + 
  # geom_bar(stat = "identity",aes(x = id2,y = inframePct, fill = "Inframe")) + 
  geom_bar(stat = "identity",aes(x = id2,y = `0AA2`, fill = "0AA")) + 
  # geom_errorbar(aes(x = id2, ymin = inframePct - inframeSe, ymax = inframePct + inframeSe), 
  #               color = "black") + 
  geom_errorbar(aes(x = id2, ymin = `0AA2` - `0AASe`, ymax =  `0AA2` + `0AASe`), 
                color = "black") +
  labs(x = "", y = "CStrand 0AA/Inframe", fill = "Type") + 
  scale_fill_manual(values = c("Inframe" = "#06a0c9", "0AA" = "#06a0c9")) + 
  scale_y_reverse(limits = c(100, 0)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                hjust = 0, 
                                                vjust = 0.5), 
                     panel.grid  = element_blank())
pdf("~/Nutstore Files/Tobin/Merged1NT/132_inframe_0aa_barplot_for_BStrand_Del1_reorder_v3.pdf", width = 20, height = 10)
ggpubr::ggarrange(p1, p2, nrow = 2, common.legend = T)
dev.off()



plot_data <- hist(sgCmp_bstrand$inframePct, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
plot_data$color <- "#0A502866"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028"
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p1 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 80, y = 8, label = anno_info, data = data.frame(a=1)) +
  xlab("Inframe Frequency") + ylab("Counts") +
  scale_y_continuous(breaks = 1 : 26, limits = c(0, 26)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100)) +
  theme_bw() +
  theme(plot.margin =  unit(c(0,0,0,0.2), "cm"), 
        panel.grid = element_blank())
# dev.off()

plot_data <- hist(sgCmp_cstrand$inframePct, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
plot_data$color <- "#0A502866"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028"
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p2 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 80, y = -5, label = anno_info, data = data.frame(a=1)) +
  scale_y_reverse(expand = c(0,0), breaks = 1 : 26, limits = c(26, 0)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), position = "top") +
  ylab("Counts") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        plot.margin =  unit(c(-1,0,1,0.2), "cm"))

pdf("~/Nutstore Files/Tobin/Merged1NT/132_inframe_histgram_for_BStrand_Del1_width5.pdf", width = 6, height = 6)
ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T)
dev.off()




plot_data <- hist(sgCmp_bstrand$`0AA2`, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
plot_data$color <- "#0A502866"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028"
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p1 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 80, y = 8, label = anno_info, data = data.frame(a=1)) +
  xlab("Inframe Frequency") + ylab("Counts") +
  scale_y_continuous(breaks = 1 : 22, limits = c(0, 22)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100)) +
  theme_bw() +
  theme(plot.margin =  unit(c(0,0,0,0.2), "cm"), 
        panel.grid = element_blank())
# dev.off()

plot_data <- hist(sgCmp_cstrand$`0AA2`, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
plot_data$color <- "#0A502866"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028"
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p2 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 80, y = -5, label = anno_info, data = data.frame(a=1)) +
  scale_y_reverse(expand = c(0,0), breaks = 1 : 22, limits = c(22, 0)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), position = "top") +
  ylab("Counts") + xlab("") + 
theme_bw() + 
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        plot.margin =  unit(c(-1,0,1,0.2), "cm"))

pdf("~/Nutstore Files/Tobin/Merged1NT/132_0aa_inframe_histgram_for_BStrand_Del1_width5.pdf", width = 6, height = 6)
ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T)
dev.off()
