setwd('~/indelphi_res/')
indelphi_118_res <- list()
indelphi_file <- read.csv("~/118_for_indelphi.csv")
for(file in indelphi_file$id2){
  indelphi_118_res[[file]] <- read.table(paste0(file, ".csv"), sep = ",", header = T)
}
indelphi_118_res <- lapply(indelphi_118_res, function(x){
  x$Length[x$Category == "del"] <- -x$Length[x$Category == "del"]
  res <- lapply(split(x$Predicted.frequency, x$Length), sum)
  res <- data.frame(indel_size = names(res), fq = unlist(res))
  res$Pct <- res$fq
  res <- res[gtools::mixedorder(res$indel_size),]
  res
})

setwd('~/forecast_result/')
forecast_118_res <- list()
forcast_file <- read.table("~/118_for_forecast.txt")
for(file in forcast_file$V1){
  forecast_118_res[[file]] <- read.table(paste0(file, "_predictedindelsummary.txt"), 
                                         sep = "\t", header = T)
}
forecast_118_res <- lapply(forecast_118_res, function(x){
  edit <- unlist(lapply(x[,1], function(y){
    unlist(strsplit(y, "[_]"))[1]
  }))
  edit <- data.frame(type = str_sub(edit, 1, 1), len = str_sub(edit, 2), counts = x[,3])
  edit$len <- as.numeric(edit$len)
  edit$len[edit$type == "D"] <- -edit$len[edit$type == "D"]
  res <- lapply(split(edit$counts, edit$len), sum)
  res <- data.frame(indel_size = names(res), fq = unlist(res))
  res$Pct <- res$fq / sum(res$fq) * 100
  res <- res[gtools::mixedorder(res$indel_size),]
  res
})

load("~/data/project/ear_project/gene_therapy_ll/Previews/total_used_seq_split.rda")
total_ids <- str_remove(total_used_seq_tissue$result, "[-].*")
total_used_seq_tissue$id <- total_ids
inframe_result_indelphi <- lapply(total_ids, function(id){
  tmp <- indelphi_118_res[[id]]
  info <- total_used_seq_tissue[total_used_seq_tissue$id == id,]
  total_reads <- sum(tmp$fq)
  
  if(info$mutType == "d"){
    tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 1, ]
    tmp$inframe <- as.numeric(tmp$indel_size) - 1
  } else {
    tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 2, ]
    tmp$inframe <- as.numeric(tmp$indel_size) + 1
  }

  
  tmp <- lapply(split(tmp$fq, tmp$inframe), sum)
  tmp <- data.frame(indel = as.numeric(names(tmp)), fq = unlist(tmp))
  tmp$pct <- tmp$fq / total_reads * 100
  res <- tmp
  res <- res[res$pct != 0,]
  res$aa <- as.numeric(res$indel) / 3
  res
})
names(inframe_result_indelphi) <- total_ids
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


inframe_result_indelphi <- lapply(total_ids, function(id){
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
names(inframe_result_indelphi) <- total_ids
inframe_result_indelphi <- data.frame(do.call(rbind, inframe_result_indelphi))



inframe_result_forecast <- lapply(total_ids, function(id){
  tmp <- forecast_118_res[[id]]
  info <- total_used_seq_tissue[total_used_seq_tissue$id == id,]
  total_reads <- sum(tmp$fq)
  
  if(info$mutType == "d"){
    tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 1, ]
    tmp$inframe <- as.numeric(tmp$indel_size) - 1
  } else {
    tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 2, ]
    tmp$inframe <- as.numeric(tmp$indel_size) + 1
  }
  
  tmp <- lapply(split(tmp$fq, tmp$inframe), sum)
  tmp <- data.frame(indel = as.numeric(names(tmp)), fq = unlist(tmp))
  tmp$pct <- tmp$fq / total_reads * 100
  res <- tmp
  res <- res[res$pct != 0,]
  res$aa <- as.numeric(res$indel) / 3
  res
})
names(inframe_result_forecast) <- total_ids
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


inframe_result_forecast <- lapply(total_ids, function(id){
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
names(inframe_result_forecast) <- total_ids
inframe_result_forecast <- data.frame(do.call(rbind, inframe_result_forecast))


print(load("~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda"))
inframe_result_processed$id2 <- str_remove(inframe_result_processed$id, "[-].*")
cell_inframe_result_processed$id2 <- str_remove(cell_inframe_result_processed$id, "[-].*")



total_used_seq_tissue$tissueInframe <- 0
total_used_seq_tissue$cellInframe <- 0
total_used_seq_tissue$forecastInframe <- 0
total_used_seq_tissue$indelphiInframe <- 0
total_used_seq_tissue$tissue0AA <- 0
total_used_seq_tissue$cell0AA<- 0
total_used_seq_tissue$forecast0AA <- 0
total_used_seq_tissue$indelphi0AA <- 0

for(i in 1 : nrow(total_used_seq_tissue)){
  id <- total_used_seq_tissue$id[i]
  total_used_seq_tissue$tissueInframe[i] <- sum(inframe_result_processed$pct[inframe_result_processed$id2 == id])
  total_used_seq_tissue$cellInframe[i] <- sum(cell_inframe_result_processed$pct[cell_inframe_result_processed$id2 == id])
  total_used_seq_tissue$forecastInframe[i] <- sum(inframe_result_forecast$pct[inframe_result_forecast$id == id])
  total_used_seq_tissue$indelphiInframe[i] <- sum(inframe_result_indelphi$pct[inframe_result_indelphi$id == id])
  total_used_seq_tissue$tissue0AA[i] <- sum(inframe_result_processed$pct2[inframe_result_processed$id2 == id & 
                                                                            inframe_result_processed$aa == "0"])
  total_used_seq_tissue$cell0AA[i] <- sum(cell_inframe_result_processed$pct2[cell_inframe_result_processed$id2 == id & 
                                                                               cell_inframe_result_processed$aa == "0"])
  total_used_seq_tissue$forecast0AA[i] <- sum(inframe_result_forecast$pct2[inframe_result_forecast$id == id & 
                                                                             inframe_result_forecast$aa == "0"])
  total_used_seq_tissue$indelphi0AA[i] <- sum(inframe_result_indelphi$pct2[inframe_result_indelphi$id == id & 
                                                                             inframe_result_indelphi$aa == "0"])

}
total_used_seq_tissue[is.na(total_used_seq_tissue)] <- 0
slope_cal <-lapply(split(total_used_seq_tissue[,c("tissueInframe", "tissue0AA", 
                                       "cellInframe", "cell0AA")], total_used_seq_tissue$id), function(x){
  suby <- abs(x[1,1]) + abs(x[1,3])
  subx <- abs(x[1,2]) - abs(x[1,4])
  if(suby == 0){
    return(0)
  }
  if(subx == 0 & suby != 0){
    return(99999)
  }
  
  return(abs(suby / subx))
})
slope_cal <- data.frame(id = names(slope_cal), slope = unlist(slope_cal))
plot_data <- total_used_seq_tissue
plot_data <- merge(plot_data, slope_cal, by="id")
plot_data$large60 <- plot_data$slope > sqrt(3)


slope_cal <-lapply(split(total_used_seq_tissue[,c("tissueInframe", "tissue0AA", 
                                                  "forecastInframe", "forecast0AA")], total_used_seq_tissue$id), function(x){
                                                    suby <- abs(x[1,1]) + abs(x[1,3])
                                                    subx <- abs(x[1,2]) - abs(x[1,4])
                                                    if(suby == 0){
                                                      return(0)
                                                    }
                                                    if(subx == 0 & suby != 0){
                                                      return(99999)
                                                    }
                                                    
                                                    return(abs(suby / subx))
                                                  })
slope_cal <- data.frame(id = names(slope_cal), slope2 = unlist(slope_cal))
plot_data <- merge(plot_data, slope_cal, by="id")
plot_data$large60_2 <- plot_data$slope2 > sqrt(3)


slope_cal <-lapply(split(total_used_seq_tissue[,c("tissueInframe", "tissue0AA", 
                                                  "indelphiInframe", "indelphi0AA")], total_used_seq_tissue$id), function(x){
                                                    suby <- abs(x[1,1]) + abs(x[1,3])
                                                    subx <- abs(x[1,2]) - abs(x[1,4])
                                                    if(suby == 0){
                                                      return(0)
                                                    }
                                                    if(subx == 0 & suby != 0){
                                                      return(99999)
                                                    }
                                                    
                                                    return(abs(suby / subx))
                                                  })
slope_cal <- data.frame(id = names(slope_cal), slope3 = unlist(slope_cal))

plot_data <- merge(plot_data, slope_cal, by="id")
plot_data$large60_3 <- plot_data$slope3 > sqrt(3)

plot_data$cellalpha <- stats::approxfun(c(sqrt(3), 1/ tan(15 / 180 * pi),15), c(0.1,0.7,1))(plot_data$slope)
plot_data$cellalpha[plot_data$slope < sqrt(3)] <- 0.1
plot_data$cellalpha[plot_data$slope > 15] <- 1

plot_data$tissuealpha <- stats::approxfun(c(sqrt(3), 1/ tan(15 / 180 * pi),15), c(0.1,0.7,1))(plot_data$slope)
plot_data$tissuealpha[plot_data$slope < sqrt(3)] <- 0.1
plot_data$tissuealpha[plot_data$slope > 15] <- 1
plot_data$forecastalpha <- stats::approxfun(c(sqrt(3), 1/ tan(15 / 180 * pi),15), c(0.1,0.7,1))(plot_data$slope2)
plot_data$forecastalpha[plot_data$slope2 < sqrt(3)] <- 0.1
plot_data$forecastalpha[plot_data$slope2 > 15] <- 1
plot_data$indelphialpha <- stats::approxfun(c(sqrt(3), 1/ tan(15 / 180 * pi),15), c(0.1,0.7,1))(plot_data$slope3)
plot_data$indelphialpha[plot_data$slope3 < sqrt(3)] <- 0.1
plot_data$indelphialpha[plot_data$slope3 > 15] <- 1

plot_data_alpha <- plot_data[,c("id", "tissuealpha", "cellalpha", "forecastalpha", "indelphialpha")]
plot_data_alpha <- reshape2::melt(plot_data_alpha)
plot_data_alpha$variable <- str_remove(plot_data_alpha$variable, "alpha")
plot_data_alpha$tmp <- paste0(plot_data_alpha$id, plot_data_alpha$variable)
colnames(plot_data_alpha)[3] <- "alpha"

plot_data_inframe <- plot_data[,c("id", "tissueInframe", "cellInframe", "forecastInframe", "indelphiInframe")]
plot_data_inframe <- reshape2::melt(plot_data_inframe)
plot_data_inframe$variable <- str_remove(plot_data_inframe$variable, "Inframe")
plot_data_inframe$tmp <- paste0(plot_data_inframe$id, plot_data_inframe$variable)
colnames(plot_data_inframe)[3] <- "Inframe"
plot_data_0aa <- plot_data[,c("id", "tissue0AA", "cell0AA", "forecast0AA", "indelphi0AA")]
plot_data_0aa <- reshape2::melt(plot_data_0aa)
plot_data_0aa$variable <- str_remove(plot_data_0aa$variable, "0AA")
plot_data_0aa$tmp <- paste0(plot_data_0aa$id, plot_data_0aa$variable)
colnames(plot_data_0aa)[3] <- "0AA"
plot_data2 <- merge(plot_data_inframe, plot_data_0aa[c(-1,-2)], by="tmp")
plot_data2 <- merge(plot_data2, plot_data_alpha[c(-1,-2)], by="tmp")

plot_data2$Inframe[plot_data2$variable == "cell"] <- -plot_data2$Inframe[plot_data2$variable == "cell"]
plot_data2$Inframe[plot_data2$variable == "indelphi"] <- plot_data2$Inframe[plot_data2$variable == "indelphi"] + 100
plot_data2$Inframe[plot_data2$variable == "forecast"] <- plot_data2$Inframe[plot_data2$variable == "forecast"] + 200

###
# plot_data2$Inframe[plot_data2$variable == "tissue"] <- -plot_data2$Inframe[plot_data2$variable == "tissue"]
# plot_data2$Inframe[plot_data2$variable == "indelphi"] <- -plot_data2$Inframe[plot_data2$variable == "indelphi"] - 100
# plot_data2$Inframe[plot_data2$variable == "forecast"] <- -plot_data2$Inframe[plot_data2$variable == "forecast"] - 200
###
plot_data_line <- plot_data2[order(plot_data2$id),]
plot_data2$type <- "Del1"
plot_data2$type[plot_data2$id %in% total_used_seq_tissue$id[total_used_seq_tissue$mutType == "i"]] <- "Ins1"
plot_data2$type[plot_data2$variable == "tissue"] <- paste("Tissue", plot_data2$type[plot_data2$variable == "tissue"], sep = "-")
plot_data2$type[plot_data2$variable != "tissue"] <- plot_data2$variable[plot_data2$variable != "tissue"]

# stat <- plot_data2 %>%
#   wilcox_test(Inframe ~ variable , paired = T)
# stat <- add_xy_position(stat, x = "type")
# stat <- stat[stat$group1 == "target",]
# stat$y.position <- c(106, 118, 130)

# pdf("~/Nutstore Files/Tobin/Previous/118_inframe_0aa_pointPlot_add_predict_v3.pdf", width = 12, height = 4)
# ggplot() + 
#   geom_line(data = plot_data_line ,aes(y = abs(`0AA`), x = Inframe, 
#                                    group = id, alpha = alpha), linewidth = 0.1,orientation = "x") + 
#   geom_point(data = plot_data2 ,aes(y = abs(`0AA`), x = Inframe, 
#                                     color = type), alpha = plot_data2$alpha) +
#   # geom_abline(slope = c(1, -1), intercept = 0) + 
#   geom_vline(xintercept = c(0,  100, 200), linewidth = 1) + 
#   labs(x = "Inframe%", y = "0AA% in Inframe") + 
#   scale_y_continuous(limits = c(0, 100), expand = c(0,1))  +
#   scale_x_continuous(limits = c(-100, 300), breaks = seq(-100, 300, 50), 
#                      labels = c(100, "Cell 50", "0", "Tissue 50", 
#                                 "100", "IndelPhi 50", "100", "ForeCasT 50", "100"), 
#                      expand = c(0,0)) + 
#   scale_color_manual(values = c("#FF6A6A", "#FF6A6A", "#FF6A6A","#B46F28FF", "#0A5028")) + 
#   theme_classic2() + theme(axis.text.x = element_text(angle = 90, 
#                                                       hjust = 1, 
#                                                       vjust = 0.5), 
#                            panel.grid  = element_blank()
#                            )
# dev.off()



plot_data3 <- plot_data2[plot_data2$type %in% c("Tissue-Del1", "cell") & 
                           plot_data2$id %in% total_used_seq_tissue$id[total_used_seq_tissue$mutType == 'd'],]
plot_data3$Inframe[plot_data3$variable == "cell"] <- -plot_data3$Inframe[plot_data3$variable == "cell"]
plot_data3$Inframe[plot_data3$variable == "tissue"] <- -plot_data3$Inframe[plot_data3$variable == "tissue"]

plot_data3 <- plot_data3[order(plot_data3$id),]
pdf("~/Nutstore Files/Tobin/Previous/87_inframe_0aa_pointPlot_tissue_cell_Del1Mut_v2.pdf", width = 5, height = 4)
ggplot(plot_data3) + 
  geom_line(aes(y = abs(`0AA`), x = Inframe, group = id), 
            linewidth = 0.1, 
            alpha = plot_data3$alpha) + 
  geom_point(aes(y = abs(`0AA`), x = Inframe, 
                 color = type), alpha = plot_data3$alpha) +
  # geom_abline(slope = c(1, -1), intercept = 0) + 
  geom_vline(xintercept = c(0), linewidth = 1) + 
  labs(x = "Inframe%", y = "0-NUM*/In-frame") + 
  scale_y_continuous(limits = c(0, 100), expand = c(0,1))  +
  scale_x_continuous(limits = c(-100, 100)) + 
  scale_color_manual(values = c("#ab1d25","#B46F28FF")) + 
  theme_classic2() + theme(axis.text.x = element_text(angle = 90, 
                                                      hjust = 1, 
                                                      vjust = 0.5), 
                           panel.grid  = element_blank()
  )
dev.off()

pdf("~/Nutstore Files/Tobin/Previous/87_inframe_0aa_Del1Mut_tissue_cell_legend_v3.pdf", width = 3, height = 7)

legend_plot <- data.frame(id = 1 : 87, 
                          color = rev(alpha("#ab1d25", sort(plot_data3$alpha[plot_data3$variable == "cell"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data$slope[plot_data$mutType == "d"] < tan(x / 180 * pi))
})) + 0.5

p2 <- ggplot(legend_plot) + geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4 ) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 87) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.3) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 87) + 0.5), 
            aes(x = x, y = 0.5, label = floor(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.3) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.3) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 87.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()


legend_plot <- data.frame(id = 1 : 87, 
                          color = rev(alpha("#B46F28FF", sort(plot_data3$alpha[plot_data3$variable != "cell"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data$slope[plot_data$mutType == "d"] < tan(x / 180 * pi))
})) + 0.5

p1 <- ggplot(legend_plot) + geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4 ) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 87) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.3) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 87) + 0.5), 
            aes(x = x, y = 0.5, label = floor(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.3) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.3) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 87.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()
ggpubr::ggarrange(p1, p2, nrow = 1)
dev.off()

plot_data4 <- plot_data2[plot_data2$type %in% c("Tissue-Ins1", "cell") & 
                           plot_data2$id %in% total_used_seq_tissue$id[total_used_seq_tissue$mutType == 'i'],]

plot_data4$Inframe[plot_data4$variable == 'tissue'] <- -plot_data4$Inframe[plot_data4$variable == 'tissue']
plot_data4$Inframe[plot_data4$variable != 'tissue'] <- -plot_data4$Inframe[plot_data4$variable != 'tissue']
plot_data4 <- plot_data4[order(plot_data4$id),]
pdf("~/Nutstore Files/Tobin/Previous/31_inframe_0aa_pointPlot_add_tissue_cell_Ins1Mut_v2.pdf", width = 5, height = 4)
ggplot() + 
  geom_line(data = plot_data4 ,aes(y = abs(`0AA`), x = Inframe, 
                                   group = id), alpha = plot_data4$alpha,
            linewidth = 0.1,orientation = "x") + 
  geom_point(data = plot_data4 ,aes(y = abs(`0AA`), x = Inframe, 
                                    color = type), alpha = plot_data4$alpha) +
  # geom_abline(slope = c(1, -1), intercept = 0) + 
  geom_vline(xintercept = c(0), linewidth = 1) + 
  labs(x = "Inframe%", y = "0-NUM*/In-frame") + 
  scale_y_continuous(limits = c(0, 100), expand = c(0,1))  +
  scale_x_continuous(limits = c(-100, 100)) + 
  scale_color_manual(values = c("#ab1d25","#0A5028")) + 
  theme_classic2() + theme(axis.text.x = element_text(angle = 90, 
                                                      hjust = 1, 
                                                      vjust = 0.5), 
                           panel.grid  = element_blank()
  )
dev.off()

pdf("~/Nutstore Files/Tobin/Previous/31_inframe_0aa_Ins1Mut-tissue_cell_legend_v3.pdf", width = 3, height = 7)

legend_plot <- data.frame(id = 1 : 31, 
                          color = rev(alpha("#ab1d25", sort(plot_data4$alpha[plot_data4$variable == "cell"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data$slope[plot_data$mutType == "i"] < tan(x / 180 * pi))
})) + 0.5
p2 <- ggplot(legend_plot) + geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4 ) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 31) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.08) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 31) + 0.5), 
            aes(x = x, y = 0.5, label = floor(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.1) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.1) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 31.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()


legend_plot <- data.frame(id = 1 : 31, 
                          color = rev(alpha("#0A5028", sort(plot_data4$alpha[plot_data4$variable != "cell"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data$slope[plot_data$mutType == "i"] < tan(x / 180 * pi))
})) + 0.5

p1 <- ggplot(legend_plot) + geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4 ) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 31) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.08) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 31) + 0.5), 
            aes(x = x, y = 0.5, label = floor(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.1) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.1) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 31.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()
ggpubr::ggarrange(p1, p2, nrow = 1)
dev.off()


plot_data3 <- plot_data2[plot_data2$type %in% c("Tissue-Del1", "indelphi") & 
                           plot_data2$id %in% total_used_seq_tissue$id[total_used_seq_tissue$mutType == 'd'],]
alpha_table <- plot_data3[plot_data3$type == "indelphi",]
rownames(alpha_table) <- alpha_table$id
plot_data3$alpha <- alpha_table[plot_data3$id, "alpha"]
plot_data3$Inframe[plot_data3$variable == 'tissue'] <- -plot_data3$Inframe[plot_data3$variable == 'tissue']
plot_data3$Inframe[plot_data3$variable != 'tissue'] <- plot_data3$Inframe[plot_data3$variable != 'tissue']- 100
plot_data3 <- plot_data3[order(plot_data3$id),]
pdf("~/Nutstore Files/Tobin/Previous/87_inframe_0aa_pointPlot_add_predict_Del1Mut_v2.pdf", width = 5, height = 4)
ggplot(plot_data3) + 
  geom_line(aes(y = abs(`0AA`), x = Inframe, group = id), 
            linewidth = 0.1, 
            alpha = plot_data3$alpha) + 
  geom_point(aes(y = abs(`0AA`), x = Inframe, 
                                    color = type), alpha = plot_data3$alpha) +
  # geom_abline(slope = c(1, -1), intercept = 0) + 
  geom_vline(xintercept = c(0), linewidth = 1) + 
  labs(x = "Inframe%", y = "0-NUM*/In-frame") + 
  scale_y_continuous(limits = c(0, 100), expand = c(0,1))  +
  scale_x_continuous(limits = c(-100, 100)) + 
  scale_color_manual(values = c("#ab1d25","#B46F28FF")) + 
  theme_classic2() + theme(axis.text.x = element_text(angle = 90, 
                                                      hjust = 1, 
                                                      vjust = 0.5), 
                           panel.grid  = element_blank()
  )
dev.off()





pdf("~/Nutstore Files/Tobin/Previous/87_inframe_0aa_Del1Mut_legend_v3.pdf", width = 3, height = 7)

legend_plot <- data.frame(id = 1 : 87, 
                          color = rev(alpha("#ab1d25", sort(plot_data3$alpha[plot_data3$variable == "indelphi"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data$slope3[plot_data$mutType == "d"] < tan(x / 180 * pi))
})) + 0.5

p2 <- ggplot(legend_plot) + geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4 ) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 87) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.3) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 87) + 0.5), 
            aes(x = x, y = 0.5, label = floor(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.3) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.3) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 87.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()


legend_plot <- data.frame(id = 1 : 87, 
                          color = rev(alpha("#B46F28FF", sort(plot_data3$alpha[plot_data3$variable != "indelphi"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data$slope3[plot_data$mutType == "d"] < tan(x / 180 * pi))
})) + 0.5

p1 <- ggplot(legend_plot) + geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4 ) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 87) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.3) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 87) + 0.5), 
            aes(x = x, y = 0.5, label = floor(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.3) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.3) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 87.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()
ggpubr::ggarrange(p1, p2, nrow = 1)
dev.off()


plot_data4 <- plot_data2[plot_data2$type %in% c("Tissue-Ins1", "indelphi") & 
                           plot_data2$id %in% total_used_seq_tissue$id[total_used_seq_tissue$mutType == 'i'],]
alpha_table <- plot_data4[plot_data4$type == "indelphi",]

rownames(alpha_table) <- alpha_table$id
plot_data4$alpha <- alpha_table[plot_data4$id, "alpha"]
plot_data4$Inframe[plot_data4$variable == 'tissue'] <- -plot_data4$Inframe[plot_data4$variable == 'tissue']
plot_data4$Inframe[plot_data4$variable != 'tissue'] <- plot_data4$Inframe[plot_data4$variable != 'tissue']- 100
plot_data4 <- plot_data4[order(plot_data4$id),]
pdf("~/Nutstore Files/Tobin/Previous/31_inframe_0aa_pointPlot_add_predict_Ins1Mut_v2.pdf", width = 5, height = 4)
ggplot() + 
  geom_line(data = plot_data4 ,aes(y = abs(`0AA`), x = Inframe, 
                                   group = id, alpha = alpha), linewidth = 0.1,orientation = "x") + 
  geom_point(data = plot_data4 ,aes(y = abs(`0AA`), x = Inframe, 
                                    color = type), alpha = plot_data4$alpha) +
  # geom_abline(slope = c(1, -1), intercept = 0) + 
  geom_vline(xintercept = c(0), linewidth = 1) + 
  labs(x = "Inframe%", y = "0-NUM*/In-frame") + 
  scale_y_continuous(limits = c(0, 100), expand = c(0,1))  +
  scale_x_continuous(limits = c(-100, 100)) + 
  scale_color_manual(values = c("#ab1d25","#0A5028")) + 
  theme_classic2() + theme(axis.text.x = element_text(angle = 90, 
                                                      hjust = 1, 
                                                      vjust = 0.5), 
                           panel.grid  = element_blank()
  )
dev.off()

pdf("~/Nutstore Files/Tobin/Previous/31_inframe_0aa_Ins1Mut_legend_v3.pdf", width = 3, height = 7)

legend_plot <- data.frame(id = 1 : 31, 
                          color = rev(alpha("#ab1d25", sort(plot_data4$alpha[plot_data4$variable == "indelphi"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data$slope3[plot_data$mutType == "i"] < tan(x / 180 * pi))
})) + 0.5
p2 <- ggplot(legend_plot) + geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4 ) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 31) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.08) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 31) + 0.5), 
            aes(x = x, y = 0.5, label = floor(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.1) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.1) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 31.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()


legend_plot <- data.frame(id = 1 : 31, 
                          color = rev(alpha("#0A5028", sort(plot_data4$alpha[plot_data4$variable != "indelphi"]))))
arctan_list <- unlist(lapply(c(0, 15, 30, 45), function(x){
  sum(1 / plot_data$slope3[plot_data$mutType == "i"] < tan(x / 180 * pi))
})) + 0.5

p1 <- ggplot(legend_plot) + geom_tile(aes(x = id, y = 1), fill = legend_plot$color,height=0.4 ) +
  geom_tile(data = data.frame(x = c(0 : 10 / 10 * 31) + 0.5), 
            aes(x = x, y = 0.7), height = 0.1, width = 0.08) + 
  geom_text(data = data.frame(x = c(0 : 10 / 10 * 31) + 0.5), 
            aes(x = x, y = 0.5, label = floor(x))) + 
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 1.15),
            height = 0.1, width = 0.1) +
  geom_tile(data = data.frame(x = arctan_list),
            aes(x = x, y = 0.85),
            height = 0.1, width = 0.1) +
  geom_text(data = data.frame(x = arctan_list, label = paste0(c(0, 15, 30, 45), "°")), 
            aes(x = x, y = 1, label = label)) + 
  geom_line(data = data.frame(x = c(0.5, 31.5), y = c(0.7 ,0.7)),
            aes(x = x, y = y)) + scale_y_continuous(limits = c(0, 1.5), expand = c(0,0)) + 
  theme_void() + coord_flip()
ggpubr::ggarrange(p1, p2, nrow = 1)
dev.off()



inframe_pval <- apply(combn(4,2), 2,function(step){
  step <- unlist(step)
  name <- paste(colnames(total_used_seq_tissue)[12 + step], collapse = "-")
  tmp <- total_used_seq_tissue[,12 + step]
  pval <- wilcox.test(tmp[,1], tmp[,2], paired = T)$p.value
  data.frame(name = name, pval = pval)
})
aa0_pval <- apply(combn(4,2), 2,function(step){
  step <- unlist(step)
  name <- paste(colnames(total_used_seq_tissue)[16 + step], collapse = "-")
  tmp <- total_used_seq_tissue[,16 + step]
  pval <- wilcox.test(tmp[,1], tmp[,2], paired = T)$p.value
  data.frame(name = name, pval = pval)
})
inframe_pval <- data.frame(rbind(data.frame(do.call(rbind, inframe_pval)), 
                                 data.frame(do.call(rbind, aa0_pval))))


openxlsx::write.xlsx(inframe_pval, file="~/Nutstore Files/Tobin/Previous/118_inframe_0aa_pointPlot_pval.xlsx", rowNames=F, colNames=T)

plot_data_slope <- plot_data[,c("id", "slope", "slope2", "slope3")]
colnames(plot_data_slope) <- c("id", "tissue_cell", "tissue_forecast", "tissue_indelphi")
for(i in 2 : 4){
  plot_data_slope[,i] <- 1 / plot_data_slope[,i]
  plot_data_slope[is.infinite(plot_data_slope[,i]),i] <- maxregion
}
maxregion <- max(unlist(plot_data_slope[,c(2:4)]))
maxy <- c()
for(i in 2 : 4){
  plot_data <- hist(plot_data_slope[,i], breaks = seq(0, maxregion, maxregion / 50))
  plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                          counts = plot_data$counts)
  maxy <- c(maxy, max(plot_data$counts))
  
}
maxy <- max(maxy)
plot_list <- list()
for(i in 2 : 4){
  plot_data <- hist(plot_data_slope[,i], breaks = seq(0, maxregion, maxregion / 50))
  plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                          counts = plot_data$counts)
# plot_data$x <- plot_data$inframePct - 2.5

  plot_list[[i - 1]] <- ggplot(plot_data) + geom_bar(aes(x = inframePct, y = counts),
                                   color = "black",
                                   stat="identity", fill = '#48D1CC') +
    geom_text(aes(x = inframePct, y = counts, label = counts)) + 
    geom_vline(xintercept = 1 / sqrt(3), linewidth = 0.1) + 
  xlab("Slope") + ylab("Counts") +
  scale_y_continuous(limits = c(0, maxy + 1)) +
  theme_bw() + ggtitle(colnames(plot_data_slope)[i]) + 
  theme(plot.margin =  unit(c(0,0,0,0.2), "cm"), 
        panel.grid = element_blank())
}
pdf("~/Nutstore Files/Tobin/Previous/118_0aa_inframe_slope_add_predict_50bn.pdf", width = 6, height =8)
ggpubr::ggarrange(plotlist = plot_list, ncol = 1)
dev.off()







plot_data <- hist(total_used_seq_tissue$indelphiInframe[total_used_seq_tissue$mutType == "d"],
                  breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)


plot_data$color <- "#B46F2877"
plot_data$color[plot_data$inframePct > 30] <-"#B46F28FF" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(total_used_seq_tissue$indelphiInframe[total_used_seq_tissue$mutType == "d"] > pct)
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")
p1 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                             color = "black",
                             stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = 5, label = anno_info, data = data.frame(a=1)) +
  xlab("Inframe Frequency") + ylab("Counts") +
  scale_y_continuous(breaks = 1 : 24, limits = c(0, 24)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), limits = c(0, 100)) +
  theme_bw() + ggtitle("IndelPhi Inframe") + 
  theme(plot.margin =  unit(c(0,0,0,0.2), "cm"), 
        panel.grid = element_blank())

plot_data <- hist(total_used_seq_tissue$forecastInframe[total_used_seq_tissue$mutType == "d"],
                  breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)


plot_data$color <- "#B46F2877"
plot_data$color[plot_data$inframePct > 30] <-"#B46F28FF" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(total_used_seq_tissue$forecastInframe[total_used_seq_tissue$mutType == "d"] > pct)
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p2 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = -5, label = anno_info, data = data.frame(a=1)) +
  scale_y_reverse(expand = c(0,0), breaks = 1 : 24, limits = c(24, 0)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), position = "top", limits = c(0, 100)) +
  ylab("Counts") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        plot.margin =  unit(c(-1,0,1,0.2), "cm"))

pdf("~/Nutstore Files/Tobin/Previous/indelphi_forecast_inframe_histgram_delMut_v1.pdf", width = 6, height = 6)
ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T)
dev.off()



plot_data <- hist(total_used_seq_tissue$indelphiInframe[total_used_seq_tissue$mutType == "i"],
                  breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)


plot_data$color <- "#0A502866"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(total_used_seq_tissue$indelphiInframe[total_used_seq_tissue$mutType == "i"] > pct)
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")
p1 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = 5, label = anno_info, data = data.frame(a=1)) +
  xlab("Inframe Frequency") + ylab("Counts") +
  scale_y_continuous(breaks = 1 : 9, limits = c(0, 9)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), limits = c(0, 100)) +
  theme_bw() + ggtitle("IndelPhi Inframe") + 
  theme(plot.margin =  unit(c(0,0,0,0.2), "cm"), 
        panel.grid = element_blank())

plot_data <- hist(total_used_seq_tissue$forecastInframe[total_used_seq_tissue$mutType == "i"],
                  breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)


plot_data$color <- "#0A502866"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(total_used_seq_tissue$forecastInframe[total_used_seq_tissue$mutType == "i"] > pct)
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p2 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = -5, label = anno_info, data = data.frame(a=1)) +
  scale_y_reverse(expand = c(0,0), breaks = 1 : 9, limits = c(9, 0)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), position = "top", limits = c(0, 100)) +
  ylab("Counts") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        plot.margin =  unit(c(-1,0,1,0.2), "cm"))

pdf("~/Nutstore Files/Tobin/Previous/indelphi_forecast_inframe_histgram_insMut_v1.pdf", width = 6, height = 6)
ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T)
dev.off()
###0AA

plot_data <- hist(total_used_seq_tissue$indelphi0AA[total_used_seq_tissue$mutType == "d"],
                  breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)


plot_data$color <- "#B46F2877"
plot_data$color[plot_data$inframePct > 30] <-"#B46F28FF" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(total_used_seq_tissue$indelphi0AA[total_used_seq_tissue$mutType == "d"] > pct)
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")
p1 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = 5, label = anno_info, data = data.frame(a=1)) +
  xlab("0AA/Inframe") + ylab("Counts") +
  scale_y_continuous(breaks = 1 : 16, limits = c(0, 16)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), limits = c(0, 100)) +
  theme_bw() + ggtitle("IndelPhi Inframe") + 
  theme(plot.margin =  unit(c(0,0,0,0.2), "cm"), 
        panel.grid = element_blank())

plot_data <- hist(total_used_seq_tissue$forecast0AA[total_used_seq_tissue$mutType == "d"],
                  breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)


plot_data$color <- "#B46F2877"
plot_data$color[plot_data$inframePct > 30] <-"#B46F28FF" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(total_used_seq_tissue$forecast0AA[total_used_seq_tissue$mutType == "d"] > pct)
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p2 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = -5, label = anno_info, data = data.frame(a=1)) +
  scale_y_reverse(expand = c(0,0), breaks = 1 : 16, limits = c(16, 0)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), position = "top", limits = c(0, 100)) +
  ylab("Counts") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        plot.margin =  unit(c(-1,0,1,0.2), "cm"))

pdf("~/Nutstore Files/Tobin/Previous/indelphi_forecast_0AA_histgram_delMut_v1.pdf", width = 6, height = 6)
ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T)
dev.off()



plot_data <- hist(total_used_seq_tissue$indelphi0AA[total_used_seq_tissue$mutType == "i"],
                  breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)


plot_data$color <- "#0A502866"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(total_used_seq_tissue$indelphi0AA[total_used_seq_tissue$mutType == "i"] > pct)
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")
p1 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = 5, label = anno_info, data = data.frame(a=1)) +
  xlab("0AA/Inframe") + ylab("Counts") +
  scale_y_continuous(breaks = 1 : 11, limits = c(0, 11)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), limits = c(0, 100)) +
  theme_bw() + ggtitle("IndelPhi Inframe") + 
  theme(plot.margin =  unit(c(0,0,0,0.2), "cm"), 
        panel.grid = element_blank())

plot_data <- hist(total_used_seq_tissue$forecast0AA[total_used_seq_tissue$mutType == "i"],
                  breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)


plot_data$color <- "#0A502866"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(total_used_seq_tissue$forecast0AA[total_used_seq_tissue$mutType == "i"] > pct)
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p2 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = -5, label = anno_info, data = data.frame(a=1)) +
  scale_y_reverse(expand = c(0,0), breaks = 1 : 11, limits = c(11, 0)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), position = "top", limits = c(0, 100)) +
  ylab("Counts") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        plot.margin =  unit(c(-1,0,1,0.2), "cm"))

pdf("~/Nutstore Files/Tobin/Previous/indelphi_forecast_0AA_histgram_insMut_v1.pdf", width = 6, height = 6)
ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T)
dev.off()

