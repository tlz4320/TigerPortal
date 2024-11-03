print(load("~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda"))
load("~/data/project/ear_project/gene_therapy_ll/Previews/total_used_seq.rda")
del_sg <- total_used_seq[total_used_seq$mutType == 'd',]
del_sg <- c(del_sg$result, paste0("m", 12 : 15, "-tissue"), paste0("m", 12 : 15))
tissue_0aa <- inframe_result_processed[inframe_result_processed$aa == "0",]
tissue_0aa <- tissue_0aa[tissue_0aa$id %in% del_sg,]
tissue_0aa$pct2[is.na(tissue_0aa$pct2)] <- 0
plot_data <- hist(tissue_0aa$pct2, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)
plot_data$color <- "#B46F2877"
plot_data$color[plot_data$inframePct > 30] <-"#B46F28FF" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")
# pdf("~/Nutstore Files/Tobin/Previous/total_cell_inframe_hist_v2.pdf", width = 6, height = 4)
p1 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                             color = "black",
                             stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = 8, label = anno_info, data = data.frame(a=1)) +
  xlab("0AA/Inframe Frequency") + ylab("Counts") +
  scale_y_continuous(breaks = 1 : 14, limits = c(0, 14)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), limits = c(0, 100)) +
  theme_bw() +
  theme(plot.margin =  unit(c(0,0,0,0.2), "cm"), 
        panel.grid = element_blank())
# dev.off()

cell_0aa <- cell_inframe_result_processed[cell_inframe_result_processed$aa == "0",]
cell_0aa <- cell_0aa[cell_0aa$id %in% del_sg,]
cell_0aa$pct2[is.na(cell_0aa$pct2)] <- 0
plot_data <- hist(cell_0aa$pct2, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)
plot_data$color <- "#B46F2877"
plot_data$color[plot_data$inframePct > 30] <-"#B46F28FF" 
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
  geom_text(x = 10, y = -5, label = anno_info, data = data.frame(a=1)) +
  scale_y_reverse(expand = c(0,0), breaks = 1 : 14, limits = c(14, 0)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), position = "top", limits = c(0, 100)) +
  ylab("Counts") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        plot.margin =  unit(c(-1,0,1,0.2), "cm"))
pdf("~/Nutstore Files/Tobin/Previous/0aa_inframe_histgram_delMut_v3_recolor.pdf", width = 6, height = 6)
ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T)
dev.off()



print(load("~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda"))
load("~/data/project/ear_project/gene_therapy_ll/Previews/total_used_seq.rda")
ins_sg <- total_used_seq[total_used_seq$mutType == 'i',]

tissue_0aa <- inframe_result_processed[inframe_result_processed$aa == "0",]
tissue_0aa <- tissue_0aa[tissue_0aa$id %in% ins_sg$result,]
tissue_0aa$pct2[is.na(tissue_0aa$pct2)] <- 0
plot_data <- hist(tissue_0aa$pct2, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)
plot_data$color <- "#0A502866"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028"
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")
# pdf("~/Nutstore Files/Tobin/Previous/total_cell_inframe_hist_v2.pdf", width = 6, height = 4)
p1 <- ggplot(plot_data) + geom_bar(aes(x = x, y = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color) +
  geom_text(x = 10, y = 8, label = anno_info, data = data.frame(a=1)) +
  xlab("0AA/Inframe Frequency") + ylab("Counts") +
  scale_y_continuous(breaks = 1 : 14, limits = c(0, 14)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), limits = c(0, 100)) +
  theme_bw() +
  theme(plot.margin =  unit(c(0,0,0,0.2), "cm"), 
        panel.grid = element_blank())
# dev.off()

cell_0aa <- cell_inframe_result_processed[cell_inframe_result_processed$aa == "0",]
cell_0aa <- cell_0aa[cell_0aa$id %in% ins_sg$result,]
cell_0aa$pct2[is.na(cell_0aa$pct2)] <- 0
plot_data <- hist(cell_0aa$pct2, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)
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
  geom_text(x = 10, y = -5, label = anno_info, data = data.frame(a=1)) +
  scale_y_reverse(expand = c(0,0), breaks = 1 : 14, limits = c(14, 0)) +
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), position = "top", limits = c(0, 100)) +
  ylab("Counts") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        plot.margin =  unit(c(-1,0,1,0.2), "cm"))
pdf("~/Nutstore Files/Tobin/Previous/0aa_inframe_histgram_insMut_v2_recolor.pdf", width = 6, height = 6)
ggpubr::ggarrange(p1, p2, ncol = 1, common.legend = T)
dev.off()
