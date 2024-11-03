#在2.1之后
setwd("~/indelphi_res/")
indelphi_result <- list()
for(file in list.files()){
  indelphi_result[[file]] <- read.table(file,sep = ",", header = T)
}
names(indelphi_result) <- str_remove(names(indelphi_result), ".csv")
indelphi_result_ins1 <- lapply(indelphi_result, function(x){
  x <- x[x$Category == "ins" & x$Length == 1,]
  x$Pct <- x$Predicted.frequency / sum(x$Predicted.frequency) * 100
  x
})
total_plot_data_sel4 <- total_plot_data_sel3
for(i in 1 : nrow(total_plot_data_sel4)){
  type <- total_plot_data_sel4$NT[i]
  type <- str_sub(type, 1, 1)
  id <- total_plot_data_sel4$sg[i]
  tmp <- indelphi_result_ins1[[id]]
  total_plot_data_sel4$Pct[i] <- tmp$Pct[tmp$Inserted.Bases == type]
}
total_plot_data_sel4 <- total_plot_data_sel4[total_plot_data_sel4$type == "-1 strand",]
total_plot_data_sel4$type <- "inDelphi"
total_plot_data_sel4 <- data.frame(rbind(total_plot_data_sel4, total_plot_data_sel3))
total_plot_data_sel4$type <- factor(total_plot_data_sel4$type, levels = c('+1 strand', '-1 strand', 'inDelphi'))
total_plot_data_sel4$shape[total_plot_data_sel4$type == "inDelphi"] <- 21 
pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_indelphi_Left_Ins1_point_link_plot_rmNN_bulge_cause.pdf", width = 6, height = 4)
ggplot(total_plot_data_sel4, aes(x = type, y = Pct, color = NT, fill = NT, group = id)) +
  geom_line(aes(x = type, y = Pct, linetype = ltp)) + 
  geom_point(size = 2, 
             shape = total_plot_data_sel4$shape, 
             color ="black") + xlab("") + ylab("Predicted Normalized(%)") +
  facet_wrap(~NT, nrow = 1) +
  scale_fill_manual(values = c("A/A" = "#2F89FC",
                               "T/T" = "#30E3CA",
                               "C/C" = "#66CD00",
                               "G/G" = "#98ABEF")) + 
  scale_color_manual(values = c("A/A" = "#2F89FC",
                                "T/T" = "#30E3CA",
                                "C/C" = "#66CD00",
                                "G/G" = "#98ABEF")) +
  theme_classic2() + theme(axis.text.x = element_blank())
dev.off()
ggplot(total_plot_data_sel4, aes(x = type, y = Pct, color = NT, fill = NT, group = id)) +
  geom_line(aes(x = type, y = Pct, linetype = ltp)) + 
  geom_point(size = 2, 
             shape = total_plot_data_sel4$shape, 
             color ="black") + xlab("") + ylab("Predicted Normalized(%)") +
  scale_fill_manual(values = c("A/A" = "#2F89FC",
                               "T/T" = "#30E3CA",
                               "C/C" = "#66CD00",
                               "G/G" = "#98ABEF")) + 
  scale_color_manual(values = c("A/A" = "#2F89FC",
                                "T/T" = "#30E3CA",
                                "C/C" = "#66CD00",
                                "G/G" = "#98ABEF")) +
  theme_classic2() + theme(axis.text.x = element_blank())

total_plot_data_sel4 <- total_plot_data_sel4[order(total_plot_data_sel4$id),]
tmp1 <- total_plot_data_sel4$Pct[total_plot_data_sel4$type == "-1 strand"]
tmp2 <- total_plot_data_sel4$Pct[total_plot_data_sel4$type == "+1 strand"]
tmp3 <- total_plot_data_sel4$Pct[total_plot_data_sel4$type == "inDelphi"]
wilcox.test(tmp1, tmp2, pair = T)
wilcox.test(tmp1, tmp3, pair = T)
