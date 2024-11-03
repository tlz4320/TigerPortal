plot_data <- openxlsx::read.xlsx("~/Documents/ear_project/gene_therapy/OTOF/相关性数据.xlsx",
                                 colNames=F)
plot_data <- reshape2::melt(plot_data, id.vars = c("X1", "X5"))
plot_data <- na.omit(plot_data)
colnames(plot_data) <- c("sgRNA", "Responders", "vf", "inframe")
library(ggpubr)

cor(plot_data$Responders, plot_data$inframe)
text_data <- data.frame(text = paste0("Pearson's r=",round(cor(plot_data$Responders, plot_data$inframe), 3), "\n",
                                              paste0("P=",scales::label_scientific()(
                                                cor.test(plot_data$Responders, plot_data$inframe)$p.value))))
plot_data_mean <- lapply(split(plot_data, plot_data$sgRNA), function(x){
  x[,4] <- mean(x[,4])
  x[1,]
})
plot_data_mean <- data.frame(do.call(rbind, plot_data_mean))
pdf("cor_of_4_sgRNA_inframe_response4.pdf", width = 4, height = 2.7)
ggplot(plot_data_mean, aes(x = as.character(round(Responders, 1)), y = inframe)) + 
  geom_point(aes(color = sgRNA)) + 
  geom_text(aes(x=2, y=45, label = text), text_data) + 
  theme_bw() +
  # geom_smooth(method = "lm", formula = y~exp(x)) + 
  # scale_y_continuous(trans = scales::log_trans(4)) + 
  # scale_x_continuous(trans = scales::exp_trans(1.1)) + 
  # theme_classic2() +
  # theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
  #                                                    ends = "last")),
  #       axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
  #                                                      ends = "last"))) + 
  ylab("Inframe(%)") + 
  xlab("Responders in electroporated cells(%)")
dev.off()


pdf("cor_of_4_sgRNA_inframe_response5.pdf", width = 4, height = 2.7)
ggplot(plot_data_mean, aes(x = Responders, y = inframe)) + 
  geom_point(aes(color = sgRNA)) + 
  # geom_text(aes(x=60, y=45, label = text), text_data) + 
  theme_bw() +
  geom_smooth(method = "lm", formula = y ~ exp(x) , se = F, data = plot_data) +
  # scale_y_continuous(trans = scales::log_trans(4)) + 
  # scale_x_continuous(trans = scales::exp_trans(1.1)) + 
  # theme_classic2() +
  # theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
  #                                                    ends = "last")),
  #       axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
  #                                                      ends = "last"))) + 
  ylab("Inframe(%)") + 
  xlab("Responders in electroporated cells\n(%)") + 
  coord_flip()
dev.off()


pdf("cor_of_4_sgRNA_inframe_response.pdf", width = 4, height = 2.7)
ggplot(plot_data, aes(x = inframe, y = Responders)) + 
  geom_point(aes(color = sgRNA)) + 
  geom_text(aes(x=60, y=45, label = text), text_data) + 
  theme_bw() +
  # geom_smooth(method = "lm", formula = y~exp(x)) + 
  # scale_y_continuous(trans = scales::log_trans(4)) + 
  # scale_x_continuous(trans = scales::exp_trans(1.1)) + 
  # theme_classic2() +
  # theme(axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
  #                                                    ends = "last")),
  #       axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
  #                                                      ends = "last"))) + 
  xlab("Inframe(%)") + 
  ylab("Responders in electroporated cells(%)")
dev.off()
