print(load("~/data/project/ear_project/gene_therapy_ll/Result/merged_indel_table_first_second.Rda"))
load("~/data/project/ear_project/gene_therapy_ll/Result/sgCmp132.rda")
total_sg_cor_each_sample <- list()
max_region <- c(-200, 30)
sgRNA_pair <- split(sgCmp132$id2, sgCmp132$id)
for(name in names(total_indel_table_first_second)){
  edit_data <- total_indel_table_first_second[[name]]
  sgids <- sgCmp132$id2
  sgids <- sgids[sgids%in% names(edit_data)]
  
  sgRNA_pair_remain <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
    x <- unlist(x)
    return(sum(x %in% names(edit_data)) == 2)
  }))]
  tmp_cor <- lapply(sgRNA_pair_remain, function(y){
    ids <- unlist(y)
    sel_sg1 <- edit_data[[ids[1]]]
    sel_sg2 <- edit_data[[ids[2]]]
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
  names(tmp_cor) <- unlist(lapply(sgRNA_pair_remain, function(x){
    paste(x, collapse = "-")
  }))
  total_sg_cor_each_sample[[name]] <- data.frame(pair = names(tmp_cor), 
                                                 cor = unlist(lapply(tmp_cor, function(x){x$estimate})),
                                                 pval = unlist(lapply(tmp_cor, function(x){x$p.value})), 
                                                 sample = name)
}
total_sg_cor_each_sample <- data.frame(do.call(rbind, total_sg_cor_each_sample))

mean_cor <- lapply(split(total_sg_cor_each_sample$cor, total_sg_cor_each_sample$pair), mean)
mean_cor <- data.frame(pair = names(mean_cor), 
                       cor = unlist(mean_cor))
mean_cor <- mean_cor[order(mean_cor$cor, decreasing = T),]
plot_data <- total_sg_cor_each_sample

plot_data$pair <- factor(plot_data$pair, levels = mean_cor$pair)
plot_data$pair2 <- factor(as.numeric(plot_data$pair), levels = 1:132)
pdf("~/Nutstore Files/Tobin/Merged1NT/sgRNA_pair_cor_each_batch_sgCmp132.pdf", width = 20, height = 4)
ggplot(plot_data, aes(x = pair, y = cor)) + 
  geom_boxplot() +
  geom_point(aes(color = sample), position = position_dodge2(width = 0.2)) + 
  geom_hline(yintercept = seq(0.1, 1, 0.1)) + xlab("") + 
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20),
        text = element_text(size=15))
dev.off()

