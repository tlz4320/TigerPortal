#在2.0_total_First_pop_peak_indel1.R

high_cor_sample <- total_pair_sg_cor[total_pair_sg_cor$pair != "Sg_9_60-Sg_15_103",]
indel1_pct <- list()
for(name in high_cor_sample$id){
  ids <- sgRNA_pair_remain[[name]]
  sel_sg1 <- total_mean[[ids[1]]]
  # sel_sg1 <- sel_sg1[sel_sg1$indel_size %in% region,]
  sel_sg2 <- total_mean[[ids[2]]]
  # sel_sg2 <- sel_sg2[sel_sg2$indel_size %in% region,]
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
  sel_sg2[,3] <- sel_sg2[,2] / sum(sel_sg2[,2])
  peak_indel <- (as.numeric(c(sel_sg1[,1], sel_sg2[,1])[which.max(c(sel_sg1[,2], sel_sg2[,2]))]))
  
  indel1_pct[[name]] <- data.frame(Percent = c(sel_sg1[which.max(sel_sg1[,2]),3],
                                               sel_sg2[which.max(sel_sg2[,2]),3]), 
                                   ID = ids)
}
indel1_pct <- data.frame(do.call(rbind, indel1_pct))
pdf("~/data/project/ear_project/gene_therapy_ll/Result/0AA_restore_result_fix.pdf", width = 20, height = 6)
ggplot(indel1_pct) + geom_bar(aes(x = ID, y = Percent*100), stat = "identity")+
  geom_text( 
    aes(x = ID, y = Percent*100,label = round(Percent, 2)* 100), vjust = -0.1) + 
  ylab("0AA restore % in indel reads")+xlab("") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        text = element_text(size=15))
dev.off()
