#在2.0_total_First_pop_peak_indel1.R
setwd("~/data/project/ear_project/gene_therapy/find1NTSim/")
tmp <- read.table("human_sg_final_anno2_addid.txt", sep="\t")
tmp2 <- read.table("threshold20000/human_sg_final_anno2_addid.txt", sep="\t")
tmp3 <- read.table("human_sg_lastNTSim_final_addid2.txt", sep = "\t")
tmp3$V14 <- lapply(1 : nrow(tmp3), function(x){
  if(str_sub(tmp3[x,2], 2, 20) == str_sub(tmp3[x, 8], 1, 19)){
    return(paste0(paste0(tmp3[x,2], "-"), " ", paste0("-",tmp3[x,8])))
  }
  return(paste0(paste0(tmp3[x,8], "-"), " ", paste0("-",tmp3[x,2])))
})
tmp4 <-unlist(c(tmp[,14], tmp2[,14], tmp3[,14]))
tmp <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[1]
}))
tmp2 <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[2]
}))
tmp4 <- data.frame(sgRNA = str_remove(c(tmp, tmp2), "[-]"), 
                   Cmp = c(tmp, tmp2))
tmp4 <- tmp4[!duplicated(tmp4$sgRNA),]
table(sgRNA$sgRNA %in% tmp4$sgRNA)
high_cor_sample <- total_pair_sg_cor[total_pair_sg_cor$cor > 0.8,]
indel1_pct <- list()
region <- c(-10: 10)
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
  sg1 <- tmp4$Cmp[tmp4$sgRNA == sgRNA$sgRNA[sgRNA$id2 == ids[1]]]
  sg2 <- tmp4$Cmp[tmp4$sgRNA == sgRNA$sgRNA[sgRNA$id2 == ids[2]]]
  use <- "sg1"
  if(peak_indel == -1){
    if(!is.na(unlist(str_match(name, "last")))){
      if(str_sub(sg1, 1, 1) == "-"){
        use <- "sg1"
      }
      else{
        use <- "sg2"
      }
    }
    else{
      if(str_sub(sg1, 1, 1) == "-" | str_sub(sg1, 20, 20) == "-"){
        use <- "sg1"
      }
      else{
        use <- "sg2"
      }
    }
  }
  else{
    if(!is.na(unlist(str_match(name, "last")))){
      if(str_sub(sg1, 1, 1) == "-"){
        use <- "sg2"
      }
      else{
        use <- "sg1"
      }
    }
    else{
      if(str_sub(sg1, 1, 1) == "-" | str_sub(sg1, 20, 20) == "-"){
        use <- "sg2"
      }
      else{
        use <- "sg1"
      }
    }
  }
  indel1_pct[[name]] <- data.frame(Percent = ifelse(use == "sg1",
                                                    sel_sg1[sel_sg1[,1] == peak_indel,3],sel_sg2[sel_sg2[,1] == peak_indel,3]), 
                                   ID = paste(ifelse(use == "sg1", ids[1], ids[2]), 
                                               "of", paste0(ids, collapse = "-"), sep = " "), peak = peak_indel, used = use)
}
indel1_pct <- data.frame(do.call(rbind, indel1_pct))
pdf("~/data/project/ear_project/gene_therapy_ll/Result/0AA_restore_result.pdf", width = 10, height = 6)
ggplot(indel1_pct) + geom_bar(aes(x = ID, y = Percent), stat = "identity")+
geom_text( 
          aes(x = ID, y = Percent,label = round(Percent, 2)), vjust = -0.1) + 
  ylab("0AA restore % in indel reads")+xlab("") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), 
        text = element_text(size=15))
dev.off()
