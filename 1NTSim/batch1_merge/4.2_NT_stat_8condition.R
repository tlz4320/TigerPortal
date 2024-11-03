#total_mean 来自2.0

indel1_sample <- total_mean[unlist(lapply(total_mean, function(x){
  x[x[,1] == 0, 2] <- 0
  abs(as.numeric(x[which.max(x[,2]),1])) == 1
}))]


sgRNA_pair_remain_indel1 <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(indel1_sample)) == 2)
}))]

total_pair_sg_cor_indel1 <- total_pair_sg_cor[total_pair_sg_cor$id %in% names(sgRNA_pair_remain_indel1),]

high_cor_sample_indel1 <- total_pair_sg_cor_indel1[total_pair_sg_cor_indel1$cor > 0.8,]

peak_indel_sample <- list()
peak_indel_sample[["Ins1"]] <- list()
peak_indel_sample[["Del1"]] <- list()
for(name in high_cor_sample_indel1$id){
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
  if(peak_indel == 1){
    peak_indel_sample[["Ins1"]][[name]] <- ids
  }
  else{
    peak_indel_sample[["Del1"]][[name]] <- ids
  }
}


####
###output result to file





ins1_peak_samples <- unlist(peak_indel_sample$Ins1)
ins1_peak_pair <- data.frame(do.call(rbind, peak_indel_sample$Ins1))
nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_diffNT_Ins1_all.txt", sep = "\t")
colnames(nt_diff) <- c("Align", "Ref", "CmpRef", "CmpWT", "ID","Percent", "Pos","Pos2", "Diff", "editNT", "changeNTs")
nt_diff <- nt_diff[nt_diff$Pos2 %in% c(-1, 1),]
nt_diff$side <- unlist(lapply(1 : nrow(nt_diff), function(i){
  cmp <- nt_diff$CmpRef[i]
  cmppair <- nt_diff$CmpWT[i]
  cmp <- unlist(strsplit(cmp, "*"))
  cmppair <- unlist(strsplit(cmppair, "*"))
  if(which(cmp == '-') == 1 | which(cmppair == '-') == 1){
    return("First")
  }
  return("Last")
}))
nt_diff <- nt_diff[order(nt_diff$Percent, decreasing = T),]
nt_diff <- nt_diff[nt_diff$ID %in% ins1_peak_samples,]
nt_diff <- data.frame(do.call(rbind, lapply(split(nt_diff, nt_diff$ID), function(x){
  x <- x[order(x$Percent, decreasing = T),]
  x[1,]
})))



paie_nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_Ins1_paired_output.txt", sep = "\t")
colnames(paie_nt_diff) <- c("Align", "Ref","CmpRef", "CmpWT","ID","Percent",  "editNT")
paie_nt_diff <- data.frame(do.call(rbind, lapply(split(paie_nt_diff, paie_nt_diff$ID), function(x){
  x <- x[order(x$Percent, decreasing = T),]
  x[1,]
})))
paie_nt_diff <- paie_nt_diff[paie_nt_diff$ID %in% ins1_peak_samples,]
rownames(paie_nt_diff) <- paie_nt_diff$ID


nt_diff_pos1_last <-  nt_diff[nt_diff$Pos2 == 1 & nt_diff$side == 'Last',]
paired_id <- unlist(lapply(nt_diff_pos1_last$ID, function(x){
  index <- c(which(ins1_peak_pair$X1 == x), which(ins1_peak_pair$X2 == x))
  setdiff(unlist(ins1_peak_pair[index, ]), x)
}))
nt_diff_pos1_last <- cbind(nt_diff_pos1_last, paie_nt_diff[paired_id,])



nt_diff_pos_1_first <-  nt_diff[nt_diff$Pos2 == -1 & nt_diff$side == 'First',]
paired_id <- unlist(lapply(nt_diff_pos_1_first$ID, function(x){
  index <- c(which(ins1_peak_pair$X1 == x), which(ins1_peak_pair$X2 == x))
  setdiff(unlist(ins1_peak_pair[index, ]), x)
}))
nt_diff_pos_1_first <- cbind(nt_diff_pos_1_first, paie_nt_diff[paired_id,])


nt_diff_pos1_first <-  nt_diff[nt_diff$Pos2 == 1 & nt_diff$side == 'First',]
paired_id <- unlist(lapply(nt_diff_pos1_first$ID, function(x){
  index <- c(which(ins1_peak_pair$X1 == x), which(ins1_peak_pair$X2 == x))
  setdiff(unlist(ins1_peak_pair[index, ]), x)
}))
nt_diff_pos1_first <- cbind(nt_diff_pos1_first, paie_nt_diff[paired_id,])

nt_diff_pos_1_last <-  nt_diff[nt_diff$Pos2 == -1 & nt_diff$side == 'Last',]












del1_peak_samples <- unlist(peak_indel_sample$Del1)
del1_peak_pair <- data.frame(do.call(rbind, peak_indel_sample$Del1))
nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_diffNT_Del1_all.txt", sep = "\t")
colnames(nt_diff) <- c("Align", "Ref", "CmpRef", "CmpWT", "ID","Percent", "Pos","Pos2", "Diff", "editNT", "changeNTs")
nt_diff <- nt_diff[nt_diff$Pos2 %in% c(-1, 1),]
nt_diff$side <- unlist(lapply(1 : nrow(nt_diff), function(i){
  cmp <- nt_diff$CmpRef[i]
  cmppair <- nt_diff$CmpWT[i]
  cmp <- unlist(strsplit(cmp, "*"))
  cmppair <- unlist(strsplit(cmppair, "*"))
  if(which(cmp == '-') == 1 | which(cmppair == '-') == 1){
    return("First")
  }
  return("Last")
}))
nt_diff <- nt_diff[order(nt_diff$Percent, decreasing = T),]
nt_diff <- nt_diff[nt_diff$ID %in% del1_peak_samples,]
nt_diff <- data.frame(do.call(rbind, lapply(split(nt_diff, nt_diff$ID), function(x){
  x <- x[order(x$Percent, decreasing = T),]
  x[1,]
})))



paie_nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_Del1_paired_output.txt", sep = "\t")
colnames(paie_nt_diff) <- c("Align", "Ref","CmpRef", "CmpWT","ID","Percent",  "editNT")
paie_nt_diff <- data.frame(do.call(rbind, lapply(split(paie_nt_diff, paie_nt_diff$ID), function(x){
  x <- x[order(x$Percent, decreasing = T),]
  x[1,]
})))
paie_nt_diff <- paie_nt_diff[paie_nt_diff$ID %in% del1_peak_samples,]
rownames(paie_nt_diff) <- paie_nt_diff$ID


nt_diff_pos1_last <-  nt_diff[nt_diff$Pos2 == 1 & nt_diff$side == 'Last',]
paired_id <- unlist(lapply(nt_diff_pos1_last$ID, function(x){
  index <- c(which(del1_peak_pair$X1 == x), which(del1_peak_pair$X2 == x))
  setdiff(unlist(del1_peak_pair[index, ]), x)
}))
nt_diff_pos1_last <- cbind(nt_diff_pos1_last, paie_nt_diff[paired_id,])



nt_diff_pos_1_first <-  nt_diff[nt_diff$Pos2 == -1 & nt_diff$side == 'First',]
paired_id <- unlist(lapply(nt_diff_pos_1_first$ID, function(x){
  index <- c(which(del1_peak_pair$X1 == x), which(del1_peak_pair$X2 == x))
  setdiff(unlist(del1_peak_pair[index, ]), x)
}))
nt_diff_pos_1_first <- cbind(nt_diff_pos_1_first, paie_nt_diff[paired_id,])


nt_diff_pos1_first <-  nt_diff[nt_diff$Pos2 == 1 & nt_diff$side == 'First',]
paired_id <- unlist(lapply(nt_diff_pos1_first$ID, function(x){
  index <- c(which(del1_peak_pair$X1 == x), which(del1_peak_pair$X2 == x))
  setdiff(unlist(del1_peak_pair[index, ]), x)
}))
nt_diff_pos1_first <- cbind(nt_diff_pos1_first, paie_nt_diff[paired_id,])

nt_diff_pos_1_last <-  nt_diff[nt_diff$Pos2 == -1 & nt_diff$side == 'Last',]








