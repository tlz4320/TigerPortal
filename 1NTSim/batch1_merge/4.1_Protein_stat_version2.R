#total_pair_sg_cor在2.0_total_First_pop_peak_indel1.R
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
sgCmp <- tmp4
rm(tmp, tmp2, tmp3, tmp4)
table(sgRNA$sgRNA %in% sgCmp$sgRNA)
Ins1_select <- list()

for(name in total_pair_sg_cor$id){
  ids <- sgRNA_pair_remain[[name]]
  sel_sg1 <- total_mean[[ids[1]]]
  # sel_sg1 <- sel_sg1[sel_sg1$indel_size %in% region,]
  sel_sg2 <- total_mean[[ids[2]]]
  # sel_sg2 <- sel_sg2[sel_sg2$indel_size %in% region,]
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
  sel_sg2[,3] <- sel_sg2[,2] / sum(sel_sg2[,2])
  peak_indel <- 1
  sg1 <- sgCmp$Cmp[sgCmp$sgRNA == sgRNA$sgRNA[sgRNA$id2 == ids[1]]]
  sg2 <- sgCmp$Cmp[sgCmp$sgRNA == sgRNA$sgRNA[sgRNA$id2 == ids[2]]]
  use <- "sg1"

  if(!is.na(unlist(str_match(name, "last")))){
    if(str_sub(sg1, 1, 1) == "-"){
      use <- "sg2"
    }else{
      use <- "sg1"
    }
  }
  else{
    if(str_sub(sg1, 1, 1) == "-" | str_sub(sg1, 21, 21) == "-"){
      use <- "sg2"
    }else{
      use <- "sg1"
    }
  }
  
  Ins1_select[[name]] <- data.frame(Percent = ifelse(use == "sg1",
                                                       sel_sg1[sel_sg1[,1] == peak_indel,3],sel_sg2[sel_sg2[,1] == peak_indel,3]), 
                                      ID = paste(ifelse(use == "sg1", ids[1], ids[2]), 
                                                 "of", paste0(ids, collapse = "-"), sep = " "), 
                                      peak = peak_indel, used = use,
                                      used_sg = ifelse(use == "sg1", ids[1], 
                                                       ids[2]))
}
Ins1_select <- data.frame(do.call(rbind, Ins1_select))


###把Indel1的另一个选出来
Del1_select <- list()
region <- c(-10: 10)
for(name in total_pair_sg_cor$id){
  ids <- sgRNA_pair_remain[[name]]
  sel_sg1 <- total_mean[[ids[1]]]
  # sel_sg1 <- sel_sg1[sel_sg1$indel_size %in% region,]
  sel_sg2 <- total_mean[[ids[2]]]
  # sel_sg2 <- sel_sg2[sel_sg2$indel_size %in% region,]
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
  sel_sg2[,3] <- sel_sg2[,2] / sum(sel_sg2[,2])
  peak_indel <- -1
  sg1 <- sgCmp$Cmp[sgCmp$sgRNA == sgRNA$sgRNA[sgRNA$id2 == ids[1]]]
  sg2 <- sgCmp$Cmp[sgCmp$sgRNA == sgRNA$sgRNA[sgRNA$id2 == ids[2]]]
  use <- "sg1"
  if(!is.na(unlist(str_match(name, "last")))){
    if(str_sub(sg1, 1, 1) == "-"){
      use <- "sg1"
    }
    else{
      use <- "sg2"
    }
  }
  else{
    if(str_sub(sg1, 1, 1) == "-" | str_sub(sg1, 21, 21) == "-"){
      use <- "sg1"
    }
    else{
      use <- "sg2"
    }
  }

  Del1_select[[name]] <- data.frame(Percent = ifelse(use == "sg1",
                                                           sel_sg1[sel_sg1[,1] == peak_indel,3],sel_sg2[sel_sg2[,1] == peak_indel,3]), 
                                          ID = paste(ifelse(use == "sg1", ids[1], ids[2]), 
                                                     "of", paste0(ids, collapse = "-"), sep = " "), 
                                          peak = peak_indel, used = use,
                                          used_sg = ifelse(use == "sg1", ids[1], 
                                                           ids[2]))
}
Del1_select <- data.frame(do.call(rbind, Del1_select))




####Read sg28 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/batch1/")
samples <- list.files(pattern = "^S")
total_edit_table <- list()
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  sgs <- list.files(pattern = "^sg")
  for(sg in sgs){
    setwd(sg)
    setwd("CRISPResso_on_nhej")
    filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
    filename <- filename[which.max(str_length(filename))]
    edit_table <- read.table(filename, 
                             sep="\t", header = T, comment.char = "")
    total_edit_table[[sample]][[sg]] <- edit_table
    setwd("../..")
  }
  setwd("..")
  rm(edit_table, sg, sgs)
}

####Read batch1_1 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
sgRNA_tmp <- read.table("sg_info.txt")
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    if(sample == "CRISPRessoPooled_on_B29" & sg %in% c("Sg_12_79", 
                                                       "Sg_12_80", 
                                                       "Sg_12_81")){
      print("rm sample")
      next
    }
    else{
      setwd(paste0("CRISPResso_on_",sg))
      filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
      filename <- filename[which.max(str_length(filename))]
      edit_table <- read.table(filename, 
                               sep="\t", header = T, comment.char = "")
      total_edit_table[[sample]][[sg]] <- edit_table
      setwd("..")
    }
    
  }
  setwd("..")
}
###Read batch1_2 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
sgRNA_tmp <- read.table("sg_info.txt")
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    if(sg %in% c("Sg_23_159", "Sg_19_133")){
      print("rm sample")
      next
    }
    if(sample == "CRISPRessoPooled_on_A99" & sg %in% c("Sg_23_157")){
      print("rm sample")
      next
    }
    if(sample == "CRISPRessoPooled_on_B99" & sg %in% c("Sg_20_137", 
                                                       "Sg_20_138", 
                                                       "Sg_20_139")){
      print("rm sample")
      next
    }
    if(sample == "CRISPRessoPooled_on_C99" & sg %in% c("Sg_19_132")){
      print("rm sample")
      next
    }
    setwd(paste0("CRISPResso_on_",sg))
    filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
    filename <- filename[which.max(str_length(filename))]
    edit_table <- read.table(filename, 
                             sep="\t", header = T, comment.char = "")
    total_edit_table[[sample]][[sg]] <- edit_table
    setwd("..")
  }
  setwd("..")
}
rm(sample, samples, sg, sgRNA_tmp, edit_table)


sgRNA <- read.xlsx("~/Nutstore Files/Tobin/First1NT/2024_1_12_integrated_design_result.xlsx")
sgRNA$id2 <- paste("Sg", sgRNA$new_id, sep = "-")
sgRNA$id2 <- str_replace_all(sgRNA$id2, "-", "_")



tmp <- unique(unlist(lapply(total_edit_table, names)))
total_edit_table_rev <- lapply(tmp, function(x){
  res <- list()
  for(name in names(total_edit_table)){
    if(x %in% names(total_edit_table[[name]])){
      res[[name]] <- total_edit_table[[name]][[x]]
    }
  }
  res
})
tmp[grep("split", tmp)] <- unlist(lapply(tmp[grep("split", tmp)], function(x){
  x <- str_remove(x, "_split")
  x <- as.numeric(str_remove(x, "sg"))
  sgRNA$id2[x]
}))
names(total_edit_table_rev) <- tmp

####get sequence after process
##首先选出Ins1的结果


total_edit_table_rev_indel1 <- total_edit_table_rev[Ins1_select$used_sg]
sgCmp_id <- merge(sgCmp, sgRNA, by="sgRNA")
sgRNA_sel <- sgRNA[sgRNA$id2 %in% Ins1_select$used_sg,]
noAA_change_seq <- list()
seqs_table <- list()
for(i in 1 : nrow(sgRNA_sel)){
  id <- sgRNA_sel$id2[i]
  edit_tables <- total_edit_table_rev[[id]]
  wt_seq <- Ins1_select$ID[Ins1_select$used_sg == id]
  wt_seq <- str_remove_all(wt_seq, id)
  wt_seq <- str_remove(str_remove(wt_seq, " of "), "-")
  wt_seq <- sgCmp_id$Cmp[sgCmp_id$id2 == wt_seq][1]
  sg <- sgCmp_id$Cmp[sgCmp_id$id2 == id][1]
  edit_table <- lapply(edit_tables, function(edit_table){
    edit_table <- edit_table[edit_table$Unedited == "False",]
    edit_table$Pct <- edit_table[,7] / sum(edit_table[,7]) * 100
    edit_table <- lapply(split(edit_table, edit_table$Aligned_Sequence), 
                         function(x){
                           if(nrow(x) == 1){
                             return(x)
                           }
                           x[,9] <- sum(x[,9])
                           return(x[1,])
                         })
    edit_table <- data.frame(do.call(rbind, edit_table))
    edit_table <- edit_table[order(edit_table$Pct, decreasing = T),]
    ##选出只有Ins1的结果
    edit_table <- edit_table[edit_table$n_deleted == 0 & 
                               edit_table$n_inserted == 1& 
                               edit_table$n_mutated <= 1,]
    edit_table$sgCmp <- sg
    edit_table$wtCmp <- wt_seq
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}
seqs_table <- list()
bad_pct <- 0
good_pct <- 0
for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  peak <- 1
  close_edit <- list()
  for(name in names(edit_table)){
    tmp <- edit_table[[name]]
    check_seq <- tmp[,2]
    if(peak < 0){
      check_seq <- tmp[,1]
    }
    poss <- unlist(lapply(check_seq, function(x){
      pos <- str_locate_all(x, "-")[[1]][,1]
      pos[which.min(abs(pos - 20))]
    }))
    close_edit[[name]] <- tmp[abs(poss - 20) <= 3,]
  }
  edit_table <- data.frame(do.call(rbind, close_edit))
  edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, edit_table$Aligned_Sequence), function(x){
    x[,7] <- mean(x[,7])
    x[,8] <- mean(x[,8])
    x[,9] <- mean(x[,9])
    x[1,]
  })))
  edit_table <- edit_table[!duplicated(edit_table$Aligned_Sequence),]
  edit_table$id <- id
  seqs_table[[id]] <- edit_table[,c(1,2, 10, 11, 12, 9)]
}

seqs_table <- data.frame(do.call(rbind, seqs_table))

ToNX::write_tb(seqs_table, 
               file="~/data/project/ear_project/gene_therapy_ll//seqs_table_Ins1_all.txt")

###接下来得到配对的Ins1的结果，用来做后面的预测
Ins1_select$pair_sg <- unlist(lapply(1:nrow(Ins1_select), function(i){
  sg <- Ins1_select$used_sg[i]
  id <- Ins1_select$ID[i]
  pair_sg <- str_remove_all(id, sg)
  str_remove(str_remove(pair_sg, " of "), "-")
}))
total_edit_table_rev_indel1 <- total_edit_table_rev[Ins1_select$pair_sg]
sgCmp_id <- merge(sgCmp, sgRNA, by="sgRNA")
sgRNA_sel <- sgRNA[sgRNA$id2 %in% Ins1_select$pair_sg,]
noAA_change_seq <- list()
seqs_table <- list()
for(i in 1 : nrow(sgRNA_sel)){
  id <- sgRNA_sel$id2[i]
  edit_tables <- total_edit_table_rev[[id]]
  wt_seq <- Ins1_select$used_sg[Ins1_select$pair_sg == id]
  wt_seq <- sgCmp_id$Cmp[sgCmp_id$id2 == wt_seq][1]
  sg <- sgCmp_id$Cmp[sgCmp_id$id2 == id][1]
  edit_table <- lapply(edit_tables, function(edit_table){
    edit_table <- edit_table[edit_table$Unedited == "False",]
    edit_table$Pct <- edit_table[,7] / sum(edit_table[,7]) * 100
    edit_table <- lapply(split(edit_table, edit_table$Aligned_Sequence), 
                         function(x){
                           if(nrow(x) == 1){
                             return(x)
                           }
                           x[,9] <- sum(x[,9])
                           return(x[1,])
                         })
    edit_table <- data.frame(do.call(rbind, edit_table))
    edit_table <- edit_table[order(edit_table$Pct, decreasing = T),]
    ##选出只有Ins1的结果
    edit_table <- edit_table[edit_table$n_deleted == 0 & 
                               edit_table$n_inserted == 1& 
                               edit_table$n_mutated <= 1,]
    edit_table$sgCmp <- sg
    edit_table$wtCmp <- wt_seq
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}
seqs_table <- list()
for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  peak <- 1
  close_edit <- list()
  for(name in names(edit_table)){
    tmp <- edit_table[[name]]
    check_seq <- tmp[,2]
    if(peak < 0){
      check_seq <- tmp[,1]
    }
    poss <- unlist(lapply(check_seq, function(x){
      pos <- str_locate_all(x, "-")[[1]][,1]
      pos[which.min(abs(pos - 20))]
    }))
    close_edit[[name]] <- tmp[abs(poss - 20) <= 3,]
  }
  edit_table <- data.frame(do.call(rbind, close_edit))
  edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, edit_table$Aligned_Sequence), function(x){
    x[,7] <- mean(x[,7])
    x[,8] <- mean(x[,8])
    x[,9] <- mean(x[,9])
    x[1,]
  })))
  edit_table <- edit_table[!duplicated(edit_table$Aligned_Sequence),]
  edit_table$id <- id
  seqs_table[[id]] <- edit_table[,c(1,2, 10, 11, 12, 9)]
}

seqs_table <- data.frame(do.call(rbind, seqs_table))

ToNX::write_tb(seqs_table, 
               file="~/data/project/ear_project/gene_therapy_ll//seqs_table_Ins1_paired.txt")



##首先选出Del1的结果


total_edit_table_rev_indel1 <- total_edit_table_rev[Del1_select$used_sg]
sgCmp_id <- merge(sgCmp, sgRNA, by="sgRNA")
sgRNA_sel <- sgRNA[sgRNA$id2 %in% Del1_select$used_sg,]
noAA_change_seq <- list()
seqs_table <- list()
for(i in 1 : nrow(sgRNA_sel)){
  id <- sgRNA_sel$id2[i]
  edit_tables <- total_edit_table_rev[[id]]
  wt_seq <- Del1_select$ID[Del1_select$used_sg == id]
  wt_seq <- str_remove_all(wt_seq, id)
  wt_seq <- str_remove(str_remove(wt_seq, " of "), "-")
  wt_seq <- sgCmp_id$Cmp[sgCmp_id$id2 == wt_seq][1]
  sg <- sgCmp_id$Cmp[sgCmp_id$id2 == id][1]
  edit_table <- lapply(edit_tables, function(edit_table){
    edit_table <- edit_table[edit_table$Unedited == "False",]
    edit_table$Pct <- edit_table[,7] / sum(edit_table[,7]) * 100
    edit_table <- lapply(split(edit_table, edit_table$Aligned_Sequence), 
                         function(x){
                           if(nrow(x) == 1){
                             return(x)
                           }
                           x[,9] <- sum(x[,9])
                           return(x[1,])
                         })
    edit_table <- data.frame(do.call(rbind, edit_table))
    edit_table <- edit_table[order(edit_table$Pct, decreasing = T),]
    #选出只有Del1的结果
    edit_table <- edit_table[edit_table$n_deleted == 1 & 
                               edit_table$n_inserted == 0 & 
                               edit_table$n_mutated <= 1,]
    edit_table$sgCmp <- sg
    edit_table$wtCmp <- wt_seq
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}
seqs_table <- list()
bad_pct <- 0
good_pct <- 0
for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  peak <- -1
  close_edit <- list()
  for(name in names(edit_table)){
    tmp <- edit_table[[name]]
    check_seq <- tmp[,2]
    if(peak < 0){
      check_seq <- tmp[,1]
    }
    poss <- unlist(lapply(check_seq, function(x){
      pos <- str_locate_all(x, "-")[[1]][,1]
      pos[which.min(abs(pos - 20))]
    }))
    close_edit[[name]] <- tmp[abs(poss - 20) <= 3,]
  }
  edit_table <- data.frame(do.call(rbind, close_edit))
  edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, edit_table$Aligned_Sequence), function(x){
    x[,7] <- mean(x[,7])
    x[,8] <- mean(x[,8])
    x[,9] <- mean(x[,9])
    x[1,]
  })))
  edit_table <- edit_table[!duplicated(edit_table$Aligned_Sequence),]
  edit_table$id <- id
  seqs_table[[id]] <- edit_table[,c(1,2, 10, 11, 12, 9)]
}

seqs_table <- data.frame(do.call(rbind, seqs_table))

ToNX::write_tb(seqs_table, 
               file="~/data/project/ear_project/gene_therapy_ll//seqs_table_Del1_all.txt")


###接下来得到配对的Del1的结果，用来做后面的预测
Del1_select$pair_sg <- unlist(lapply(1:nrow(Del1_select), function(i){
  sg <- Del1_select$used_sg[i]
  id <- Del1_select$ID[i]
  pair_sg <- str_remove_all(id, sg)
  str_remove(str_remove(pair_sg, " of "), "-")
}))
total_edit_table_rev_indel1 <- total_edit_table_rev[Del1_select$pair_sg]
sgCmp_id <- merge(sgCmp, sgRNA, by="sgRNA")
sgRNA_sel <- sgRNA[sgRNA$id2 %in% Del1_select$pair_sg,]
noAA_change_seq <- list()
seqs_table <- list()

for(i in 1 : nrow(sgRNA_sel)){
  id <- sgRNA_sel$id2[i]
  edit_tables <- total_edit_table_rev[[id]]
  wt_seq <- Del1_select$used_sg[Del1_select$pair_sg == id]
  wt_seq <- sgCmp_id$Cmp[sgCmp_id$id2 == wt_seq][1]
  sg <- sgCmp_id$Cmp[sgCmp_id$id2 == id][1]
  edit_table <- lapply(edit_tables, function(edit_table){
    edit_table <- edit_table[edit_table$Unedited == "False",]
    edit_table$Pct <- edit_table[,7] / sum(edit_table[,7]) * 100
    edit_table <- lapply(split(edit_table, edit_table$Aligned_Sequence), 
                         function(x){
                           if(nrow(x) == 1){
                             return(x)
                           }
                           x[,9] <- sum(x[,9])
                           return(x[1,])
                         })
    edit_table <- data.frame(do.call(rbind, edit_table))
    edit_table <- edit_table[order(edit_table$Pct, decreasing = T),]
    #选出只有Del1的结果
    edit_table <- edit_table[edit_table$n_deleted == 1 & 
                               edit_table$n_inserted == 0 & 
                               edit_table$n_mutated <= 1,]
    edit_table$sgCmp <- sg
    edit_table$wtCmp <- wt_seq
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}
seqs_table <- list()
for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  peak <- -1
  close_edit <- list()
  for(name in names(edit_table)){
    tmp <- edit_table[[name]]
    check_seq <- tmp[,2]
    if(peak < 0){
      check_seq <- tmp[,1]
    }
    poss <- unlist(lapply(check_seq, function(x){
      pos <- str_locate_all(x, "-")[[1]][,1]
      pos[which.min(abs(pos - 20))]
    }))
    close_edit[[name]] <- tmp[abs(poss - 20) <= 3,]
  }
  edit_table <- data.frame(do.call(rbind, close_edit))
  edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, edit_table$Aligned_Sequence), function(x){
    x[,7] <- mean(x[,7])
    x[,8] <- mean(x[,8])
    x[,9] <- mean(x[,9])
    x[1,]
  })))
  edit_table <- edit_table[!duplicated(edit_table$Aligned_Sequence),]
  edit_table$id <- id
  seqs_table[[id]] <- edit_table[,c(1,2, 10, 11, 12, 9)]
}

seqs_table <- data.frame(do.call(rbind, seqs_table))

ToNX::write_tb(seqs_table, 
               file="~/data/project/ear_project/gene_therapy_ll//seqs_table_Del1_paired.txt")





nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_diffNT_Ins1_all_new.txt", sep = "\t")
colnames(nt_diff) <- c("Align", "Ref", "CmpRef", "CmpWT", "ID","Percent", "Pos","Pos2", "Diff", "editNT", "changeNTs")


library(ggseqlogo)




##第二种做编辑产生插入后产生shift的碱基
top1_edit <-  lapply(split(nt_diff, nt_diff$ID), function(x){
  x[which.max(x$Percent),]
})
top1_edit <- data.frame(do.call(rbind, top1_edit))
editNts <- lapply(split(top1_edit, top1_edit$Pos2), function(x){
  #各种碱基需要算个加权平均
  ave_change <- mean(x$Diff)
  x <- x[x$Diff != 0,]
  tmp <- lapply(1 : nrow(x), function(i){
    data.frame(NT = unlist(strsplit(x$changeNTs[i], "*")), 
               Pct = x$Percent[i])
  })
  tmp <- data.frame(do.call(rbind, tmp))
  tmp$score <- tmp$Pct / sum(tmp$Pct) * ave_change
  tmp$Pos <- unique(x$Pos2)
  tmp
})
editNts <- data.frame(do.call(rbind, editNts))
editNts$Pos <- as.character(editNts$Pos)
plot_mat <- matrix(0, ncol = 8, nrow = 4)
rownames(plot_mat) <- c("A", "T", "C", "G")
colnames(plot_mat) <- c("-5","-4", "-3", "-2", "-1", "1", "2", "3")
for(i in  c("-5","-4", "-3", "-2", "-1", "1", "2", "3")){
  tmp <- editNts[editNts$Pos == i,]
  for(j in 1 : nrow(tmp)){
    plot_mat[tmp$NT[j], i] <- plot_mat[tmp$NT[j], i] + tmp$score[j]
  }
}

# Generate sequence logo
pdf("~/data/project/ear_project/gene_therapy_ll/Result/NT_change_Ins1_bulge.pdf", width = 6, height = 4)
ggplot() + geom_logo(data = plot_mat, method='custom', seq_type='dna') + 
  geom_vline(xintercept = 5.5, linetype="dashed") +
  scale_x_continuous(labels = c("-5","-4", "-3", "-2", "-1", "+1", "+2", "+3"), breaks = 1 : 8) +
  theme_logo() + ylab("Change NTs Counts") + xlab("Distance (target bulge site)") + ggtitle("Insert 1 Bulge") + 
  theme_classic2() + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25))
dev.off()



nt_diff <- read.table("~/data/project/ear_project/gene_therapy_ll/seqs_table_diffNT_Del1_all.txt", sep = "\t")
colnames(nt_diff) <- c("Align", "Ref", "CmpRef", "CmpWT", "ID","Percent", "Pos","Pos2", "Diff", "editNT", "changeNTs", "editBegin")


library(ggseqlogo)




##第二种做编辑产生插入后产生shift的碱基
top1_edit <-  lapply(split(nt_diff, nt_diff$ID), function(x){
  x[which.max(x$Percent),]
})
top1_edit <- data.frame(do.call(rbind, top1_edit))
editNts <- lapply(split(top1_edit, top1_edit$Pos2), function(x){
  #各种碱基需要算个加权平均
  ave_change <- mean(x$Diff)
  x <- x[x$Diff != 0,]
  tmp <- lapply(1 : nrow(x), function(i){
    data.frame(NT = unlist(strsplit(x$changeNTs[i], "*")), 
               Pct = x$Percent[i])
  })
  tmp <- data.frame(do.call(rbind, tmp))
  tmp$score <- tmp$Pct / sum(tmp$Pct) * ave_change
  tmp$Pos <- unique(x$Pos2)
  tmp
})
editNts <- data.frame(do.call(rbind, editNts))
editNts$Pos <- as.character(editNts$Pos)
plot_mat <- matrix(0, ncol = 8, nrow = 4)
rownames(plot_mat) <- c("A", "T", "C", "G")
colnames(plot_mat) <-  c("-5","-4", "-3", "-2", "-1", "1", "2", "3")
for(i in  c("-5","-4", "-3", "-2", "-1", "1", "2", "3")){
  tmp <- editNts[editNts$Pos == i,]
  for(j in 1 : nrow(tmp)){
    plot_mat[tmp$NT[j], i] <- plot_mat[tmp$NT[j], i] + tmp$score[j]
  }
}

# Generate sequence logo

pdf("~/data/project/ear_project/gene_therapy_ll/Result/NT_change_Del1_bulge.pdf", width = 6, height = 4)
ggplot() + geom_logo(data = plot_mat, method='custom', seq_type='dna') + 
  geom_vline(xintercept = 5.5, linetype="dashed") +
  scale_x_continuous(labels = c("-5","-4", "-3", "-2", "-1", "+1", "+2", "+3"), breaks = 1 : 8) +
  theme_logo() + ylab("Change NTs Counts") + xlab("Distance (target bulge site)") + ggtitle("Delete 1 Bulge") + 
  theme_classic2() + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25))
dev.off()
