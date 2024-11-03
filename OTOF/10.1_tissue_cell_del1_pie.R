seq_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/FunctioalAnalysisInput.xlsx", 3)
print(load("~/Nutstore Files/Tobin/Previous/tota_log.rda"))
library(stringr)
library(Biostrings)
source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")


total_used_seq <- lapply(total_log, function(x){
  log_name <- names(x)
  tmp <- unlist(lapply(x, function(y){
    y <- unlist(strsplit(y, "[ ]"))
    y[which(y == "-a") + 1]
  }))
  sg <- unlist(lapply(x, function(y){
    y <- unlist(strsplit(y, "[ ]"))
    y[which(y == "-g") + 1]
  }))
  data.frame(result = log_name, seq = tmp, sg = sg)
})
total_used_seq <- data.frame(do.call(rbind, total_used_seq))
table(duplicated(total_used_seq$result))
table(duplicated(paste0(total_used_seq$result,total_used_seq$seq)))
total_used_seq <- total_used_seq[!duplicated(total_used_seq$result),]
total_used_seq$seq <- str_to_lower(total_used_seq$seq)
seq_info$WT.CDS <- str_to_lower(seq_info$WT.CDS)
total_used_seq$result <- str_remove(total_used_seq$result, "CRISPResso_on_")

total_used_seq <- total_used_seq[total_used_seq$result %in% seq_info$Experiment.Number,]

table(unlist(lapply(1 : nrow(total_used_seq), function(i){
  seq <- total_used_seq$seq[i]
  cmpseq <- seq_info$WT.CDS[seq_info$Experiment.Number == total_used_seq$result[i]]
  cmpseq_mut <- seq_info$Mutant.CDS[seq_info$Experiment.Number == total_used_seq$result[i]]
  # str_locate(seq, cmpseq)[1]
  # str_length(cmpseq) - str_length(cmpseq_mut)
  # return(str_length(seq) >= str_length(cmpseq))
})))
total_used_seq$mutSeq <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  seq <- total_used_seq$seq[i]
  cmpseq <- seq_info$WT.CDS[seq_info$Experiment.Number == total_used_seq$result[i]]
  cmpseq_mut <- seq_info$Mutant.CDS[seq_info$Experiment.Number == total_used_seq$result[i]]
  str_to_lower(str_replace(seq, cmpseq, cmpseq_mut))
}))

total_used_seq$cds <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  seq <- total_used_seq$seq[i]
  seq_info$WT.CDS[seq_info$Experiment.Number == total_used_seq$result[i]]
  
}))

###是切割位点的前一位 注意一下
total_used_seq$cutSite <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  seq <- total_used_seq$seq[i]
  grna_cri <- str_to_lower(total_used_seq$sg[i])
  grna <- str_to_lower(seq_info$GRNA[seq_info$Experiment.Number == total_used_seq$result[i]])
  if(grna != grna_cri){
    print(total_used_seq$result[i])
  }
  rev <- F
  if(is.na(str_locate(seq, grna)[1])){
    rev <- T
    grna <- str_to_lower(as.character(reverseComplement(DNAString(grna))))
  }
  if(is.na(str_locate(seq, grna)[1])){
    while(is.na(str_locate(seq, grna)[1])){
      grna <- str_sub(grna, 2)
    }
  }
  pos <- str_locate(seq, grna)
  if(rev){
    return(pos[1] + 2)
  }
  else{
    return(pos[2] - 3)
  }
}))
total_used_seq$mutType <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  
  grna <- seq_info$Locus[seq_info$Experiment.Number == total_used_seq$result[i]]
  str_sub(grna, -2, -2)
}))
total_used_seq$mutNT <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  
  grna <- seq_info$Locus[seq_info$Experiment.Number == total_used_seq$result[i]]
  str_sub(grna, -1, -1)
}))
total_used_seq$mutPos <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  seq <- unlist(strsplit(total_used_seq$seq[i], "*"))
  seq_mut <- str_to_lower(unlist(strsplit(total_used_seq$mutSeq[i], "*")))
  cutPos <- total_used_seq$cutSite[i]
  type <- total_used_seq$mutType[i]
  mut <- str_to_lower(total_used_seq$mutNT[i])
  for(i in 1 : length(seq)){
    if(seq[i] != seq_mut[i]){
      #找到第一个出现差异的位点，接下来把这个差异位点尽可能往切割位点贴近
      if(type == 'd'){
        mutNT <- seq[i]
        if(mutNT != mut){
          print("WTF")
        }
        if(i  > cutPos){
          while(seq_mut[i - 1] == mutNT & i > cutPos){
            i <- i - 1
          }
        }
        return(i)
      }
      if(type == 'i'){
        mutNT <- seq_mut[i]
        if(mutNT != mut){
          print("WTF")
        }
        if(i  > cutPos){
          while(seq[i - 1] == mutNT & i > cutPos){
            i <- i - 1
          }
        }
        return(i)
      }
      
    }
  }
}))
total_used_seq$mutNT_find <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  seq <- unlist(strsplit(total_used_seq$seq[i], "*"))
  seq_mut <- str_to_lower(unlist(strsplit(total_used_seq$mutSeq[i], "*")))
  cutPos <- total_used_seq$cutSite[i]
  type <- total_used_seq$mutType[i]
  mut <- str_to_lower(total_used_seq$mutNT[i])
  for(i in 1 : length(seq)){
    if(seq[i] != seq_mut[i]){
      #找到第一个出现差异的位点，接下来把这个差异位点尽可能往切割位点贴近
      if(type == 'd'){
        mutNT <- seq[i]
      }
      if(type == 'i'){
        mutNT <- seq_mut[i]
      }
      return(mutNT)
    }
  }
}))

total_used_seq$strand <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  seq <- total_used_seq$seq[i]
  grna_cri <- str_to_lower(total_used_seq$sg[i])
  grna <- str_to_lower(seq_info$GRNA[seq_info$Experiment.Number == total_used_seq$result[i]])
  if(grna != grna_cri){
    print(total_used_seq$result[i])
  }
  rev <- F
  if(is.na(str_locate(seq, grna)[1])){
    rev <- T
    grna <- str_to_lower(as.character(reverseComplement(DNAString(grna))))
  }
  if(rev){
    return("-")
  }
  return("+")
}))

total_used_seq$mutNT_find <- str_to_upper(total_used_seq$mutNT_find)
tmp <- total_used_seq[total_used_seq$mutNT_find != total_used_seq$mutNT,]

load(file="~/Nutstore Files/Tobin/Previous//old_rep123_edit.rda")
total_old_edit_table <- list(Rep1 = old_rep1_edit, Rep2 = old_rep2_edit, Rep3 = old_rep3_edit)

total_old_edit_table_rev <- list()
for(name in names(total_old_edit_table)){
  rep <- total_old_edit_table[[name]]
  
  for(sp in names(rep)){
    if(!sp %in% names(total_old_edit_table_rev)){
      total_old_edit_table_rev[[sp]] <- list()
    }
    total_old_edit_table_rev[[sp]][[name]] <- rep[[sp]]
  }
}

names(total_old_edit_table_rev) <- str_remove(names(total_old_edit_table_rev), "CRISPResso_on_")


total_used_seq_tissue <- total_used_seq[grep("051", total_used_seq$result),]
total_used_seq_tissue <- total_used_seq_tissue[total_used_seq_tissue$mutType == 'i',]

noAA_change_seq <- list()

for(i in 1 : nrow(total_used_seq_tissue)){
  id <- total_used_seq_tissue$result[i]
  edit_tables <- total_old_edit_table_rev[[id]]
  edit_table <- lapply(edit_tables, function(edit_table){
    edit_table <- edit_table[edit_table$n_inserted + edit_table$n_deleted != 0,]
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
    edit_table <- edit_table[edit_table$n_deleted == 1 & 
                               edit_table$n_inserted == 0& 
                               edit_table$n_mutated <= 1,]
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}





seqs_table <- list()
for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  
  edit_table <- data.frame(do.call(rbind, edit_table))
  if(nrow(edit_table) == 0)
    next
  if(sum(edit_table[,8]) < 0.1)
    next
  edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, edit_table$Aligned_Sequence), function(x){
    x[,7] <- mean(x[,7])
    x[,8] <- mean(x[,8])
    x[,9] <- mean(x[,9])
    x[1,]
  })))
  edit_table$id <- id
  seqs_table[[id]] <- edit_table[,c(1,2, 8, 9, 10)]
}

seqs_table <- data.frame(do.call(rbind, seqs_table))

# ToNX::write_tb(seqs_table, 
# file="~/data/project/ear_project/gene_therapy_ll/Previews/Result/seqs_table_Ins1_all.txt")



seqs_table_split <- split(seqs_table, seqs_table$id)

seqs_table_stat <- lapply(names(seqs_table_split), function(x){
  seq_table <- seqs_table_split[[x]]
  seq_table <- seq_table[order(seq_table$Pct,decreasing = T),]
  dis <- c(total_used_seq_tissue$mutPos - total_used_seq_tissue$cutSite)[total_used_seq_tissue$result == x]
  cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
  seq_table$dis <- ifelse(dis <= 0, dis - 1, dis)
  mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
  strand <- total_used_seq_tissue$strand[total_used_seq_tissue$result == x]
  mutNT <- total_used_seq_tissue$mutNT[total_used_seq_tissue$result == x]
  #insert NT will be cap 
  mutseq <- str_to_lower(total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x])
  sg <- str_to_lower(total_used_seq_tissue$sg[total_used_seq_tissue$result == x])
  wtseq <-  unlist(strsplit(total_used_seq_tissue$seq[total_used_seq_tissue$result == x], "*"))
  seq_table$delPos <- unlist(lapply(1 : nrow(seq_table), function(i){
    editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
    refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
    editpos <- which(editseq == '-')
    editpos <- editpos[which.min(abs(editpos - 20))]
  }))
  
  seq_table$diff <- unlist(lapply(1 : nrow(seq_table), function(i){
    editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
    refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
    editpos <- which(editseq == '-')
    editpos <- editpos[which.min(abs(editpos - 20))]
    editNT <- str_to_lower(refseq[editpos])
    # ref <- seq_table$Reference_Sequence[i]
    #算出和cutsite的距离
    editpos <- editpos - 20
    ref <- sg
    if(cutsite + 20 > length(wtseq)){
      ref <- str_sub(ref, 1, -(cutsite + 21 - length(wtseq)))
    }
    if(cutsite < 20){
      ref <- str_sub(ref, 20 - cutsite)
    }
    if(strand == '-'){
      ref <- as.character(reverseComplement(DNAString(ref)))
      # ref <- paste(c(wtseq[(cutsite - 5) : (cutsite - 3)], ref), collapse = "")
      #遇到了一个坑爹的情况 就是插入这个碱基刚好和sgRNA形成了一个shift这样下面match的时候是会考虑突变而非insert而出现错误
      #因此需要把这个碱基插回去
      if(mutpos > (cutsite - 2) & mutpos <= (cutsite + str_length(sg) - 3)){
        ref <- unlist(strsplit(ref, "*"))
        ref <- paste(c(ref[ 1 : (mutpos - (cutsite - 2))], mutNT, ref[ (mutpos - (cutsite - 2) + 1) : length(ref)]), collapse = "")
      }
    }else{
      # ref <- paste(c(ref, wtseq[(cutsite + 4) : (cutsite + 6)]), collapse = "")
      if(mutpos > (cutsite - str_length(sg) + 4) & mutpos <= (cutsite + 3)){
        ref <- unlist(strsplit(ref, "*"))
        ref <- paste(c(ref[ 1 : (mutpos - (cutsite - str_length(sg) + 4))], mutNT, ref[ (mutpos - (cutsite - str_length(sg) + 4) + 1) : length(ref)]), collapse = "")
      }
    }
    findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
                                        subject = mutseq, max.mismatch = 3, with.indels = T)
    mutseqs <- unlist(strsplit(mutseq, "*"))
    if(strand == "+"){
      deletepos <- editpos + findRef@ranges@start + findRef@ranges@width - 3
      mutseqs <- mutseqs[-deletepos]
    }else{
      deletepos <- editpos + findRef@ranges@start + 2
      mutseqs <- mutseqs[-deletepos]
    }
    return(sum(mutseqs != wtseq))
  }))
  
  seq_table
})
seqs_table_stat <- data.frame(do.call(rbind, seqs_table_stat))



###计算Del 左右的统计
wt_seqs <- lapply(unique(seqs_table$id), function(id){
  edit_table <- total_old_edit_table_rev[[id]][[1]]
  edit_table <- edit_table[edit_table$Unedited == "True",]
  getWtSeq(edit_table)
})
names(wt_seqs) <- unique(seqs_table$id)


del1_pct_table <- lapply(unique(seqs_table$id), function(id){
  ins1_table <- seqs_table[seqs_table$id == id,]
  wt_seq <- wt_seqs[[id]]
  tmp <- getStat2(ins1_table, F)
  tmp$id <- id
  tmp
})
names(del1_pct_table) <- unique(seqs_table$id)
del1_pct_table_pos <- lapply(del1_pct_table, function(tmp){
  tmp <- lapply(split(tmp, tmp$pos), function(x){
    x$Pct <- sum(x$Pct)
    x[1,]
  })
  tmp <- do.call(rbind, tmp)
}) 

besides_nt_stat <- list()
nt_mats <- lapply(del1_pct_table, function(tmp){
  plot_mat <- matrix(0, ncol = 17, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(20, 21),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  tmp_left <- tmp[tmp$pos <= 20,]
  tmp_right <- tmp[tmp$pos > 20,]
  
  #首先把15-20的算一下，15等于是-6～-5之间插入
  #20则是-1～1之间插入
  #之后算21-23的  21也是-1～1之间插入  所以得分开算一下
  if(nrow(tmp_left) > 0){
    for(i in 1 : nrow(tmp_left)){
      mat_pos <- (tmp_left$pos[i] - 14) * 2
      plot_mat[tmp_left$NT[i], mat_pos] <- plot_mat[tmp_left$NT[i], mat_pos] + tmp_left$Pct[i]
    }
  }
  
  if(nrow(tmp_right) > 0){
    for(i in 1 : nrow(tmp_right)){
      mat_pos <- (tmp_right$pos[i] - 14) * 2 - 2
      plot_mat[tmp_right$NT[i], mat_pos] <- plot_mat[tmp_right$NT[i], mat_pos] + tmp_right$Pct[i]
    }
  }
  
  plot_mat
  
})
for(id in names(nt_mats)){
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- wt_seqs[[id]]
  
  #还是过滤1一下的结果
  i <- 12
  leftnt <- tmp_seq_mat[20]
  rightnt <- tmp_seq_mat[21]
  
  strand <- total_used_seq_tissue[str_remove(total_used_seq_tissue$result, "[-].*") ==
                                    str_remove(id,"[-].*"), "strand"]
  if(strand == "-"){
    leftnt <- as.character(Biostrings::reverseComplement(DNAString(leftnt)))
    rightnt <- as.character(Biostrings::reverseComplement(DNAString(rightnt)))
    tmp_nt <- leftnt
    leftnt <- rightnt
    rightnt <- tmp_nt
    rownames(tmp_pct_mat) <- as.character(Biostrings::reverseComplement(DNAStringSet(rownames(tmp_pct_mat))))
  }
  if(!leftnt %in% names(besides_nt_stat)){
    besides_nt_stat[[leftnt]] <- list()
  }
  tmp_besides_nt_stat <- besides_nt_stat[[leftnt]]
  if(!rightnt %in% names(tmp_besides_nt_stat)){
    tmp_nt_stat <- c(0,0,0,0)
    names(tmp_nt_stat) <- c("A", "T", "C", "G")
    tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
  }
  tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
  
  for(nt in rownames(tmp_pct_mat)){
    tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_pct_mat[nt, i]
  }
  tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
  besides_nt_stat[[leftnt]] <- tmp_besides_nt_stat
  
  
}
output_table <- lapply(names(besides_nt_stat), function(nt){
  x <- besides_nt_stat[[nt]]
  tmp <- lapply(x, function(y){
    y <- unlist(y)
    y / sum(y) * 100
  })
  tmp <- data.frame(do.call(rbind, tmp))
  tmp$rightNT <- rownames(tmp)
  tmp$leftNT <- nt
  tmp
})
output_table <- data.frame(do.call(rbind, output_table))
output_table <- output_table[,c(6, 5, 1,2,3,4)]
openxlsx::write.xlsx(output_table, 
                     file="~/Nutstore Files/Tobin/Previous/tissue_del1_left_right_nt_pie_data_fixerror.xlsx", 
                     rowNames=F, colNames=T)


###下面是细胞的
sample_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)
sample_info <- data.frame(cellid = sample_info$`Cell-w-W`, 
                         tissueid = sample_info$`Tissue-w-W(Mix)`,
                         mutation = sample_info$GRNAa)
noAA_change_seq <- list()

for(i in 1 : nrow(total_used_seq_tissue)){
  id <- total_used_seq_tissue$result[i]
  id <- sample_info$cellid[sample_info$tissueid == id]
  edit_tables <- total_old_edit_table_rev[[id]]
  edit_table <- lapply(edit_tables, function(edit_table){
    edit_table <- edit_table[edit_table$n_inserted + edit_table$n_deleted != 0,]
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
    edit_table <- edit_table[edit_table$n_deleted == 1 & 
                               edit_table$n_inserted == 0& 
                               edit_table$n_mutated <= 1,]
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}

seqs_table <- list()
for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  
  edit_table <- data.frame(do.call(rbind, edit_table))
  if(nrow(edit_table) == 0)
    next
  if(sum(edit_table[,8]) < 0.1)
    next
  edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, edit_table$Aligned_Sequence), function(x){
    x[,7] <- mean(x[,7])
    x[,8] <- mean(x[,8])
    x[,9] <- mean(x[,9])
    x[1,]
  })))
  edit_table$id <- id
  seqs_table[[id]] <- edit_table[,c(1,2, 8, 9, 10)]
}

seqs_table <- data.frame(do.call(rbind, seqs_table))

# ToNX::write_tb(seqs_table, 
# file="~/data/project/ear_project/gene_therapy_ll/Previews/Result/seqs_table_Ins1_all.txt")



seqs_table_split <- split(seqs_table, seqs_table$id)


####接下来要计算来自组织的Left+Right预测Insert的结果
wt_seqs <- lapply(unique(seqs_table$id), function(id){
  edit_table <- total_old_edit_table_rev[[id]][[1]]
  edit_table <- edit_table[edit_table$Unedited == "True",]
  if(nrow(edit_table) == 0){
    edit_table <- total_old_edit_table_rev[[id]][[1]]
  }
  getWtSeq(edit_table)
})
names(wt_seqs) <- unique(seqs_table$id)


del1_pct_table <- lapply(unique(seqs_table$id), function(id){
  ins1_table <- seqs_table[seqs_table$id == id,]
  wt_seq <- wt_seqs[[id]]
  tmp <- getStat2(ins1_table, F)
  tmp$id <- id
  tmp
})
names(del1_pct_table) <- unique(seqs_table$id)
del1_pct_table_pos <- lapply(del1_pct_table, function(tmp){
  tmp <- lapply(split(tmp, tmp$pos), function(x){
    x$Pct <- sum(x$Pct)
    x[1,]
  })
  tmp <- do.call(rbind, tmp)
}) 

besides_nt_stat <- list()
nt_mats <- lapply(del1_pct_table, function(tmp){
  plot_mat <- matrix(0, ncol = 17, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(20, 21),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  tmp_left <- tmp[tmp$pos <= 20,]
  tmp_right <- tmp[tmp$pos > 20,]
  
  #首先把15-20的算一下，15等于是-6～-5之间插入
  #20则是-1～1之间插入
  #之后算21-23的  21也是-1～1之间插入  所以得分开算一下
  if(nrow(tmp_left) > 0){
    for(i in 1 : nrow(tmp_left)){
      mat_pos <- (tmp_left$pos[i] - 14) * 2
      plot_mat[tmp_left$NT[i], mat_pos] <- plot_mat[tmp_left$NT[i], mat_pos] + tmp_left$Pct[i]
    }
  }
  
  if(nrow(tmp_right) > 0){
    for(i in 1 : nrow(tmp_right)){
      mat_pos <- (tmp_right$pos[i] - 14) * 2 - 2
      plot_mat[tmp_right$NT[i], mat_pos] <- plot_mat[tmp_right$NT[i], mat_pos] + tmp_right$Pct[i]
    }
  }
  
  plot_mat
  
})
for(id in names(nt_mats)){
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- wt_seqs[[id]]
  
  #还是过滤1一下的结果
  i <- 12
  leftnt <- tmp_seq_mat[20]
  rightnt <- tmp_seq_mat[21]
  strand <- total_used_seq_tissue[str_remove(total_used_seq_tissue$result, "[-].*") ==
                                    str_remove(id,"[-].*"), "strand"]
  if(strand == "-"){
    leftnt <- as.character(Biostrings::reverseComplement(DNAString(leftnt)))
    rightnt <- as.character(Biostrings::reverseComplement(DNAString(rightnt)))
    tmp_nt <- leftnt
    leftnt <- rightnt
    rightnt <- tmp_nt
    rownames(tmp_pct_mat) <- as.character(Biostrings::reverseComplement(DNAStringSet(rownames(tmp_pct_mat))))
  }
  if(!leftnt %in% names(besides_nt_stat)){
    besides_nt_stat[[leftnt]] <- list()
  }
  tmp_besides_nt_stat <- besides_nt_stat[[leftnt]]
  if(!rightnt %in% names(tmp_besides_nt_stat)){
    tmp_nt_stat <- c(0,0,0,0)
    names(tmp_nt_stat) <- c("A", "T", "C", "G")
    tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
  }
  tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
  
  for(nt in rownames(tmp_pct_mat)){
    tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_pct_mat[nt, i]
  }
  tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
  besides_nt_stat[[leftnt]] <- tmp_besides_nt_stat
  
  
}
output_table <- lapply(names(besides_nt_stat), function(nt){
  x <- besides_nt_stat[[nt]]
  tmp <- lapply(x, function(y){
    y <- unlist(y)
    y / sum(y) * 100
  })
  tmp <- data.frame(do.call(rbind, tmp))
  tmp$rightNT <- rownames(tmp)
  tmp$leftNT <- nt
  tmp
})
output_table <- data.frame(do.call(rbind, output_table))
output_table <- output_table[,c(6, 5, 1,2,3,4)]
openxlsx::write.xlsx(output_table, 
                     file="~/Nutstore Files/Tobin/Previous/cell_del1_left_right_nt_pie_data_fixerror.xlsx", 
                     rowNames=F, colNames=T)






tissue_data <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/tissue_del1_left_right_nt_pie_data.xlsx")
cell_data <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/cell_del1_left_right_nt_pie_data.xlsx")

fake_plot_mat <- matrix(runif(8*4), nrow = 8, ncol = 4)
rownames(fake_plot_mat) <- c("A","A_cell", "T","T_cell",
                             "C","C_cell", "G", "G_cell")
colnames(fake_plot_mat) <- c("A", "T", "C", "G")


col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(col) <-  c("A", 'T','C','G')
ht <- Heatmap(fake_plot_mat,
              name = "Name",
              column_split = NULL,
              column_title=NULL,
              row_title = NULL,
              row_split = rep(c("", " ", "  ", "   "), c(2,2,2,2)),
              show_column_names = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              top_annotation = HeatmapAnnotation("colname" = anno_empty(height = unit(2, "cm"), border = F)),
              left_annotation = rowAnnotation("rowname" = anno_empty(height = unit(2, "cm"), border = F)),
              right_annotation = NULL,
              show_row_names = F,
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                leftnt <- rownames(fake_plot_mat)[i]
                rightnt <- colnames(fake_plot_mat)[j]
                leftnt_real <- str_remove(leftnt, "[_].*")
                if(leftnt == leftnt_real){
                  pct <- unlist(tissue_data[tissue_data$leftNT == leftnt_real & 
                                              tissue_data$rightNT == rightnt, c(3:6)])
                }
                else{
                  pct <- unlist(cell_data[cell_data$leftNT == leftnt_real & 
                                            cell_data$rightNT == rightnt, c(3:6)])
                }
                if(length(pct) != 0){
                  ToNX::grid.pie(pct,x = x,y = y, w= w, h=h, 
                                 colors = col,size.scales = 0.3,
                                 init.angle = 90)
                }
                
                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "DeleteNT", type = "grid", pch = 16,
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)
pdf("~/Nutstore Files/Tobin/Previous/tissue_cell_left_right_nt_pie_chart_groupPct_insMut.pdf", 
    width = 6, height = 10)
draw(ht, annotation_legend_list = lgd_list)
library(gridtext)
decorate_annotation("colname", slice = 1, {
  tg <- richtext_grob(gt_render(c("A", "T", "C", "G")), 
                      rot = 0, 
                      x = unit(1 : 4 / 4 - 0.5 / 4, "npc"),
                      y=unit(0, "npc"), hjust = 0.5, 
                      gp = gpar(fontsize = 50, 
                                col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
  grid.draw(tg)
  invisible(tg)
})
for(i in 1 : 4){
  decorate_annotation("rowname", slice = i, {
    tg <- richtext_grob(gt_render(c("A", "T", "C", "G")[i]), 
                        rot = 0, 
                        x = unit(0.8, "npc"),
                        y=unit(0.5, "npc"), hjust = 0.5, 
                        gp = gpar(fontsize = 50, 
                                  col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")[i]))
    grid.draw(tg)
    invisible(tg)
    grid.segments(1.5, 0.1, 1.5, 0.9, gp = gpar(lwd = 4))
  })
}

dev.off()
