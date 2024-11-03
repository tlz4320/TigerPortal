###需要统计inframe之中0AA的Identity变化的情况，
###之前在Volcano统计的比较简单 这次单独都拉出来进行比较
library(Biostrings)
library(stringr)
source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
seq_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/FunctioalAnalysisInput.xlsx", 3)
setwd("~/data/project/ear_project/gene_therapy_ll/Previews/")
reps <- list.files(pattern = "^Rep.*data")
total_log <- list()
for(rep in reps){
  setwd(rep)
  if(!rep %in% names(total_log)){
    total_log[[rep]] <- list()
  }
  files <- list.files(pattern = "CRISPResso")
  for(file in files){
    tmp <- read.table(paste0(file,"/CRISPResso_RUNNING_LOG.txt"), sep="\t", fill = T)
    total_log[[rep]][[file]] <- unlist(tmp[3,1])
  }
  setwd('..')
}




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
  str_replace(seq, cmpseq, cmpseq_mut)
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
  seq_mut <- unlist(strsplit(total_used_seq$mutSeq[i], "*"))
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
  seq_mut <- unlist(strsplit(total_used_seq$mutSeq[i], "*"))
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

load(file="Result/old_rep123_edit.rda")
print(load("~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data_newer_newer.rda"))

otof_tissue_sg12_15_rep1 <- otof_tissue_sg12_15_resote[grep("_1",names(otof_tissue_sg12_15_resote))]
otof_tissue_sg12_15_rep2 <- otof_tissue_sg12_15_resote[grep("_2",names(otof_tissue_sg12_15_resote))]
otof_tissue_sg12_15_rep3 <- otof_tissue_sg12_15_resote[grep("_3",names(otof_tissue_sg12_15_resote))]
otof_tissue_sg12_15_rep4 <- otof_tissue_sg12_15_resote[grep("_4",names(otof_tissue_sg12_15_resote))]
names(otof_tissue_sg12_15_rep1) <- paste0(str_to_lower(str_remove(names(otof_tissue_sg12_15_rep1), "[_].*")), '-tissue')
names(otof_tissue_sg12_15_rep2) <- paste0(str_to_lower(str_remove(names(otof_tissue_sg12_15_rep2), "[_].*")), '-tissue')
names(otof_tissue_sg12_15_rep3) <- paste0(str_to_lower(str_remove(names(otof_tissue_sg12_15_rep3), "[_].*")), '-tissue')
names(otof_tissue_sg12_15_rep4) <- paste0(str_to_lower(str_remove(names(otof_tissue_sg12_15_rep4), "[_].*")), '-tissue')



for(name in names(otof_tissue_sg12_15_rep1)){
  old_rep1_edit[[name]] <- otof_tissue_sg12_15_rep1[[name]]
}
for(name in names(otof_tissue_sg12_15_rep2)){
  old_rep2_edit[[name]] <- otof_tissue_sg12_15_rep2[[name]]
}
for(name in names(otof_tissue_sg12_15_rep3)){
  old_rep3_edit[[name]] <- otof_tissue_sg12_15_rep3[[name]]
}
old_rep4 <- otof_tissue_sg12_15_rep4
total_old_edit_table <- list(Rep1 = old_rep1_edit, Rep2 = old_rep2_edit, Rep3 = old_rep3_edit, 
                             Rep4 = old_rep4)

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
total_used_seq_tissue <- total_used_seq_tissue[total_used_seq_tissue$mutType == 'd',]
otof_sg12_15_info <- total_used_seq_tissue[rep(which(total_used_seq_tissue$result == "105-051"), 4),]
otof_sg12_15_info$result <- paste0("m", 12:15, "-tissue")
otof_sg12_15_info$sg <- str_to_upper(c("gcccactgccgttcggg", "tgcccactgccgttcgg",
                                       "gtgcccactgccgttcg", "cgtgcccactgccgttc"))
otof_sg12_15_info$cutSite <- c(41 : 44)
otof_sg12_15_info$strand <- "-"
total_used_seq_tissue <- data.frame(rbind(total_used_seq_tissue, otof_sg12_15_info))

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
    edit_table <- edit_table[edit_table$n_deleted == 0 & 
                               edit_table$n_inserted == 1& 
                               edit_table$n_mutated <= 1,]
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}




for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  remain_table <- list()
  for(rep in names(edit_table)){
    if(nrow(edit_table[[rep]]) == 0){
      next
    }
    if(sum(edit_table[[rep]]$Pct) < 1){
      next
    }
    remain_table[[rep]] <- edit_table[[rep]]
  }
  if(length(remain_table) == 0){
    print(id)
  }
  noAA_change_seq[[id]] <- remain_table
}



####接下来要计算来自组织的Left+Right预测Insert的结果
wt_seqs <- lapply(unique(seqs_table$id), function(id){
  edit_table <- total_old_edit_table_rev[[id]][[1]]
  edit_table <- edit_table[edit_table$Unedited == "True",]
  getWtSeq(edit_table)
})
names(wt_seqs) <- unique(seqs_table$id)


ins1_pct_table <- lapply(unique(seqs_table$id), function(id){
  ins1_table <- seqs_table[seqs_table$id == id,]
  wt_seq <- wt_seqs[[id]]
  tmp <- getStat2(ins1_table, T)
  tmp$id <- id
  tmp
})
names(ins1_pct_table) <- unique(seqs_table$id)
ins1_pct_table_pos <- lapply(ins1_pct_table, function(tmp){
  tmp <- lapply(split(tmp, tmp$pos), function(x){
    x$Pct <- sum(x$Pct)
    x[1,]
  })
  tmp <- do.call(rbind, tmp)
}) 

seqs_table_stat <- lapply(names(noAA_change_seq), function(x){
  seq_tables <- noAA_change_seq[[x]]
  lapply(seq_tables, function(seq_table){
    seq_table <- seq_table[order(seq_table$Pct,decreasing = T),]
    dis <- c(total_used_seq_tissue$mutPos - total_used_seq_tissue$cutSite)[total_used_seq_tissue$result == x]
    cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
    seq_table$dis <- ifelse(dis <= 0, dis - 1, dis)
    mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
    
    strand <- total_used_seq_tissue$strand[total_used_seq_tissue$result == x]
    #insert NT will be cap 
    mutseq <- str_to_lower(total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x])
    sg <- str_to_lower(total_used_seq_tissue$sg[total_used_seq_tissue$result == x])
    wtseq <-  unlist(strsplit(total_used_seq_tissue$seq[total_used_seq_tissue$result == x], "*"))
    
    seq_table$NT <- unlist(lapply(1 : nrow(seq_table), function(i){
      editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
      refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
      editpos <- which(refseq == '-')
      editpos <- editpos[which.min(abs(editpos - 20))]
      str_to_lower(editseq[editpos])
    }))
    seq_table$insPos <- unlist(lapply(1 : nrow(seq_table), function(i){
      editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
      refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
      editpos <- which(refseq == '-')
      editpos[which.min(abs(editpos - 20))]
    }))
    seq_table$diff <- unlist(lapply(1 : nrow(seq_table), function(i){
      editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
      refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
      editpos <- which(refseq == '-')
      editpos <- editpos[which.min(abs(editpos - 20))]
      editNT <- str_to_lower(editseq[editpos])
      ref <- seq_table$Reference_Sequence[i]
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
      } 
      findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
                                          subject = mutseq, max.mismatch = 3, with.indels = T)
      mutseqs <- unlist(strsplit(mutseq, "*"))
      if(strand == "+"){
        insertpos <- editpos + findRef@ranges@start + findRef@ranges@width - 4
        mutseqs <- c(mutseqs[1 : (insertpos - 1)], editNT, mutseqs[insertpos : length(mutseqs)])
      }else{
        insertpos <- editpos + findRef@ranges@start + 2
        mutseqs <- c(mutseqs[1 : (insertpos - 1)], editNT, mutseqs[insertpos : length(mutseqs)])
      }
      return(sum(mutseqs != wtseq))
    }))
    seq_table
  })
})
###改成了分开样本处理，接下来是要先根据位置+NT进行一个合并用于后续计算Se

names(seqs_table_stat) <- names(noAA_change_seq)
seqs_table_stat <- seqs_table_stat[names(seqs_table_stat) != "50-051"]
###789三列进行加和
seqs_table_merge <- lapply(seqs_table_stat, function(seq_tables){
  lapply(seq_tables, function(seq_table){
    pos_nt <- paste(seq_table$insPos, seq_table$NT, by="-")
    seq_table$groupPct <- seq_table$Pct / sum(seq_table$Pct) * 100
    seq_table_split <- lapply(split(seq_table, pos_nt), function(sp){
      sp[,7] <- sum(sp[,7])
      sp[,8] <- sum(sp[,8])
      sp[,9] <- sum(sp[,9])
      sp$groupPct <- sum(sp$groupPct)
      sp[1,]
    })
    data.frame(do.call(rbind, seq_table_split))
  })
})

seqs_table_mean <- lapply(names(seqs_table_merge), function(id){
  seq_tables <- seqs_table_merge[[id]]
  counts <- length(seq_tables)
  seq_table <- data.frame(do.call(rbind, seq_tables))
  pos_nt <- paste(seq_table$insPos, seq_table$NT, by="-")
  seq_table_split <- lapply(split(seq_table, pos_nt), function(sp){
    sp[,7] <- sum(sp[,7]) / counts
    sp[,8] <- sum(sp[,8])/ counts
    sp[,9] <- sum(sp[,9])/ counts
    sp$groupPct <- sum(sp$groupPct)/ counts
    sp[1,]
  })
  seq_table <- data.frame(do.call(rbind, seq_table_split))
  seq_table$id <- id
  seq_table
})
seqs_table_mean <- data.frame(do.call(rbind, seqs_table_mean))
openxlsx::write.xlsx(seqs_table_mean, file="~/Nutstore Files/Tobin/Previous/total_ins1_0AA_stat.xlsx",rowNames=F, colNames=T)
###然后基于diff进行合并后 计算se
seqs_table_se <- lapply(names(seqs_table_merge), function(id){
  seq_tables <- seqs_table_merge[[id]]
  counts <- length(seq_tables)
  seq_tables <- lapply(seq_tables, function(seq_table){
    seq_table_split <- lapply(split(seq_table, seq_table$diff), function(sp){
      sp[,7] <- sum(sp[,7])
      sp[,8] <- sum(sp[,8])
      sp[,9] <- sum(sp[,9])
      sp$groupPct <- sum(sp$groupPct)
      sp[1,]
    })
    data.frame(do.call(rbind, seq_table_split))
  })
  
  seq_table <- data.frame(do.call(rbind, seq_tables))
  
  seq_table$se <- NULL
  seq_table$groupSe <- NULL
  seq_table_split <- lapply(split(seq_table, seq_table$diff), function(sp){
    se <- sp[,9]
    groupse <- sp$groupPct
    if(length(se) < counts){
      se <- unlist(c(se, rep(0, counts - length(se))))
      groupse <- unlist(c(groupse, rep(0, counts - length(se))))
    }
    sp$se <- plotrix::std.error(unlist(se))
    sp$groupSe <- plotrix::std.error(unlist(groupse))
    sp[,7] <- sum(sp[,7]) / counts
    sp[,8] <- sum(sp[,8])/ counts
    sp[,9] <- sum(sp[,9])/ counts
    sp$groupPct <- sum(sp$groupPct)/ counts
    sp[1,]
  })
  seq_table <- data.frame(do.call(rbind, seq_table_split))
  seq_table$id <- id
  seq_table
})
names(seqs_table_se) <- names(seqs_table_merge)
seqs_table_se <- data.frame(do.call(rbind, seqs_table_se))
openxlsx::write.xlsx(seqs_table_se, file="~/Nutstore Files/Tobin/Previous/total_ins1_0AA_stat_Se.xlsx",rowNames=F, colNames=T)

seqs_table_mean <- lapply(names(seqs_table_merge), function(id){
  seq_tables <- seqs_table_merge[[id]]
  counts <- length(seq_tables)
  seq_tables <- lapply(seq_tables, function(seq_table){
    if(any(seq_table$diff > 6)){
      seq_table_6 <- seq_table[seq_table$diff > 6,]
      
      seq_table_split <- lapply(split(seq_table_6, seq_table_6$NT), function(sp){
        sp[,7] <- sum(sp[,7])
        sp[,8] <- sum(sp[,8])
        sp[,9] <- sum(sp[,9])
        sp$groupPct <- sum(sp$groupPct)
        sp[1,]
      })
      seq_table_6 <- data.frame(do.call(rbind, seq_table_split))
      #单独把大于6的拿出来
      #原因是因为本身我们关心6范围内的数据
      #另外如果把全范围画出来会非常乱
      seq_table_6$insPos <- 9999
      seq_table_6$diff <- 7
      seq_table <- data.frame(rbind(seq_table_6, seq_table[seq_table$diff <= 6,]))
    }
    seq_table
  })
  
  seq_table <- data.frame(do.call(rbind, seq_tables))
  pos_nt <- paste(seq_table$insPos, seq_table$NT, by="-")
  seq_table_split <- lapply(split(seq_table, pos_nt), function(sp){
    sp[,7] <- sum(sp[,7]) / counts
    sp[,8] <- sum(sp[,8])/ counts
    sp[,9] <- sum(sp[,9])/ counts
    sp$groupPct <- sum(sp$groupPct)/ counts
    sp[1,]
  })
  seq_table <- data.frame(do.call(rbind, seq_table_split))
  seq_table$id <- id
  seq_table
})
seqs_table_mean <- data.frame(do.call(rbind, seqs_table_mean))
tmp <- seqs_table_mean
tmp$diff[tmp$diff == 7] <- ">6"
openxlsx::write.xlsx(tmp, file="~/Nutstore Files/Tobin/Previous/total_ins1_0AA_stat_merge_large6.xlsx",rowNames=F, colNames=T)



seqs_table_se <- lapply(names(seqs_table_merge), function(id){
  seq_tables <- seqs_table_merge[[id]]
  counts <- length(seq_tables)
  seq_tables <- lapply(seq_tables, function(seq_table){
    seq_table$diff[seq_table$diff > 6] <- 7
    seq_table_split <- lapply(split(seq_table, seq_table$diff), function(sp){
      sp[,7] <- sum(sp[,7])
      sp[,8] <- sum(sp[,8])
      sp[,9] <- sum(sp[,9])
      sp$groupPct <- sum(sp$groupPct)
      sp[1,]
    })
    data.frame(do.call(rbind, seq_table_split))
  })
  
  seq_table <- data.frame(do.call(rbind, seq_tables))
  
  seq_table$se <- NULL
  seq_table$groupSe <- NULL
  seq_table_split <- lapply(split(seq_table, seq_table$diff), function(sp){
    se <- sp[,9]
    groupse <- sp$groupPct
    if(length(se) < counts){
      se <- unlist(c(se, rep(0, counts - length(se))))
      groupse <- unlist(c(groupse, rep(0, counts - length(se))))
    }
    sp$se <- plotrix::std.error(unlist(se))
    sp$groupSe <- plotrix::std.error(unlist(groupse))
    sp[,7] <- sum(sp[,7]) / counts
    sp[,8] <- sum(sp[,8])/ counts
    sp[,9] <- sum(sp[,9])/ counts
    sp$groupPct <- sum(sp$groupPct)/ counts
    sp[1,]
  })
  seq_table <- data.frame(do.call(rbind, seq_table_split))
  seq_table$id <- id
  seq_table
})
names(seqs_table_se) <- names(seqs_table_merge)
seqs_table_se <- data.frame(do.call(rbind, seqs_table_se))
tmp <- seqs_table_se
tmp$diff[tmp$diff == 7] <- ">6"
openxlsx::write.xlsx(tmp, file="~/Nutstore Files/Tobin/Previous/total_ins1_0AA_stat_Se_merge_large6.xlsx",rowNames=F, colNames=T)


load("~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda")
sp_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)
sp_info <- split_name <- data.frame(
  tissueid = c(sp_info$`Tissue-w-W(Mix)`,
               paste0("m", 12:15, "-tissue")),
  mutation = c(sp_info$GRNAa, paste0("OTOF-1236dC-g", 3:6)))



seqs_table_mean <- merge(sp_info, seqs_table_mean, by.x="tissueid", by.y="id")
seqs_table_se <- merge(sp_info, seqs_table_se, by.x="tissueid", by.y="id")
cols = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(cols) <- c("A", 'T','C','G')
cols2 <- c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")
names(cols2) <- c("-2", "-1", "0", "1", "2", "Others")

plot_list <- list()

for(id in sp_info$tissueid){
  if(!id %in% seqs_table_mean$tissueid)
    next
  plot_data <- inframe_result_processed[inframe_result_processed$id == id,]
  if(nrow(plot_data) != 6){
    append <- data.frame("indel" = 0,
                         "fq" = 0,
                         "pct" = 0,
                         "aa" = c("-2", "-1", "0", "1", "2", "Others"), 
                         "se" = 0,
                         "pct2" = 0,
                         "id" = id)
    append <- append[!append$aa %in% plot_data$aa,]
    plot_data <- data.frame(rbind(append, plot_data))
  }
  plot_data$aa  <- factor(plot_data$aa, 
                          levels = c("-2", "-1", "0", "1", "2", "Others"))
  inframe_plot <- ggplot(plot_data, aes(x = aa, y = pct2, fill = aa)) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymin=pct2-se, ymax=pct2+se), width=.2, color = "black") + 
    ylab("Percent in Indel") + xlab("number of AA change") + 
    scale_x_discrete(breaks = c("-2", "-1", "0", "1", "2", "Others")) + 
    scale_fill_manual(values = cols2) + 
    theme_bw() + theme(panel.grid = element_blank())
  plot_data_0aa <- seqs_table_mean[seqs_table_mean$tissueid == id,]
  plot_data_se <- seqs_table_se[seqs_table_se$tissueid == id,]
  if(!all(0 : 7 %in% plot_data_0aa$diff)){
    needappend <- setdiff(c(0 : 7), plot_data_0aa$diff)
    append <- lapply(needappend, function(d){
      tmp <- plot_data_0aa[1,]
      tmp$Pct <- tmp$groupPct <- 0
      tmp$diff <- d
      tmp
    })
    append <- data.frame(do.call(rbind, append))
    plot_data_0aa <- data.frame(rbind(append, plot_data_0aa))
    append <- lapply(needappend, function(d){
      tmp <- plot_data_se[1,]
      tmp$Pct <- tmp$groupPct <- 0
      tmp$diff <- d
      tmp$se <- 0
      tmp
    })
    append <- data.frame(do.call(rbind, append))
    plot_data_se <- data.frame(rbind(append, plot_data_se))
    
    
  }
  plot_data_0aa <- plot_data_0aa[plot_data_0aa$diff %in% 0 : 7,]
  plot_data_se <- plot_data_se[plot_data_se$diff %in% 0 : 7,]
  
  plot_data_0aa$diff <- factor(plot_data_0aa$diff, levels = c(0 : 7))
  plot_data_se$diff <- factor(plot_data_se$diff, levels = c(0 : 7))
  plot_data_0aa$NT <- stringr::str_to_upper(plot_data_0aa$NT)
  noaa_plot <- ggplot(plot_data_0aa, aes(x = diff, y = groupPct, fill = NT)) + 
    geom_bar(stat = "identity", position = "stack") + 
    geom_errorbar(data = plot_data_se, 
                  aes(x = diff, y = groupPct, 
                      ymin=groupPct-se, ymax=groupPct+se), width=.2, color = "black") + 
    ylab("Percent in 0AA") + xlab("number of NT identity change") + 
    scale_x_discrete(breaks = c(0:7), labels = c(0:6, ">6")) + 
    scale_fill_manual(values = cols) + 
    theme_bw() + theme(panel.grid = element_blank())
  dis <- unique(plot_data_0aa$dis)
  gene_name <- unique(plot_data_0aa$mutation)
  plot_list[[id]] <- (inframe_plot + ggtitle(gene_name)) * 
    (noaa_plot + ggtitle(paste0("Cut site distance:", dis)))
}



pdf("~/Nutstore Files/Tobin/Previous/total_ins1_tissue_inframe_0aa_stat.pdf", width = 8, height = 4)
for(plot in plot_list){
  print(plot)
}
dev.off()





