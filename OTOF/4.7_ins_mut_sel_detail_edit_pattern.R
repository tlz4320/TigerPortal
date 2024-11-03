###这个脚本依赖于4.4的统计结果
total_sample_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)
total_sample_info <- split_name <- data.frame(
  tissueid = total_sample_info$`Tissue-w-W(Mix)`,
  mutation = total_sample_info$GRNAa)
total_sample_info <- merge(total_sample_info, total_used_seq_tissue, by.x="tissueid", by.y="result")
total_sample_info <- total_sample_info[total_sample_info$tissueid %in% names(seqs_table_split),]
arround_len <- 6




edit_stat_of_distance <- lapply(names(seqs_table_split), function(x){
  seq_table <- seqs_table_split[[x]]
  seq_table <- seq_table[order(seq_table$Pct,decreasing = T),]
  mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
  cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
  dis <- c(total_used_seq_tissue$cutSite - total_used_seq_tissue$mutPos)[total_used_seq_tissue$result == x]
  
  seq_table$dis <- dis
  # mutseq <- total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x]
  wtseq <-  total_used_seq_tissue$seq[total_used_seq_tissue$result == x]
  seq_table$dis2mut <- unlist(lapply(1 : nrow(seq_table), function(i){
    editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
    refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
    editpos <- which(editseq == '-')
    editpos <- editpos[which.min(abs(editpos - 20))]
    editNT <- str_to_lower(refseq[editpos])
    ref <- seq_table$Aligned_Sequence[i]
    if(cutsite + 20 > str_length(wtseq)){
      ref <- str_sub(ref, 1, -(cutsite + 21 - str_length(wtseq)))
    }
    if(cutsite < 20){
      ref <- str_sub(ref, 20 - cutsite)
    }
    ref <- str_remove_all(ref, "-")
    findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
                                        subject = wtseq, max.mismatch = 5, with.indels = T)
    #找不到的情况
    if(length(findRef@ranges@start) == 0){
      print(x)
      print(i)
      return(-99)
    }
    res <- (findRef@ranges@start[1] + editpos - 1) - mutpos
    res
  }))
  seq_table$editNT <- unlist(lapply(1 : nrow(seq_table), function(i){
    editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
    refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
    editpos <- which(editseq == '-')
    editpos <- editpos[which.min(abs(editpos - 20))]
    editNT <- str_to_lower(refseq[editpos])
    editNT
  }))
  seq_table
})
edit_stat_of_distance <- data.frame(do.call(rbind, edit_stat_of_distance))
table(edit_stat_of_distance$dis2mut)



edit_stat_of_distance <- merge(edit_stat_of_distance, total_sample_info, by.x="id", by.y = 'tissueid')
edit_stat_of_distance$gene <- unlist(lapply(edit_stat_of_distance$mutation, function(x){
  unlist(strsplit(x, "[-]"))[2]
}))
edit_stat_of_distance$loca <- unlist(lapply(edit_stat_of_distance$mutation, function(x){
  unlist(strsplit(x, "[-]"))[3]
}))
edit_stat_of_distance$sg <- unlist(lapply(edit_stat_of_distance$mutation, function(x){
  unlist(strsplit(x, "[-]"))[4]
}))


edit_stat_of_distance$geneFactor <- factor(edit_stat_of_distance$gene, levels = c(
  "Myo7a","USH1C","PCDH15","SANS", "USH2A","OTOF"
))
edit_stat_of_distance <- edit_stat_of_distance[order(edit_stat_of_distance$geneFactor, edit_stat_of_distance$loca,edit_stat_of_distance$sg, edit_stat_of_distance$Pct),]

plot_mat <- matrix(0, nrow = arround_len * 2, ncol=length(unique(edit_stat_of_distance$id)))
rownames(plot_mat) <- as.character(c(((-arround_len): (-1)) , 1 : arround_len))
colnames(plot_mat) <- unique(edit_stat_of_distance$id)
for(i in 1 : nrow(edit_stat_of_distance)){
  dis <- edit_stat_of_distance$dis2mut[i]
  if(dis >= 0){
    dis <- dis + 1
  }
  dis <- as.character(dis)
  if(!dis %in% rownames(plot_mat))
    next
  id <- edit_stat_of_distance$id[i]
  plot_mat[dis, id] <- plot_mat[dis, id] + edit_stat_of_distance$Pct[i]
}
plot_mat <- t(plot_mat)
dis_mat <- matrix(0, nrow = arround_len * 2, ncol=length(unique(edit_stat_of_distance$id)))
rownames(dis_mat) <- as.character(c(((-arround_len): (-1)) , 1 : arround_len))
colnames(dis_mat) <- unique(edit_stat_of_distance$id)
sp_info <- edit_stat_of_distance[,c("id", "dis", "sg", "gene", "loca", "geneFactor", "strand")]
sp_info <- sp_info[!duplicated(sp_info$id),]
for(i in 1 : nrow(sp_info)){
  dis <- sp_info$dis[i] + 1
  if(dis >= 0){
    dis <- dis + 1
  }
  dis <- as.character(dis)
  if(!dis %in% rownames(dis_mat))
    next
  id <- sp_info$id[i]
  if(sp_info$strand[i] == "-"){
    dis_mat[dis, id] <- -1
  }
  else{
    dis_mat[dis, id] <- 1
  }
  
}
dis_mat <- t(dis_mat)
seq_mat <- matrix('', nrow = arround_len * 2, ncol = length(unique(edit_stat_of_distance$id)))
rownames(seq_mat) <- as.character(c(((-arround_len): (-1)) , 1 : arround_len))
colnames(seq_mat) <- unique(edit_stat_of_distance$id)

seq_arround_mut <- lapply(1 : nrow(total_sample_info), function(i){
  seq <- unlist(strsplit(total_sample_info$seq[i], "*"))
  mutpos <- total_sample_info$mutPos[i]
  str_to_upper(seq[(mutpos - arround_len) : (mutpos + arround_len - 1) ])
})
for(i in 1 : length(seq_arround_mut)){
  seq_mat[,total_sample_info$tissueid[i]] <- seq_arround_mut[[i]]
}
seq_mat <- t(seq_mat)





deleteNT_pct <- lapply(split(edit_stat_of_distance,edit_stat_of_distance$id), function(x){
  tmp <- lapply(split(x$Pct, x$editNT), sum)
  tmp <- data.frame(editNT = names(tmp), Pct = unlist(tmp))
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  tmp$id <- unique(x$id)
  tmp
})
deleteNT_pct <- data.frame(do.call(rbind, deleteNT_pct))
deleteNT_pct <- deleteNT_pct[deleteNT_pct$editNT != "n",]
insert_mat <- matrix(0, ncol = 4, nrow = length(unique(edit_stat_of_distance$id)))
colnames(insert_mat) <- c("a", 't','c','g')
rownames(insert_mat) <- unique(edit_stat_of_distance$id)
for(i in 1 : nrow(deleteNT_pct)){
  insert_mat[deleteNT_pct$id[i], deleteNT_pct$editNT[i]] <- insert_mat[deleteNT_pct$id[i], deleteNT_pct$editNT[i]] + 
    deleteNT_pct$Pct[i]
}


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

indel_stat <- lapply(total_old_edit_table_rev, function(x){
  indel <- unlist(lapply(x, function(y){
    indel <- y$n_inserted + y$n_deleted
    sum(y$X.Reads.1[indel != 0])
  }))
  se <- plotrix::std.error(indel)
  data.frame(indel = mean(indel), se = se)
})
indel_stat <- data.frame(do.call(rbind, indel_stat))
indel_stat$id <- rownames(indel_stat)
indel_stat$se[is.na(indel_stat$se)] <- 0
indel_stat <- indel_stat[unique(edit_stat_of_distance$id),]

del1_pct <- lapply(total_old_edit_table_rev, function(x){
  indel <- unlist(lapply(x, function(y){
    indel <- y$n_inserted + y$n_deleted
    sum(y$X.Reads.1[y$n_deleted == 1]) / sum(y$X.Reads.1[indel != 0]) * 100
  }))
  se <- plotrix::std.error(indel)
  data.frame(indel = mean(indel), se = se)
})
del1_pct <- data.frame(do.call(rbind, del1_pct))
del1_pct$id <- rownames(del1_pct)
del1_pct$se[is.na(del1_pct$se)] <- 0






lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16, 
         legend_gp = gpar(fill = 2:5))
)

outside_sel <- abs(rowSums(dis_mat))
dis_mat <- dis_mat[outside_sel == 1,]
plot_mat <- plot_mat[outside_sel == 1,]
seq_mat <- seq_mat[outside_sel == 1,]
sp_info <- sp_info[sp_info$id %in% rownames(plot_mat),]
insert_mat <- insert_mat[outside_sel == 1,]
del1_pct <- del1_pct[sp_info$id,]
indel_stat <- indel_stat[sp_info$id,]





orf_table <- lapply(1 : nrow(total_sample_info), function(i){
  cds <- total_sample_info$cds[i]
  seq <- total_sample_info$seq[i]
  cds_loc <- str_locate(seq, cds)[1,1]
  cds_loc_last <- str_locate(seq, cds)[1,2]
  seq <- unlist(strsplit(total_sample_info$seq[i], "*"))
  mutpos <- total_sample_info$mutPos[i]
  mutpos_in_cds <- mutpos - cds_loc + 1
  
  
  cds <- unlist(strsplit(cds, "*"))
  if(mutpos_in_cds < 6){
    cds <- c(seq[(cds_loc - 6) : (cds_loc - 1)],cds)
    mutpos_in_cds <- mutpos_in_cds + 6
  }
  if(length(cds) - mutpos_in_cds < 6){
    cds <- c(cds, seq[(cds_loc_last + 1) : (cds_loc_last + 6)])
  }
  tmp_res <- c()
  res <- list()
  first <- T
  insertNT <- 7
  for(i in 1 : (mutpos_in_cds + arround_len)){
    if(i %% 3 == 1){
      if(length(tmp_res) != 0){
        res[[length(res) + 1]] <- paste0(tmp_res, collapse = "")
      }
      tmp_res <- c()
    }
    if(i >= mutpos_in_cds - arround_len){
      if(first){
        if(i %% 3 == 2){
          tmp_res <- c(tmp_res, cds[i - 1])
          insertNT <- 8
        }
        if(i %% 3 == 0){
          tmp_res <- c(tmp_res, cds[c(i - 2, i - 1)])
          insertNT <- 9
        }
        first <- F
      }
      
      tmp_res <- c(tmp_res, cds[i])
    }
  }
  if(length(tmp_res) != 0){
    if(length(tmp_res) != 3){
      if(length(tmp_res) == 2){
        tmp_res <- c(tmp_res, cds[i + 1])
      }else{
        tmp_res <- c(tmp_res, cds[c(i + 1, i + 2)])
      }
    }
    res[[length(res) + 1]] <- paste0(tmp_res, collapse = "")
  }
  c(paste0(unlist(res), collapse = ""), insertNT)
})
names(orf_table) <- total_sample_info$tissueid
orf_table <- data.frame(do.call(rbind, orf_table))
orf_table$id <- rownames(orf_table)
orf_table <- merge(orf_table, total_sample_info[,c(1,2)], by.x="id", by.y="tissueid")
rownames(orf_table) <- orf_table$id


sel_for_detail_plot <- sp_info[sp_info$id %in% c("20-051", "62-051", "64-051","104-051"),]
orf_table <- orf_table[sel_for_detail_plot$id,]

orf_table$bestEdit <- unlist(lapply(orf_table$id, function(x){
  seq_table <- seqs_table_split[[x]]
  seq_table <- seq_table[order(seq_table$Pct,decreasing = T),]
  dis <- c(total_used_seq_tissue$mutPos - total_used_seq_tissue$cutSite)[total_used_seq_tissue$result == x]
  mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
  cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
  seq_table$dis <- ifelse(dis <= 0, dis - 1, dis)
  seq <- total_used_seq$seq[total_used_seq$result == x]
  cmpseq <- seq_info$WT.CDS[seq_info$Experiment.Number == x]
  cmpseq_mut <- seq_info$Mutant.CDS[seq_info$Experiment.Number == x]
  strand <- total_used_seq_tissue$strand[total_used_seq_tissue$result == x]
  sg <- str_to_lower(total_used_seq_tissue$sg[total_used_seq_tissue$result == x])
  mutNT <- total_used_seq_tissue$mutNT[total_used_seq_tissue$result == x]
  
  mutseq <- str_to_lower(total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x])
  wtseq <-  unlist(strsplit(total_used_seq_tissue$seq[total_used_seq_tissue$result == x], "*"))
  i <- 1
  editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
  refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
  editpos <- which(editseq == '-')
  editpos <- editpos[which.min(abs(editpos - 20))]
  editNT <- str_to_lower(refseq[editpos])
  # ref <- seq_table$Reference_Sequence[i]
  #算出和cutsite的距离
  editpos <- editpos - 20
  ref <- sg
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
                                      subject = mutseq, max.mismatch = 1, with.indels = T, min.mismatch = 0)
  mutseqs <- unlist(strsplit(mutseq, "*"))
  mutseqs[mutpos] <- str_to_upper(mutseqs[mutpos])
  if(strand == "+"){
    deletepos <- editpos + findRef@ranges@start + findRef@ranges@width - 4
    mutseqs <- mutseqs[-deletepos]
  }else{
    deletepos <- editpos + findRef@ranges@start + 2
    mutseqs <- mutseqs[-deletepos]
  }
  mutseqs <- paste0(mutseqs, collapse = "")
  # editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
  # refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
  # editpos <- which(editseq == '-')
  # editpos <- editpos[which.min(abs(editpos - 20))]
  # editNT <- str_to_lower(refseq[editpos])
  # ref <- seq_table$Reference_Sequence[i]
  # 
  # if(cutsite + 20 > length(wtseq)){
  #   ref <- str_sub(ref, 1, -(cutsite + 21 - length(wtseq)))
  # }
  # if(cutsite < 20){
  #   ref <- str_sub(ref, 20 - cutsite)
  # }
  # ref <- str_remove_all(ref, "-")
  # findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
  #                                     subject = total_used_seq_tissue$seq[total_used_seq_tissue$result == x], max.mismatch = 4, with.indels = T)
  # deletepos <- findRef@ranges@start + editpos - 1
  # findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
  #                                     subject = mutseq, max.mismatch = 4, with.indels = T)
  # #找不到的情况
  # if(length(findRef@ranges@start) == 0){
  #   return(-1)
  # }
  # mutseqs <- unlist(strsplit(mutseq, "*"))
  # if(cutsite + 6 < mutpos){
  #   matchpos <- findRef@ranges@start
  #   mutseqs <- c(mutseqs[1 : (matchpos + editpos - 2)], mutseqs[(matchpos + editpos) : length(mutseqs)])
  #   mutseqs <- paste0(mutseqs, collapse = "")
  # }else{
  #   matchpos <- length(mutseqs) - (findRef@ranges@start + findRef@ranges@width - 1) + 1
  #   editpos <- str_length(ref) - editpos + 1
  #   mutseqs <- rev(mutseqs)
  #   mutseqs <- c(mutseqs[1 : (matchpos + editpos - 1)], mutseqs[(matchpos + editpos + 1) : length(mutseqs)])
  #   mutseqs <- rev(mutseqs)
  #   mutseqs <- paste0(mutseqs, collapse = "")
  # }
  orf <- orf_table$X1[orf_table$id == x]
  findRef <- Biostrings::matchPattern(pattern= str_to_lower(orf),
                                      subject = mutseqs, max.mismatch = 2, with.indels = T)
  if(length(findRef@ranges@start) != 1){
    print(x)
  }
  str_sub(mutseqs, findRef@ranges@start, findRef@ranges@start + findRef@ranges@width - 1)
}))
orf_table$editCut <- unlist(lapply(orf_table$id, function(x){
  seq_table <- seqs_table_split[[x]]
  seq_table <- seq_table[order(seq_table$Pct,decreasing = T),]
  dis <- c(total_used_seq_tissue$mutPos - total_used_seq_tissue$cutSite)[total_used_seq_tissue$result == x]
  mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
  cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
  seq_table$dis <- ifelse(dis <= 0, dis - 1, dis)
  seq <- total_used_seq$seq[total_used_seq$result == x]
  cmpseq <- seq_info$WT.CDS[seq_info$Experiment.Number == x]
  cmpseq_mut <- seq_info$Mutant.CDS[seq_info$Experiment.Number == x]
  strand <- total_used_seq_tissue$strand[total_used_seq_tissue$result == x]
  sg <- str_to_lower(total_used_seq_tissue$sg[total_used_seq_tissue$result == x])
  mutNT <- total_used_seq_tissue$mutNT[total_used_seq_tissue$result == x]
  
  mutseq <- str_to_lower(total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x])
  wtseq <-  unlist(strsplit(total_used_seq_tissue$seq[total_used_seq_tissue$result == x], "*"))
  i <- 1
  editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
  refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
  editpos <- which(editseq == '-')
  editpos <- editpos[which.min(abs(editpos - 20))]
  editNT <- str_to_lower(refseq[editpos])
  # ref <- seq_table$Reference_Sequence[i]
  #算出和cutsite的距离
  editpos <- editpos - 20
  # if(strand == "+" & editpos == 0){
  #   return(0)
  # }
  # if(strand == "-" & editpos == 1){
  #   return(0)
  # }
  if(editpos <= 0){
    return(-1)
  }
  return(1)
}))
orf_table$mutCut <- unlist(lapply(orf_table$id, function(x){
  seq_table <- seqs_table_split[[x]]
  seq_table <- seq_table[order(seq_table$Pct,decreasing = T),]
  dis <- c(total_used_seq_tissue$mutPos - total_used_seq_tissue$cutSite)[total_used_seq_tissue$result == x]
  mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
  cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
  seq_table$dis <- ifelse(dis <= 0, dis - 1, dis)
  seq <- total_used_seq$seq[total_used_seq$result == x]
  cmpseq <- seq_info$WT.CDS[seq_info$Experiment.Number == x]
  cmpseq_mut <- seq_info$Mutant.CDS[seq_info$Experiment.Number == x]
  strand <- total_used_seq_tissue$strand[total_used_seq_tissue$result == x]
  sg <- str_to_lower(total_used_seq_tissue$sg[total_used_seq_tissue$result == x])
  mutNT <- total_used_seq_tissue$mutNT[total_used_seq_tissue$result == x]
  
  mutseq <- str_to_lower(total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x])
  wtseq <-  unlist(strsplit(total_used_seq_tissue$seq[total_used_seq_tissue$result == x], "*"))
  i <- 1
  editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
  refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
  editpos <- which(editseq == '-')
  editpos <- editpos[which.min(abs(editpos - 20))]
  editNT <- str_to_lower(refseq[editpos])
  # ref <- seq_table$Reference_Sequence[i]
  #算出和cutsite的距离
  editpos <- editpos - 20
  ref <- sg
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
                                      subject = mutseq, max.mismatch = 1, with.indels = T, min.mismatch = 0)
  mutseqs <- unlist(strsplit(mutseq, "*"))
  mutseqs[mutpos] <- str_to_upper(mutseqs[mutpos])
  if(strand == "+"){
    cutpos <- findRef@ranges@start + findRef@ranges@width - 4

  }else{
    cutpos <- findRef@ranges@start + 2
  }
  mutseqs <- paste0(mutseqs, collapse = "")
  orf <- orf_table$X1[orf_table$id == x]
  findRef <- Biostrings::matchPattern(pattern= str_to_lower(orf),
                                      subject = mutseqs, max.mismatch = 2, with.indels = T)
  cutpos - findRef@ranges@start + 1
}))
orf_table$mutSeq <- unlist(lapply(orf_table$id, function(x){

  mutseq <- total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x]
  mutPos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
  
  mutseq <- unlist(strsplit(mutseq,"*"))
  mutseq[mutPos] <- str_to_upper(mutseq[mutPos])
  
  mutseq <- paste(mutseq, collapse = "")
  seq <- orf_table$X1[orf_table$id == x]
  findRef <- Biostrings::matchPattern(pattern= str_to_lower(seq),
                                      subject = mutseq, max.mismatch = 2, with.indels = T)
  
  str_sub(mutseq, findRef@ranges@start, findRef@ranges@start + findRef@ranges@width - 1)
}))


edit_mat <- lapply(1 : nrow(orf_table), function(i){
  seqs <- unlist(strsplit(orf_table$bestEdit[i], "*"))
  
  mutseq <- str_sub(orf_table$mutSeq[i], 1, 15)
  pos <- which(seqs %in% c("A", "T", 'C','G'))
  if(length(pos) == 0)
    pos <- -1
  mutpos <- which(unlist(strsplit(mutseq, "*"))  %in% c("A", "T", 'C','G'))
  
  data.frame(seq = c(orf_table$X1[i],mutseq,orf_table$bestEdit[i]), 
             type = c("WT","Mut", "Top Edit"), 
             pos = c(orf_table$X2[i], mutpos, pos), 
             cutpos = c(orf_table$X2[i], orf_table$mutCut[i], orf_table$mutCut[i] + orf_table$editCut[i]),
             mutation = orf_table$mutation[i],
             id = orf_table$id[i])
  
})

edit_mat <- data.frame(do.call(rbind, edit_mat))
tmp <- merge(edit_mat,sel_for_detail_plot, by = 'id')
rownames(tmp) <- paste0(tmp$id, tmp$type)
tmp <- tmp[paste0(edit_mat$id, edit_mat$type),]
edit_mat <- tmp
# edit_mat <- edit_mat[order(edit_mat$geneFactor),]
edit_seq_mat <- matrix("", nrow = nrow(edit_mat), ncol = 15)
rownames(edit_seq_mat) <- rownames(edit_mat)
for(i in 1 : nrow(edit_mat)){
  edit_seq_mat[i,] <- unlist(strsplit(str_to_upper(edit_mat$seq[i]), "*"))
}
fake_plot_mat <- matrix(runif(nrow(edit_mat) * 15), nrow = nrow(edit_mat), ncol = 15)

edit_mat$pos <- as.numeric(edit_mat$pos)
edit_mat$cutpos <- as.numeric(edit_mat$cutpos)
left_anno <- 
  rowAnnotation(gene = anno_empty(border = F, width = unit(20, "mm")),
                "line1" = 
                  anno_empty(border = FALSE, width = unit(2, "mm")),
                loc = anno_empty(border = F, width = unit(20, "mm")),
                "line2" = 
                  anno_empty(border = FALSE, width = unit(2, "mm")),
                sg = anno_empty(border = F, width = unit(20, "mm"))
  )
pdf("Result/Del1_edit_between(+-1)_result_v2.pdf", width = 10, height = 6)
ht <- Heatmap(fake_plot_mat, 
              name = "EditDis",
              column_split = rep(c(" ", "  ", "   ", "    ", "     "), c(3,3,3,3,3)),
              column_title=NULL,
              row_title = NULL,
              row_split = edit_mat$geneFactor,
              show_column_names = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 24),
              column_names_side = "top",
              left_annotation = left_anno,
              right_annotation = NULL,
              show_row_names = F,
              cell_fun = function(j, i, x, y, w, h, fill) {
                # print(i)
                
                # transform the matrix into a long vector
                nt = pindex(edit_seq_mat, i, j)
                strand = edit_mat$strand[i]
              
                pos = edit_mat$pos[i]
                if(pos == -1){
                  cutpos = edit_mat$pos[i - 2] + edit_mat$dis[i] + 1
                }else{
                  cutpos = pos + edit_mat$dis[i] + 1
                }
                type = edit_mat$type[i]
                
                dis = edit_mat$dis[i]
                strand = edit_mat$strand[i]
                if(edit_mat$dis[i] >= 1 & type != "WT"){
                  cutpos = cutpos + 1
                }
                if(type == 'Mut'){
                  if(strand == "+" & dis < 0){
                    cutpos = cutpos + 1
                  }
                  if(strand == '-' & dis > 1){
                    cutpos = cutpos - 1
                  }
                }
                if(type != "WT"){
                  cutpos = edit_mat$cutpos[i] + 1
                }
                # `j` here is also a vector with the same length of `v`
                
                col2 = "white"
                grid.rect(x, y, w, h, gp = gpar(fill = col2, col = col2))
                
                # if(dis== 1){
                #   grid.points(x, y, gp = gpar(col = "red"))
                # }
                x_from <- x - unit(as.numeric(w) / 2, "npc")
                x_to <- x + unit(as.numeric(w) / 2, "npc")
                y_from <- y - unit(as.numeric(h) / 2, "npc")
                y_to <- y + unit(as.numeric(h) / 2, "npc")
                if(j == pos && type != "WT"){
                  col = c(2 : 5)
                  names(col) <- c("A", 'T','C','G')
                  grid.text(nt, x, y, 
                            gp = gpar(fontsize = 20, col = col[nt]))
                }else{
                  grid.text(nt, x, y, 
                            gp = gpar(fontsize = 18))
                }
                if(j == cutpos){
                  if(strand == "+"){
                    
                    grid.polygon(x = c(x_from - unit(as.numeric(w) / 8, "npc"), 
                                       x_from, 
                                       x_from + unit(as.numeric(w) / 8, "npc")),
                                 y = c(y + unit(as.numeric(h) / 3, "npc"), 
                                       y - unit(as.numeric(h) / 3, "npc"), 
                                       y + unit(as.numeric(h) / 3, "npc")), 
                                 gp = gpar(col = "red", fill = "red"))
                    
                  }
                  if(strand == "-"){
                    
                    grid.polygon(x = c(x_from - unit(as.numeric(w) / 8, "npc"), 
                                       x_from, 
                                       x_from + unit(as.numeric(w) / 8, "npc")),
                                 y = c(y - unit(as.numeric(h) / 3, "npc"), 
                                       y + unit(as.numeric(h) / 3, "npc"), 
                                       y - unit(as.numeric(h) / 3, "npc")), 
                                 gp = gpar(col = "red", fill = "red"))
                    
                  }
                }
                
                
              }, layer_fun = function(j, i, x, y, w, h, fill) {
                ind_mat = restore_matrix(j, i, x, y)
                for(ir in seq_len(nrow(ind_mat))) {
                  ind1 = ind_mat[ir, 1] # previous column
                  ind2 = ind_mat[ir, ncol(ind_mat)]
                  grid.segments(x[ind1], y[ind1] - unit(as.numeric(h[ind1]) / 3, "npc"), x[ind2], y[ind2] - unit(as.numeric(h[ind2]) / 3, "npc"),
                                gp = gpar(col = "green", lwd = 2))
                }
                
                
              },col=c("white", "blue"),show_heatmap_legend = F,
              border = F)
ht <- draw(ht, annotation_legend_list = lgd_list)
lens <- unique(edit_mat$gene)
pos_line1 <- 0
##因为是从下到上画的 注意需要反过来画
edit_mat  <- edit_mat[nrow(edit_mat) : 1,]
edit_mat$sgname <- edit_mat$type
edit_mat$sgname[!edit_mat$sgname %in% c("Mut","WT")] <- paste0("MT+", edit_mat$sg[!edit_mat$sgname %in% c("Mut","WT")])
for(i in 1:length(unique(lens))) {
  len <- sum(edit_mat$gene == lens[i])
  decorate_annotation("gene", slice = i, {
    tg <- richtext_grob(gt_render(lens[i]), 
                        rot = 0, 
                        x = unit(0.5, "npc"),
                        y=unit(0.5, "npc"), hjust = 0.5)
    grid.draw(tg)
    invisible(tg)
  })
  locs_count <- length(unique(edit_mat$loca[edit_mat$gene == lens[i]]))
  locs <- unique(edit_mat$loca[edit_mat$gene == lens[i]])
  sg_count <- sum(edit_mat$gene == lens[i])
  if(locs_count == 1){
    decorate_annotation("line1", slice = i, {
      grid.lines(y = unit(c(0.5, 0.5), "npc"),
                 x=unit(c(0, 1), "npc"))
    })
  }else{
    decorate_annotation("line1", slice = i, {
      grid.lines(y = unit(c(1 / sg_count / 2, 
                            1 / sg_count / 2), "npc"),
                 x=unit(c(0, 1), "npc"))
      grid.lines(y = unit(c(1 / sg_count / 2, 
                            1 - 1 / sg_count / 2), "npc"),
                 x=unit(c(0, 0), "npc"))
      grid.lines(y = unit(c(1 - 1 / sg_count / 2,
                            1 - 1 / sg_count / 2), "npc"),
                 x=unit(c(0, 1), "npc"))
    })
  }
  
  pos <- 0
  
  for(l in locs){
    
    ll <- sum(edit_mat$loca[edit_mat$gene == lens[i]] == l)
    pos <- pos + 1 / sg_count * ll / 2
    decorate_annotation("loc", slice = i, {
      tg <- richtext_grob(gt_render(l), 
                          rot = 0, 
                          x = unit(0.5, "npc"),
                          y=unit(pos, "npc"), hjust = 0.5)
      grid.draw(tg)
      invisible(tg)
    })
    if(ll == 1){
      decorate_annotation("line2", slice = i, {
        grid.lines(y = unit(c(pos, pos), "npc"),
                   x=unit(c(0, 1), "npc"))
      })
    }
    else{
      decorate_annotation("line2", slice = i, {
        grid.lines(y = unit(c(pos - 1 / sg_count * ll / 2 + 1 / sg_count / 2, pos - 1 / sg_count * ll / 2 + 1 / sg_count / 2), "npc"),
                   x=unit(c(0, 1), "npc"))
        grid.lines(y = unit(c(pos - 1 / sg_count * ll / 2 + 1 / sg_count / 2, 
                              pos + 1 / sg_count * ll / 2 - 1 / sg_count / 2), "npc"),
                   x=unit(c(0, 0), "npc"))
        grid.lines(y = unit(c(pos + 1 / sg_count * ll / 2 - 1 / sg_count / 2,
                              pos + 1 / sg_count * ll / 2 - 1 / sg_count / 2), "npc"),
                   x=unit(c(0, 1), "npc"))
      })
    }
    
    
    
    pos <- pos + 1 / sg_count * ll / 2
    
    
    
  }
  sgs <- edit_mat$sgname[edit_mat$gene == lens[i]]
  pos <- 1 / sg_count / 2
  for(s in sgs){
    decorate_annotation("sg", slice = i, {
      tg <- richtext_grob(gt_render(s), 
                          rot = 0, 
                          x = unit(0.5, "npc"),
                          y=unit(pos, "npc"), hjust = 0.5)
      grid.draw(tg)
      invisible(tg)
    })
    pos <- pos + 1 / sg_count
  }
  od <- row_order(ht)
}



dev.off()









