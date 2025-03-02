###这个脚本依赖于4.1的统计结果
total_sample_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)
total_sample_info <- data.frame(
  tissueid = c(total_sample_info$`Tissue-w-W(Mix)`,
               paste0("m", 12:15, "-tissue")),
  mutation = c(total_sample_info$GRNAa, paste0("w-OTOF-1236dC-g", 3 : 6)))
total_sample_info <- merge(total_sample_info, total_used_seq_tissue, by.x="tissueid", by.y="result")
total_sample_info <- total_sample_info[total_sample_info$tissueid %in% names(seqs_table_split),]
arround_len <- 15




edit_stat_of_distance <- lapply(names(seqs_table_split), function(x){
  seq_table <- seqs_table_split[[x]]
  seq_table <- seq_table[order(seq_table$Pct,decreasing = T),]
  mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
  cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
  dis <- c(total_used_seq_tissue$cutSite - total_used_seq_tissue$mutPos)[total_used_seq_tissue$result == x]
  
  seq_table$dis <- ifelse(dis >= 0, dis + 1, dis)
  # mutseq <- total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x]
  wtseq <-  total_used_seq_tissue$seq[total_used_seq_tissue$result == x]
  seq_table$dis2mut <- unlist(lapply(1 : nrow(seq_table), function(i){
    editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
    refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
    editpos <- which(refseq == '-')
    editpos <- editpos[which.min(abs(editpos - 20))]
    editNT <- str_to_lower(editseq[editpos])
    ref <- seq_table$Reference_Sequence[i]
    if(cutsite + 20 > str_length(wtseq)){
      ref <- str_sub(ref, 1, -(cutsite + 21 - str_length(wtseq)))
    }
    if(cutsite < 20){
      ref <- str_sub(ref, 20 - cutsite)
    }
    ref <- str_remove_all(ref, "-")
    findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
                                        subject = wtseq, max.mismatch = 4, with.indels = T)
    #找不到的情况
    if(length(findRef@ranges@start) == 0){
      return(-99)
    }
    res <- (findRef@ranges@start[1] + editpos - 1) - mutpos
    res
  }))
  seq_table$editNT <- unlist(lapply(1 : nrow(seq_table), function(i){
    editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
    refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
    editpos <- which(refseq == '-')
    editpos <- editpos[which.min(abs(editpos - 20))]
    editNT <- str_to_lower(editseq[editpos])
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
  "Myo7a","USH1C","CDH23","PCDH15","USH2A","VLGR1","whirlin","TMC1","OTOF"
))
edit_stat_of_distance <- edit_stat_of_distance[order(edit_stat_of_distance$geneFactor, edit_stat_of_distance$loca,edit_stat_of_distance$sg, edit_stat_of_distance$Pct),]

plot_mat <- matrix(0, nrow = arround_len * 2 + 1, ncol=length(unique(edit_stat_of_distance$id)))
rownames(plot_mat) <- as.character(c(-arround_len : arround_len))
colnames(plot_mat) <- unique(edit_stat_of_distance$id)
for(i in 1 : nrow(edit_stat_of_distance)){
  dis <- as.character(edit_stat_of_distance$dis2mut[i])
  if(!dis %in% rownames(plot_mat))
    next
  id <- edit_stat_of_distance$id[i]
  plot_mat[dis, id] <- plot_mat[dis, id] + edit_stat_of_distance$Pct[i]
}
plot_mat <- t(plot_mat)
dis_mat <- matrix(0, nrow = arround_len * 2 + 1, ncol=length(unique(edit_stat_of_distance$id)))
rownames(dis_mat) <- as.character(c(-arround_len:arround_len))
colnames(dis_mat) <- unique(edit_stat_of_distance$id)
sp_info <- edit_stat_of_distance[,c("id", "dis", "sg", "gene", "loca", "geneFactor", "strand")]
sp_info <- sp_info[!duplicated(sp_info$id),]
for(i in 1 : nrow(sp_info)){
  
  if(sp_info$dis[i] < 0){
    dis <- as.character(sp_info$dis[i] + 1)
  }else{
    dis <- as.character(sp_info$dis[i])
  }
  if(!dis %in% as.character((-arround_len) : arround_len))
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
seq_mat <- matrix('', nrow = arround_len * 2 + 1, ncol = length(unique(edit_stat_of_distance$id)))
rownames(seq_mat) <- as.character(c(-arround_len:arround_len))
colnames(seq_mat) <- unique(edit_stat_of_distance$id)

seq_arround_mut <- lapply(1 : nrow(total_sample_info), function(i){
  seq <- unlist(strsplit(total_sample_info$seq[i], "*"))
  mutpos <- total_sample_info$mutPos[i]
  str_to_upper(seq[(mutpos - arround_len) : (mutpos + arround_len) ])
})
for(i in 1 : length(seq_arround_mut)){
  seq_mat[,total_sample_info$tissueid[i]] <- seq_arround_mut[[i]]
}
seq_mat <- t(seq_mat)




insertNT_pct <- lapply(split(edit_stat_of_distance,edit_stat_of_distance$id), function(x){
  tmp <- lapply(split(x$Pct, x$editNT), sum)
  tmp <- data.frame(editNT = names(tmp), Pct = unlist(tmp))
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  tmp$id <- unique(x$id)
  tmp
})
insertNT_pct <- data.frame(do.call(rbind, insertNT_pct))
insertNT_pct <- insertNT_pct[insertNT_pct$editNT != "n",]
insert_mat <- matrix(0, ncol = 4, nrow = length(unique(edit_stat_of_distance$id)))
colnames(insert_mat) <- c("a", 't','c','g')
rownames(insert_mat) <- unique(edit_stat_of_distance$id)
for(i in 1 : nrow(insertNT_pct)){
  insert_mat[insertNT_pct$id[i], insertNT_pct$editNT[i]] <- insert_mat[insertNT_pct$id[i], insertNT_pct$editNT[i]] + 
    insertNT_pct$Pct[i]
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




ins1_pct <- lapply(total_old_edit_table_rev, function(x){
  indel <- unlist(lapply(x, function(y){
    indel <- y$n_inserted + y$n_deleted
    sum(y$X.Reads.1[y$n_inserted == 1]) / sum(y$X.Reads.1[indel != 0]) * 100
  }))
  se <- plotrix::std.error(indel)
  data.frame(indel = mean(indel), se = se)
})
ins1_pct <- data.frame(do.call(rbind, ins1_pct))
ins1_pct$id <- rownames(ins1_pct)
ins1_pct$se[is.na(ins1_pct$se)] <- 0




lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16, 
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)

outside_sel <- abs(rowSums(dis_mat))
dis_mat <- dis_mat[outside_sel == 1,]
plot_mat <- plot_mat[outside_sel == 1,]
seq_mat <- seq_mat[outside_sel == 1,]
sp_info <- sp_info[sp_info$id %in% rownames(plot_mat),]
insert_mat <- insert_mat[outside_sel == 1,]
ins1_pct <- ins1_pct[sp_info$id,]
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
  if(mutpos_in_cds < 15){
    cds <- c(seq[(cds_loc - 15) : (cds_loc - 1)],cds)
    mutpos_in_cds <- mutpos_in_cds + 15
  }
  if(length(cds) - mutpos_in_cds < 15){
    cds <- c(cds, seq[(cds_loc_last + 1) : (cds_loc_last + 15)])
  }
  tmp_res <- c()
  res <- list()
  first <- T
  delNT <- 16
  for(i in 1 : (mutpos_in_cds + arround_len + 3)){
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
          delNT <- 17
        }
        if(i %% 3 == 0){
          tmp_res <- c(tmp_res, cds[c(i - 2, i - 1)])
          delNT <- 18
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
  c(paste0(unlist(res), collapse = ""), delNT)
})
names(orf_table) <- total_sample_info$tissueid
orf_table <- data.frame(do.call(rbind, orf_table))
orf_table$id <- rownames(orf_table)
orf_table <- merge(orf_table, total_sample_info[,c(1,2)], by.x="id", by.y="tissueid")
rownames(orf_table) <- orf_table$id


sel_for_detail_plot <- sp_info[(sp_info$dis == 1 & sp_info$strand == "+") | 
                                 (sp_info$dis == -1 & sp_info$strand == "-"),]
####计算+1 -2 -5的编辑产物
#首先从5.2脚本的结果inframe_result_processed确认是否存在这三种编辑
sel_top3_indel <- lapply(sel_for_detail_plot$id, function(x){
  tmp <- inframe_result_processed[inframe_result_processed$id == x,]
  tmp <- tmp[order(tmp$pct, decreasing = T),]
  tmp <- tmp[tmp$indel %in% c(0,-3,-6),]
  tmp$indel <- tmp$indel + 1
  tmp
})
names(sel_top3_indel) <- sel_for_detail_plot$id
#确定都有 然后开始挑选这三种做预测

seqs_table_stat <- lapply(sel_for_detail_plot$id, function(x){
  indels <- sel_top3_indel[[x]]$indel
  res <- lapply(indels, function(indel_len) {
    cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
    mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
    strand <- total_used_seq_tissue$strand[total_used_seq_tissue$result == x]
    #insert NT will be cap 
    mutNT <- total_used_seq_tissue$mutNT[total_used_seq_tissue$result == x]
    mutseq <- str_to_lower(total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x])
    sg <- str_to_lower(total_used_seq_tissue$sg[total_used_seq_tissue$result == x])
    wtseq <-  unlist(strsplit(total_used_seq_tissue$seq[total_used_seq_tissue$result == x], "*"))
    edit_tables <- total_old_edit_table_rev[[x]]
    edit_table <- lapply(edit_tables, function(edit_table){
      edit_table <- edit_table[edit_table$n_inserted + edit_table$n_deleted != 0,]
      inframe_sum <- edit_table[(edit_table$n_inserted - edit_table$n_deleted) %% 3 == 1,]
      edit_table$Pct <- edit_table[,7] / sum(edit_table[,7]) * 100
      edit_table$inframePct <- edit_table[,7] / sum(inframe_sum[,7]) * 100
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
      edit_table <- edit_table[(edit_table$n_inserted - 
                                   edit_table$n_deleted) == indel_len & 
                                 edit_table$n_mutated <= 1,]

      edit_table
    })
    edit_table <- data.frame(do.call(rbind, edit_table))
    if(nrow(edit_table) == 0)
      return(data.frame(id = x, indel = indel_len , 
                            seq = "", Pct = 0, 
                            inframePct = 0, pctSe = 0, 
                        inframeSe = 0))
    edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, 
                                                         edit_table$Aligned_Sequence), function(x){
      x[,7] <- sum(x[,7]) / length(edit_tables)
      x[,8] <- sum(x[,8]) / length(edit_tables)
      x$pctse <- plotrix::std.error(unlist(x[,9]))
      x$inframepctse <- plotrix::std.error(unlist(x[,10]))
      x[,9] <- sum(x[,9]) / length(edit_tables)
      x[,10] <- sum(x[,10]) / length(edit_tables)
      x[1,]
    })))
    edit_table$id <- x
    edit_table <- edit_table[,c(1,2, 8, 9, 10, 11, 12, 13)]
    edit_table <- edit_table[order(edit_table$X.Reads, decreasing = T),]
    data.frame(id = x, indel = indel_len , 
               seq = edit_table$Aligned_Sequence[1], Pct = edit_table$Pct[1], 
               inframePct = edit_table$inframePct[1], pctSe = edit_table$pctse[1], 
               inframeSe = edit_table$inframepctse[1])
  })
  res <- data.frame(do.call(rbind, res))
})

sel_info <- total_used_seq_tissue[total_used_seq_tissue$result %in% 
                                                 sel_for_detail_plot$id,]

seqs_table_stat <- data.frame(do.call(rbind, seqs_table_stat))
seqs_table_stat <- seqs_table_stat[order(seqs_table_stat$id, seqs_table_stat$indel),]
seqs_table_stat <- seqs_table_stat[seqs_table_stat$seq != "",]
seqs_table_stat <- seqs_table_stat[seqs_table_stat$id != "22-051",]
seqs_table_stat$predict_seq <- c("CCCGAG------------GAACGGCAGTGGGCACGGT",
                                 "CCCGAGGGCGTGCCC--GAACGGCAGTGGGCACGGTT", 
                                 "CCCGAGGGCGTGCCCCCCGAACGGCAGTGGGCACGGT", 
                                 "CCCCTGAGGCCTCG---------AGAGCCTGGGGA",
                                 "CCCCTGAGGCCTCGG--CGGAAAGAGCCTGGGGA",
                                 "CCCCTGAGGCCTCGGCCCCGGAAAGAGCCTGGGG", 
                                 "TGG--------TCAGGGAGTGGAGCCGGCAGACCA", 
                                 "TGGGTAGATGATCA--AGTGGAGCCGGCAGACCAC", 
                                 "TGGGTAGATGATCAGGGAGTGGAGCCGGCAGACCA",
                                 "GACAACCCTC-----TCAGCAAGCCGGCCTACTT",
                                 "GACAACCCTCCCACC--AGCAAGCCGGCCTACTT",
                                 "GACAACCCTCCCACCTTCAGCAAGCCGGCCTACT",
                                 "CAGTACTCTGTCCTC-----GGCTATTACGACAATA", 
                                 "CAGTACTCTGTCCTCT--ATGGCTATTACGACAATA",
                                 "CAGTACTCTGTCCTCTTTTATGGCTATTACGACAAT")

seqs_table_stat$predict_seq <- str_sub(seqs_table_stat$predict_seq, 1, 2 * arround_len)
orf_table <- orf_table[sel_for_detail_plot$id,]

orf_table$bestEdit <- unlist(lapply(orf_table$id, function(x){
  seq_table <- seqs_table_split[[x]]
  seq_table <- seq_table[order(seq_table$Pct,decreasing = T),]
  dis <- c(total_used_seq_tissue$mutPos - total_used_seq_tissue$cutSite)[total_used_seq_tissue$result == x]
  mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
  cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
  seq_table$dis <- ifelse(dis <= 0, dis - 1, dis)
  mutseq <- total_used_seq_tissue$mutSeq[total_used_seq_tissue$result == x]
  wtseq <-  unlist(strsplit(total_used_seq_tissue$seq[total_used_seq_tissue$result == x], "*"))
  i <- 1
  
  editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
  refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
  editpos <- which(refseq == '-')
  editpos <- editpos[which.min(abs(editpos - 20))]
  editNT <- str_to_lower(editseq[editpos])
  ref <- seq_table$Reference_Sequence[i]
  
  if(cutsite + 20 > length(wtseq)){
    ref <- str_sub(ref, 1, -(cutsite + 21 - length(wtseq)))
  }
  if(cutsite < 20){
    ref <- str_sub(ref, 20 - cutsite)
  }
  ref <- str_remove_all(ref, "-")
  findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
                                      subject = total_used_seq_tissue$seq[total_used_seq_tissue$result == x], 
                                      max.mismatch = 4, with.indels = T)
  insertpos <- findRef@ranges@start + editpos - 1
  findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
                                      subject = mutseq, max.mismatch = 4, with.indels = T)
  #找不到的情况
  if(length(findRef@ranges@start) == 0){
    return(-1)
  }
  mutseqs <- unlist(strsplit(mutseq, "*"))
  
  if(insertpos <= mutpos ){
    matchpos <- findRef@ranges@start
    mutseqs <- c(mutseqs[1 : (matchpos + editpos - 2)], str_to_upper(editNT), mutseqs[(matchpos + editpos - 1) : length(mutseqs)])
    mutseqs <- paste0(mutseqs, collapse = "")
  }else{
    matchpos <- length(mutseqs) - (findRef@ranges@start + findRef@ranges@width - 1) + 1
    editpos <- str_length(ref) - editpos + 1
    mutseqs <- rev(mutseqs)
    mutseqs <- c(mutseqs[1 : (matchpos + editpos - 1)], str_to_upper(editNT), mutseqs[(matchpos + editpos) : length(mutseqs)])
    mutseqs <- rev(mutseqs)
    mutseqs <- paste0(mutseqs, collapse = "")
  }
  orf <- orf_table$X1[orf_table$id == x]
  findRef <- Biostrings::matchPattern(pattern= str_to_lower(orf),
                                      subject = mutseqs, max.mismatch = 2, with.indels = T)
  if(length(findRef@ranges@start) != 1){
    print("Bad")
  }
  str_sub(mutseqs, findRef@ranges@start, findRef@ranges@start + findRef@ranges@width - 1)
}))


orf_table <- orf_table[orf_table$id != "22-051",]
edit_mat <- lapply(1 : nrow(orf_table), function(i){
  seqs <- unlist(strsplit(orf_table$bestEdit[i], "*"))
  mutseq <- seqs[!seqs %in% c("A","T",'C','G')]
  mutseq <- mutseq[1:30]
  seqs <- seqs[1:30]
  pos <- which(seqs %in% c("A", "T", 'C','G'))
  tmp_inframe <- seqs_table_stat[seqs_table_stat$id == orf_table$id[i],]
  tmp_inframe <- tmp_inframe[order(tmp_inframe$indel, decreasing = T),]
  data.frame(seq = c(str_sub(orf_table$X1[i], 1, 30),
                     paste(mutseq, collapse = ""), 
                     tmp_inframe$predict_seq), 
             type = c("WT","Mut", tmp_inframe$indel), 
             pos = c(orf_table$X2[i],pos, pos, rep(0, nrow(tmp_inframe) - 1)), 
             mutation = orf_table$mutation[i],
             id = orf_table$id[i])
  
})

edit_mat <- data.frame(do.call(rbind, edit_mat))
tmp <- merge(edit_mat,sel_for_detail_plot, by = 'id')
rownames(tmp) <- paste0(tmp$id, tmp$type)
tmp <- tmp[paste0(edit_mat$id, edit_mat$type),]
edit_mat <- tmp
edit_seq_mat <- matrix("", nrow = nrow(edit_mat), ncol = 2 * arround_len)
for(i in 1 : nrow(edit_mat)){
  edit_seq_mat[i,] <- unlist(strsplit(str_to_upper(edit_mat$seq[i]), "*"))
}

edit_mat$pos <- as.numeric(edit_mat$pos)

edit_mat_sel <- edit_mat[edit_mat$id != "22-051",]
fake_plot_mat <- matrix(runif(nrow(edit_mat_sel) * 2 * arround_len), nrow = nrow(edit_mat_sel), ncol = 2 * arround_len)
edit_seq_mat_sel <- edit_seq_mat[edit_mat$id != "22-051",]
library(ComplexHeatmap)
library(circlize)
library(gridtext)
left_anno <- 
  rowAnnotation(gene = anno_empty(border = F, width = unit(20, "mm")),
                "line1" = 
                  anno_empty(border = FALSE, width = unit(2, "mm")),
                loc = anno_empty(border = F, width = unit(20, "mm")),
                "line2" = 
                  anno_empty(border = FALSE, width = unit(2, "mm")),
                sg = anno_empty(border = F, width = unit(20, "mm"))
  )
insert_mat_sel <- insert_mat[rownames(insert_mat) %in% edit_mat_sel$id,]
insert_mat_sel_append <- matrix(0, nrow = nrow(insert_mat_sel) * 3, ncol = ncol(insert_mat_sel))
colnames(insert_mat_sel_append) <- colnames(insert_mat_sel)
for(i in 1 : nrow(insert_mat_sel)){
  insert_mat_sel_append[i * 3 - 2,] <- insert_mat_sel[i,]
}
indel_stat_sel <- indel_stat[indel_stat$id %in% edit_mat_sel$id,]
indel_stat_sel_append <- data.frame(do.call(rbind, lapply(split(indel_stat_sel, indel_stat_sel$id), function(x){
  x <- data.frame(indel = c(x$indel, 0, 0), se = c(x$se, 0, 0), id = x$id)
})))
indel_stat_sel_append$id <- factor(indel_stat_sel_append$id, levels = unique(edit_mat_sel$id))
indel_stat_sel_append <- indel_stat_sel_append[order(indel_stat_sel_append$id),]


inframe_pct_sel_append <- data.frame(do.call(rbind, lapply(indel_stat_sel$id, function(x){
  tmp_inframe <- seqs_table_stat[seqs_table_stat$id == x,]
  tmp_inframe <- tmp_inframe[order(tmp_inframe$indel, decreasing = T),]
  data.frame(indel = c(0, 0, tmp_inframe$inframePct), 
  se = c(0, 0, tmp_inframe$inframeSe), 
  id = x)
})))
inframe_pct_sel_append$id <- factor(inframe_pct_sel_append$id, levels = unique(edit_mat_sel$id))
inframe_pct_sel_append <- inframe_pct_sel_append[order(inframe_pct_sel_append$id),]
inframe_pct_sel_append$se[is.na(inframe_pct_sel_append$se)] <- 0
right_anno <- rowAnnotation(
  "% in Inframe" = anno_barplot(x = inframe_pct_sel_append$indel, 
                                     gp = gpar(lwd = 0.1, fill = "#BBBBBB"),width = unit(2, "cm"), 
                                     bar_width = 0.5, ylim = c(0, max(inframe_pct_sel_append$indel + 
                                                                        inframe_pct_sel_append$se + 3))))

pdf("~/Nutstore Files/Tobin/Previous/Ins1_edit_between(+-1)_result_v8.pdf", width = 11, height = 9)
ht <- Heatmap(fake_plot_mat, 
              name = "EditDis",
              #column_split = rep(unlist(lapply(1 : 10, 
              #                                 function(x){
              #                                   paste(rep(" ", x), collapse = "")})),
              #                   rep(3, 10)),
              column_title=NULL,
              row_title = NULL,
              row_split = edit_mat_sel$geneFactor,
              show_column_names = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 24),
              column_names_side = "top",
              left_annotation = left_anno,
              right_annotation = right_anno,
              show_row_names = F,
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                # transform the matrix into a long vector
                nt = pindex(edit_seq_mat_sel, i, j)
                strand = edit_mat_sel$strand[i]
                pos = edit_mat_sel$pos[i]
                type = edit_mat_sel$type[i]
                cutpos = pos
                if(edit_mat_sel$dis[i] == 1){
                  cutpos = cutpos + 1
                }
                if(strand == "+" & type== "Mut"){
                  cutpos = cutpos - 1
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
                if(j == pos & type != "Mut"){
                  col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                  names(col) <- c("A", 'T','C','G')
                  grid.text(nt, x, y, 
                            gp = gpar(fontsize = 20, col = col[nt]))
                }else{
                  grid.text(nt, x, y, 
                            gp = gpar(fontsize = 18))
                }
                if(j == cutpos & pos != 0 & type != "1"){
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
                
                
              }, 
              # layer_fun = function(j, i, x, y, w, h, fill) {
              # ind_mat = restore_matrix(j, i, x, y)
              # for(ir in seq_len(nrow(ind_mat))) {
              #   ind1 = ind_mat[ir, 1] # previous column
              #   ind2 = ind_mat[ir, ncol(ind_mat)]
              #   grid.segments(x[ind1], y[ind1] - unit(as.numeric(h[ind1]) / 3, "npc"), x[ind2], y[ind2] - unit(as.numeric(h[ind2]) / 3, "npc"),
              #                 gp = gpar(col = "green", lwd = 2))
              # }
              # },
              col=c("white", "white"),show_heatmap_legend = F,
              border = F)
ht <- draw(ht, annotation_legend_list = lgd_list)
lens <- unique(edit_mat_sel$gene)
pos_line1 <- 0
##因为是从下到上画的 注意需要反过来画
edit_mat_sel  <- edit_mat_sel[nrow(edit_mat_sel) : 1,]
edit_mat_sel$sgname <- edit_mat_sel$type
edit_mat_sel$sgname[!edit_mat_sel$sgname %in% c("Mut","WT")] <- 
  paste0("Predict of", edit_mat_sel$type[!edit_mat_sel$sgname %in% c("Mut","WT")])
for(i in 1:length(unique(lens))) {
  len <- sum(edit_mat_sel$gene == lens[i])
  decorate_annotation("gene", slice = i, {
    tg <- richtext_grob(gt_render(lens[i]), 
                        rot = 0, 
                        x = unit(0.5, "npc"),
                        y=unit(0.5, "npc"), hjust = 0.5)
    grid.draw(tg)
    invisible(tg)
  })
  locs_count <- length(unique(edit_mat_sel$loca[edit_mat_sel$gene == lens[i]]))
  locs <- unique(edit_mat_sel$loca[edit_mat_sel$gene == lens[i]])
  sg_count <- sum(edit_mat_sel$gene == lens[i])
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
    
    ll <- sum(edit_mat_sel$loca[edit_mat_sel$gene == lens[i]] == l)
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
  sgs <- edit_mat_sel$sgname[edit_mat_sel$gene == lens[i]]
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
  decorate_annotation("% in Inframe", slice = i, {
    od = rev(od[[i]])
    grid.segments(inframe_pct_sel_append$indel[od] - inframe_pct_sel_append$se[od],
                  seq_along(od), 
                  inframe_pct_sel_append$indel[od] + inframe_pct_sel_append$se[od],
                  seq_along(od), default.units = "native")
    grid.segments(inframe_pct_sel_append$indel[od] - inframe_pct_sel_append$se[od], 
                  seq_along(od) - 0.1, 
                  inframe_pct_sel_append$indel[od] - inframe_pct_sel_append$se[od], 
                  seq_along(od) + 0.1, default.units = "native")
    grid.segments(inframe_pct_sel_append$indel[od] + inframe_pct_sel_append$se[od], 
                  seq_along(od) - 0.1, 
                  inframe_pct_sel_append$indel[od] + inframe_pct_sel_append$se[od], 
                  seq_along(od) + 0.1, default.units = "native")
    
  })
}



dev.off()









