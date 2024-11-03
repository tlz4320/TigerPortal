###先跑5.2再这个脚本依赖于4.4的统计结果
total_sample_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)
total_sample_info <- split_name <- data.frame(
  tissueid = total_sample_info$`Tissue-w-W(Mix)`,
  mutation = total_sample_info$GRNAa)
total_sample_info <- merge(total_sample_info, total_used_seq_tissue, by.x="tissueid", by.y="result")
total_sample_info <- total_sample_info[total_sample_info$tissueid %in% names(seqs_table_split),]
arround_len <- 15




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


####计算+1 -2 -5的编辑产物
#首先从5.2脚本的结果inframe_result_processed确认是否存在这三种编辑
sel_top3_indel <- lapply(sel_for_detail_plot$id, function(x){
  tmp <- inframe_result_processed[inframe_result_processed$id == x,]
  tmp <- tmp[order(tmp$pct, decreasing = T),]
  tmp <- tmp[tmp$indel %in% c(3, 0,-3),]
  tmp$indel <- tmp$indel - 1
  tmp
})
#确定都有 然后开始挑选这三种做预测

seqs_table_stat <- lapply(sel_for_detail_plot$id, function(x){
  res <- lapply(c(2, -1, -4), function(indel_len) {
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
      inframe_sum <- edit_table[(edit_table$n_inserted - edit_table$n_deleted) %% 3 == 2,]
      edit_table$Pct <- edit_table[,7] / sum(edit_table[,7]) * 100
      edit_table$inframePct <- edit_table[,7] / sum(inframe_sum[,7]) * 100
      edit_table <- lapply(split(edit_table, edit_table$Aligned_Sequence), 
                           function(x){
                             if(nrow(x) == 1){
                               return(x)
                             }
                             x[,9] <- sum(x[,9])
                             x[,10] <- sum(x[,10])
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
               seq = edit_table$Aligned_Sequence, Pct = edit_table$Pct, 
               inframePct = edit_table$inframePct, pctSe = edit_table$pctse, 
               inframeSe = edit_table$inframepctse)
  })
  res <- data.frame(do.call(rbind, res))
  res <- res[order(res$inframePct, decreasing = T),]
  res[1:4,]
})

sel_info <- total_used_seq_tissue[total_used_seq_tissue$result %in% 
                                    sel_for_detail_plot$id,]

seqs_table_stat <- data.frame(do.call(rbind, seqs_table_stat))
seqs_table_stat <- seqs_table_stat[order(seqs_table_stat$id, seqs_table_stat$indel),]
seqs_table_stat <- seqs_table_stat[seqs_table_stat$seq != "",]


orf_table <- orf_table[sel_for_detail_plot$id,]


edit_mat <- data.frame(seq = c(
  "cagctgacaccccgacgctccaggtgtgtc", 
  "cagctgacaccccggacgctccaggtgtgt",
  "CAGCTGACACCCCG-ACGCTCCAGGTGTGT",
  "CAGCTGACACCCCGG-CGCTCCAGGTGTGT",
  "CAGCTGACACCC----CGCTCCAGGTGTGT",
  "ctcaaccccaagctggtgggcaagctgaag", 
  "ctcaaccccaagctgggtgggcaagctgaa",
  "CTCAACCCCAAGC-GGGTGGGCAAGCTGAA",
  "CTCAACCCAAAGCT-GGTGGGCAAGCTGAA",
  "CTCAACCCCAAGCTTTGGGTGGGCAAGCTG",
  "ctcaaccccaagctggtgggcaagctgaag",
  "ctcaaccccaagctgggtgggcaagctgaa",
  "CTCAACCCCAAGC-GGGTGGGCAAGCTGAA",
  "CTCAACCCCAAGCT-GGTGGGCAAGCTGAA",
  "CTCAACCCCAAGCCCTGGGTGGGCAAGCTG",
  "cacaaggccaacgagacggatgaggacgac",
  "cacaaggccaacgaggacggatgaggacga",
  "CACAAGGCCAACG-GGACGGATGAGGACGA", 
  "CACAAGGCCAACGA-GACGGATGAGGACGA", 
  "CACAAGGCCAAC-AGGACGGATGAGGACGA"), 
  type = c("WT", "Mut", -1, -1, -4,
           "WT", "Mut", -1, -1, 2,
           "WT", "Mut", -1, -1, 2,
           "WT", "Mut",-1, -1, -1), 
  pos = c(0, 15, 0, 0, 0, 
          0, 15, 0, 0, 15,
          0, 15, 0, 0, 14, 
          0, 15, 0, 0, 0), 
  cutpos = c(15, 15, 0, 0, 0, 
             14, 15, 0, 0, 0, 
             13, 14, 0, 0, 0,
             13, 14, 0, 0, 0), 
  pct = c(0, 0, 46.373520 , 12.214824, 11.813294,
          0 ,0, 37.381085, 23.791478, 18.389368, 
          0 ,0 ,34.385988, 25.156141, 6.572421, 
          0, 0, 35.117649, 14.055234, 11.474134), 
  se = c(0,0, 0.7449722, 1.9568140, 0.5304478,
         0 ,0 ,3.6298122, 3.5305311, 1.5347223, 
         0, 0, 1.1518608, 3.9609222, 0,
         0 ,0, 1.7693998, 5.4771518, 1.9428667),
  id = rep(c("20-051", "62-051", "64-051", "104-051"), rep(5,4))
)



edit_mat <- edit_mat[edit_mat$id %in% c("20-051", "62-051"),]
tmp <- merge(edit_mat,sel_for_detail_plot, by = 'id')
# rownames(tmp) <- paste0(tmp$id, tmp$type)
# tmp <- tmp[paste0(edit_mat$id, edit_mat$type),]

edit_mat <- tmp
edit_mat <- edit_mat[order(edit_mat$geneFactor),]
edit_mat$seq <- str_to_upper(edit_mat$seq)
# edit_mat <- edit_mat[order(edit_mat$geneFactor),]
edit_seq_mat <- matrix("", nrow = nrow(edit_mat),ncol = 30)
rownames(edit_seq_mat) <- rownames(edit_mat)
for(i in 1 : nrow(edit_mat)){
  edit_seq_mat[i,] <- unlist(strsplit(str_to_upper(edit_mat$seq[i]), "*"))
}
fake_plot_mat <- matrix(runif(nrow(edit_mat) * 30), nrow = nrow(edit_mat), ncol = 30)

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

inframe_data <- edit_mat
right_anno <- rowAnnotation(
  "% in Inframe" = anno_barplot(x = inframe_data$pct, 
                                gp = gpar(lwd = 0.1, fill = "#BBBBBB"),width = unit(2, "cm"), 
                                bar_width = 0.5, ylim = c(0, max(inframe_data$pct + 
                                                                   inframe_data$se + 3))))
pdf("~/Nutstore Files/Tobin/Previous/Del1_edit_between(+-1)_result_v4.pdf", width = 11, height = 5)
ht <- Heatmap(fake_plot_mat, 
              name = "EditDis",
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
              right_annotation = right_anno,
              show_row_names = F,
              cell_fun = function(j, i, x, y, w, h, fill) {
                # print(i)
                
                # transform the matrix into a long vector
                nt = pindex(edit_seq_mat, i, j)
                strand = edit_mat$strand[i]
                
                pos = edit_mat$pos[i]
                cutpos = edit_mat$cutpos[i] + 1
                type = edit_mat$type[i]
                
                dis = edit_mat$dis[i]
                strand = edit_mat$strand[i]
                col2 = "white"
                grid.rect(x, y, w, h, gp = gpar(fill = col2, col = col2))
                
                # if(dis== 1){
                #   grid.points(x, y, gp = gpar(col = "red"))
                # }
                x_from <- x - unit(as.numeric(w) / 2, "npc")
                x_to <- x + unit(as.numeric(w) / 2, "npc")
                y_from <- y - unit(as.numeric(h) / 2, "npc")
                y_to <- y + unit(as.numeric(h) / 2, "npc")
                if(j == pos & type == "Mut"){
                  col = c(2 : 5)
                  names(col) <- c("A", 'T','C','G')
                  grid.text(nt, x, y, 
                            gp = gpar(fontsize = 20, col = col[nt]))

                } else if(type == 2 & j %in% c(pos, pos + 1)){
                  col = c(2 : 5)
                  names(col) <- c("A", 'T','C','G')
                  grid.text(nt, x, y, 
                            gp = gpar(fontsize = 20, col = col[nt]))
                }else{
                  grid.text(nt, x, y, 
                            gp = gpar(fontsize = 18))
                }
                if(j == cutpos && type %in% c("WT", "Mut")){
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
                
                
              },col=c("white", "white"),show_heatmap_legend = F,
              border = F)
ht <- draw(ht, annotation_legend_list = lgd_list)

pos_line1 <- 0
lens <- unique(edit_mat$gene)
##因为是从下到上画的 注意需要反过来画
edit_mat  <- edit_mat[nrow(edit_mat) : 1,]

edit_mat$sgname <- edit_mat$type
edit_mat$sgname[!edit_mat$sgname %in% c("Mut","WT")] <- 
  paste0("Predict of", edit_mat$type[!edit_mat$sgname %in% c("Mut","WT")])
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
  decorate_annotation("% in Inframe", slice = i, {
    od = rev(od[[i]])
    grid.segments(inframe_data$pct[od] - inframe_data$se[od],
                  seq_along(od), 
                  inframe_data$pct[od] + inframe_data$se[od],
                  seq_along(od), default.units = "native")
    grid.segments(inframe_data$pct[od] - inframe_data$se[od], 
                  seq_along(od) - 0.1, 
                  inframe_data$pct[od] - inframe_data$se[od], 
                  seq_along(od) + 0.1, default.units = "native")
    grid.segments(inframe_data$pct[od] + inframe_data$se[od], 
                  seq_along(od) - 0.1, 
                  inframe_data$pct[od] + inframe_data$se[od], 
                  seq_along(od) + 0.1, default.units = "native")
    
  })
}



dev.off()









