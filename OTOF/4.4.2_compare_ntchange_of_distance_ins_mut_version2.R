###这个脚本依赖于5.2 4.4的统计结果
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
  seq_table$delPos <- unlist(lapply(1 : nrow(seq_table), function(i){
    editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
    refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
    editpos <- which(editseq == '-')
    editpos[which.min(abs(editpos - 20))]
  }))
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
edit_stat_of_distance <- edit_stat_of_distance[edit_stat_of_distance$delPos %in% c(20,21),]


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



# color_func = colorRamp2(c(0, 100), c("white", "#1432EB"))
left_anno <- 
  rowAnnotation(gene = anno_empty(border = F, width = unit(20, "mm")),
                loc = anno_empty(border = F, width = unit(20, "mm")),
                sg = anno_empty(border = F, width = unit(20, "mm"))
  )


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
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)

wt_0aa <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/wt_of_0aa_in_Del1.xlsx")

outside_sel <- abs(rowSums(dis_mat))
outside_sel[!names(outside_sel) %in% wt_0aa$id] <- 0
sp_info <- sp_info[sp_info$id %in% rownames(plot_mat),]
insert_mat <- insert_mat[outside_sel == 1,]
wt_0aa <- wt_0aa[wt_0aa$id %in% rownames(insert_mat),]
rownames(sp_info) <- sp_info$id
sp_info <- sp_info[wt_0aa$id,]
insert_mat <- insert_mat[sp_info$id,]
del1_pct <- del1_pct[sp_info$id,]
indel_stat <- indel_stat[sp_info$id,]
plot_mat <- plot_mat[sp_info$id,]
dis_mat <- dis_mat[sp_info$id,]
seq_mat <- seq_mat[sp_info$id,]


cols = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(cols) <- c("a", "t", "c", "g")
load("~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda")
inframe_spec <- reshape2::dcast(inframe_result_processed, id ~ aa, value.var = "pct2")
rownames(inframe_spec) <- inframe_spec[,1]
inframe_spec <- inframe_spec[,-1]
inframe_spec[is.na(inframe_spec)] <- 0
inframe_spec <- inframe_spec[sp_info$id, c("-2", "-1", "0", "1", "2", "Others")]

inframe_ratio <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/previous_tissue_inframe_in_indel.xlsx")

rownames(inframe_ratio) <- inframe_ratio$result
inframe_ratio <- inframe_ratio[sp_info$id,]

cols2 <- c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")
names(cols2) <- colnames(inframe_spec)
inframe_ratio$inframeSe[is.na(inframe_ratio$inframeSe)] <- 0




right_anno <- rowAnnotation(InsertNT = anno_barplot2(insert_mat, 
                                                     gp = gpar(fill = cols, lwd = 0.1), 
                                                     bar_width = 0.5, 
                                                     width = unit(2, "cm")), 
                            emp = anno_empty(width = unit(2, "mm"), border = F),
                            "Inframe Spec" = anno_barplot2(x = inframe_spec, width = unit(2, "cm"), 
                                                           gp = gpar(lwd = 0.1, fill = c(cols2)),
                                                           bar_width = 0.5, ylim = ceiling(c(0, 100))), 
                            emp2 = anno_empty(width = unit(2, "mm"), border = F),
                            "Inframe in Indel%" = anno_barplot(x = inframe_ratio$inframePct, 
                                                               gp = gpar(lwd = 0.1, fill = "#BBBBBB"),
                                                               width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 100)), 
                            "WT in 0AA%" = anno_barplot(x = wt_0aa$wtPct, 
                                                        gp = gpar(lwd = 0.1, fill = "#BBBBBB"),
                                                        width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 100)))
total_sample_info$cds_start <- unlist(lapply(1 : nrow(total_sample_info), function(index){
  str_locate(total_sample_info$seq[index], total_sample_info$cds[index])[1,1]
}))
total_sample_info$mutcds <- unlist(lapply(1 : nrow(total_sample_info), function(i){
  start <- total_sample_info$cds_start[i]
  mut <- total_sample_info$mutPos[i]
  while(start > mut){
    start = start - 3
  }
  while(T){
    if(start <= mut & (start + 2) >= mut){
      if(start == mut){
        return(0)
      }
      if((start + 1) == mut){
        return(1)
      }
      return(2)
    }
    start = start + 3
  }
  
}))
rownames(total_sample_info) <- total_sample_info$tissueid
sp_info$mutcds <- total_sample_info[sp_info$id, "mutcds"]

pdf("~/Nutstore Files/Tobin/Previous/Del1_edit_distance_add_WT_0aa.pdf", width = 11, height = 2)
ht <- Heatmap(plot_mat, 
              name = "EditDis",
              column_split = c(rep(" ", arround_len), rep("   ", arround_len)),
              column_title=NULL,
              row_title = NULL,
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
                v = pindex(plot_mat, i, j) 
                dis = pindex(dis_mat, i, j)
                nt = pindex(seq_mat, i, j)
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
                grid.polyline(x = c(x_from, x_from, x_to , x_to), 
                              y = c(y_from,y_to ,y_from, y_to), 
                              gp = gpar(col = "#D0D0D0"), 
                              id = rep(1 : (length(x) * 2), 2))
                
                # if(j == 7){
                #   col = c(2 : 5)
                #   names(col) <- c("A", 'T','C','G')
                #   grid.text(nt, x, y, 
                #             gp = gpar(fontsize = 20, col = col[nt]))
                # }else{
                #   grid.text(nt, x, y, 
                #             gp = gpar(fontsize = 18))
                # }
                grid.text(nt, x, y, 
                          gp = gpar(fontsize = 18))
                if(dis == 1){
                  
                  grid.polygon(x = c(x_from - unit(as.numeric(w) / 8, "npc"), 
                                     x_from, 
                                     x_from + unit(as.numeric(w) / 8, "npc")),
                               y = c(y + unit(as.numeric(h) / 3, "npc"), 
                                     y - unit(as.numeric(h) / 3, "npc"), 
                                     y + unit(as.numeric(h) / 3, "npc")), 
                               gp = gpar(col = "red", fill = "red"))
                  
                }
                if(dis == -1){
                  
                  grid.polygon(x = c(x_from - unit(as.numeric(w) / 8, "npc"), 
                                     x_from, 
                                     x_from + unit(as.numeric(w) / 8, "npc")),
                               y = c(y - unit(as.numeric(h) / 3, "npc"), 
                                     y + unit(as.numeric(h) / 3, "npc"), 
                                     y - unit(as.numeric(h) / 3, "npc")), 
                               gp = gpar(col = "red", fill = "red"))
                  
                }
                
              }, col=c("white", "blue"),show_heatmap_legend = F,
              border = F)
ht <- draw(ht, annotation_legend_list = lgd_list)
lens <- unique(sp_info$gene)
rownames(sp_info) <- sp_info$id
sp_info <- sp_info[rownames(plot_mat),]
pos_line1 <- 0
##因为是从下到上画的 注意需要反过来画
sp_info  <- sp_info[nrow(sp_info) : 1,]
decorate_heatmap_body("EditDis", {
  cell_width <- 1 / 6
  cell_height <- 1 / nrow(sp_info)
  for(i in 1 : nrow(sp_info)){
    mutcds <- sp_info$mutcds[i]
    grid.segments(cell_width * (6 - mutcds) - cell_width / 50 , cell_height * (i - 1) + cell_height / 20, 
                  cell_width * (9 - mutcds) + cell_width / 50, cell_height * (i - 1) + cell_height / 20,
                  gp = gpar(col = "green", lwd = 2))
  }
  
  
})

decorate_annotation("gene", slice = 1, {
  tg <- richtext_grob(gt_render(sp_info$gene), 
                      rot = 0, 
                      x = unit(0.5, "npc"),
                      y=unit((1 : nrow(sp_info)) / nrow(sp_info) - 
                               1 / nrow(sp_info) / 2, "npc"), hjust = 0.5)
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("loc", slice = 1, {
  tg <- richtext_grob(gt_render(sp_info$loc), 
                      rot = 0, 
                      x = unit(0.5, "npc"),
                      y=unit((1 : nrow(sp_info)) / nrow(sp_info) - 
                               1 / nrow(sp_info) / 2, "npc"), hjust = 0.5)
  grid.draw(tg)
  invisible(tg)
})
decorate_annotation("sg", slice = 1, {
  tg <- richtext_grob(gt_render(sp_info$sg), 
                      rot = 0, 
                      x = unit(0.5, "npc"),
                      y=unit((1 : nrow(sp_info)) / nrow(sp_info) - 
                               1 / nrow(sp_info) / 2, "npc"), hjust = 0.5)
  grid.draw(tg)
  invisible(tg)
})
od <- row_order(ht)
decorate_annotation("Inframe in Indel%", slice = 1, {
  od = rev(od)
  grid.segments(inframe_ratio$inframePct[od] - inframe_ratio$inframeSe[od],seq_along(od), 
                inframe_ratio$inframePct[od] + inframe_ratio$inframeSe[od],seq_along(od), default.units = "native")
  grid.segments(inframe_ratio$inframePct[od] - inframe_ratio$inframeSe[od], 
                seq_along(od) - 0.1, inframe_ratio$inframePct[od] - inframe_ratio$inframeSe[od], 
                seq_along(od) + 0.1, default.units = "native")
  grid.segments(inframe_ratio$inframePct[od] + inframe_ratio$inframeSe[od], 
                seq_along(od) - 0.1, inframe_ratio$inframePct[od] + inframe_ratio$inframeSe[od], 
                seq_along(od) + 0.1, default.units = "native")
})
decorate_annotation("WT in 0AA%", slice = 1, {
  od = rev(od)
  grid.segments(wt_0aa$wtPct[od] - wt_0aa$se[od],seq_along(od), 
                wt_0aa$wtPct[od] + wt_0aa$se[od],seq_along(od), default.units = "native")
  grid.segments(wt_0aa$wtPct[od] - wt_0aa$se[od], 
                seq_along(od) - 0.1, wt_0aa$wtPct[od] - wt_0aa$se[od], 
                seq_along(od) + 0.1, default.units = "native")
  grid.segments(wt_0aa$wtPct[od] + wt_0aa$se[od], 
                seq_along(od) - 0.1, wt_0aa$wtPct[od] + wt_0aa$se[od], 
                seq_along(od) + 0.1, default.units = "native")
})


dev.off()

