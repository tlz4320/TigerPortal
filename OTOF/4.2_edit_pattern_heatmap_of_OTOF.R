###这个脚本依赖4.1的seqs_table_split等结果
otof_sample_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)
otof_sample_info <- split_name <- data.frame(
  tissueid = otof_sample_info$`Tissue-w-W(Mix)`,
  mutation = otof_sample_info$GRNAa)
otof_sample_info <- merge(otof_sample_info, total_used_seq_tissue, by.x="tissueid", by.y="result")
otof_sample_info <- otof_sample_info[grep("OTOF", otof_sample_info$mutation),]
otof_sample_info <- otof_sample_info[otof_sample_info$tissueid %in% names(seqs_table_split),]
otof_edit_result <- seqs_table_split[otof_sample_info$tissueid]
arround_len <- 6

edit_stat_of_distance <- lapply(names(otof_edit_result), function(x){
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



edit_stat_of_distance <- merge(edit_stat_of_distance, otof_sample_info, by.x="id", by.y = 'tissueid')
edit_stat_of_distance$gene <- unlist(lapply(edit_stat_of_distance$mutation, function(x){
  unlist(strsplit(x, "[-]"))[2]
}))
edit_stat_of_distance$loca <- unlist(lapply(edit_stat_of_distance$mutation, function(x){
  unlist(strsplit(x, "[-]"))[3]
}))
edit_stat_of_distance$sg <- unlist(lapply(edit_stat_of_distance$mutation, function(x){
  unlist(strsplit(x, "[-]"))[4]
}))
edit_stat_of_distance <- edit_stat_of_distance[order(edit_stat_of_distance$gene, edit_stat_of_distance$loca,edit_stat_of_distance$sg, edit_stat_of_distance$Pct),]

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
sp_info <- edit_stat_of_distance[,c("id", "dis", "sg", "gene", "loca")]
sp_info <- sp_info[!duplicated(sp_info$id),]
for(i in 1 : nrow(sp_info)){
  if(sp_info$dis[i] < 0){
    dis <- as.character(sp_info$dis[i] + 1)
  }else{
    dis <- as.character(sp_info$dis[i])
  }
  id <- sp_info$id[i]
  dis_mat[dis, id] <- 1
}
dis_mat <- t(dis_mat)
seq_mat <- matrix('', nrow = arround_len * 2 + 1, ncol = length(unique(edit_stat_of_distance$id)))
rownames(seq_mat) <- as.character(c(-arround_len:arround_len))
colnames(seq_mat) <- unique(edit_stat_of_distance$id)

seq_arround_mut <- lapply(1 : nrow(otof_sample_info), function(i){
  seq <- unlist(strsplit(otof_sample_info$seq[i], "*"))
  mutpos <- otof_sample_info$mutPos[i]
  str_to_upper(seq[(mutpos - arround_len) : (mutpos + arround_len) ])
})
for(i in 1 : length(seq_arround_mut)){
  seq_mat[,otof_sample_info$tissueid[i]] <- seq_arround_mut[[i]]
}
seq_mat <- t(seq_mat)




color_func = colorRamp2(c(0, 100), c("white", "#1432EB"))
left_anno <- 
  rowAnnotation(gene = anno_empty(border = F, width = unit(20, "mm")),
                "line1" = 
                  anno_empty(border = FALSE, width = unit(2, "mm")),
                    loc = anno_empty(border = F, width = unit(20, "mm")),
                    "line2" = 
                      anno_empty(border = FALSE, width = unit(2, "mm")),
                    sg = anno_empty(border = F, width = unit(20, "mm"))
  )

insertNT_pct <- lapply(split(edit_stat_of_distance,edit_stat_of_distance$id), function(x){
  tmp <- lapply(split(x$Pct, x$editNT), sum)
  tmp <- data.frame(editNT = names(tmp), Pct = unlist(tmp))
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  tmp$id <- unique(x$id)
  tmp
})
insertNT_pct <- data.frame(do.call(rbind, insertNT_pct))
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
    sum(y$X.Reads.1[indel == 1]) / sum(y$X.Reads.1[indel != 0]) * 100
  }))
  se <- plotrix::std.error(indel)
  data.frame(indel = mean(indel), se = se)
})
ins1_pct <- data.frame(do.call(rbind, ins1_pct))
ins1_pct$id <- rownames(ins1_pct)
ins1_pct$se[is.na(ins1_pct$se)] <- 0
ins1_pct <- ins1_pct[unique(edit_stat_of_distance$id),]



right_anno <- rowAnnotation(InsertNT = anno_barplot(insert_mat, gp = gpar(fill = 2 : 5), 
                                                  bar_width = 0.5, 
                                                  width = unit(2, "cm")), 
                            emp = anno_empty(width = unit(2, "mm"), border = F),
                            "Indel%" = anno_barplot(x = indel_stat$indel, width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 8)), 
                            emp2 = anno_empty(width = unit(2, "mm"), border = F),
                            "Ins1 in Indel%" = anno_barplot(x = ins1_pct$indel, width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 100)))

lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16, 
         legend_gp = gpar(fill = 2:5))
)



pdf("Previews/Result/Ins1_edit_distance_heatmap_v4.pdf", width = 12, height = 4.5)

ht <- Heatmap(plot_mat, 
        name = "EditDis",
        column_split = c(rep(" ", arround_len),"  ", rep("   ", arround_len)),
        column_title=NULL,
        row_split = NULL,
        show_column_names = T, 
        cluster_rows = F, 
        cluster_columns = F, 
        column_names_rot = 0,
        column_names_centered = T,
        column_names_side = "top",
        column_order = as.character(-6 : 6),
        left_annotation = left_anno,
        right_annotation = right_anno,
        show_row_names = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          
          # transform the matrix into a long vector
          v = pindex(plot_mat, i, j) 
          dis = pindex(dis_mat, i, j)
          nt = pindex(seq_mat, i, j)
          # `j` here is also a vector with the same length of `v`
          col = color_func(v)
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

          if(j == 7){
            col = c(2 : 5)
            names(col) <- c("A", 'T','C','G')
            grid.text(nt, x, y, 
                      gp = gpar(fontsize = 20, col = col[nt]))
          }else{
            grid.text(nt, x, y, 
                      gp = gpar(fontsize = 18))
          }
          if(dis == 1){
     
              grid.polygon(x = c(x_from - unit(as.numeric(w) / 4, "npc"), 
                                 x_from, 
                                 x_from + unit(as.numeric(w) / 4, "npc")),
                            y = c(y + unit(as.numeric(h) / 4, "npc"), 
                                  y - unit(as.numeric(h) / 4, "npc"), 
                                  y + unit(as.numeric(h) / 4, "npc")), 
                            gp = gpar(col = "red"))
            
          }
          
          
        }, col=c("white", "blue"),show_heatmap_legend = F,
        border = F)
draw(ht, annotation_legend_list = lgd_list)
lens <- unique(sp_info$gene)
rownames(sp_info) <- sp_info$id
sp_info <- sp_info[rownames(plot_mat),]
pos_line1 <- 0
##因为是从下到上画的 注意需要反过来画
sp_info  <- sp_info[nrow(sp_info) : 1,]
for(i in 1:length(unique(lens))) {
  len <- sum(sp_info$gene == lens[i])
  decorate_annotation("gene", slice = i, {
    tg <- richtext_grob(gt_render(lens[i]), 
                        rot = 0, 
                        x = unit(0.5, "npc"),
                        y=unit(0.5, "npc"), hjust = 0.5)
    grid.draw(tg)
    invisible(tg)
  })
  locs_count <- length(unique(sp_info$loca[sp_info$gene == lens[i]]))
  locs <- unique(sp_info$loca[sp_info$gene == lens[i]])
  sg_count <- sum(sp_info$gene == lens[i])
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
    
    ll <- sum(sp_info$loca[sp_info$gene == lens[i]] == l)
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
  sgs <- sp_info$sg[sp_info$gene == lens[i]]
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
}
decorate_annotation("Indel%", {
  od = 6 : 1
  #pushViewport(viewport(xscale = c(0.5, length(od) + 0.5), yscale = c(0, 1)))
  grid.segments(indel_stat$indel[od] - indel_stat$se[od],seq_along(od), indel_stat$indel[od] + indel_stat$se[od],seq_along(od), default.units = "native")
  grid.segments(indel_stat$indel[od] - indel_stat$se[od], seq_along(od) - 0.1, indel_stat$indel[od] - indel_stat$se[od], seq_along(od) + 0.1, default.units = "native")
  grid.segments(indel_stat$indel[od] + indel_stat$se[od], seq_along(od) - 0.1, indel_stat$indel[od] + indel_stat$se[od], seq_along(od) + 0.1, default.units = "native")
  # grid.function(function(y) {list(y=y, x=0.2)},
  #               range="y", gp=gpar(lty=2) )
  #grid.yaxis()
  #popViewport()
})
decorate_annotation("Ins1 in Indel%", {
  od = 6 : 1
  #pushViewport(viewport(xscale = c(0.5, length(od) + 0.5), yscale = c(0, 1)))
  grid.segments(ins1_pct$indel[od] - ins1_pct$se[od],seq_along(od), ins1_pct$indel[od] + ins1_pct$se[od],seq_along(od), default.units = "native")
  grid.segments(ins1_pct$indel[od] - ins1_pct$se[od], seq_along(od) - 0.1, ins1_pct$indel[od] - ins1_pct$se[od], seq_along(od) + 0.1, default.units = "native")
  grid.segments(ins1_pct$indel[od] + ins1_pct$se[od], seq_along(od) - 0.1, ins1_pct$indel[od] + ins1_pct$se[od], seq_along(od) + 0.1, default.units = "native")
  # grid.function(function(y) {list(y=y, x=0.2)},
  #               range="y", gp=gpar(lty=2) )
  #grid.yaxis()
  #popViewport()
})

dev.off()

