#这个文件是统计了所有Del1的蛋白质修复的结果

library(Biostrings)
library(stringr)
source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
seq_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/FunctioalAnalysisInput.xlsx", 3)
setwd("~/data/project/ear_project/gene_therapy_ll/Previews/")
reps <- list.files(pattern = "^Rep.*data")
load("~/data/project/ear_project/gene_therapy_ll/Previews/total_used_seq.rda")


load(file="Result/old_rep123_edit.rda")
print(load("~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data_newer_newer.rda"))
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
    ##选出Inframe的结果
    edit_table <- edit_table[(edit_table$n_inserted - edit_table$n_deleted) %% 3 == 2 & 
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
  seqs_table[[id]] <- edit_table
}

seqs_table <- data.frame(do.call(rbind, seqs_table))
seqs_table <- seqs_table[seqs_table$n_deleted == 0 | seqs_table$n_inserted == 0,]
seqs_table <- seqs_table[abs(seqs_table$n_inserted - seqs_table$n_deleted + 1) / 3 <= 1,]
# need_sg <- openxlsx::read.xlsx("~/data/project/ear_project/gene_therapy_ll/Previews/volcano_+-3_sel_sg.xlsx")
# seqs_table <- seqs_table[seqs_table$id %in% need_sg$id,]
seqs_table_split <- split(seqs_table, seqs_table$id)
seqs_table_stat <- lapply(names(seqs_table_split), function(x){
  seq_table <- seqs_table_split[[x]]
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
  mutNT <- total_used_seq_tissue$mutNT[total_used_seq_tissue$result == x]
  
  if(strand == '-'){
    sg <- as.character(reverseComplement(DNAString(sg)))
  } 
  if(strand == '-'){
    if(mutpos > (cutsite - 2) & mutpos <= (cutsite + str_length(sg) - 3)){
      sg <- unlist(strsplit(sg, "*"))
      sg <- paste(c(sg[ 1 : (mutpos - (cutsite - 2))], mutNT, sg[ (mutpos - (cutsite - 2) + 1) : length(sg)]), collapse = "")
    }
  }else{
    if(mutpos > (cutsite - str_length(sg) + 4) & mutpos <= (cutsite + 3)){
      sg <- unlist(strsplit(sg, "*"))
      sg <- paste(c(sg[ 1 : (mutpos - (cutsite - length(sg) + 4))], mutNT, sg[ (mutpos - (cutsite - length(sg) + 4) + 1) : length(sg)]), collapse = "")
    }
  }
  seq_table$predict <- unlist(lapply(1 : nrow(seq_table), function(i){
    editseq <- unlist(strsplit(seq_table$Aligned_Sequence[i], "*"))
    refseq <- unlist(strsplit(seq_table$Reference_Sequence[i], "*"))
    findRef <- Biostrings::matchPattern(pattern= str_to_lower(sg),
                                        subject = mutseq, max.mismatch = 3, with.indels = T)
    mutseqs <- unlist(strsplit(mutseq, "*"))
    hasDelete <- 0
    for(editpos in 10 : 30){
      if(editseq[editpos] == "-" | refseq[editpos] == "-"){
        isInsert <- (refseq[editpos] == "-")
        editNT <- str_to_lower(editseq[editpos])
        if(strand == "+"){
          insertpos <- editpos  - 20  + findRef@ranges@start + findRef@ranges@width - 4 + hasDelete
        }else{
          insertpos <- editpos  - 20  + findRef@ranges@start + 2 + hasDelete
        }
        if(isInsert){
          mutseqs <- c(mutseqs[1 : (insertpos - 1)], editNT, mutseqs[insertpos : length(mutseqs)])
        } else {
          mutseqs <- c(mutseqs[1 : (insertpos - 1)], mutseqs[(insertpos + 1) : length(mutseqs)])
          hasDelete <- hasDelete - 1
        }
      }
    }
    return(paste(mutseqs, collapse = ""))
  }))
  seq_table
})
seqs_table_stat <- data.frame(do.call(rbind, seqs_table_stat))
total_used_seq_tissue$cds_start <- unlist(lapply(1 : nrow(total_used_seq_tissue), function(index){
  str_locate(total_used_seq_tissue$seq[index], total_used_seq_tissue$cds[index])[1,1]
}))
seqs_table_stat$aa <- (seqs_table_stat$n_inserted - seqs_table_stat$n_deleted + 1) / 3
total_used_seq_tissue$realAA <- unlist(lapply(1 : nrow(total_used_seq_tissue), function(i){
  paste(seqinr::translate(seqinr::s2c(total_used_seq_tissue$cds[i]), ambiguous = T), collapse = "")
}))
seqs_table_stat$predictAA <- unlist(lapply(1 : nrow(seqs_table_stat), function(i){
  cds <- seqs_table_stat$predict[i]
  realcds <- total_used_seq_tissue$cds[total_used_seq_tissue$result == seqs_table_stat$id[i]]
  cds <- str_sub(cds, total_used_seq_tissue$cds_start[total_used_seq_tissue$result == seqs_table_stat$id[i]])
  paste(seqinr::translate(seqinr::s2c(cds), ambiguous = T), collapse = "")
}))

seqs_table_stat_split <- split(seqs_table_stat, seqs_table_stat$id)
seqs_table_stat_split <- lapply(seqs_table_stat_split, function(x){
  x$Pct2 <- x$Pct / sum(x$Pct) * 100
  x
})
seqs_table_stat_split <- do.call(rbind, seqs_table_stat_split)
seqs_table_stat_split$aadiff <- 0

# seqs_table_stat_split <- seqs_table_stat_split[seqs_table_stat_split$id %in% need_sg$id,]
library(parallel)
par_result <- mclapply(1 : nrow(seqs_table_stat_split), function(i){
  predictAA <- seqs_table_stat_split$predictAA[i]
  realAA <- total_used_seq_tissue$realAA[total_used_seq_tissue$result == seqs_table_stat_split$id[i]]
  aa <- seqs_table_stat_split$aa[i]
  if(str_length(predictAA) + abs(aa) < str_length(realAA)){
    print(i)
    print("Why?")
  }
  if(aa < 0){
    predictAA <- str_sub(predictAA, 1, str_length(realAA) - 1)
  }
  if(aa > 0){
    predictAA <- str_sub(predictAA, 1, str_length(realAA) + 1)
  }
  if(aa == 0){
    predictAA <- str_sub(predictAA, 1, str_length(realAA))
  }
  comp <- Biostrings::pairwiseAlignment(predictAA, realAA)
  pt <- Biostrings::alignedPattern(comp)
  sj <- Biostrings::alignedSubject(comp)
  pt <- unlist(strsplit(as.character(pt), "*"))
  sj <- unlist(strsplit(as.character(sj), "*"))
  sum(pt != sj) - abs(aa)
}, mc.cores = 10)
#这里太慢了
seqs_table_stat_split$aadiff <- unlist(par_result)
# save(seqs_table_stat_split,total_used_seq_tissue, file="~/Nutstore Files/Tobin/Previous/AA_stat_ins1Mut.rda")

seqs_table_stat_split$isCut <- unlist(lapply(1 : nrow(seqs_table_stat_split), function(x){
  isInsert <- seqs_table_stat_split$n_inserted[x] != 0
  if(isInsert){
    pos <- unlist(strsplit(seqs_table_stat_split$Reference_Sequence[x], "*"))
    pos <- which(pos == "-")
    pos <- pos[which.min(abs(pos - 20))]
  } else {
    pos <- unlist(strsplit(seqs_table_stat_split$Aligned_Sequence[x], "*"))
    pos <- which(pos == "-")
    pos <- pos[which.min(abs(pos - 20))]
  }
  pos %in% c(20, 21)
  
}))

table(seqs_table_stat_split$aadiff)
tmp <- seqs_table_stat_split[seqs_table_stat_split$n_inserted == 4,]
seqs_table_stat_split_sel <- seqs_table_stat_split[seqs_table_stat_split$aadiff < 2,]

seqs_table_stat_sel <- split(seqs_table_stat_split_sel, seqs_table_stat_split_sel$id)
need_sg <- openxlsx::read.xlsx("~/data/project/ear_project/gene_therapy_ll/Previews/volcano_+-3_sel_sg_Ins1Mut.xlsx")
need_sg <- need_sg[order(need_sg$diff),]
seqs_table_stat_sel <- seqs_table_stat_sel[need_sg$id]
seqs_table_stat_sel_mat <- lapply(seqs_table_stat_sel, function(x){
  x <- x[x$isCut,]
  result <- matrix(0, nrow = 2, ncol = 3)
  colnames(result) <- c("-1", "0", "1")
  rownames(result) <- c("0", "1")
  for(i in 1 : nrow(x)){
    result[as.character(x$aadiff[i]), as.character(x$aa[i])] <- 
      result[as.character(x$aadiff[i]), as.character(x$aa[i])] + 
      x$Pct2[i]
  }
  result
})
has_wt <- unlist(lapply(seqs_table_stat_sel_mat, function(x){
  return(x["0","0"] > 1)
}))
need_sg <- need_sg[has_wt,]
seqs_table_stat_sel_mat <- seqs_table_stat_sel_mat[has_wt]
names(seqs_table_stat_sel_mat) == need_sg$id
plot_mat <- do.call(rbind, seqs_table_stat_sel_mat)
rowSplit <- unlist(lapply(1 : nrow(need_sg), function(i){
  rep(paste(rep(" ", i), collapse = ""), 2)
}))

left_anno <- rowAnnotation(
  sgname1 = anno_empty(width = unit(5, "cm"), border = F),
  left1 = anno_empty(width = unit(2, "cm"), border = F), 
  show_annotation_name = F)

pdf("~/Nutstore Files/Tobin/Previous/AA_counts_ID_heatmap_Ins1Mut.pdf", width = 6, height = 4)
Heatmap(plot_mat, cluster_columns = F, cluster_rows = F, 
        show_row_names = F, show_column_names = F,
        row_order = 1 : nrow(plot_mat),
        cell_fun = function(j, i, x, y, w, h, fill){
          val <- pindex(plot_mat, i, j)
          grid.text(round(val, 2), x, y)
        },
        left_annotation = left_anno,
        top_annotation = HeatmapAnnotation(AA = anno_empty(height = unit(1, "cm"), border = F)),
        row_split = rowSplit, row_title = NULL, col = c("white", "red"))

for(i in 1 : nrow(need_sg)){
  decorate_annotation("sgname1",slice = i, {
    
    tg <- textGrob(gt_render(need_sg$mutation[i]), 
                   rot = 0, 
                   x = unit(0.4, "npc"),
                   y=unit(0.5, "npc"), hjust = 0.5, 
                   gp = gpar(fontsize = 10, col = "black"))
    grid.draw(tg)
    invisible(tg)
  })
  decorate_annotation("left1",slice = i, {
    
    tg <- textGrob(gt_render(paste0("Identity Change:", c(1, 0))), 
                   rot = 0, 
                   x = unit(0.2, "npc"),
                   y=unit(c(0.25, 0.75), "npc"), hjust = 0.5, 
                   gp = gpar(fontsize = 10, col = "black"))
    grid.draw(tg)
    invisible(tg)
  })
}
decorate_annotation("AA",slice = 1, {
  
  tg <- textGrob(gt_render(c(-1, 0, 1)), 
                 rot = 0, 
                 x = unit(c(1/ 6, 1 / 2, 5/6), "npc"),
                 y=unit(0.5, "npc"), hjust = 0.5, 
                 gp = gpar(fontsize = 15, col = "black"))
  grid.draw(tg)
  invisible(tg)
})
dev.off()


merged_sample_aa <- lapply(seqs_table_stat_sel_mat, function(x){
  counts <- sum(colSums(x) * c(1, 0, 1)) / 100
  ids <- sum(rowSums(x) * c(0, 1)) / 100
  data.frame(counts = counts, ids = ids)
})
merged_sample_aa <- data.frame(do.call(rbind, merged_sample_aa))
merged_sample_aa$id <- need_sg$mutation
merged_sample_aa <- merged_sample_aa[order(merged_sample_aa$ids),]
merged_sample_aa$color <- "black"
merged_sample_aa$color[1:3] <- 'red'
pdf("~/Nutstore Files/Tobin/Previous/AA_counts_ID_one_mat_Ins1Mut.pdf", width = 7.5, height = 6)
ggplot(merged_sample_aa, aes(x = counts, y = ids)) + 
  geom_point(color = merged_sample_aa$color) + 
  geom_text_repel(aes(label = id), min.segment.length = unit(0, 'lines')) + 
  xlab("Number Change") + ylab("Identity Change") + 
  theme_bw()
dev.off()



plot_data2 <- lapply(seqs_table_stat_sel_mat, function(x){
  ids <- sum(rowSums(x) * c(0, 1)) / 100
  x <- unlist(x[,"0"])
  x <- x / sum(x) * 100
  
  data.frame(counts = x[1] / 100, ids = ids)
})
plot_data2 <- data.frame(do.call(rbind, plot_data2))
plot_data2$id <- need_sg$mutation
plot_data2 <- plot_data2[order(plot_data2$ids),]
plot_data2$color <- "black"
plot_data2$color[1:3] <- 'red'
pdf("~/Nutstore Files/Tobin/Previous/AA_counts_ID_one_mat2_Ins1Mut.pdf", width = 7.5, height = 6)
ggplot(plot_data2, aes(x = counts, y = ids)) + 
  geom_point(color = plot_data2$color) + 
  geom_text_repel(aes(label = id), min.segment.length = unit(0, 'lines')) + 
  xlab("WT% In Del1") + ylab("Total Identity Change") + 
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1)) + 
  theme_bw()
dev.off()





print(load("~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda"))
tissue_0aa <- inframe_result_processed[inframe_result_processed$aa == "0",]

tissue_inframe <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/previous_tissue_inframe_in_indel.xlsx")
tissue_inframe <- tissue_inframe[,c("result", "inframePct")]


rownames(plot_data2) <- plot_data2$id
rownames(need_sg) <- need_sg$mutation
plot_data2$id2 <- need_sg[plot_data2$id, "id"]
plot_data3 <- merge(plot_data2, tissue_inframe, by.x="id2", by.y="result")
plot_data3 <- merge(plot_data3, tissue_0aa, by.x="id2", by.y="id")

#把不能预测的删掉
plot_data3 <- plot_data3[!plot_data3$id2 %in% c("m12-tissue", "22-051", "37-051", "65-05", "33-051"),]

pdf("~/Nutstore Files/Tobin/Previous/AA_counts_ID_one_mat_Ins1_v2.pdf", width = 7.5, height = 6)
ggplot(plot_data3, aes(x = counts, y = inframePct)) + 
  geom_point(aes(fill = pct2), shape = 21, size = 3, strock = 1) + 
  scale_color_gradientn(colours = c("blue",  "red"),limits = c(0,100),
                        name = "Inframe%") + theme_bw() + 
  ggnewscale::new_scale_colour() + 
  scale_fill_gradientn(colours = c("#FFEBCD", "#F09B28"), name = "0AA/Inframe") + 
  geom_text_repel(aes(label = id), min.segment.length = unit(0, 'lines')) + 
  xlab("WT% In Del1") + ylab("Inframe%") + 
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,100)) + 
  theme_bw()
dev.off()
