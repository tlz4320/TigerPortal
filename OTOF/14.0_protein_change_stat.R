#这个文件是统计了所有Ins1的蛋白质修复的结果

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
    ##选出Inframe的结果
    edit_table <- edit_table[(edit_table$n_inserted - edit_table$n_deleted) %% 3 == 1 & 
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
seqs_table <- seqs_table[abs(seqs_table$n_inserted - seqs_table$n_deleted - 1) / 3 <= 1,]
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
  if(strand == '-'){
    sg <- as.character(reverseComplement(DNAString(sg)))
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
seqs_table_stat$aa <- (seqs_table_stat$n_inserted - seqs_table_stat$n_deleted - 1) / 3
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
# save(seqs_table_stat_split,total_used_seq_tissue, file="~/Nutstore Files/Tobin/Previous/AA_stat.rda")

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
need_sg <- openxlsx::read.xlsx("~/data/project/ear_project/gene_therapy_ll/Previews/volcano_+-3_sel_sg.xlsx")
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

pdf("~/Nutstore Files/Tobin/Previous/AA_counts_ID_heatmap_v2.pdf", width = 6, height = 8)
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
pdf("~/Nutstore Files/Tobin/Previous/AA_counts_ID_one_mat_v2.pdf", width = 7.5, height = 6)
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

pdf("~/Nutstore Files/Tobin/Previous/AA_counts_ID_one_mat_v5.pdf", width = 7.5, height = 6)
ggplot(plot_data3, aes(x = counts, y = inframePct)) + 
  geom_point(aes(fill = pct2), shape = 21, size = 3, strock = 1) + 
  scale_color_gradientn(colours = c("blue",  "red"),limits = c(0,100),
                        name = "Inframe%") + theme_bw() + 
  ggnewscale::new_scale_colour() + 
  scale_fill_gradientn(colours = c("#FFEBCD", "#F09B28"), name = "0AA/Inframe") + 
  geom_text_repel(aes(label = id), min.segment.length = unit(0, 'lines')) + 
  xlab("WT% In Ins1") + ylab("Inframe%") + 
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,100)) + 
  theme_bw()
dev.off()
