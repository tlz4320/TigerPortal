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

seq_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/FunctioalAnalysisInput.xlsx", 3)
print(load("~/Nutstore Files/Tobin/Previous/tota_log.rda"))
library(stringr)
library(Biostrings)

total_used_seq_tissue <- total_used_seq[grep("051", total_used_seq$result),]
otof_sg12_15_info <- total_used_seq_tissue[rep(which(total_used_seq_tissue$result == "105-051"), 4),]
otof_sg12_15_info$result <- paste0("m", 12:15, "-tissue")
otof_sg12_15_info$sg <- str_to_upper(c("gcccactgccgttcggg", "tgcccactgccgttcgg",
                                       "gtgcccactgccgttcg", "cgtgcccactgccgttc"))
otof_sg12_15_info$cutSite <- c(41 : 44)
otof_sg12_15_info$strand <- "-"
total_used_seq_tissue <- data.frame(rbind(total_used_seq_tissue, otof_sg12_15_info))

load(file="Result/old_rep123_edit.rda")
print(load("~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data_newer_newer.rda"))


rep1_new_otof <- otof_cell_sg12_15_edit[seq(1,8, 2)]
rep2_new_otof <- otof_cell_sg12_15_edit[seq(2,8, 2)]
for(name in names(rep1_new_otof)){
  name2 <- str_remove(name, "[-].*")
  old_rep1_edit[[name2]] <- rep1_new_otof[[name]]
}
for(name in names(rep2_new_otof)){
  name2 <- str_remove(name, "[-].*")
  old_rep2_edit[[name2]] <- rep2_new_otof[[name]]
}

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

tmp <- lapply(total_used_seq_tissue$result, function(id){
  info <- total_used_seq_tissue[total_used_seq_tissue$result == id,]
  tables <- total_old_edit_table_rev[[id]]
  res <- lapply(tables, function(tmp){
    indel <- tmp$n_deleted + tmp$n_inserted
    tmp <- tmp[indel != 0,]
    if(nrow(tmp) == 0){
      return(data.frame(fq = 0, pct = 0))
    }
    tmp$inframe <- tmp$n_inserted - tmp$n_deleted
    total_reads <- sum(tmp[,7])
    if(info$mutType == "d"){
      tmp <- tmp[tmp$inframe %% 3 == 1, ]
      tmp$inframe <- tmp$inframe - 1
    } else {
      tmp <- tmp[tmp$inframe %% 3 == 2, ]
      tmp$inframe <- tmp$inframe + 1
    }
    if(nrow(tmp) == 0){
      return(data.frame(fq = 0, pct = 0))
    }
    return(data.frame(fq = sum(tmp[,7]), pct = sum(tmp[,7]) / total_reads) * 100)
  })
  res <- data.frame(do.call(rbind, res))
  res <- res[res$pct != 0,]
  if(nrow(res) == 0){
    return(data.frame(fq = 0, pct = 0, se = 0))
  }
  res$se <- plotrix::std.error(res$pct)
  res$fq <- sum(res$fq) / length(tables)
  res$pct <- sum(res$pct) / length(tables)
  res[1,]
})
tmp <- data.frame(do.call(rbind, tmp))
total_used_seq_tissue$inframePct <- tmp$pct
total_used_seq_tissue$inframeSe <- tmp$se
plot_data <- hist(total_used_seq_tissue$inframePct, breaks = seq(from=0, to=93, by=3))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
plot_data$color <- "#FAB0B0"
plot_data$color[plot_data$inframePct > 30] <-"#FA1E1E" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 1.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c("Total Counts:118", paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")
pdf("~/Nutstore Files/Tobin/Previous/total_inframe_hist_v2.pdf", width = 6, height = 4)
ggplot(plot_data) + geom_bar(aes(x = x, y = counts), 
                             color = "black",
                             stat="identity", fill = plot_data$color) + 
  geom_text(x = 80, y = 8, label = anno_info, data = data.frame(a=1)) + 
  xlab("Inframe Frequency") + ylab("Counts") + 
  scale_y_continuous(breaks = 1 : 10) + 
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90)) + 
  theme_bw() + 
  theme(panel.grid = element_blank())
dev.off()

openxlsx::write.xlsx(total_used_seq_tissue, file="~/Nutstore Files/Tobin/Previous/previous_tissue_inframe_in_indel.xlsx", rowNames=F, colNames=T)


#还有细胞的数据  之前没统计 重新统计一次
total_used_seq_cell <- total_used_seq[grep("001", total_used_seq$result),]
otof_sg12_15_info <- total_used_seq_cell[rep(which(total_used_seq_cell$result == "105-001"), 4),]
otof_sg12_15_info$result <- paste0("m", 12:15)
otof_sg12_15_info$sg <- str_to_upper(c("gcccactgccgttcggg", "tgcccactgccgttcgg",
                                       "gtgcccactgccgttcg", "cgtgcccactgccgttc"))
otof_sg12_15_info$cutSite <- c(41 : 44)
otof_sg12_15_info$strand <- "-"
total_used_seq_cell <- data.frame(rbind(total_used_seq_cell, otof_sg12_15_info))
tmp <- lapply(total_used_seq_cell$result, function(id){
  info <- total_used_seq_cell[total_used_seq_cell$result == id,]
  tables <- total_old_edit_table_rev[[id]]
  res <- lapply(tables, function(tmp){
    indel <- tmp$n_deleted + tmp$n_inserted
    tmp <- tmp[indel != 0,]
    if(nrow(tmp) == 0){
      return(data.frame(fq = 0, pct = 0))
    }
    tmp$inframe <- tmp$n_inserted - tmp$n_deleted
    total_reads <- sum(tmp[,7])
    if(info$mutType == "d"){
      tmp <- tmp[tmp$inframe %% 3 == 1, ]
      tmp$inframe <- tmp$inframe - 1
    } else {
      tmp <- tmp[tmp$inframe %% 3 == 2, ]
      tmp$inframe <- tmp$inframe + 1
    }
    if(nrow(tmp) == 0){
      return(data.frame(fq = 0, pct = 0))
    }
    return(data.frame(fq = sum(tmp[,7]), pct = sum(tmp[,7]) / total_reads) * 100)
  })
  res <- data.frame(do.call(rbind, res))
  res <- res[res$pct != 0,]
  if(nrow(res) == 0){
    return(data.frame(fq = 0, pct = 0, se = 0))
  }
  res$se <- plotrix::std.error(res$pct)
  res$fq <- sum(res$fq) / length(tables)
  res$pct <- sum(res$pct) / length(tables)
  res[1,]
})
tmp <- data.frame(do.call(rbind, tmp))
total_used_seq_cell$inframePct <- tmp$pct
total_used_seq_cell$inframeSe <- tmp$se
openxlsx::write.xlsx(total_used_seq_cell, file="~/Nutstore Files/Tobin/Previous/previous_cell_inframe_in_indel.xlsx", rowNames=F, colNames=T)
plot_data <- hist(total_used_seq_cell$inframePct, breaks = seq(from=0, to=93, by=3))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
plot_data$color <- "#FAB0B0"
plot_data$color[plot_data$inframePct > 30] <-"#FA1E1E" 
#plot_data$color[plot_data$inframePct > 50] <-"#FA1E1E" 
plot_data$x <- plot_data$inframePct - 1.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c("Total Counts:118", paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")
pdf("~/Nutstore Files/Tobin/Previous/total_cell_inframe_hist_v2.pdf", width = 6, height = 4)
ggplot(plot_data) + geom_bar(aes(x = x, y = counts), 
                             color = "black",
                             stat="identity", fill = plot_data$color) + 
  geom_text(x = 80, y = 8, label = anno_info, data = data.frame(a=1)) + 
  xlab("Inframe Frequency") + ylab("Counts") + 
  scale_y_continuous(breaks = 1 : 10) + 
  scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90)) + 
  theme_bw() + 
  theme(panel.grid = element_blank())
dev.off()

inframe_result <- lapply(total_used_seq_tissue$result, function(id){
  info <- total_used_seq_tissue[total_used_seq_tissue$result == id,]
  tables <- total_old_edit_table_rev[[id]]
  res <- lapply(tables, function(tmp){
    indel <- tmp$n_deleted + tmp$n_inserted
    tmp <- tmp[indel != 0,]
    if(nrow(tmp) == 0){
      return(data.frame(indel = 0, fq = 0, pct = 0))
    }
    tmp$inframe <- tmp$n_inserted - tmp$n_deleted
    total_reads <- sum(tmp[,7])
    if(info$mutType == "d"){
      tmp <- tmp[tmp$inframe %% 3 == 1, ]
      tmp$inframe <- tmp$inframe - 1
    } else {
      tmp <- tmp[tmp$inframe %% 3 == 2, ]
      tmp$inframe <- tmp$inframe + 1
    }
    if(nrow(tmp) == 0){
      return(data.frame(indel = 0, fq = 0, pct = 0))
    }
    tmp <- lapply(split(tmp[,7], tmp$inframe), sum)
    tmp <- data.frame(indel = as.numeric(names(tmp)), fq = unlist(tmp))
    tmp$pct <- tmp$fq / total_reads * 100
    tmp
  })
  res <- data.frame(do.call(rbind, res))
  res <- res[res$pct != 0,]
  res$aa <- as.numeric(res$indel) / 3
  res <- lapply(split(res, res$aa), function(tmp){
    tmp[,2] <- sum(tmp[,2]) / length(tables)
    tmp[,3] <- sum(tmp[,3]) / length(tables)
    tmp[1,]
  })
  res <- data.frame(do.call(rbind, res))
  res
})
names(inframe_result) <- total_used_seq_tissue$result
inframe_result_processed <- lapply(inframe_result, function(x){
  if(nrow(x) == 0)
    return(x)
  x <- x[order(abs(x$aa)),]
  x2 <- x
  x2$aa[abs(x2$aa) > 2] <- "Others"
  x2$fq[abs(x$aa) > 2] <- sum(x2$fq[abs(x$aa) > 2])
  x2$pct[abs(x$aa) > 2] <- sum(x2$pct[abs(x$aa) > 2])
  if(any(abs(x$aa) > 2)){
    x <- x2[c(which(abs(x$aa) <= 2), which(abs(x$aa) > 2)[1]),] 
  }
  x
})


inframe_result_processed <- lapply(total_used_seq_tissue$result, function(id){
  tmp <- inframe_result_processed[[id]]
  if(nrow(tmp) == 0){
    tmp <- data.frame(indel=0,
                      fq=0,
                      pct=0,
                      aa=0)
  }
  tmp$pct <- tmp$pct / sum(tmp$pct) * 100
  tmp$id <- id
  tmp
})
names(inframe_result_processed) <- total_used_seq_tissue$result
inframe_result_processed <- data.frame(do.call(rbind, inframe_result_processed))
ids <- gtools::mixedsort(unique(inframe_result_processed$id))
inframe_result_processed$id <- factor(inframe_result_processed$id, levels = rev(ids))
plot_data <- inframe_result_processed
colnames(plot_data)[4] <- "Number of AAChange"
pdf("~/Nutstore Files/Tobin/Previous/total_tissue_inframe_aa_change.pdf", width = 6, height = 12)
ggplot(plot_data) + 
  geom_bar(aes(x = pct, y = id, fill = `Number of AAChange`), 
           stat = "identity", position = "stack") + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

# pdf("~/Nutstore Files/Tobin/Previous/total_inframe_aa_change_pie.pdf", width = 12, height = 12)
# ggplot(data=plot_data, aes(x = 0, y = pct, fill=`Number of AAChange`)) + 
#   geom_col(position = "fill") + 
#   facet_wrap(~id, ncol = 10) +
#   coord_polar(theta = "y") +
#   theme_void()
# dev.off()





cell_inframe_result <- lapply(total_used_seq_cell$result, function(id){
  info <- total_used_seq_cell[total_used_seq_cell$result == id,]
  tables <- total_old_edit_table_rev[[id]]
  res <- lapply(tables, function(tmp){
    indel <- tmp$n_deleted + tmp$n_inserted
    tmp <- tmp[indel != 0,]
    if(nrow(tmp) == 0){
      return(data.frame(indel = 0, fq = 0, pct = 0))
    }
    tmp$inframe <- tmp$n_inserted - tmp$n_deleted
    total_reads <- sum(tmp[,7])
    if(info$mutType == "d"){
      tmp <- tmp[tmp$inframe %% 3 == 1, ]
      tmp$inframe <- tmp$inframe - 1
    } else {
      tmp <- tmp[tmp$inframe %% 3 == 2, ]
      tmp$inframe <- tmp$inframe + 1
    }
    if(nrow(tmp) == 0){
      return(data.frame(indel = 0, fq = 0, pct = 0))
    }
    tmp <- lapply(split(tmp[,7], tmp$inframe), sum)
    tmp <- data.frame(indel = as.numeric(names(tmp)), fq = unlist(tmp))
    tmp$pct <- tmp$fq / total_reads * 100
    tmp
  })
  res <- data.frame(do.call(rbind, res))
  res <- res[res$pct != 0,]
  res$aa <- as.numeric(res$indel) / 3
  res <- lapply(split(res, res$aa), function(tmp){
    tmp[,2] <- sum(tmp[,2]) / length(tables)
    tmp[,3] <- sum(tmp[,3]) / length(tables)
    tmp[1,]
  })
  res <- data.frame(do.call(rbind, res))
  res
})
names(cell_inframe_result) <- total_used_seq_cell$result
cell_inframe_result_processed <- lapply(cell_inframe_result, function(x){
  if(nrow(x) == 0)
    return(x)
  x <- x[order(abs(x$aa)),]
  x2 <- x
  x2$aa[abs(x2$aa) > 2] <- "Others"
  x2$fq[abs(x$aa) > 2] <- sum(x2$fq[abs(x$aa) > 2])
  x2$pct[abs(x$aa) > 2] <- sum(x2$pct[abs(x$aa) > 2])
  if(any(abs(x$aa) > 2)){
    x <- x2[c(which(abs(x$aa) <= 2), which(abs(x$aa) > 2)[1]),] 
  }
  x
})


cell_inframe_result_processed <- lapply(total_used_seq_cell$result, function(id){
  tmp <- cell_inframe_result_processed[[id]]
  if(nrow(tmp) == 0){
    tmp <- data.frame(indel=0,
                      fq=0,
                      pct=0,
                      aa=0)
  }
  tmp$pct <- tmp$pct / sum(tmp$pct) * 100
  tmp$id <- id
  tmp
})
names(cell_inframe_result_processed) <- total_used_seq_cell$result
cell_inframe_result_processed <- data.frame(do.call(rbind, cell_inframe_result_processed))
ids <- gtools::mixedsort(unique(cell_inframe_result_processed$id))
cell_inframe_result_processed$id <- factor(cell_inframe_result_processed$id, levels = rev(ids))
plot_data <- cell_inframe_result_processed
colnames(plot_data)[4] <- "Number of AAChange"

pdf("~/Nutstore Files/Tobin/Previous/total_cell_inframe_aa_change.pdf", width = 6, height = 12)
ggplot(plot_data) + 
  geom_bar(aes(x = pct, y = id, fill = `Number of AAChange`), 
           stat = "identity", position = "stack") + 
  theme_bw() + theme(panel.grid = element_blank())
dev.off()

# save(cell_inframe_result_processed, inframe_result_processed, file="~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda")




plot_data_tissue <- inframe_result_processed
plot_data_tissue$type <- 'tissue'

plot_data_cell <- cell_inframe_result_processed
plot_data_cell$type <- 'cell'
plot_data_tissue <- plot_data_tissue[plot_data_tissue$aa == 0,]
plot_data_cell <- plot_data_cell[plot_data_cell$aa == 0,]
plot_data <- data.frame(rbind(plot_data_tissue, plot_data_cell))
plot_data <- plot_data[!is.na(plot_data$pct),]
pdf("~/Nutstore Files/Tobin/Previous/0AA_Inframe_boxplot_tissue_cell.pdf", width = 6, height = 4)
ggplot(plot_data, aes(x= type, y = pct)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.35),
             shape= 21,
             color = "black", fill = "gray") + 
  stat_summary(fun = "median", colour = "black", size = 4,
               geom = "text", 
               aes(label = paste0("Median=",round(after_stat(y), 2))),
               vjust = c(-9,15), hjust = c(1.5, -0.5)) +
  xlab("") + ylab("0AA Inframe%") +
  theme_bw()
dev.off()
