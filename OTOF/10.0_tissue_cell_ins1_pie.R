library(Biostrings)
library(stringr)
source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
seq_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/FunctioalAnalysisInput.xlsx", 3)
setwd("~/data/project/ear_project/gene_therapy_ll/Previews/")
reps <- list.files(pattern = "^Rep.*data")
load("~/data/project/ear_project/gene_therapy_ll/Previews/total_used_seq.rda")
tmp <- total_used_seq[total_used_seq$mutNT_find != total_used_seq$mutNT,]

load(file="Result/old_rep123_edit.rda")

print(load("~/data/project/ear_project/gene_therapy_ll/Previews/Result/Rep123_indel_stat.rda"))
sample_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)

names(old_tissue) <- str_remove(names(old_tissue), "CRISPResso_on_")
names(old_cell) <- str_remove(names(old_cell), "CRISPResso_on_")
old_cell_order <- old_cell[names(old_cell) %in% sample_info$`Cell-w-W`]
split_name <- data.frame(cellid = sample_info$`Cell-w-W`, 
                         tissueid = sample_info$`Tissue-w-W(Mix)`,
                         mutation = sample_info$GRNAa)
split_name$mutation <- unlist(lapply(split_name$mutation, function(x){
  x <- unlist(strsplit(x, "[-]"))
  paste(x[-length(x)], collapse = "-")
}))

split_name_cell <- split_name[split_name$cellid %in% names(old_cell),]
split_name_cell <- na.omit(split_name_cell)

split_name_tissue <- split_name[split_name$tissueid %in% names(old_tissue),]
split_name_tissue <- na.omit(split_name_tissue)



split_name_cell$mutation <- unlist(lapply(split_name_cell$mutation, function(x){
  x <- unlist(str_split(x, "[-]"))
  x[2] <- stringr::str_to_title(x[2])
  paste(x, collapse = "-")
}))
split_name_cell$type <- unlist(lapply(split_name_cell$mutation, function(x){
  str_sub(x, str_length(x) - 1,  str_length(x) - 1)
}))


split_name_tissue$mutation <- unlist(lapply(split_name_tissue$mutation, function(x){
  x <- unlist(str_split(x, "[-]"))
  x[2] <- stringr::str_to_title(x[2])
  paste(x, collapse = "-")
}))
split_name_tissue$type <- unlist(lapply(split_name_tissue$mutation, function(x){
  str_sub(x, str_length(x) - 1,  str_length(x) - 1)
}))

split_name_tissue <- data.frame(rbind(split_name_tissue, 
                                      data.frame("cellid" = paste0("m", 12 : 15, "-tissue"),
                                                 "tissueid" = paste0("m", 12 : 15, "-tissue"),
                                                 "mutation" = 'w-Otof-1236dC', 
                                                 "type" = 'd')))
split_name_cell <- data.frame(rbind(split_name_cell, 
                                    data.frame("cellid" = paste0("m", 12 : 15),
                                               "tissueid" = paste0("m", 12 : 15),
                                               "mutation" = 'w-Otof-1236dC', 
                                               "type" = 'd')))

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

besides_nt_stat <- list()
nt_mats <- lapply(ins1_pct_table, function(tmp){
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
  strand <- total_used_seq_tissue[total_used_seq_tissue$result == id, "strand"]
  ###前面全部都是按照正链统计的，但是发现还有一些负链的，需要反过来算
  if(strand == "-"){
    leftnt <- as.character(Biostrings::reverseComplement(DNAString(leftnt)))
    rightnt <- as.character(Biostrings::reverseComplement(DNAString(rightnt)))
    tmp_nt <- leftnt
    leftnt <- rightnt
    rightnt <- tmp_nt
    rownames(tmp_pct_mat) <- as.character(Biostrings::reverseComplement(DNAStringSet(rownames(tmp_pct_mat))))
  }
  if(leftnt == "T" & rightnt == "A"){
    print(id)
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
                     file="~/Nutstore Files/Tobin/Previous/tissue_ins1_left_right_nt_pie_data.xlsx", 
                     rowNames=F, colNames=T)


fake_plot_mat <- matrix(runif(4*4), nrow = 4, ncol = 4)
rownames(fake_plot_mat) <- c("A", "T", "C", "G")
colnames(fake_plot_mat) <- c("A", "T", "C", "G")
pdf("~/Nutstore Files/Tobin/Previous/tissue_left_right_nt_pie_chart_groupPct.pdf", width = 6, height = 5)
col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(col) <-  c("A", 'T','C','G')
ht <- Heatmap(fake_plot_mat,
              name = "Name",
              column_split = NULL,
              column_title=NULL,
              row_title = NULL,
              row_split = NULL,
              show_column_names = T,
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 24, col = col),
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 24, col = col),
              column_names_side = "top", top_annotation = NULL,
              left_annotation = NULL,
              right_annotation = NULL,
              show_row_names = T,
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {

                leftnt <- rownames(fake_plot_mat)[i]
                rightnt <- colnames(fake_plot_mat)[j]
                pct <- besides_nt_stat[[leftnt]][[rightnt]]
                if(!is.null(pct)){
                  labels <- names(pct)
                  pct <- unlist(pct)
                  col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                  names(col) <- c("A", 'T','C','G')
                  pct <- c(0, cumsum(pct)/sum(pct))
                  dx <- diff(pct)
                  nx <- length(dx)
                  twopi <- 2 * pi
                  init.angle <- 0
                  t2xy <- function(t) {
                    t2p <- twopi * t + init.angle * pi/180
                    list(x = 0.8 * cos(t2p), y = 0.8 * sin(t2p))
                  }
                  maxx <- max(unlist(lapply(1 : nx, function(i){
                    n <- max(10, floor(200 * dx[i]))
                    P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                    max(abs(P$x))
                  })))
                  maxy <- max(unlist(lapply(1 : nx, function(i){
                    n <- max(10, floor(200 * dx[i]))
                    P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                    max(abs(P$y))
                  })))
                  facy <- as.numeric(h) / maxy * 0.3
                  facx <- as.numeric(w) / maxx * 0.3
                  for (i in 1L:nx) {
                    n <- max(10, floor(200 * dx[i]))
                    P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))

                    grid.polygon(unit(as.numeric(x) + c(P$x * facx, 0), units = "npc"),
                                 unit(as.numeric(y) + c(P$y * facy, 0), units = "npc"),
                                 gp = gpar(fill = col[labels[i]]))
                    # print(P)
                    # polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = 45,
                    #         border = border[i], col = col[i], lty = lty[i])
                  }
                }


              }, col = c("white", "white"), show_heatmap_legend = F,

              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16,
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)
draw(ht, annotation_legend_list = lgd_list)
dev.off()



###计算细胞的结果
noAA_change_seq <- list()

for(i in 1 : nrow(total_used_seq_tissue)){
  id <- total_used_seq_tissue$result[i]
  if(id %in% split_name_cell$tissueid){
    id <- split_name_cell$cellid[split_name_cell$tissueid == id]
  } else {
    id <- str_remove(id, "-tissue")
  }
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

besides_nt_stat <- list()
nt_mats <- lapply(ins1_pct_table, function(tmp){
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
  
  if(leftnt == "T" & rightnt == "A"){
    print(id)
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
                     file="~/Nutstore Files/Tobin/Previous/cell_ins1_left_right_nt_pie_data.xlsx", 
                     rowNames=F, colNames=T)


fake_plot_mat <- matrix(runif(4*4), nrow = 4, ncol = 4)
rownames(fake_plot_mat) <- c("A", "T", "C", "G")
colnames(fake_plot_mat) <- c("A", "T", "C", "G")
pdf("~/Nutstore Files/Tobin/Previous/cell_left_right_nt_pie_chart_groupPct.pdf", width = 6, height = 5)
col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(col) <-  c("A", 'T','C','G')
ht <- Heatmap(fake_plot_mat,
              name = "Name",
              column_split = NULL,
              column_title=NULL,
              row_title = NULL,
              row_split = NULL,
              show_column_names = T,
              cluster_rows = F,
              cluster_columns = F,
              row_names_gp = gpar(fontsize = 24, col = col),
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 24, col = col),
              column_names_side = "top", top_annotation = NULL,
              left_annotation = NULL,
              right_annotation = NULL,
              show_row_names = T,
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                leftnt <- rownames(fake_plot_mat)[i]
                rightnt <- colnames(fake_plot_mat)[j]
                pct <- besides_nt_stat[[leftnt]][[rightnt]]
                if(!is.null(pct)){
                  labels <- names(pct)
                  pct <- unlist(pct)
                  col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                  names(col) <- c("A", 'T','C','G')
                  pct <- c(0, cumsum(pct)/sum(pct))
                  dx <- diff(pct)
                  nx <- length(dx)
                  twopi <- 2 * pi
                  init.angle <- 0
                  t2xy <- function(t) {
                    t2p <- twopi * t + init.angle * pi/180
                    list(x = 0.8 * cos(t2p), y = 0.8 * sin(t2p))
                  }
                  maxx <- max(unlist(lapply(1 : nx, function(i){
                    n <- max(10, floor(200 * dx[i]))
                    P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                    max(abs(P$x))
                  })))
                  maxy <- max(unlist(lapply(1 : nx, function(i){
                    n <- max(10, floor(200 * dx[i]))
                    P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                    max(abs(P$y))
                  })))
                  facy <- as.numeric(h) / maxy * 0.3
                  facx <- as.numeric(w) / maxx * 0.3
                  for (i in 1L:nx) {
                    n <- max(10, floor(200 * dx[i]))
                    P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                    
                    grid.polygon(unit(as.numeric(x) + c(P$x * facx, 0), units = "npc"),
                                 unit(as.numeric(y) + c(P$y * facy, 0), units = "npc"),
                                 gp = gpar(fill = col[labels[i]]))
                    # print(P)
                    # polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = 45,
                    #         border = border[i], col = col[i], lty = lty[i])
                  }
                }
                
                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16,
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)
draw(ht, annotation_legend_list = lgd_list)
dev.off()





tissue_data <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/tissue_ins1_left_right_nt_pie_data.xlsx")
cell_data <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/cell_ins1_left_right_nt_pie_data.xlsx")

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
                           init.angle = 0)
                }
                
                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16,
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)
pdf("~/Nutstore Files/Tobin/Previous/tissue_cell_left_right_nt_pie_chart_groupPct_fixerror.pdf", 
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

