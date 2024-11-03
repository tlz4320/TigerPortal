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
nt_counts_mats <- lapply(ins1_pct_table, function(tmp){
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

nt_mats <- lapply(ins1_pct_table, function(tmp){
  plot_mat <- matrix(0, ncol = 17, nrow = 4)
  rownames(plot_mat) <- c("A", "T", "C", "G")
  tmp <- tmp[tmp$pos %in% c(15 : 23),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
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
for(id in names(nt_counts_mats)){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- wt_seqs[[id]]

    #还是过滤1一下的结果
  i <- 12
  leftnt <- tmp_seq_mat[20]
  rightnt <- tmp_seq_mat[21]
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
  
  for(nt in rownames(tmp_counts_mat)){
    tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_counts_mat[nt, i]
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
# openxlsx::write.xlsx(output_table, 
#                      file="~/Nutstore Files/Tobin/Previous/tissue_ins1_left_right_nt_pie_data.xlsx", 
#                      rowNames=F, colNames=T)


fake_plot_mat <- matrix(runif(4*4), nrow = 4, ncol = 4)
rownames(fake_plot_mat) <- c("A", "T", "C", "G")
colnames(fake_plot_mat) <- c("A", "T", "C", "G")
# pdf("~/Nutstore Files/Tobin/Previous/tissue_left_right_nt_pie_chart_groupPct.pdf", width = 6, height = 5)
# col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
# names(col) <-  c("A", 'T','C','G')
# ht <- Heatmap(fake_plot_mat, 
#               name = "Name",
#               column_split = NULL,
#               column_title=NULL,
#               row_title = NULL,
#               row_split = NULL,
#               show_column_names = T, 
#               cluster_rows = F, 
#               cluster_columns = F, 
#               row_names_gp = gpar(fontsize = 24, col = col),
#               column_names_rot = 0,
#               column_names_centered = T,
#               column_names_gp = gpar(fontsize = 24, col = col),
#               column_names_side = "top", top_annotation = NULL,
#               left_annotation = NULL,
#               right_annotation = NULL,
#               show_row_names = T,
#               row_names_side = "left",
#               cell_fun = function(j, i, x, y, w, h, fill) {
#                 
#                 leftnt <- rownames(fake_plot_mat)[i]
#                 rightnt <- colnames(fake_plot_mat)[j]
#                 pct <- besides_nt_stat[[leftnt]][[rightnt]]
#                 if(!is.null(pct)){
#                   labels <- names(pct)
#                   pct <- unlist(pct)
#                   col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
#                   names(col) <- c("A", 'T','C','G')
#                   pct <- c(0, cumsum(pct)/sum(pct))
#                   dx <- diff(pct)
#                   nx <- length(dx)
#                   twopi <- 2 * pi
#                   init.angle <- 0
#                   t2xy <- function(t) {
#                     t2p <- twopi * t + init.angle * pi/180
#                     list(x = 0.8 * cos(t2p), y = 0.8 * sin(t2p))
#                   }
#                   maxx <- max(unlist(lapply(1 : nx, function(i){
#                     n <- max(10, floor(200 * dx[i]))
#                     P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
#                     max(abs(P$x))
#                   })))
#                   maxy <- max(unlist(lapply(1 : nx, function(i){
#                     n <- max(10, floor(200 * dx[i]))
#                     P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
#                     max(abs(P$y))
#                   })))
#                   facy <- as.numeric(h) / maxy * 0.3
#                   facx <- as.numeric(w) / maxx * 0.3
#                   for (i in 1L:nx) {
#                     n <- max(10, floor(200 * dx[i]))
#                     P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
#                     
#                     grid.polygon(unit(as.numeric(x) + c(P$x * facx, 0), units = "npc"), 
#                                  unit(as.numeric(y) + c(P$y * facy, 0), units = "npc"), 
#                                  gp = gpar(fill = col[labels[i]]))
#                     # print(P)
#                     # polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = 45,
#                     #         border = border[i], col = col[i], lty = lty[i])
#                   }
#                 }
# 
#                 
#               }, col = c("white", "white"), show_heatmap_legend = F,
#               
#               border = F)
# lgd_list = list(
#   Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16, 
#          legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
# )
# draw(ht, annotation_legend_list = lgd_list)
# dev.off()



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
    
    # if(cutsite + 20 > length(wtseq)){
    #   ref <- str_sub(ref, 1, -(cutsite + 21 - length(wtseq)))
    # }
    # if(cutsite < 20){
    #   ref <- str_sub(ref, 20 - cutsite)
    # }
    
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
    
    
    # ref <- str_remove_all(ref, "-")
    # findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
    #                                     subject = total_used_seq_tissue$seq[total_used_seq_tissue$result == x], max.mismatch = 4, with.indels = T)
    # insertpos <- findRef@ranges@start + editpos - 1
    # findRef <- Biostrings::matchPattern(pattern= str_to_lower(ref),
    #                                     subject = mutseq, max.mismatch = 4, with.indels = T)
    # #找不到的情况
    # if(length(findRef@ranges@start) == 0){
    #   return(-1)
    # }
    # mutseqs <- unlist(strsplit(mutseq, "*"))
    # if(insertpos <= mutpos ){
    #   matchpos <- findRef@ranges@start
    #   mutseqs <- c(mutseqs[1 : (matchpos + editpos - 2)], str_to_upper(editNT), mutseqs[(matchpos + editpos - 1) : length(mutseqs)])
    #   if(length(mutseqs) != length(wtseq)){
    #     print(x)
    #     print(i)
    #   }
    #   return(sum(mutseqs != wtseq))
    # }else{
    #   matchpos <- length(mutseqs) - (findRef@ranges@start + findRef@ranges@width - 1) + 1
    #   editpos <- str_length(ref) - editpos + 1
    #   mutseqs <- rev(mutseqs)
    #   mutseqs <- c(mutseqs[1 : (matchpos + editpos - 1)], editNT, mutseqs[(matchpos + editpos) : length(mutseqs)])
    #   mutseqs <- rev(mutseqs)
    #   if(length(mutseqs) != length(wtseq)){
    #     print(x)
    #     print(i)
    #   }
    #   return(sum(mutseqs != wtseq))
    # }
  }))
  seq_table
})
seqs_table_stat <- data.frame(do.call(rbind, seqs_table_stat))



table(seqs_table_stat$diff == -1,seqs_table_stat$id )

sp_info <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/USH_Loci_DataSum.xlsx", 2)
sp_info <- split_name <- data.frame(
                                    tissueid = c(sp_info$`Tissue-w-W(Mix)`,
                                                 paste0("m", 12:15, "-tissue")),
                                    mutation = c(sp_info$GRNAa, paste0("OTOF-1236dC-g", 3:6)))

#选出delete在1,2,3号碱基或者在4号碱基且会改变左侧碱基类型的认为是不可预测的
unpredictable_sample <- total_used_seq_tissue[(total_used_seq_tissue$mutPos > total_used_seq_tissue$cutSite & 
                                                total_used_seq_tissue$mutPos <= 3 + total_used_seq_tissue$cutSite) | 
                                                total_used_seq_tissue$mutPos == total_used_seq_tissue$cutSite & 
                                                unlist(lapply(total_used_seq_tissue$sg, function(x){
                                                  x <- unlist(strsplit(x, "*"))
                                                  x[length(x) - 3] != x[length(x) - 4]
                                                })),]
ins1_lr_nt_stat <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Merged1NT/merged_ins1_left_right_nt_pie_data.xlsx")

unpredictable_sample$diff <- unlist(lapply(1: nrow(unpredictable_sample), function(i){
  wt_seq <- unpredictable_sample$seq[i]
  wt_seq <- str_to_upper(unlist(strsplit(wt_seq, "*")))
  
  mutpos <- unpredictable_sample$mutPos[i]
  cutpos <- unpredictable_sample$cutSite[i]
  sg <- str_to_upper(unlist(strsplit(total_used_seq_tissue$sg, "*")))
  del_pos <- mutpos - cutpos
  sg <- sg[-(length(sg) - 3 + del_pos)]
  leftnt <- sg[length(sg) - 3]
  rightnt <- sg[length(sg) - 2]
  
  pred <- ins1_lr_nt_stat[ins1_lr_nt_stat$leftNT == leftnt&ins1_lr_nt_stat$rightNT == rightnt,c(-1,-2)]
  pred <- unlist(pred)
  diff <- 0
  for(j in 1 : 4){
    pred_seq <- wt_seq[-mutpos]
    pred_seq <- c(pred_seq[1 : (cutpos - 1)], names(pred)[j], pred_seq[cutpos : length(pred_seq)])

    diff <- diff + sum(pred_seq != wt_seq) * pred[j] / 100
  }
  diff
}))

plot_data <- lapply(split(seqs_table_stat, seqs_table_stat$id), function(x){
  x <- x[x$insPos %in% c(20,21),]
  top1_ins1 <- lapply(split(x$Pct, x$insPos), sum)
  top1_ins1 <- max(unlist(top1_ins1))
  data.frame(diff = sum((x$Pct / sum(x$Pct)) * x$diff), dis = x$dis[1], id = x$id[1], 
             top1Ins1Pct = top1_ins1)
})

plot_data <- data.frame(do.call(rbind, plot_data))
plot_data <- merge(plot_data, sp_info, by.x="id", by.y="tissueid")
rownames(plot_data) <- plot_data$id
unpredictable_sample_sel <- unpredictable_sample[unpredictable_sample$result %in% plot_data$id,]
plot_data[unpredictable_sample_sel$result, ]$diff <- unpredictable_sample_sel$diff

plot_data <- plot_data[order(plot_data$dis),]
plot_data$color <- "Other"
plot_data$color[plot_data$id %in% c("105-051")] <- "OTOF-1236dC-g1"
plot_data$label <- NA
plot_data <- plot_data[order(plot_data$diff),]
plot_data$label[1:5] <- plot_data$mutation[1:5]


library(ggtext)

plot_data$dis[plot_data$dis > 0] <- plot_data$dis[plot_data$dis > 0]  - 1
library(ggrepel)
plot_data$ng <- 2
plot_data$ng[c(2:4)] <- -1

# tmp <- plot_data[plot_data$diff < 3 & abs(plot_data$dis) <= 3,]
# openxlsx::write.xlsx(tmp, file="~/data/project/ear_project/gene_therapy_ll/Previews/volcano_+-3_sel_sg.xlsx")
pdf("~/Nutstore Files/Tobin/Previous/volcano_plot_of_distance_ntdiff_del1Mut_v6_b_r_fix.pdf", width = 9, height = 6)
ggplot(plot_data) + geom_point(aes(x = dis, y =diff, color = color, 
                                   fill = top1Ins1Pct), 
                               position = position_jitter(seed = 1), 
                               shape = 21, size = 3) + 
  geom_text_repel(aes(x = dis, y = diff, label = label, color = color), 
                  position = position_jitter(seed = 1), 
                  min.segment.length = unit(0, 'lines'), 
                  hjust = plot_data$ng) + 
  xlab("Distance cleavage-del (NT)")+
  ylab("Number Change (NT)") + 
  scale_fill_gradientn(colours = c("blue","white", "red")) + 
  scale_x_continuous(breaks = (-13 : 13), labels = c(-13 : -1, 1 : 14)) + 
  scale_y_continuous(breaks = seq(0, 15, 2)) + 
  scale_color_manual(values = c("black", "red")) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()



pdf("~/Nutstore Files/Tobin/Previous/volcano_plot_of_distance_ntdiff_del1Mut_v6_r.pdf", width = 9, height = 6)
ggplot(plot_data) + geom_point(aes(x = dis, y =diff, color = color, 
                                   fill = top1Ins1Pct), 
                               position = position_jitter(seed = 1), 
                               shape = 21, size = 3) + 
  geom_text_repel(aes(x = dis, y = diff, label = label, color = color), 
                  position = position_jitter(seed = 1), 
                  min.segment.length = unit(0, 'lines'), 
                  hjust = plot_data$ng) + 
  xlab("Distance cleavage-del (NT)")+
  ylab("Number Change (NT)") + 
  scale_fill_gradientn(colours = c("white", "red")) + 
  scale_x_continuous(breaks = (-13 : 13), labels = c(-13 : -1, 1 : 14)) + 
  scale_y_continuous(breaks = seq(0, 15, 2)) + 
  scale_color_manual(values = c("black", "red")) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
table(plot_data$dis < 6 & plot_data$dis > -5)

print(load("~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda"))
tissue_0aa <- inframe_result_processed[inframe_result_processed$aa == "0",]

tissue_inframe <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/previous_tissue_inframe_in_indel.xlsx")
tissue_inframe <- tissue_inframe[,c("result", "inframePct")]
plot_data2 <- merge(plot_data, tissue_inframe,by.x="id", by.y="result")
plot_data2 <- merge(plot_data2, tissue_0aa, by="id")


pdf("~/Nutstore Files/Tobin/Previous/volcano_plot_of_distance_ntdiff_del1Mut_v7_fix.pdf", width = 9, height = 6)
ggplot(plot_data2) + geom_point(aes(x = dis, y =diff, color = inframePct, 
                                   fill = pct2), 
                               position = position_jitter(seed = 1), 
                               shape = 21, size = 3, stroke = 1) + 
  scale_color_gradientn(colours = c("blue",  "red"),limits = c(0,100),
                        name = "Inframe%") + theme_bw() + 
  ggnewscale::new_scale_colour() + 
  geom_text_repel(aes(x = dis, y = diff, label = label, color = color), 
                  position = position_jitter(seed = 1), 
                  min.segment.length = unit(0, 'lines'), 
                  hjust = plot_data$ng) + 
  xlab("Distance cleavage-del (NT)")+
  ylab("Number Change (NT)") + 
  scale_fill_gradientn(colours = c("#FFEBCD", "#F09B28"), name = "0AA/Inframe") + 
  scale_x_continuous(breaks = (-13 : 13), labels = c(-13 : -1, 1 : 14)) + 
  scale_y_continuous(breaks = seq(0, 15, 2)) + 
  scale_color_manual(values = c("black", "red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

plot_data_sub <- plot_data[plot_data$dis %in% as.character(-6 : 6),]
#为了把0去掉
plot_data_sub$dis[plot_data_sub$dis > 0] <- plot_data_sub$dis[plot_data_sub$dis > 0]  - 1
pdf("~/data/project/ear_project/gene_therapy_ll/Previews/Result/volcano_plot_of_distance_ntdiff_(-6)-6.pdf", width = 6, height = 4)
ggplot(plot_data_sub) + geom_point(aes(x = dis, y =diff, color = color), 
                               position = position_jitter(seed = 1), shape = 1, size = 3) + 
  geom_text_repel(aes(x = dis, y = diff, label = label, color = color), 
                  position = position_jitter(seed = 1), 
                  min.segment.length = unit(0, 'lines'), max.iter = 1000000, max.time = 100000, max.overlaps = 10) + 
  scale_x_continuous(breaks = (-6 : 5), labels = c(-6 : -1, 1 : 6)) + 
  scale_color_manual(values = c("black", "red")) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()




edit_stat_of_distance <- lapply(names(seqs_table_split), function(x){
  seq_table <- seqs_table_split[[x]]
  seq_table <- seq_table[order(seq_table$Pct,decreasing = T),]
  mutpos <- total_used_seq_tissue$mutPos[total_used_seq_tissue$result == x]
  cutsite <- total_used_seq_tissue$cutSite[total_used_seq_tissue$result == x]
  dis <- c(total_used_seq_tissue$mutPos - total_used_seq_tissue$cutSite)[total_used_seq_tissue$result == x]
  
  seq_table$dis <- ifelse(dis <= 0, dis - 1, dis)
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
      return(-1)
    }
    res <- mutpos - (findRef@ranges@start[1] + editpos - 2)
    if(res <= 0){
      res <- res - 1
    }
    res
  }))
  seq_table
})
edit_stat_of_distance <- data.frame(do.call(rbind, edit_stat_of_distance))
table(edit_stat_of_distance$dis2mut)



edit_stat_of_distance <- merge(edit_stat_of_distance, sp_info, by.x="id", by.y = 'tissueid')
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

plot_mat <- matrix(0, nrow = 40, ncol=length(unique(edit_stat_of_distance$id)))
rownames(plot_mat) <- as.character(c(-20 : -1, 1 : 20))
colnames(plot_mat) <- unique(edit_stat_of_distance$id)
for(i in 1 : nrow(edit_stat_of_distance)){
  dis <- as.character(edit_stat_of_distance$dis2mut[i])
  if(!dis %in% rownames(plot_mat))
    next
  id <- edit_stat_of_distance$id[i]
  plot_mat[dis, id] <- plot_mat[dis, id] + edit_stat_of_distance$Pct[i]
}




# tmp <- edit_stat_of_distance[abs(edit_stat_of_distance$dis2mut) > 15,]
# sum(edit_stat_of_distance$Pct)
# sum(plot_mat)
# unlist(lapply(split(edit_stat_of_distance$Pct, edit_stat_of_distance$mutation), sum))
# View(plot_mat)
sp_info <- edit_stat_of_distance[,c(1,6,7,8,9,10, 11)]
sp_info <- sp_info[!duplicated(sp_info$id),]
sp_info$id == colnames(plot_mat)


dis_mat <- matrix(0, nrow = 40, ncol=length(unique(edit_stat_of_distance$id)))
rownames(dis_mat) <- as.character(c(-20 : -1, 1 : 20))
colnames(dis_mat) <- unique(edit_stat_of_distance$id)
for(i in 1 : nrow(sp_info)){
  dis <- as.character(sp_info$dis[i])
  id <- sp_info$id[i]
  dis_mat[dis, id] <- 1
}

color_func = colorRamp2(c(0, 100), c("white", "#1432EB"))
heatanno <- 
  HeatmapAnnotation(gene = anno_empty(border = F, height = unit(20, "mm")),
                    loc = anno_empty(border = F, height = unit(20, "mm")),
                    "line" = 
                      anno_empty(border = FALSE, height = unit(2, "mm")),
                    sg = anno_empty(border = F, height = unit(20, "mm"))
  )
pdf("Previews/Result/Ins1_edit_distance_heatmap.pdf", width = 20, height = 12)
Heatmap(plot_mat, 
              name = "EditDis",
              column_split = sp_info$gene,
              column_title=NULL,
              row_split = c(rep(" ", 20), rep("  ", 20)),
              show_column_names = F, 
              cluster_rows = F, cluster_columns = F, 
              column_title_rot = 45,
              column_names_side = "top",
              top_annotation = heatanno,
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                # transform the matrix into a long vector
                v = pindex(plot_mat, i, j) 
                dis = pindex(dis_mat, i, j)
                # `j` here is also a vector with the same length of `v`
                col = color_func(v)
                col2 = "white"
                grid.rect(x, y, w, h, gp = gpar(fill = col, col = col2))
                if(dis== 1){
                  grid.points(x, y, gp = gpar(col = "red"))
                }
                x_from <- x - unit(as.numeric(w) / 2, "npc")
                x_to <- x + unit(as.numeric(w) / 2, "npc")
                y_from <- y - unit(as.numeric(h) / 2, "npc")
                y_to <- y + unit(as.numeric(h) / 2, "npc")
                grid.polyline(x = c(x_from, x_to, x_from, x_to), 
                              y = c(y_from,y_from,y_to , y_to), 
                              gp = gpar(col = "#D0D0D0"), 
                              id = rep(1 : (length(x) * 2), 2))
                grid.text(sprintf("%.2f", v), x, y, 
                          gp = gpar(fontsize = 3))
                
                
              }, col=c("white", "blue"),show_heatmap_legend = T,
              border = F)
# 
# draw(ht, annotation_legend_list = pd)

# for(i in 1 : length(unique(split_col))){
#   for(j in 1 : ifelse(is.null(restore), 2, 3)){
#     if(j != 1){
#       decorate_heatmap_body(heatmap = title, row_slice = j, column_slice = i,
#                             {
#                               grid.lines(x = unit(c(0, 1), "npc"),
#                                          y=unit(c(1, 1), "npc"), 
#                                          gp = gpar(lwd = 2))
#                             })
#     }
#     decorate_heatmap_body(heatmap = title, row_slice = j, column_slice = i,
#                           {
#                             grid.lines(x = unit(c(0, 0), "npc"),
#                                        y=unit(c(0, 1), "npc"), 
#                                        gp = gpar(lwd = 1))
#                           })
#     decorate_heatmap_body(heatmap = title, row_slice = j, column_slice = i,
#                           {
#                             grid.lines(x = unit(c(1, 1), "npc"),
#                                        y=unit(c(0, 1), "npc"), 
#                                        gp = gpar(lwd = 1))
#                           })
#     # }
#   }
#   
# }
  lens <- unique(sp_info$gene)
  
  
  
  for(i in 1:length(unique(lens))) {
    len <- sum(sp_info$gene == lens[i])
    decorate_annotation("gene", slice = i, {
      tg <- richtext_grob(gt_render(lens[i]), 
                          rot = 90, 
                          x = unit(0.5, "npc"),
                          y=unit(0, "npc"), hjust = 0)
      grid.draw(tg)
      invisible(tg)
    })
    locs_count <- length(unique(sp_info$loca[sp_info$gene == lens[i]]))
    locs <- unique(sp_info$loca[sp_info$gene == lens[i]])
    sg_count <- sum(sp_info$gene == lens[i])
    pos <- 0
    for(l in locs){
      
      ll <- sum(sp_info$loca[sp_info$gene == lens[i]] == l)
      pos <- pos + 1 / sg_count * ll / 2
      decorate_annotation("loc", slice = i, {
        tg <- richtext_grob(gt_render(l), 
                            rot = 90, 
                            x = unit(pos, "npc"),
                            y=unit(0, "npc"), hjust = 0)
        grid.draw(tg)
        invisible(tg)
      })
      if(ll == 1){
        decorate_annotation("line", slice = i, {
          grid.lines(x = unit(c(pos, pos), "npc"),
                     y=unit(c(0, 1), "npc"))
        })
      }
      else{
        decorate_annotation("line", slice = i, {
          grid.lines(x = unit(c(pos - 1 / sg_count * ll / 2 + 1 / sg_count / 2, pos - 1 / sg_count * ll / 2 + 1 / sg_count / 2), "npc"),
                     y=unit(c(0, 1), "npc"))
          grid.lines(x = unit(c(pos - 1 / sg_count * ll / 2 + 1 / sg_count / 2, 
                                pos + 1 / sg_count * ll / 2 - 1 / sg_count / 2), "npc"),
                     y=unit(c(1, 1), "npc"))
          grid.lines(x = unit(c(pos + 1 / sg_count * ll / 2 - 1 / sg_count / 2,
                                pos + 1 / sg_count * ll / 2 - 1 / sg_count / 2), "npc"),
                     y=unit(c(0, 1), "npc"))
        })
      }
      
      
      
      pos <- pos + 1 / sg_count * ll / 2
      
      
      
    }
    sgs <- sp_info$sg[sp_info$gene == lens[i]]
    pos <- 1 / sg_count / 2
    for(s in sgs){
      decorate_annotation("sg", slice = i, {
        tg <- richtext_grob(gt_render(s), 
                            rot = 90, 
                            x = unit(pos, "npc"),
                            y=unit(0, "npc"), hjust = 0)
        grid.draw(tg)
        invisible(tg)
      })
      pos <- pos + 1 / sg_count
    }
  }
dev.off()
