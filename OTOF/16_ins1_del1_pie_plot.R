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
    ##选出只有Del1的结果
    edit_table <- edit_table[edit_table$n_deleted == 1 & 
                               edit_table$n_inserted == 0& 
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


wt_seqs <- lapply(unique(seqs_table$id), function(id){
  edit_table <- total_old_edit_table_rev[[id]][[1]]
  edit_table <- edit_table[edit_table$Unedited == "True",]
  getWtSeq(edit_table)
})
names(wt_seqs) <- unique(seqs_table$id)
del1_pct_table <- lapply(unique(seqs_table$id), function(id){
  del1_table <- seqs_table[seqs_table$id == id,]
  wt_seq <- wt_seqs[[id]]
  tmp <- getStat2(del1_table, F)
  tmp$id <- id
  tmp
})
names(del1_pct_table) <- unique(seqs_table$id)


noAA_change_seq_cell <- list()

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
    ##选出只有Del1的结果
    edit_table <- edit_table[edit_table$n_deleted == 1 & 
                               edit_table$n_inserted == 0& 
                               edit_table$n_mutated <= 1,]
    edit_table
  })
  noAA_change_seq_cell[[id]] <- edit_table
  
}

seqs_table <- list()
for(id in names(noAA_change_seq_cell)){
  edit_table <- noAA_change_seq_cell[[id]]
  
  edit_table <- data.frame(do.call(rbind, edit_table))
  if(nrow(edit_table) == 0){
    print(id)
    next
  }
    
  if(sum(edit_table[,8]) < 0.1){
    print(id)
    next
  }
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


seqs_table_split <- split(seqs_table, seqs_table$id)


del1_pct_table_cell <- lapply(unique(seqs_table$id), function(id){
  ins1_table <- seqs_table[seqs_table$id == id,]
  wt_seq <- wt_seqs[[id]]
  tmp <- getStat2(ins1_table, F)
  tmp$id <- id
  tmp
})
names(del1_pct_table_cell) <- unique(seqs_table$id)



total_plot_data <- list()
common_id <- names(del1_pct_table)[str_remove(names(del1_pct_table), "[-].*") %in% 
                                     str_remove(names(del1_pct_table_cell), "[-].*")]
sameLR <- list()
for(name in common_id){
  leftnt <- wt_seqs[[name]][20]
  rightnt <- wt_seqs[[name]][21]
  if(leftnt == rightnt){
    sameLR[[length(sameLR) + 1]] <- name
  }
  strand <- total_used_seq_tissue$strand[total_used_seq_tissue$result == name]
  tmp <- del1_pct_table[[name]]
  
  id <- name
  if(id %in% split_name_cell$tissueid){
    id <- split_name_cell$cellid[split_name_cell$tissueid == id]
  } else {
    id <- str_remove(id, "-tissue")
  }
  tmp2 <- del1_pct_table_cell[[id]]
  tmp <- tmp[tmp$pos %in% c(20, 21),]
  tmp <- tmp[tmp$NT %in% c("A", "T", "C", "G"),]
  tmp$Pct <- tmp$Pct / sum(tmp$Pct) * 100
  tmp2 <- tmp2[tmp2$pos %in% c(20, 21),]
  tmp2 <- tmp2[tmp2$NT %in% c("A", "T", "C", "G"),]
  tmp2$Pct <- tmp2$Pct / sum(tmp2$Pct) * 100
  if(strand == "-"){
    leftnt <- as.character(Biostrings::reverseComplement(DNAString(leftnt)))
    rightnt <- as.character(Biostrings::reverseComplement(DNAString(rightnt)))
    tmp_nt <- leftnt
    leftnt <- rightnt
    rightnt <- tmp_nt
    tmp$NT <- as.character(Biostrings::reverseComplement(DNAStringSet(tmp$NT)))
    tmp2$NT <- as.character(Biostrings::reverseComplement(DNAStringSet(tmp2$NT)))
  }
  
  
  diff <- sum(tmp$Pct[tmp$NT == leftnt]) / sum(tmp$Pct)
  
  
  
  diff_cell <- sum(tmp2$Pct[tmp2$NT == leftnt]) / sum(tmp2$Pct)
  
  total_plot_data[[id]] <- data.frame(diff = c(diff, diff_cell), 
                                      type=c("Tissue", "Cell"),
                                      id = name, 
                                      NT = leftnt)
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data$diff[is.na(total_plot_data$diff)] <- 0




###后面是读取预测的结果
setwd('~/indelphi_res//')
indelphi_118_res <- list()
for(file in list.files()){
  indelphi_118_res[[file]] <- read.table(file,sep = ",", header = T)
}
names(indelphi_118_res) <- str_remove(names(indelphi_118_res), ".csv")
indelphi_res_del1 <- lapply(indelphi_118_res, function(x){
  x <- x[x$Category == "del" & x$Length == 1,]
  x$Pct <- x$Predicted.frequency / sum(x$Predicted.frequency) * 100
  colnames(x)[4] <- "NT"
  x
})
# save(indelphi_118_res, file="~/data/project/ear_project/gene_therapy_ll/Result/118_indelphi.rda")

setwd('~/forecast_result//')
forcast_file <- read.table("~/118_for_forecast.txt")
forecast_118_res_del1 <- list()
total_used_seq_tissue$id <- str_remove(total_used_seq_tissue$result, "[-].*")
for(file in total_used_seq_tissue$id){
  tmp1 <- read.table(paste0(file, "_predictedindelsummary.txt"), 
                     sep = "\t", header = T)
  colnames(tmp1) <- c("V1", "NT", "Count")
  tmp2 <- read.table(paste0(file, "_predictedreads.txt"), 
                     sep = "\t", header = F)
  sg <- str_to_upper(total_used_seq_tissue$sg[total_used_seq_tissue$id == file][1])
  if(is.null(sg) | any(is.na(str_locate(tmp2[1,2], sg)))){
    print(file)
  }
  cutpos <- unlist(str_locate(tmp2[1,2], sg)[1]) + str_length(sg) - 3
  tmp2 <- merge(tmp2, tmp1, by.y='V1', by.x = "V3")
  type <- unlist(lapply(tmp2[,1], function(x){
    unlist(strsplit(x, "[_]"))[1]
  }))
  tmp2 <- tmp2[type == "D1",]
  sg <- str_to_upper(total_used_seq_tissue$sg[total_used_seq_tissue$result == paste0(file, "-051")][1])
  sg <- unlist(strsplit(sg, "*"))
  leftnt <- sg[length(sg) - 3]
  rightnt <- sg[length(sg) - 2]
  if(leftnt == rightnt){
    tmp2$NT <- leftnt
  } else {
    pos <- as.numeric(str_sub(tmp2[,1], -1, -1))
    tmp2$NT[which.min(pos)] <- leftnt
    tmp2$NT[which.max(pos)] <- rightnt
  }

  
  forecast_118_res_del1[[file]] <- tmp2
}

predict_data <- list()
for(name in common_id){
  sg <- str_to_upper(total_used_seq_tissue$sg[total_used_seq_tissue$result == name][1])
  if(is.null(sg)){
    print(name)
  }
  leftnt <- unlist(strsplit(sg, "*"))
  leftnt <- leftnt[length(leftnt) - 3]
  name <- str_remove(name, "[-].*")
  
  tmp <- indelphi_res_del1[[name]]
  
  indelphi_diff <- sum(tmp$Pct[tmp$Genotype.position == 0]) / sum(tmp$Pct)
  
  tmp <- forecast_118_res_del1[[name]]
  forecast_diff <- sum(tmp$Count[tmp$NT == leftnt]) / sum(tmp$Count)
  predict_data[[name]] <- data.frame(diff = c(forecast_diff,indelphi_diff), 
                                     type=c("FORECasT", "inDelphi"), 
                                     id = name, 
                                     NT = leftnt)
}
predict_data <- data.frame(do.call(rbind, predict_data))
predict_data$diff[is.na(predict_data$diff)] <- 0

total_plot_data2 <- data.frame(rbind(total_plot_data, predict_data))
total_plot_data2$type <- factor(total_plot_data2$type, levels = c( "Tissue","Cell", "inDelphi", "FORECasT"))


tmp <- total_plot_data2
tmp$type <- as.numeric(tmp$type)
tmp <- beeswarm::beeswarm(diff ~ type, tmp, corral = "none", corralWidth = 0.2, spacing = 0.4)
total_plot_data2$beex <- 0
for(i in 1 : 4){
  total_plot_data2[as.numeric(total_plot_data2$type) == i, "beex"] <- tmp$x[tmp$x.orig == i]
}
total_plot_data2$id <- str_remove(total_plot_data2$id, "[-].*")
total_plot_data2 <- total_plot_data2[order(total_plot_data2$id),]
total_plot_data2 <- total_plot_data2[!total_plot_data2$id %in% str_remove(unlist(sameLR), "[-].*"),]

pdf("~/Nutstore Files/Tobin/Previous/ins1Mut_link_plot_add_predict.pdf", width = 6, height = 4)
library(ggpubr)
ggplot(total_plot_data2, aes(x = type, y = diff)) + 
  geom_violin(aes(fill = type), 
              alpha = 0.5, 
              scale = "width", 
              adjust = 1, 
              trim = TRUE) + 
  ggnewscale::new_scale("fill") +
  geom_point(aes(x = beex, fill=NT),
             size = 2, 
             shape = 24,
             color = "black") + 
  geom_line(aes(x = beex, y = diff, color = NT, 
                group = id), 
            linewidth=0.1) + 
  stat_compare_means(comparisons = list(c("Tissue", "Cell"), 
                                        c("Cell", "inDelphi"), 
                                        c("Cell", "FORECasT"),
                                        c("Tissue", "inDelphi"), 
                                        c("Tissue", "FORECasT")), paired = T) +
  scale_fill_manual(values = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")) + 
  scale_color_manual(values = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")) + 
  xlab("") + ylab("Normalized Delete Left NT(%)") +
  theme_bw()
dev.off()
total_plot_data3 <- total_plot_data2[total_plot_data2$type != "inDelphi",]
total_plot_data3$type <- as.character(total_plot_data3$type)
total_plot_data3$type <- factor(total_plot_data3$type, levels = c( "Tissue","Cell", "FORECasT"))
total_plot_data3$beex[total_plot_data3$type == "FORECasT"] <- total_plot_data3$beex[total_plot_data3$type == "FORECasT"] - 1
pdf("~/Nutstore Files/Tobin/Previous/ins1Mut_link_plot_add_predict_rmIndelphi.pdf", width = 6, height = 4)

ggplot(total_plot_data3, aes(x = type, y = diff)) + 
  geom_violin(aes(fill = type), 
              alpha = 0.5, 
              scale = "width", 
              adjust = 1, 
              trim = TRUE) + 
  ggnewscale::new_scale("fill") +
  geom_point(aes(x = beex, fill=NT),
             size = 2, 
             shape = 24,
             color = "black") + 
  geom_line(aes(x = beex, y = diff, color = NT, 
                group = id), 
            linewidth=0.1) + 
  stat_compare_means(comparisons = list(c("Tissue", "Cell"), 
                                        c("Cell", "FORECasT"),
                                        c("Tissue", "FORECasT")), paired = T) +
  scale_fill_manual(values = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")) + 
  scale_color_manual(values = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")) + 
  xlab("") + ylab("Normalized Delete Left NT(%)") +
  theme_bw()
dev.off()


total_plot_data3 <- total_plot_data2[total_plot_data2$type != "FORECasT",]
total_plot_data3$type <- as.character(total_plot_data3$type)
total_plot_data3$type <- factor(total_plot_data3$type, levels = c( "Tissue","Cell", "inDelphi"))
pdf("~/Nutstore Files/Tobin/Previous/ins1Mut_link_plot_add_predict_rmForecast.pdf", width = 5.3, height = 4)

ggplot(total_plot_data3, aes(x = type, y = diff)) + 
  geom_violin(aes(fill = type), 
              alpha = 0.5, 
              scale = "width", 
              adjust = 1, 
              trim = TRUE) + 
  ggnewscale::new_scale("fill") +
  geom_point(aes(x = beex, fill=NT),
             size = 2, 
             shape = 24,
             color = "black") + 
  geom_line(aes(x = beex, y = diff, color = NT, 
                group = id), 
            linewidth=0.1) + 
  stat_compare_means(comparisons = list(c("Tissue", "Cell"), 
                                        c("Cell", "inDelphi"),
                                        c("Tissue", "inDelphi")), paired = T) +
  scale_fill_manual(values = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")) + 
  scale_color_manual(values = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")) + 
  xlab("") + ylab("Normalized Delete Left NT(%)") +
  theme_bw()
dev.off()


