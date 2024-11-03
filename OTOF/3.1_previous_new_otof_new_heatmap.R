###load 之前114个sgRNA的结果
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

###读取新测的4个sgRNA结果 然后合并到之前的结果里面去
print(load("~/data/project/ear_project/gene_therapy_ll/Otof_cell_tissue_new_data_newer.rda"))
rep1_new_otof <- otof_cell_sg12_15[seq(1,8, 2)]
rep2_new_otof <- otof_cell_sg12_15[seq(2,8, 2)]
for(name in names(rep1_new_otof)){
  name2 <- str_remove(name, "[-].*")
  old_rep1[[name2]] <- rep1_new_otof[[name]]
}
for(name in names(rep2_new_otof)){
  name2 <- str_remove(name, "[-].*")
  old_rep2[[name2]] <- rep2_new_otof[[name]]
}
split_name_cell <- data.frame(rbind(split_name_cell, 
                                    data.frame("cellid" = paste0("m", 12 : 15),
                                               "tissueid" = paste0("m", 12 : 15),
                                               "mutation" = 'w-Otof-1236dC', 
                                               "type" = 'd')))
otof_tissue_sg7_15_rep1 <- otof_tissue_sg7_15[c(9,10,11)]
otof_tissue_sg7_15_rep2 <- otof_tissue_sg7_15[c(2:5)]
names(otof_tissue_sg7_15_rep1) <- c("m12-tissue", "m13-tissue", "m14-tissue")
names(otof_tissue_sg7_15_rep2) <- c("m12-tissue", "m13-tissue", "m14-tissue", "m15-tissue")


for(name in names(otof_tissue_sg7_15_rep1)){
  old_rep1[[name]] <- otof_tissue_sg7_15_rep1[[name]]
}
for(name in names(otof_tissue_sg7_15_rep2)){
  old_rep2[[name]] <- otof_tissue_sg7_15_rep2[[name]]
}
split_name_tissue <- data.frame(rbind(split_name_tissue, 
                                      data.frame("cellid" = paste0("m", 12 : 15, "-tissue"),
                                                 "tissueid" = paste0("m", 12 : 15, "-tissue"),
                                                 "mutation" = 'w-Otof-1236dC', 
                                                 "type" = 'd')))

###因为加入了几个新的样本 所以重新开始计算均值
old_rep123 <- list()
total_samples <- unique(c(names(old_rep1), names(old_rep2), names(old_rep3)))
old_result <- list(Rep1 = old_rep1, Rep2 = old_rep2, Rep3 = old_rep3)

for(name in total_samples){
  
  sample_count <- sum(unlist(lapply(old_result, function(x){
    name %in% names(x)
  })))
  tmp_rep123 <- list()
  for(sp in names(old_result)){
    if(name %in% names(old_result[[sp]])){
      tmp_rep123[[sp]] <- old_result[[sp]][[name]]
    }
  }
  tmp_rep123 <- data.frame(do.call(rbind, tmp_rep123))
  tmp_rep123 <- lapply(split(tmp_rep123$fq, tmp_rep123$indel_size), sum)
  tmp_rep123 <- data.frame(indel_size = as.numeric(names(tmp_rep123)), fq = unlist(tmp_rep123) / sample_count)
  tmp_rep123 <- tmp_rep123[order(tmp_rep123$indel_size),]
  old_rep123[[name]] <- tmp_rep123
}
old_tissue <- old_rep123[str_remove(names(old_rep123), "CRISPResso_on_") %in% split_name_tissue$tissueid]
old_cell <- old_rep123[str_remove(names(old_rep123), "CRISPResso_on_") %in% split_name_cell$cellid]

###接下来根据是Insert突变还是Delete突变来计算Del1与Ins1的编辑比例

names(old_tissue) <- str_remove(names(old_tissue), "CRISPResso_on_")
names(old_cell) <- str_remove(names(old_cell), "CRISPResso_on_")
# old_result <- list(Rep1 = old_rep1, Rep2 = old_rep2, Rep3 = old_rep3)
# total_sample <- unique(c(split_name_tissue$tissueid, split_name_cell$cellid))


# old_sample_noAA <- lapply(old_result, function(x){
#   
#   names(x) <- str_remove(names(x), "CRISPResso_on_")
#   tmp_sample <- total_sample[total_sample %in% names(x) ]
#   tmp <- lapply(tmp_sample, function(y){
#     if(y %in% split_name_cell$cellid){
#       type <- split_name_cell$type[split_name_cell$cellid == y]
#     }else{
#       type <- split_name_tissue$type[split_name_tissue$tissueid == y]
#     }
#     stat <- x[[y]]
#     stat <- stat[stat[,1] != 0,]
#     if(type == 'd'){
#       
#       pct <- stat[,2] / sum(stat[,2]) * 100
#       return(pct[stat[,1] == 1])
#     }
#     if(type == 'i'){
#       pct <- stat[,2] / sum(stat[,2]) * 100
#       return(pct[stat[,1] == -1])
#     }
#   })
#   names(tmp) <- tmp_sample
#   tmp
# })
# noAA <- matrix(NA, nrow = 3, ncol = length(total_sample))
# rownames(noAA) <- c("Rep1", "Rep2", "Rep3")
# colnames(noAA) <- total_sample
# for(rep in names(old_sample_noAA)){
#   rep_res <- old_sample_noAA[[rep]]
#   for(sample in names(rep_res)){
#     noAA[rep, sample] <- rep_res[[sample]]
#   }
# }
split_name_cell <- data.frame(do.call(rbind, 
                                      lapply(split(split_name_cell, 
                                                   split_name_cell$mutation), 
                                             function(x){
                                               x$id <- 1 : nrow(x)
                                               x
                                             })))
split_name_tissue <- data.frame(do.call(rbind, 
                                      lapply(split(split_name_tissue, 
                                                   split_name_tissue$mutation), 
                                             function(x){
                                               x$id <- 1 : nrow(x)
                                               x
                                             })))
insert_mut_cell <- split_name_cell[split_name_cell$type == 'i',]
delete_mut_cell <- split_name_cell[split_name_cell$type == 'd',]
insert_mut_tissue <- split_name_tissue[split_name_tissue$type == 'i',]
delete_mut_tissue <- split_name_tissue[split_name_tissue$type == 'd',]

# insert_cell_noAA <- noAA[,insert_mut_cell$cellid]
# insert_mut_cell$Pct <- unlist(colMeans(insert_cell_noAA, na.rm = T))
# delete_cell_noAA <- noAA[,delete_mut_cell$cellid]
# delete_mut_cell$Pct <- unlist(colMeans(delete_cell_noAA, na.rm = T))

insert_mut_cell$Pct <- unlist(lapply(insert_mut_cell$cellid, function(x){
  stat <- old_cell[[x]]
  stat[stat[,1] == 0,2] <- 0
  pct <- stat[,2] / sum(stat[,2]) * 100
  return(pct[stat[,1] == -1])
}))
delete_mut_cell$Pct <- unlist(lapply(delete_mut_cell$cellid, function(x){
  stat <- old_cell[[x]]
  stat[stat[,1] == 0,2] <- 0
  pct <- stat[,2] / sum(stat[,2]) * 100
  return(pct[stat[,1] == 1])
}))

# insert_tissue_noAA <- noAA[,insert_mut_tissue$tissueid]
# insert_mut_tissue$Pct <- unlist(colMeans(insert_tissue_noAA, na.rm = T))
# delete_tissue_noAA <- noAA[,delete_mut_tissue$tissueid]
# delete_mut_tissue$Pct <- unlist(colMeans(delete_tissue_noAA, na.rm = T))
# plotStat_merge_v2(insertOne_mean_reorder, noAA = insertOne_noAA,  title = "insert")
# plotStat_merge_v2(deleteOne_mean_reorder, noAA = deleteOne_noAA,  title = "delete")

insert_mut_tissue$Pct <- unlist(lapply(insert_mut_tissue$tissueid, function(x){
  stat <- old_tissue[[x]]
  stat[stat[,1] == 0,2] <- 0
  pct <- stat[,2] / sum(stat[,2]) * 100
  return(pct[stat[,1] == -1])
}))
delete_mut_tissue$Pct <- unlist(lapply(delete_mut_tissue$tissueid, function(x){
  stat <- old_tissue[[x]]
  stat[stat[,1] == 0,2] <- 0
  pct <- stat[,2] / sum(stat[,2]) * 100
  return(pct[stat[,1] == 1])
}))

insert_mut_cell_reorder <- insert_mut_cell[order(insert_mut_cell$Pct, decreasing = T),]
delete_mut_cell_reorder <- delete_mut_cell[order(delete_mut_cell$Pct, decreasing = T),]
insert_mut_tissue_reorder <- insert_mut_tissue[order(insert_mut_tissue$Pct, decreasing = T),]
delete_mut_tissue_reorder <- delete_mut_tissue[order(delete_mut_tissue$Pct, decreasing = T),]



old_cell_insert <- old_cell[insert_mut_cell_reorder$cellid]
names(old_cell_insert) <- paste(insert_mut_cell_reorder$mutation, 
                                insert_mut_cell_reorder$id, sep = "-")

old_cell_delete <- old_cell[delete_mut_cell_reorder$cellid]
names(old_cell_delete) <- paste(delete_mut_cell_reorder$mutation, 
                                delete_mut_cell_reorder$id, sep = "-")


old_tissue_insert <- old_tissue[insert_mut_tissue_reorder$tissueid]
names(old_tissue_insert) <- paste(insert_mut_tissue_reorder$mutation, 
                                  insert_mut_tissue_reorder$id, sep = "-")
old_tissue_delete <- old_tissue[delete_mut_tissue_reorder$tissueid]
names(old_tissue_delete) <- paste(delete_mut_tissue_reorder$mutation, 
                                  delete_mut_tissue_reorder$id, sep = "-")
# tmp <- unlist(lapply(insert_mut_tissue_reorder$mutation, function(x){unlist(strsplit(x, "[-]"))[2]}))
# tmp <- unlist(lapply(delete_mut_tissue_reorder$mutation, function(x){unlist(strsplit(x, "[-]"))[2]}))
# 
# edit_stat_of_distance$geneFactor <- factor(edit_stat_of_distance$gene, levels = c(
#   "Myo7a","USH1C","CDH23","PCDH15","USH2A","VLGR1","whirlin","TMC1","OTOF"
# ))

pdf("~/data/project/ear_project/gene_therapy_ll/Previews/Result/indel_stat_heatmap_tissue_insert_new.pdf", width = 10, height = 12)
plotStat(old_tissue_insert, region = -18 : 5, title = "Tissue_insert", showname = T,
         shownumber = T)
dev.off()
pdf("~/data/project/ear_project/gene_therapy_ll/Previews/Result/indel_stat_heatmap_tissue_delete_new.pdf", width = 30, height = 12)

plotStat(old_tissue_delete, region = -18 : 5, title = "Tissue_delete", showname = T,
         shownumber = T)
dev.off()
pdf("~/data/project/ear_project/gene_therapy_ll/Previews/Result/indel_stat_heatmap_cell_insert_new.pdf", width = 10, height = 12)
plotStat(old_cell_insert, region = -18 : 5, title = "Cell_insert", showname = T,
         shownumber = T)
dev.off()
pdf("~/data/project/ear_project/gene_therapy_ll/Previews/Result/indel_stat_heatmap_cell_delete_new.pdf", width = 30, height = 12)

plotStat(old_cell_delete, region = -18 : 5, title = "Cell_delete", showname = T,
         shownumber = T)
dev.off()
