print(load("~/data/project/ear_project/gene_therapy_ll/Result/indel_stat_total.rda"))
names(sg28) <- str_remove(names(sg28), "_.*")
sgs <- paste0("sg", c(1 : 7, 10, 11, 12, 14 : 17, 19:21,24,25, 27))

sg28_single_sel <- sg28[sgs]
sg28_single_sel <- lapply(sg28_single_sel, function(x){
  samples <- str_remove(names(x), "./")
  names(x) <- samples
  multi_samples <- samples[grep("^B", samples)]
  single_samples <- x[!samples %in% multi_samples]
  names(single_samples) <- c("A", "B", "C")
  single_samples[unlist(lapply(single_samples, function(y){
    sum(y[,2]) != 0
  }))]
})
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/First1NT/2024_1_12_integrated_design_result.xlsx")
sgRNA$id2 <- paste("Sg", sgRNA$new_id, sep = "-")
sgRNA$id2 <- str_replace_all(sgRNA$id2, "-", "_")
sg28 <- read.xlsx("~/data/project/ear_project/gene_therapy/find1NTSim/sg28_redesign_fixbad_fix_sgprimer.xlsx")
sg28$order <- 1 : 28
table(sg28$sgRNA == sgRNA$sgRNA[1:28])
for(i in 1 : length(sg28_single_sel)){
  n <- as.numeric(str_remove(names(sg28_single_sel)[i], "sg"))
  names(sg28_single_sel)[i] <- sgRNA$id2[n]
}



setwd("~/data/project/ear_project/gene_therapy_ll/batch2/")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/indel_stat2.rda")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599/indel_stat.rda")


batch2 <- batch2_1
for(name in names(batch2_2)){
  batch2[[name]] <- batch2_2[[name]]
}
total_sg <- unique(unlist(lapply(batch2, names)))
batch2_rev <- lapply(total_sg, function(sg){
  res <- list()
  for(name in names(batch2)){
    n <- str_remove(name, "CRISPRessoPooled_on_")
    n <- str_remove(n, "[0-9]+")
    if(sg %in% names(batch2[[name]])){
      res[[n]] <- batch2[[name]][[sg]]
    }
  }
  res
})
names(batch2_rev) <- total_sg
total_data <- batch2_rev
for(name in names(sg28_single_sel)){
  total_data[[name]] <- sg28_single_sel[[name]]
}


max_region <- c(-200, 30)
sg28_mean <- lapply(sg28_single_sel, function(x){
  if(length(x) == 1){
    return(x[[1]])
  }
  fq <- unlist(rowMeans(data.frame(do.call(cbind, lapply(x, function(y){y[,2]})))))
  x[[1]][,2] <- fq
  x[[1]]
})
batch2_total_mean <- batch2_mean
for(name in names(batch2_mean2)){
  batch2_total_mean[[name]] <- batch2_mean2[[name]]
}
total_mean <- batch2_total_mean
for(name in names(sg28_mean)){
  total_mean[[name]] <- sg28_mean[[name]]
}

####

indel1_sample <- total_mean[unlist(lapply(total_mean, function(x){
  x[x[,1] == 0, 2] <- 0
  abs(as.numeric(x[which.max(x[,2]),1])) == 1
}))]

sgRNA_pair <- split(sgRNA$id2, sgRNA$id)
sgRNA_pair_remain <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(indel1_sample)) == 2)
}))]
remain_sample <- unlist(sgRNA_pair_remain)
sgRNA$pos <- unlist(lapply(sgRNA$id, function(x){
  if(!is.na(str_match(x, "last"))){
    return(1)
  }
  unlist(strsplit(x, "[_-]"))[2]
}))
remain_sample <- sgRNA[sgRNA$id2 %in% remain_sample,]
remain_sample <- remain_sample[order(remain_sample$pos, remain_sample$id),]

total_mean_remain <- total_mean[remain_sample$id2]





indel1_pct <- list()
for(name in names(sgRNA_pair_remain)){
  ids <- sgRNA_pair_remain[[name]]
  sel_sg1 <- total_mean[[ids[1]]]
  # sel_sg1 <- sel_sg1[sel_sg1$indel_size %in% region,]
  sel_sg2 <- total_mean[[ids[2]]]
  # sel_sg2 <- sel_sg2[sel_sg2$indel_size %in% region,]
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
  sel_sg2[,3] <- sel_sg2[,2] / sum(sel_sg2[,2])
  sg1_se <- plotrix::std.error(unlist(lapply(total_data[[ids[1]]], function(sel_sg1){
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
    sel_sg1[which.max(sel_sg1[,2]),3]
  })))
  sg2_se <- plotrix::std.error(unlist(lapply(total_data[[ids[2]]], function(sel_sg2){
    sel_sg2[sel_sg2[,1] == 0, 2] <- 0
    sel_sg2[,3] <- sel_sg2[,2] / sum(sel_sg2[,2])
    sel_sg2[which.max(sel_sg2[,2]),3]
  })))
  
  indel1_pct[[name]] <- data.frame(Percent = c(sel_sg1[which.max(sel_sg1[,2]),3],
                                               sel_sg2[which.max(sel_sg2[,2]),3]), 
                                   ID = ids, 
                                   se = c(sg1_se, sg2_se))
}
indel1_pct <- data.frame(do.call(rbind, indel1_pct))





###去掉那些不匹配的样本
sgRNA_pair <- split(sgRNA$id2, sgRNA$id)
sgRNA_pair_rm <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(total_mean)) != 2)
}))]
sgRNA_pair_remain <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(total_mean)) == 2)
}))]
###计算均值的相关性，因为本次数据每次样本是分开测得 不好AA BB CC一一对比了
total_pair_sg_cor <- lapply(sgRNA_pair_remain, function(y){
  ids <- unlist(y)
  sel_sg1 <- total_mean[[ids[1]]]
  sel_sg2 <- total_mean[[ids[2]]]
  max_del <- min(c(sel_sg1[,1], sel_sg2[,1]))
  max_del <- max(c(max_del, max_region[1]))
  max_in <- max(c(sel_sg1[,1], sel_sg2[,1]))
  max_in <- min(c(max_in, max_region[2]))
  region <- max_del : max_in
  tmp_region <- region[!region %in% sel_sg1[,1]]
  if(length(tmp_region) != 0){
    sel_sg1 <- data.frame(rbind(sel_sg1, data.frame(indel_size = tmp_region, fq = 0)))
  }
  tmp_region <- region[!region %in% sel_sg2[,1]]
  if(length(tmp_region) != 0){
    sel_sg2 <- data.frame(rbind(sel_sg2, data.frame(indel_size = tmp_region, fq = 0)))
  }
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  common_indel <- intersect(sel_sg1[,1], sel_sg2[,1])
  sel_sg1 <- sel_sg1[sel_sg1[,1] %in% common_indel,]
  sel_sg2 <- sel_sg2[sel_sg2[,1] %in% common_indel,]
  sel_sg1 <- sel_sg1[order(sel_sg1[,1]),]
  sel_sg2 <- sel_sg2[order(sel_sg2[,1]),]
  # rm0 <- sel_sg1[,2] == 0 & sel_sg2[,2] == 0
  #cor.test(sel_sg1[!rm0,2], sel_sg2[!rm0,2])
  cor.test(sel_sg1[,2], sel_sg2[,2])
})
names(total_pair_sg_cor) <- unlist(lapply(sgRNA_pair_remain, function(x){
  paste(x, collapse = "-")
}))
total_pair_sg_cor <- data.frame(pair = names(total_pair_sg_cor), 
                                cor = unlist(lapply(total_pair_sg_cor, function(x){x$estimate})),
                                pval = unlist(lapply(total_pair_sg_cor, function(x){x$p.value})))

sg_name <- data.frame(id = names(sgRNA_pair_remain), 
                      pair = unlist(lapply(sgRNA_pair_remain, function(x){
                        paste(x, collapse = "-")
                      })))
sg_name$pos <- unlist(lapply(sg_name$id, function(x){
  if(!is.na(str_match(x, "last"))){
    return(1)
  }
  unlist(strsplit(x, "[_-]"))[2]
}))
total_pair_sg_cor <- merge(total_pair_sg_cor, sg_name, by = "pair")

total_pair_sg_cor$pos <- as.integer(total_pair_sg_cor$pos)


sample_order <- apply(total_pair_sg_cor, 1, function(x){
  x <- unlist(x)
  sg <- unlist(strsplit(x[1], "[-]"))
  data.frame(sg = sg, cor = as.numeric(x[2]), pair = x[1], id = x[4], pos = x[5])
})
sample_order <- data.frame(do.call(rbind, sample_order))

sample_order$indel <- unlist(lapply(unlist(sample_order$pair), function(sg){
  ids <- unlist(strsplit(sg, "[-]"))
  sel_sg1 <- total_mean[[ids[1]]]
  sel_sg2 <- total_mean[[ids[2]]]
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
  sel_sg2[,3] <- sel_sg2[,2] / sum(sel_sg2[,2])
  peak_indel <- (as.numeric(c(sel_sg1[,1], sel_sg2[,1])[which.max(c(sel_sg1[,2], sel_sg2[,2]))]))
}))
sample_order$common <- unlist(lapply(unlist(sample_order$pair), function(sg){
  ids <- unlist(strsplit(sg, "[-]"))
  sel_sg1 <- total_mean[[ids[1]]]
  sel_sg2 <- total_mean[[ids[2]]]
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  sel_sg1 <- sel_sg1[sel_sg1[,1] %in% c(-10:10),]
  sel_sg2 <- sel_sg2[sel_sg2[,1] %in% c(-10:10),]
  sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
  sel_sg2[,3] <- sel_sg2[,2] / sum(sel_sg2[,2])
  sel_sg1[which.max(sel_sg1[,2]),1] == sel_sg2[which.max(sel_sg2[,2]),1]
}))
sample_order <- sample_order[order(sample_order$common,
                                   -abs(sample_order$indel), 
                                   sample_order$indel,
                                   sample_order$cor, 
                                   sample_order$id, 
                                   decreasing = T),]




###把cmp信息拿出来看看取出insert一个碱基与delete一个碱基的区分开来
tmp <- read.table("~/data/project/ear_project/gene_therapy/find1NTSim/human_sg_final_anno2_addid.txt", sep="\t")
tmp2 <- read.table("~/data/project/ear_project/gene_therapy/find1NTSim/threshold20000/human_sg_final_anno2_addid.txt", sep="\t")
tmp3 <- read.table("~/data/project/ear_project/gene_therapy/find1NTSim/human_sg_lastNTSim_final_addid2.txt", sep = "\t")
tmp3$V14 <- lapply(1 : nrow(tmp3), function(x){
  if(str_sub(tmp3[x,2], 2, 20) == str_sub(tmp3[x, 8], 1, 19)){
    return(paste0(paste0(tmp3[x,2], "-"), " ", paste0("-",tmp3[x,8])))
  }
  return(paste0(paste0(tmp3[x,8], "-"), " ", paste0("-",tmp3[x,2])))
})
tmp4 <-unlist(c(tmp[,14], tmp2[,14], tmp3[,14]))
tmp <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[1]
}))
tmp2 <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[2]
}))
tmp4 <- data.frame(sgRNA = str_remove(c(tmp, tmp2), "[-]"), 
                   Cmp = c(tmp, tmp2))
tmp4 <- tmp4[!duplicated(tmp4$sgRNA),]
sgCmp <- tmp4
rm(tmp, tmp2, tmp3, tmp4)

sgCmp <- merge(sgCmp, sgRNA, by="sgRNA")
sgCmp$test <- unlist(lapply(1 : nrow(sgCmp), function(i){
  cmp <- sgCmp$Cmp[i]
  cmp_pair <- sgCmp$Cmp[sgCmp$id == sgCmp$id[i]]
  cmp_pair <- cmp_pair[!cmp_pair %in% cmp]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  if(!(this_pos == 1 | pair_pos == 1)){
    return("Last")
  }
  return("First")
}))
sgCmp$isInsert <- unlist(lapply(1 : nrow(sgCmp), function(i){
  cmp <- sgCmp$Cmp[i]
  cmp_pair <- sgCmp$Cmp[sgCmp$id == sgCmp$id[i]]
  cmp_pair <- cmp_pair[!cmp_pair %in% cmp]
  cmp <- unlist(strsplit(cmp, "*"))
  cmp_pair <- unlist(strsplit(cmp_pair, "*"))
  this_pos <- which(cmp == '-')
  pair_pos <- which(cmp_pair == '-')
  #在第一位有一个插入情况下，说明后面可能是last或者非last的情况
  if(this_pos == 1 | pair_pos == 1){
    
    #说明是last位有插入 且在这个sgRNA插入
    if(pair_pos == length(cmp)){
      return("Insert")
    }
    #说明是last位有插入 且不在这个sgRNA插入
    if(this_pos == length(cmp)){
      return("Delete")
    }

    #说明不是last位有插入，那就看谁在第一位有-那就是后面有插入
    if(this_pos == 1){
      return("Insert")
    }else{
      return("Delete")
    }
    
  } else {
    #说明是最后一个碱基出现了shift，那就看最早的-是出现在谁那边谁就有delete
    if(this_pos < pair_pos){
      return("Delete")
    } else{
      return("Insert")
    }
  }
}))
rownames(sgCmp) <- sgCmp$id2

sample_order$isInsert <- sgCmp[sample_order$sg, "isInsert"]
sample_order$test <- sgCmp[sample_order$sg, "test"]

sample_order_insertOne <- sample_order[sample_order$isInsert == "Insert",]
sample_order_deleteOne <- sample_order[sample_order$isInsert != "Insert",]

total_mean_reorder <- total_mean[sample_order$sg]

insertOne_mean_reorder <- total_mean[sample_order_insertOne$sg]
deleteOne_mean_reorder <- total_mean[sample_order_deleteOne$sg]



indel1_pct <- list()
for(name in names(total_mean_reorder)){
  sg1_se <- unlist(lapply(total_data[[name]], function(sel_sg1){
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    type <- sample_order$isInsert[sample_order$sg == name]
    
    sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
    if(type == "Insert"){
      return(sel_sg1[sel_sg1[,1] == -1,3])
    } else{
      return(sel_sg1[sel_sg1[,1] == 1,3])
    }
    
  }))
  
  indel1_pct[[name]] <- data.frame(Percent = c(sg1_se), 
                                   ID = name, sample = names(sg1_se))
}
indel1_pct <- data.frame(do.call(rbind, indel1_pct))
noAA <- matrix(NA, nrow = 3, ncol = length(total_mean_reorder))
colnames(noAA) <- names(total_mean_reorder)
rownames(noAA) <- c("A", "B", "C")
for(i in 1 : nrow(indel1_pct)){
  noAA[indel1_pct$sample[i], indel1_pct$ID[i]] <- indel1_pct$Percent[i]
}

insertOne_noAA <- noAA[,sample_order_insertOne$sg]
deleteOne_noAA <- noAA[,sample_order_deleteOne$sg]

# pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_indel1.pdf", width = 30, height = 5)
# plotStat_merge(total_mean_reorder, noAA = noAA,  title = "tmp")
# dev.off()


pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_indel1_v4.pdf", width = 20, height = 5)
plotStat_merge_v2(insertOne_mean_reorder, noAA = insertOne_noAA,  title = "insert")
plotStat_merge_v2(deleteOne_mean_reorder, noAA = deleteOne_noAA,  title = "delete")

dev.off()



###按照indel1的比例进行排序

sample_order$indel1Pct <- unlist(lapply(sample_order$sg, function(name){
  sg1_se <- lapply(total_data[[name]], function(sel_sg1){
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    type <- sample_order$isInsert[sample_order$sg == name]
    
    sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
    sel_sg1[abs(sel_sg1[,1]) == 1, ]
  })
  sg1_se <- data.frame(do.call(rbind, sg1_se))
  sg1_se <- lapply(split(sg1_se$V3, sg1_se$indel_size), mean)
  sg1_se <- data.frame(indel = names(sg1_se), Pct = unlist(sg1_se))
  sg1_se$Pct[sg1_se$indel == '1'] - sg1_se$Pct[sg1_se$indel == '-1']
}))

sample_order$del1Pct <- unlist(lapply(sample_order$sg, function(name){
  sg1_se <- lapply(total_data[[name]], function(sel_sg1){
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    type <- sample_order$isInsert[sample_order$sg == name]
    
    sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
    sel_sg1[sel_sg1[,1] == -1, ]
  })
  sg1_se <- data.frame(do.call(rbind, sg1_se))
  mean(sg1_se$V3)
}))

sample_order$del4Pct <- unlist(lapply(sample_order$sg, function(name){
  sg1_se <- lapply(total_data[[name]], function(sel_sg1){
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    type <- sample_order$isInsert[sample_order$sg == name]
    
    sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
    sel_sg1[sel_sg1[,1] == -4, ]
  })
  sg1_se <- data.frame(do.call(rbind, sg1_se))
  mean(sg1_se$V3)
}))

sample_order$Ins2Pct <- unlist(lapply(sample_order$sg, function(name){
  sg1_se <- lapply(total_data[[name]], function(sel_sg1){
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    type <- sample_order$isInsert[sample_order$sg == name]
    
    sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
    sel_sg1[sel_sg1[,1] == 2, ]
  })
  sg1_se <- data.frame(do.call(rbind, sg1_se))
  mean(sg1_se$V3)
}))

sample_order$Ins1Pct <- unlist(lapply(sample_order$sg, function(name){
  sg1_se <- lapply(total_data[[name]], function(sel_sg1){
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    type <- sample_order$isInsert[sample_order$sg == name]
    
    sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
    sel_sg1[sel_sg1[,1] == 1, ]
  })
  sg1_se <- data.frame(do.call(rbind, sg1_se))
  mean(sg1_se$V3)
}))

sample_order$Ins4Pct <- unlist(lapply(sample_order$sg, function(name){
  sg1_se <- lapply(total_data[[name]], function(sel_sg1){
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    type <- sample_order$isInsert[sample_order$sg == name]
    
    sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
    sel_sg1[sel_sg1[,1] == 4, ]
  })
  sg1_se <- data.frame(do.call(rbind, sg1_se))
  mean(sg1_se$V3)
}))

sample_order$Del2Pct <- unlist(lapply(sample_order$sg, function(name){
  sg1_se <- lapply(total_data[[name]], function(sel_sg1){
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    type <- sample_order$isInsert[sample_order$sg == name]
    
    sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
    sel_sg1[sel_sg1[,1] == -2, ]
  })
  sg1_se <- data.frame(do.call(rbind, sg1_se))
  mean(sg1_se$V3)
}))

sample_order2 <- sample_order
#取平均的结果不好
# for(i in seq(1, nrow(sample_order2), 2)){
  # sample_order2$indel1Pct[i] <- sample_order2$indel1Pct[i + 1] <-
    # (sample_order2$indel1Pct[i] + sample_order2$indel1Pct[i + 1]) / 2
# }
#以insert为主
# for(i in seq(1, nrow(sample_order2), 2)){
# 
#   sample_order2$indel1Pct[i] <- sample_order2$indel1Pct[i + 1] <-
#    sample_order2$indel1Pct[c(i, (i + 1))][sample_order2$isInsert[c(i, (i + 1))] == "Insert"]
# }
# 以delete为主
# for(i in seq(1, nrow(sample_order2), 2)){
# 
#   sample_order2$indel1Pct[i] <- sample_order2$indel1Pct[i + 1] <-
#    sample_order2$indel1Pct[c(i, (i + 1))][sample_order2$isInsert[c(i, (i + 1))] == "Delete"]
# }
# sample_order2$indel1Pct <- sample_order2$indel1Pct * ifelse(sample_order2$common, 2, 0.1)
# sample_order2 <- sample_order2[order(sample_order2$indel1Pct,
#                                      sample_order2$pair,
# 
#                                      decreasing = T),]

sample_order_insertOne <- sample_order2[sample_order2$isInsert == "Insert",]
sample_order_insertOne <- sample_order_insertOne[order(sample_order_insertOne$del1Pct, sample_order_insertOne$del4Pct, sample_order_insertOne$Ins2Pct, decreasing = T),]
sample_order_deleteOne <- sample_order2[sample_order2$isInsert != "Insert",]
sample_order_deleteOne <- sample_order_deleteOne[order(sample_order_deleteOne$Ins1Pct, sample_order_deleteOne$Ins4Pct, sample_order_deleteOne$Del2Pct, decreasing = T),]

# sample_order_insertOne <- sample_order_insertOne[order(sample_order_insertOne$indel1Pct, decreasing = T), ]
# sample_order_deleteOne <- sample_order_deleteOne[order(sample_order_deleteOne$indel1Pct, decreasing = T), ]


insertOne_mean_reorder <- total_mean[sample_order_insertOne$sg]
deleteOne_mean_reorder <- total_mean[sample_order_deleteOne$sg]

insertOne_noAA <- noAA[,sample_order_insertOne$sg]
deleteOne_noAA <- noAA[,sample_order_deleteOne$sg]

# insertOne_mean_reorder <- insertOne_mean_reorder[rev(1:73)]
# insertOne_noAA <- insertOne_noAA[,rev(1:73)]
pdf("~/data/project/ear_project/gene_therapy_ll/Result/merged_indel1_v9.pdf", width = 20, height = 5)
plotStat_merge_v3(insertOne_mean_reorder, noAA = insertOne_noAA,  title = "Insert Bulge", plotOut = F)
plotStat_merge_v3(deleteOne_mean_reorder, noAA = deleteOne_noAA,  title = "Delete Bulge", plotOut = F)

# dev.off()



 ####unused

test1 <- total_edit_pattern$C$Sg_7_44
test1 <- test1[test1$indel_size != 0,]
test1 <- test1[order(test1$fq ,decreasing = T),]
test2 <- total_edit_pattern$C$Sg_21_145
test2 <- test2[test2$indel_size != 0,]
test2 <- test2[order(test2$fq, decreasing = T),]


plot_data <- total_pair_sg_cor
table(unlist(lapply(total_pair_sg_cor$pair, function(x){unlist(strsplit(x, "[-]")) %in% sgCmp$id2})))
plot_data$last <- unlist(lapply(total_pair_sg_cor$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))
  if(sum(x %in% sgCmp$id2[sgCmp$test == "First"]) != 0){
    return("First")
  } else {
    return("Last")
  }
}))
pdf("~/data/project/ear_project/gene_therapy_ll/Result/cor_of_differentNT_position.pdf", width = 6, height = 4)
ggplot(plot_data, aes(x = last, y = cor)) +
  geom_boxplot(outlier.size = 0, aes(color = last)) + 
  geom_point(aes(color = last), position = position_dodge2(width = 0.7)) +
  theme_bw() + 
  ggpubr::stat_compare_means(comparisons = list(c("First", "Last")), method = "wilcox") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        text = element_text(size=15))
dev.off()

