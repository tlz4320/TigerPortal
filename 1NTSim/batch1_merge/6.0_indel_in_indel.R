#在2.0_total_First_pop_peak_indel1.R
setwd("~/data/project/ear_project/gene_therapy/find1NTSim/")
tmp <- read.table("human_sg_final_anno2_addid.txt", sep="\t")
tmp2 <- read.table("threshold20000/human_sg_final_anno2_addid.txt", sep="\t")
tmp3 <- read.table("human_sg_lastNTSim_final_addid2.txt", sep = "\t")
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
table(sgRNA$sgRNA %in% sgCmp$sgRNA)
high_cor_sample <- total_pair_sg_cor[total_pair_sg_cor$cor > 0.8,]
not_prefer_select <- list()
region <- c(-10: 10)
for(name in high_cor_sample$id){
  ids <- sgRNA_pair_remain[[name]]
  sel_sg1 <- total_mean[[ids[1]]]
  # sel_sg1 <- sel_sg1[sel_sg1$indel_size %in% region,]
  sel_sg2 <- total_mean[[ids[2]]]
  # sel_sg2 <- sel_sg2[sel_sg2$indel_size %in% region,]
  sel_sg1[sel_sg1[,1] == 0, 2] <- 0
  sel_sg2[sel_sg2[,1] == 0, 2] <- 0
  sel_sg1[,3] <- sel_sg1[,2] / sum(sel_sg1[,2])
  sel_sg2[,3] <- sel_sg2[,2] / sum(sel_sg2[,2])
  peak_indel <- (as.numeric(c(sel_sg1[,1], sel_sg2[,1])[which.max(c(sel_sg1[,2], sel_sg2[,2]))]))
  sg1 <- sgCmp$Cmp[sgCmp$sgRNA == sgRNA$sgRNA[sgRNA$id2 == ids[1]]]
  sg2 <- sgCmp$Cmp[sgCmp$sgRNA == sgRNA$sgRNA[sgRNA$id2 == ids[2]]]
  use <- "sg1"
  if(peak_indel < 0){
    if(!is.na(unlist(str_match(name, "last")))){
      if(str_sub(sg1, 1, 1) == "-"){
        use <- "sg1"
      }
      else{
        use <- "sg2"
      }
    }
    else{
      if(str_sub(sg1, 1, 1) == "-" | str_sub(sg1, 21, 21) == "-"){
        use <- "sg1"
      }
      else{
        use <- "sg2"
      }
    }
  }  else{
    if(!is.na(unlist(str_match(name, "last")))){
      if(str_sub(sg1, 1, 1) == "-"){
        use <- "sg2"
      }else{
        use <- "sg1"
      }
    }
    else{
      if(str_sub(sg1, 1, 1) == "-" | str_sub(sg1, 21, 21) == "-"){
        use <- "sg2"
      }else{
        use <- "sg1"
      }
    }
  }
  if(use == "sg1"){
    use <- "sg2"
  }
  else{
    use <- "sg1"
  }
  not_prefer_select[[name]] <- data.frame(Percent = ifelse(use == "sg1",
                                                       sel_sg1[sel_sg1[,1] == peak_indel,3],sel_sg2[sel_sg2[,1] == peak_indel,3]), 
                                      ID = paste(ifelse(use == "sg1", ids[1], ids[2]), 
                                                 "of", paste0(ids, collapse = "-"), sep = " "), 
                                      peak = peak_indel, used = use,
                                      used_sg = ifelse(use == "sg1", ids[1], 
                                                       ids[2]))
}
not_prefer_select <- data.frame(do.call(rbind, not_prefer_select))

####get sequence after process
sgCmp_id <- merge(sgCmp, sgRNA, by="sgRNA")
sgRNA_sel <- sgRNA[sgRNA$id2 %in% not_prefer_select$used_sg,]
doAA_change_seq <- list()
for(i in 1 : nrow(sgRNA_sel)){
  id <- sgRNA_sel$id2[i]
  edit_tables <- total_edit_table_rev[[id]]
  wt_seq <- not_prefer_select$ID[not_prefer_select$used_sg == id]
  wt_seq <- str_remove_all(wt_seq, id)
  wt_seq <- str_remove(str_remove(wt_seq, " of "), "-")
  wt_seq <- sgCmp_id$Cmp[sgCmp_id$id2 == wt_seq][1]
  sg <- sgCmp_id$Cmp[sgCmp_id$id2 == id][1]
  if(not_prefer_select$peak[not_prefer_select$used_sg == id] < 0){
    edit_table <- lapply(edit_tables, function(edit_table){
      edit_table <- edit_table[edit_table$Unedited == "False",]
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
      edit_table <- edit_table[edit_table$n_deleted > 0 & 
                                 edit_table$n_deleted <= 5 & 
                                 edit_table$n_inserted == 0 & 
                                 edit_table$n_mutated <= 1,]
      edit_table$sgCmp <- sg
      edit_table$wtCmp <- wt_seq
      edit_table
    })
  }
  else{
    edit_table <- lapply(edit_tables, function(edit_table){
      edit_table <- edit_table[edit_table$Unedited == "False",]
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
      edit_table <- edit_table[edit_table$n_deleted == 0 & 
                                 edit_table$n_inserted > 0&
                                 edit_table$n_inserted <= 5&
                                 edit_table$n_mutated <= 1,]
      edit_table$sgCmp <- sg
      edit_table$wtCmp <- wt_seq
      edit_table
    })
  }
  doAA_change_seq[[id]] <- edit_table
  
}
seqs_table2 <- list()
bad_pct <- 0
good_pct <- 0
for(id in names(doAA_change_seq)){
  edit_table <- doAA_change_seq[[id]]
  peak <- not_prefer_select$peak[not_prefer_select$used_sg == id]
  close_edit <- list()
  for(name in names(edit_table)){
    tmp <- edit_table[[name]]
    check_seq <- tmp[,2]
    if(peak < 0){
      check_seq <- tmp[,1]
    }
    poss <- unlist(lapply(check_seq, function(x){
      pos <- str_locate_all(x, "-")[[1]][,1]
      pos[which.min(abs(pos - 20))]
    }))
    close_edit[[name]] <- tmp[abs(poss - 20) <= 3,]
  }
  edit_table <- data.frame(do.call(rbind, close_edit))
  edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, edit_table$Aligned_Sequence), function(x){
    x[,7] <- mean(x[,7])
    x[,8] <- mean(x[,8])
    x[,9] <- mean(x[,9])
    x[1,]
  })))
  edit_table <- edit_table[!duplicated(edit_table$Aligned_Sequence),]
  edit_table$id <- id
  seqs_table2[[id]] <- edit_table[,c(1,2, 10, 11, 12, 9)]
}
seqs_table_out2 <- lapply(seqs_table2, function(x){
  x <- x[order(x[,6], decreasing = T),]
  x[1:3,]
})
seqs_table2 <- data.frame(do.call(rbind, seqs_table2))
seqs_table_out2 <- data.frame(do.call(rbind, seqs_table_out2))
##在4.0脚本里面也输出了一个seqs_table_out_less  但是注意得改一下 默认用的是indel1

ToNX::write_tb(seqs_table_out_less, file="~/data/project/ear_project/gene_therapy_ll/indel_in_indel_pairOne.txt")
ToNX::write_tb(seqs_table_out2, file="indel_in_indel.txt")



setwd("~/data/project/ear_project/gene_therapy_ll/")
result <- read.table("indel_in_indel_res.txt")
colnames(result) <- c("Edited", "Reference", "sgRNA", "paired_sgRNA", "name", "Percent", "Bulge_Position", "Type", "ChangedNT", "Position(to Cut Site)")
result <- merge(result, not_prefer_select[,c("peak","used_sg")], by.x="name", by.y = "used_sg")

pair_result <- read.table("indel_in_indel_pairOne_res.txt")
colnames(pair_result) <- c("Edited ", "Reference ", "sgRNA ", "paired_sgRNA ", "name ", "Percent ", "Bulge_Position ", "Type ", "ChangedNT ", "Position(to Cut Site) ")

total_result <- list()
for(id in high_cor_sample_sel$pair){
  ids <- unlist(strsplit(id, "[-]"))
  res1 <- result[result$name %in% ids,]
  res2 <- pair_result[pair_result$name %in% ids,]
  res <- data.frame(cbind(res1, res2))
  res$pair <- id
  res <- res[,c(ncol(res), 1 : (ncol(res) - 1))]
  total_result[[id]] <- res
}
total_result <- data.frame(do.call(rbind, total_result))






openxlsx::write.xlsx(list(insert = total_result[total_result$peak > 0,], 
                          del = total_result[total_result$peak < 0,]), file="indel_in_indel_res2.xlsx", rowNames=F, colNames=T)







####把高相关性的结果具体的编辑结果拷贝到一个文件夹去
####暂时不用了  先放着
high_cor_sample_sel <- high_cor_sample[unlist(lapply(high_cor_sample$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))
  sum(x %in% result$name) > 0
})),]
high_cor_sample_sel$pos <- unlist(lapply(high_cor_sample_sel$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))
  unique(result$Bulge_Position[result$name %in% x])
}))
openxlsx::write.xlsx(high_cor_sample_sel, file="~/data/project/ear_project/gene_therapy_ll/high_cor/sample_info.xlsx", rowNames=F, colNames=T)
lapply(high_cor_sample$pair, function(x){
  x <- unlist(strsplit(x, "[-]"))
  lapply(x, function(y){
    yy <- unlist(strsplit(y, "[_]"))
    if(as.numeric(yy[3]) < 29){
      setwd("~/data/project/ear_project/gene_therapy_ll/batch1/")
      samples <- list.files(pattern = "Sg")
      for(i in 1 : length(samples)){
        setwd(samples[i])
        if(file.exists(paste0("sg", yy[3], "_split"))){
          file.copy(paste0("sg", yy[3], "_split/CRISPResso_on_nhej"), 
                    paste0("~/data/project/ear_project/gene_therapy_ll/high_cor/Rep", i, "/"),
                    recursive = T)
          file.rename(paste0("~/data/project/ear_project/gene_therapy_ll/high_cor/Rep", i, "/CRISPResso_on_nhej"), 
                      paste0("~/data/project/ear_project/gene_therapy_ll/high_cor/Rep", i, "/", y))
        }
        setwd("..")
      }
    }else{
      setwd("~/data/project/ear_project/gene_therapy_ll/batch2/total_output/")
      samples <- c("A", "B", "C")
      for(i in 1 : length(samples)){
        setwd(samples[i])
        if(file.exists(paste0("CRISPResso_on_", y))){
          file.copy(paste0("CRISPResso_on_", y), 
                    paste0("~/data/project/ear_project/gene_therapy_ll/high_cor/Rep", i, "/"),
                    recursive = T)
          file.rename(paste0("~/data/project/ear_project/gene_therapy_ll/high_cor/Rep", i, "/CRISPResso_on_", y), 
                      paste0("~/data/project/ear_project/gene_therapy_ll/high_cor/Rep", i, "/", y))
        }
        setwd("..")
      }
    }
  })
})

