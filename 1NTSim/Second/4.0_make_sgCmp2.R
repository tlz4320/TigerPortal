setwd("~/data/project/ear_project/gene_therapy/find1NTSim/")
tmp <- read.table("human_sg_final_anno2_addid.txt", sep="\t")
tmp2 <- read.table("threshold30000/human_sg_final_anno_threshold30000.txt", sep="\t")
tmp2$V16 <- lapply(1 : nrow(tmp2), function(x){
  return(paste0(tmp2$V8[x], " ", tmp2$V15[x]))
})
tmp3 <- read.table("human_sg_lastNTSim_final_addid2.txt", sep = "\t")
tmp3$V14 <- lapply(1 : nrow(tmp3), function(x){
  if(str_sub(tmp3[x,2], 2, 20) == str_sub(tmp3[x, 8], 1, 19)){
    return(paste0(paste0(tmp3[x,2], "-"), " ", paste0("-",tmp3[x,8])))
  }
  return(paste0(paste0(tmp3[x,8], "-"), " ", paste0("-",tmp3[x,2])))
})
tmp4 <-unlist(c(tmp[,14], tmp2[,16], tmp3[,14]))
tmp <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[1]
}))
tmp2 <- unlist(lapply(tmp4, function(x){
  unlist(strsplit(x, "[ ]"))[2]
}))
tmp4 <- data.frame(sgRNA = str_remove(c(tmp, tmp2), "[-]"), 
                   Cmp = c(tmp, tmp2))
tmp4 <- tmp4[!duplicated(tmp4$sgRNA),]

sgCmp2 <- tmp4
rm(tmp, tmp2, tmp3, tmp4)
sgRNA2 <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Second1NT/second_batch_design_total160_replace_bad_fix.xlsx")
sgCmp2 <- merge(sgCmp2, sgRNA2, by="sgRNA")
sgCmp2$id2 <- paste0("Sg_",str_replace(sgCmp2$ID_new, "-", "_"))
sgCmp2$pos <- unlist(lapply(sgCmp2$ID, function(x){
  unlist(strsplit(x, "[_]"))[2]
}))
sgCmp2$test <- unlist(lapply(1 : nrow(sgCmp2), function(i){
  cmp <- sgCmp2$Cmp[i]
  cmp_pair <- sgCmp2$Cmp[sgCmp2$ID == sgCmp2$ID[i]]
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
sgCmp2$isInsert <- unlist(lapply(1 : nrow(sgCmp2), function(i){
  cmp <- sgCmp2$Cmp[i]
  cmp_pair <- sgCmp2$Cmp[sgCmp2$ID == sgCmp2$ID[i]]
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
sgCmp2 <- sgCmp2[order(sgCmp2$ID_new),]
save(sgCmp2, file="~/data/project/ear_project/gene_therapy_ll/Result/second_sgCmp.rda")


bulge_pos <- lapply(split(sgCmp2, sgCmp2$ID), function(tmp){
  ins_one <- tmp$Cmp[tmp$isInsert == 'Insert']
  del_one <- tmp$Cmp[tmp$isInsert != 'Insert']
  ins_one_seq <- unlist(strsplit(ins_one, "*"))
  del_one_seq <- unlist(strsplit(del_one, "*"))
  ins_nt <- ins_one_seq[del_one_seq == "-"]
  c(ins_one, del_one, ins_nt)

})
bulge_pos <- data.frame(do.call(rbind, bulge_pos))
bulge_pos$id <- rownames(bulge_pos)
ToNX::write_tb(bulge_pos, "~/data/project/ear_project/gene_therapy_ll/Second_sgRNA_cmp.txt")
