###

setwd("~/data/project/ear_project/gene_therapy_ll/Second/Batch1/")
print(load("../second_sg.rda"))
print(load("indel_stat.rda"))

###去掉那些不匹配的样本
second_sg_pair <- split(second_sg$ID2, second_sg$ID)
second_sg_pair_rm <- second_sg_pair[unlist(lapply(second_sg_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(second_batch1_indel_mean)) != 2)
}))]
second_sg_pair_remain <- second_sg_pair[unlist(lapply(second_sg_pair, function(x){
  x <- unlist(x)
  return(sum(x %in% names(second_batch1_indel_mean)) == 2)
}))]
second_batch1_pair_cor <- lapply(second_sg_pair_remain, function(y){
  IDs <- unlist(y)
  sel_sg1 <- second_batch1_indel_mean[[IDs[1]]]
  sel_sg2 <- second_batch1_indel_mean[[IDs[2]]]
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
names(second_batch1_pair_cor) <- unlist(lapply(second_sg_pair_remain, function(x){
  paste(x, collapse = "-")
}))
second_batch1_pair_cor <- data.frame(pair = names(second_batch1_pair_cor), 
                                cor = unlist(lapply(second_batch1_pair_cor, function(x){x$estimate})),
                                pval = unlist(lapply(second_batch1_pair_cor, function(x){x$p.value})))

sg_name <- data.frame(ID = names(second_sg_pair_remain), 
                      pair = unlist(lapply(second_sg_pair_remain, function(x){
                        paste(x, collapse = "-")
                      })))
sg_name$pos <- unlist(lapply(sg_name$ID, function(x){
  if(!is.na(str_match(x, "last"))){
    return(1)
  }
  unlist(strsplit(x, "[_-]"))[2]
}))
second_batch1_pair_cor <- merge(second_batch1_pair_cor, sg_name, by = "pair")



names(second_batch1_indel) <- c("A", "B", "C")
total_sg_cor_each_sample <- list()


second_batch1_pair_cor_sample <- list()

for(name in names(second_batch1_indel)){
  edit_data <- second_batch1_indel[[name]]
  second_sg_pair_remain <- second_sg_pair[unlist(lapply(second_sg_pair, function(x){
    x <- unlist(x)
    return(sum(x %in% names(second_batch1_indel_mean)) == 2)
  }))]
  tmp_cor <- lapply(second_sg_pair_remain, function(y){
    ids <- unlist(y)
    sel_sg1 <- edit_data[[ids[1]]]
    sel_sg2 <- edit_data[[ids[2]]]
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
    cor.test(sel_sg1[,2], sel_sg2[,2])
  })
  names(tmp_cor) <- unlist(lapply(second_sg_pair_remain, function(x){
    paste(x, collapse = "-")
  }))
  second_batch1_pair_cor_sample[[name]] <- data.frame(pair = names(tmp_cor), 
                                                 cor = unlist(lapply(tmp_cor, function(x){x$estimate})),
                                                 pval = unlist(lapply(tmp_cor, function(x){x$p.value})), 
                                                 sample = name)
}



# save(second_batch1_pair_cor, second_batch1_pair_cor_sample, file = "processed_cor_result.rda")














