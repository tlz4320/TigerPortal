print(load("~/data/project/ear_project/gene_therapy_ll/batch1/Result/indel_stat_total.rda"))
names(sg28) <- str_remove(names(sg28), "_.*")
sgs <- paste0("sg", c(1 : 28))

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
print(load("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/indel_stat2.rda"))
print(load("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599//indel_stat.rda"))
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


load("~/data/project/ear_project/gene_therapy_ll/Result/total_indel_table_second.rda")

second_indel_mean <- lapply(names(total_indel_table_second[[1]]), function(x){
  max_del <- min(unlist(lapply(total_indel_table_second, function(y){
    if(!x %in% names(y)){
      return(0)
    }
    sel_sg <- y[[x]]
    min(sel_sg[,1])
  })))
  max_in <- max(unlist(lapply(total_indel_table_second, function(y){
    if(!x %in% names(y)){
      return(0)
    }
    sel_sg <- y[[x]]
    max(sel_sg[,1])
  })))
  res <- lapply(total_indel_table_second, function(y){
    sel_sg <- y[[x]]
    if(is.null(sel_sg))
      return(sel_sg)
    region <- max_del : max_in
    tmp_region <- region[!region %in% sel_sg[,1]]
    if(length(tmp_region) != 0){
      sel_sg <- data.frame(rbind(sel_sg, data.frame(indel_size = tmp_region, fq = 0)))
    }
    sel_sg <- sel_sg[sel_sg[,1] %in% region,]
    sel_sg <- sel_sg[order(sel_sg[,1]),]
    sel_sg[sel_sg[,1] == 0, 2] <- 0 
    sel_sg
  })
  res[!unlist(lapply(res, is.null))]
})
names(second_indel_mean) <- names(total_indel_table_second[[1]])
second_indel_mean <- lapply(second_indel_mean, function(x){
  if(length(x) == 1){
    return(x[[1]])
  }
  fq <- unlist(rowMeans(data.frame(do.call(cbind, lapply(x, function(y){y[,2]})))))
  x[[1]][,2] <- fq
  x[[1]]
})


total_mean_first_second <- total_mean
for(name in names(second_indel_mean)){
  total_mean_first_second[[name]] <- second_indel_mean[[name]]
}

total_indel_table_first_second <- total_indel_table_second
names(total_indel_table_first_second) <- c("A","B","C")
for(batch in 1 : 3){
  sps <- names(batch2_1[[batch]])
  for(sp in sps){
    total_indel_table_first_second[[batch]][[sp]] <- batch2_1[[batch]][[sp]]
  }
  sps <- names(batch2_2[[batch]])
  for(sp in sps){
    total_indel_table_first_second[[batch]][[sp]] <- batch2_2[[batch]][[sp]]
  }
}
for(sp in names(sg28_single_sel)){
  for(batch in names(sg28_single_sel[[sp]])){
    total_indel_table_first_second[[batch]][[sp]] <- sg28_single_sel[[sp]][[batch]]
  }
}

tmp <- unique(unlist(lapply(total_indel_table_first_second, names)))
total_indel_table_first_second_rev <- lapply(tmp, function(x){
  res <- list()
  for(name in names(total_indel_table_first_second)){
    if(x %in% names(total_indel_table_first_second[[name]])){
      res[[name]] <- total_indel_table_first_second[[name]][[x]]
    }
  }
  res
})

names(total_indel_table_first_second_rev) <- tmp
save(total_indel_table_first_second,total_mean_first_second, total_indel_table_first_second_rev, 
     file="~/data/project/ear_project/gene_therapy_ll/Result/merged_indel_table_first_second.Rda")


