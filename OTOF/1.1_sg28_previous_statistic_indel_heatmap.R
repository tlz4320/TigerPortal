library(stringr)
setwd("~/data/project/ear_project/gene_therapy_ll/Previews/old_293T_data/CRISPResso_report/")
samples <- list.files(pattern = "^CRI")
total_insert_stat <- list()
region <- -18 : 18
for(sample in samples){
  setwd(sample)
  total_insert_stat[[sample]] <- list()
  mapping_rate <- read.table("indel_histogram.txt", 
                             sep="\t", header = T)
  total_insert_stat[[sample]] <- mapping_rate
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
old_293t <- total_insert_stat

setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep1-replace-by0829data/")
samples <- list.files(pattern = "^CRI")
total_insert_stat <- list()
for(sample in samples){
  setwd(sample)
  total_insert_stat[[sample]] <- list()
  mapping_rate <- read.table("Indel_histogram.txt", 
                             sep="\t", header = T)
  total_insert_stat[[sample]] <- mapping_rate
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
old_rep1 <- total_insert_stat

setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep2-replace-by0829data/")
samples <- list.files(pattern = "^CRI")
total_insert_stat <- list()
for(sample in samples){
  setwd(sample)
  total_insert_stat[[sample]] <- list()
  mapping_rate <- read.table("Indel_histogram.txt", 
                             sep="\t", header = T)
  total_insert_stat[[sample]] <- mapping_rate
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
old_rep2 <- total_insert_stat

setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep3-replace-by0829data/")
samples <- list.files(pattern = "^CRI")
total_insert_stat <- list()
for(sample in samples){
  setwd(sample)
  total_insert_stat[[sample]] <- list()
  mapping_rate <- read.table("Indel_histogram.txt", 
                             sep="\t", header = T)
  total_insert_stat[[sample]] <- mapping_rate
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples)
old_rep3 <- total_insert_stat

table(names(old_rep3) %in% names(old_rep1))
old_rep1_2 <- list()

for(name in names(old_rep1)){
  tmp_rep1 <- old_rep1[[name]]
  if(!name %in% names(old_rep2)){
    old_rep1_2[[name]] <- tmp_rep1
    next
  }
  tmp_rep2 <- old_rep2[[name]]
  region <- min(c(tmp_rep1[,1], tmp_rep2[,1])) : 
    max(c(tmp_rep1[,1], tmp_rep2[,1]))
  tmp_region <- region[!region %in% tmp_rep1[,1]]
  if(length(tmp_region) != 0){
    tmp_rep1 <- data.frame(rbind(tmp_rep1, data.frame(indel_size = tmp_region, fq = 0)))
  }
  tmp_region <- region[!region %in% tmp_rep2[,1]]
  if(length(tmp_region) != 0){
    tmp_rep2 <- data.frame(rbind(tmp_rep2, data.frame(indel_size = tmp_region, fq = 0)))
  }
  tmp_rep1 <- tmp_rep1[order(tmp_rep1[,1]),]
  tmp_rep2 <- tmp_rep2[order(tmp_rep2[,1]),]
  tmp_rep1[,2] <- (tmp_rep1[,2] + tmp_rep2[,2]) / 2
  old_rep1_2[[name]] <- tmp_rep1
}


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



old_tissue <- old_rep123[str_ends(names(old_rep123),"051")]
old_cell <- old_rep123[str_ends(names(old_rep123),"001")]

setwd("~/data/project/ear_project/gene_therapy_ll")
samples <- list.dirs(recursive = F)
samples <- samples[grep("(S|B1)", samples)]
total_insert_stat <- list()
region <- -18 : 18
for(sample in samples){
  setwd(sample)
  sgs <- list.files(pattern = "^sg")
  total_insert_stat[[sample]] <- list()
  for(sg in sgs){
    setwd(sg)
    mapping_rate <- read.table("CRISPResso_on_nhej/Indel_histogram.txt", 
                               sep="\t", header = T)
    total_insert_stat[[sample]][[sg]] <- mapping_rate
    setwd("..")
  }
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples, sg, sgs)
sg28 <- total_insert_stat

sg28 <- lapply(names(sg28[[1]]), function(x){
  max_del <- min(unlist(lapply(sg28, function(y){
    if(!x %in% names(y)){
      return(0)
    }
    sel_sg <- y[[x]]
    min(sel_sg[,1])
  })))
  
  max_in <- max(unlist(lapply(sg28, function(y){
    if(!x %in% names(y)){
      return(0)
    }
    sel_sg <- y[[x]]
    max(sel_sg[,1])
  })))
  res <- lapply(sg28, function(y){
    sel_sg <- y[[x]]
    region <- max_del : max_in
    tmp_region <- region[!region %in% sel_sg[,1]]
    if(length(tmp_region) != 0){
      sel_sg <- data.frame(rbind(sel_sg, data.frame(indel_size = tmp_region, fq = 0)))
    }
    sel_sg <- sel_sg[sel_sg[,1] %in% region,]
    sel_sg <- sel_sg[order(sel_sg[,1]),]
    sel_sg
  })
  res
})
names(sg28) <- names(total_insert_stat[[1]])

sg28_multi <- lapply(sg28, function(x){
  samples <- str_remove(names(x), "./")
  names(x) <- samples
  multi_samples <- samples[grep("^B", samples)]
  indel_size <- x[[multi_samples[1]]][,1] 
  multi_samples <- data.frame(do.call(rbind,lapply(multi_samples, function(y){
    x[[y]][,2]
  })))
  multi_samples <- unlist(colMeans(multi_samples))
  data.frame(indel_size = indel_size, fq = multi_samples)
})

sg28_single <- lapply(sg28, function(x){
  samples <- str_remove(names(x), "./")
  names(x) <- samples
  multi_samples <- samples[grep("^B", samples)]
  single_samples <- samples[!samples %in% multi_samples]
  indel_size <- x[[single_samples[1]]][,1] 
  single_samples <- data.frame(do.call(rbind,lapply(single_samples, function(y){
    x[[y]][,2]
  })))
  single_samples <- unlist(colMeans(single_samples))
  data.frame(indel_size = indel_size, fq = single_samples)
})

save(old_293t, old_rep1, old_rep1_2, old_tissue, old_cell, old_rep123,
     old_rep2, sg28,sg28_multi, sg28_single, old_rep3,
     file="~/data/project/ear_project/gene_therapy_ll/batch1/Result/indel_stat_total.rda")
save(old_rep1, old_tissue, old_cell, old_rep123,
     old_rep2, old_rep3,
     file="~/data/project/ear_project/gene_therapy_ll/Previews//Result/Rep123_indel_stat.rda")


# names(sg28_multi) <- str_remove(names(sg28_multi), "_.*")
# names(sg28_single) <- str_remove(names(sg28_single), "_.*")
# sgs <- paste0("sg", c(1 : 7, 10, 11, 12, 14 : 17, 19:21,24,25, 27))
# sg28_multi_sel <- sg28_multi[sgs]
# sg28_single_sel <- sg28_single[sgs]
# names(old_tissue) <- str_remove(names(old_tissue), "CRISPResso_on_")
# names(old_cell) <- str_remove(names(old_cell), "CRISPResso_on_")
# pdf("Result/indel_stat_heatmap_sg28.pdf", width = 6, height = 4)
# plotStat(sg28_single_sel, region = -18 : 5, title = "Sg28 Single", showname = T, shownumber = T)
# plotStat(sg28_multi_sel, region = -18 : 5, title = "Sg28 Multi", showname = T, shownumber = T)
# dev.off()
# pdf("Result/indel_stat_heatmap_293T.pdf", width = 15, height = 6)
# plotStat(old_293t, region = -18 : 5, title = "293T", showname = F,
#          shownumber = T)
# dev.off()
# pdf("Result/indel_stat_heatmap_old.pdf", width = 30, height = 12)
# plotStat(old_cell, region = -18 : 5, title = "Cell", showname = T,
#          split_col = split_name_cell$mutation,
#          shownumber = T)
# plotStat(old_tissue, region = -18 : 5, title = "Tissue", showname = T,
#          shownumber = T)
# dev.off()
