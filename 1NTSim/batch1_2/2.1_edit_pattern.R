library(stringr)
setwd("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
total_insert_stat <- list()
sgRNA_tmp <- read.table("sg_info.txt")
for(sample in samples){
  setwd(sample)
  total_insert_stat[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    if(sample == "CRISPRessoPooled_on_B29" & sg %in% c("Sg_12_79", 
                                                        "Sg_12_80", 
                                                        "Sg_12_81")){
      print("rm sample")
      next
    }
    else{
      setwd(paste0("CRISPResso_on_",sg))
      edit_info <- read.table("Indel_histogram.txt", 
                              sep="\t", header = T)
      total_insert_stat[[sample]][[sg]] <- edit_info
      setwd("..")
    }

  }
  setwd("..")
  rm(edit_info)
}
rm(sample, samples, sg)
batch2_1 <- total_insert_stat
batch2_mean <- lapply(names(total_insert_stat[[1]]), function(x){
  max_del <- min(unlist(lapply(total_insert_stat, function(y){
    if(!x %in% names(y)){
      return(0)
    }
    sel_sg <- y[[x]]
    min(sel_sg[,1])
  })))
  max_in <- max(unlist(lapply(total_insert_stat, function(y){
    if(!x %in% names(y)){
      return(0)
    }
    sel_sg <- y[[x]]
    max(sel_sg[,1])
  })))
  res <- lapply(total_insert_stat, function(y){
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
names(batch2_mean) <- names(total_insert_stat[[1]])
batch2_mean <- lapply(batch2_mean, function(x){
  if(length(x) == 1){
    return(x[[1]])
  }
  fq <- unlist(rowMeans(data.frame(do.call(cbind, lapply(x, function(y){y[,2]})))))
  x[[1]][,2] <- fq
  x[[1]]
})
rm(total_insert_stat)
# save(batch2_1, batch2_mean, file="indel_stat.rda")
batch2_each <- list()
col_split <- c()
for(name in names(batch2_1[[1]])){
  for(name2 in names(batch2_1)){
    if(name %in% names(batch2_1[[name2]])){
      n <- paste(name, name2, sep="_")
      batch2_each[[n]] <- batch2_1[[name2]][[name]]
      col_split <- c(col_split, name)
    }
  }
}
rm(name, name2, n)



pdf("../Result/Batch2_edit_mean2.pdf", width = 20, height = 8)
plotStat(batch2_mean, region = -18 : 5, title = "Batch2", showname = T,
         shownumber = T)
dev.off()
pdf("../Result/Batch2_edit_each2.pdf", width = 20, height = 8)
plotStat(batch2_each, region = -18 : 5, title = "Batch2", showname = F,
         shownumber = T, split_col = col_split)
dev.off()


plot_data <- lapply(names(batch2_1[[1]]), function(x){
  region <- c(-19: 19)
  
  res <- lapply(batch2_1, function(y){
    sel_sg <- y[[x]]
    if(is.null(sel_sg))
      return(sel_sg)
    long_del <- sum(sel_sg[sel_sg[,1] < min(region), 2])
    long_ins <- sum(sel_sg[sel_sg[,1] > max(region), 2])
    long_res <- data.frame(indel_size = c("<=-20", ">=20"), fq = c(long_del, long_ins))
    tmp_region <- region[!region %in% sel_sg[,1]]
    sel_sg <- sel_sg[sel_sg[,1] %in% region,]
    if(length(tmp_region) != 0){
      sel_sg <- data.frame(rbind(sel_sg, data.frame(indel_size = tmp_region, fq = 0)))
      sel_sg <- sel_sg[order(sel_sg[,1]),]
    }
    sel_sg <- data.frame(rbind(sel_sg, long_res))
    sel_sg[sel_sg[,1] == 0, 2] <- 0 
    sel_sg
  })
  res[!unlist(lapply(res, is.null))]
})
names(plot_data) <- names(batch2_1[[1]])




pdf("../Result/Edit_pattern_each_sample2.pdf", width = 10, height = 6)
for(name in names(plot_data)){
  tmp_data <- plot_data[[name]]
  max_pct <- max(unlist(lapply(tmp_data, function(y){
      max(y[,2] / sum(y[,2]))
  })))
  for(sample in names(tmp_data)){
    tmp_data[[sample]]$sample <- str_remove(sample, "CRISPRessoPooled_on_")
    tmp_data[[sample]]$pct <- tmp_data[[sample]][, 2] / sum(tmp_data[[sample]][, 2])
  }
  x <- data.frame(do.call(rbind, tmp_data))
  x$indel_size <- factor(x$indel_size, levels = c("<=-20", -19 : 19, ">=20"))
  print(ggplot(x, aes(x = indel_size, y = pct * 100, fill = sample)) + 
    geom_bar(stat="identity", position = "dodge") + 
    scale_y_continuous(limits = c(0, max_pct * 100)) + 
    ylab("% in indel reads") + xlab("") + 
    theme_bw() + ggtitle(name))
}
dev.off()
