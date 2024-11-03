setwd("~/data/project/ear_project/gene_therapy_ll/Second/Batch1/240318-A00599B/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
second_batch1_indel <- list()
print(load("../../second_sg.rda"))
second_sg <- second_sg[second_sg$ID_new %in% c(paste("24", 179:196, sep = "-")),]
second_sg <- second_sg[!second_sg$ID_new %in% c("24-193", "24-194"),]
for(sample in samples){
  setwd(sample)
  second_batch1_indel[[sample]] <- list()
  for(sg in second_sg$ID2){
    
    setwd(paste0("CRISPResso_on_",sg))
    edit_info <- read.table("Indel_histogram.txt", 
                            sep="\t", header = T)
    second_batch1_indel[[sample]][[sg]] <- edit_info
    setwd("..")
    
    
  }
  setwd("..")
  rm(edit_info)
}
rm(sample, samples, sg)

second_batch1_indel_mean <- lapply(names(second_batch1_indel[[1]]), function(x){
  max_del <- min(unlist(lapply(second_batch1_indel, function(y){
    if(!x %in% names(y)){
      return(0)
    }
    sel_sg <- y[[x]]
    min(sel_sg[,1])
  })))
  max_in <- max(unlist(lapply(second_batch1_indel, function(y){
    if(!x %in% names(y)){
      return(0)
    }
    sel_sg <- y[[x]]
    max(sel_sg[,1])
  })))
  res <- lapply(second_batch1_indel, function(y){
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
names(second_batch1_indel_mean) <- names(second_batch1_indel[[1]])
second_batch1_indel_mean <- lapply(second_batch1_indel_mean, function(x){
  if(length(x) == 1){
    return(x[[1]])
  }
  fq <- unlist(rowMeans(data.frame(do.call(cbind, lapply(x, function(y){y[,2]})))))
  x[[1]][,2] <- fq
  x[[1]]
})


second_batch1_edit_table <- list()

samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
print(load("../../second_sg.rda"))
second_sg <- second_sg[second_sg$ID_new %in% c(paste("24", 179:196, sep = "-")),]
second_sg <- second_sg[!second_sg$ID_new %in% c("24-193", "24-194"),]
for(sample in samples){
  setwd(sample)
  second_batch1_edit_table[[sample]] <- list()
  for(sg in second_sg$ID2){
    setwd(paste0("CRISPResso_on_",sg))
    filename <- list.files(pattern = "Alleles_frequency_table_around.*txt$")
    filename <- filename[which.max(str_length(filename))]
    edit_table <- read.table(filename, 
                             sep="\t", header = T, comment.char = "")
    second_batch1_edit_table[[sample]][[sg]] <- edit_table
    setwd("..")
  }
  setwd("..")
}


save(second_batch1_indel, second_batch1_indel_mean,second_batch1_edit_table, file="../indel_stat.rda")
