load("~/data/project/ear_project/gene_therapy_ll/Result/sgCmp132.rda")
total_samples <- list()
for(batch in names(total_indel_table_first_second)){
  tmp <- total_indel_table_first_second[[batch]]
  total_samples[[batch]] <- list() 
  for(sg in names(tmp)){
    if(sg %in% sgCmp132$id2){
      total_samples[[batch]][[sg]] <- sg
    }
    
  }
  total_samples[[batch]] <- unlist(total_samples[[batch]])
}
sgCmp <- sgCmp[gtools::mixedorder(sgCmp$id2),]
total_fastq <- read.table("~/data/project/ear_project/gene_therapy_ll/total_fastq.txt")
total_fastq <- total_fastq[-grep("noBC", total_fastq$V1),]
id <- unlist(lapply(total_fastq, function(x){
  if(str_starts(x, "batch1")){
    x <- unlist(strsplit(x, "[_]"))
    x <- x[length(x)]
    x <- str_remove(str_remove(x, "Sg"), ".fq")
    x <- as.numeric(x)
    return(sgCmp$id2[x])
  }
  x <- unlist(strsplit(x, "[/]"))
  x <- x[length(x)]
  x <- str_remove(str_remove(x, "AMPL_"), ".fastq.gz")
  return(x)
}))
total_fastq <- data.frame(fastq = total_fastq, id = id)
total_fastq <- total_fastq[total_fastq$id %in% sgCmp132$id2,]
total_fastq <- total_fastq[-grep("B1-B4",total_fastq$fastq),]
tmp <- data.frame(table(total_fastq$id))
tmp <- tmp[order(tmp$Freq, decreasing = T),]
for(i in tmp$Var1[tmp$Freq == 6]){
  if(i == "Sg_7_49"){
    rmid <- intersect(which(total_fastq$id == i), grep("240131-A00599", total_fastq$fastq))
    total_fastq <- total_fastq[-rmid,]
  }else{
    #其中179-196用的是batch1的 161到178用的是Batch2的 所以还得分开
    sel_179_196 <- sgCmp[179 : 196,]
    if(i %in% sel_179_196$id2){
      rmid <- intersect(which(total_fastq$id == i), grep("Batch2", total_fastq$fastq))
    } else{
      rmid <- intersect(which(total_fastq$id == i), grep("Batch1", total_fastq$fastq))
      
    }

    total_fastq <- total_fastq[-rmid,]
  }
}
tmp <- data.frame(table(total_fastq$id))
tmp <- tmp[order(tmp$Freq, decreasing = T),]
table(sgCmp132$id2 %in% total_fastq$id)
total_fastq$batch <- unlist(lapply(total_fastq$fastq, function(x){
  if(!is.na(unlist(str_match(x, "Sg1-28-1")))){
    return("Rep1")
  }
  if(!is.na(unlist(str_match(x, "Sg1-28-2")))){
    return("Rep2")
  }
  if(!is.na(unlist(str_match(x, "Sg1-28-3")))){
    return("Rep3")
  }
  if(!is.na(unlist(str_match(x, "CRISPRessoPooled_on_A")))){
    return("Rep1")
  }
  if(!is.na(unlist(str_match(x, "CRISPRessoPooled_on_B")))){
    return("Rep2")
  }
  if(!is.na(unlist(str_match(x, "CRISPRessoPooled_on_C")))){
    return("Rep3")
  }
})) 
  
setwd("~/data/project/ear_project/gene_therapy_ll/")
names(total_samples) <- c("Rep1", "Rep2", "Rep3")
for(batch in names(total_samples)){
  sgs <- total_samples[[batch]]
  for(sg in sgs){
    path <- total_fastq$fastq[total_fastq$id == sg & total_fastq$batch == batch]
    if(length(path) != 1){
      print(sg)
    }
    #因为前28是没有压缩的 重新压缩一下
    if(sg %in% sgCmp$id2[1:28]){
      file.copy(path, paste0("fastq/", batch, "/", sg, ".fastq"), overwrite = T)
    }else{
      file.copy(path, paste0("fastq/", batch, "/", sg, ".fastq.gz"), overwrite = T)
    }
    
  }
  
  
  
}


old_stat <- list()
setwd("Previews/Rep123_used_file/renamed/")
for(file in c("Rep1", "Rep2", "Rep3")){
  old_stat[[file]] <- read.table(paste0(file, "/", file,"_depth.txt"))
}
old_stat <- unlist(old_stat)

bulge_stat <- list()
setwd("~/data/project/ear_project/gene_therapy_ll/fastq/")
for(file in c("Rep1", "Rep2", "Rep3")){
  bulge_stat[[file]] <- read.table(paste0(file, "/", file,"_depth.txt"))
}
bulge_stat <- unlist(bulge_stat)
all_file_stat <- c(old_stat, bulge_stat)
mean(all_file_stat) / 4
mean(old_stat)/4
