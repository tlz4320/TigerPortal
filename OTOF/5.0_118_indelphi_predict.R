setwd("~/data/project/ear_project/gene_therapy_ll/Previews/")
batchs <- list.files(pattern = "Rep[1-3]-")
total_result <- list()
for(batch in batchs){
  setwd(batch)
  files <- list.files()
  for(file in files){
    tmp <- read.table(paste0(file,"/CRISPResso_RUNNING_LOG.txt"), sep="\t", fill = T)
    total_result[[file]] <- unlist(tmp[3,1])
  }
  setwd("..")
}


total_used_seq <- lapply(total_result, function(y){

  y <- unlist(strsplit(y, "[ ]"))
  seq <- y[which(y == "-a") + 1]
  sg <- y[which(y == "-g") + 1]

  data.frame(seq = seq, sg = sg)
})
total_used_seq <- data.frame(do.call(rbind, total_used_seq))
total_used_seq$sample <- names(total_result)
total_used_seq$sample <- str_remove(total_used_seq$sample, "CRISPResso_on_")
total_used_seq$id <- unlist(lapply(total_used_seq$sample, function(x){
  unlist(strsplit(x, "[-]"))[1]
}))
total_used_seq <- total_used_seq[str_ends(total_used_seq$sample, "001") | 
                                   str_ends(total_used_seq$sample, "051"),]
total_used_seq <- total_used_seq[order(total_used_seq$id),]
# table(unlist(lapply(split(total_used_seq$seq, total_used_seq$id), function(x){
#   x[1] == x[2]
# })))
# table(unlist(lapply(split(total_used_seq$sg, total_used_seq$id), function(x){
#   x[1] == x[2]
# })))
total_used_seq <- total_used_seq[duplicated(total_used_seq$id),]

total_used_seq$seq2 <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  seq <- str_to_upper(total_used_seq$seq[i])
  sg <- str_to_upper(total_used_seq$sg[i])
  if(str_detect(seq, sg)){
    return(seq)
  }
  seq <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAString(seq)
  ))
  if(str_detect(seq, sg)){
    return(seq)
  }
  return("??")
} ))
total_used_seq$cutpos <- unlist(lapply(1 : nrow(total_used_seq), function(i){
  seq <- str_to_upper(total_used_seq$seq2[i])
  sg <- str_to_upper(total_used_seq$sg[i])
  res <- str_locate(seq, sg)
  if(nrow(res) != 1){
    print("Something error")
  }
  res[1,1] + str_length(sg) - 4
}))
another_4_seq <- data.frame(id2 = c("m12", "m13", "m14", "m15"), 
                            genomeSeq = "AGTCCCTCTGCTCGGTAAATTTTCACATAGAACCGTGCCCACTGCCGTTCGGGGGGCACGCCCTCGGGGAGCAGCAAGTTCCTATGAGCACAGGCAGCAT",
                            cutpos = c(50, 49, 48, 47))
indelphi_file <- data.frame(id2 = total_used_seq$id, 
                            genomeSeq = total_used_seq$seq2, 
                            cutpos = total_used_seq$cutpos)
indelphi_file <- data.frame(rbind(indelphi_file, another_4_seq))
write.csv(indelphi_file, file='~/118_for_indelphi.csv', row.names = F)

forcast_file <- data.frame(ID = indelphi_file$id2, 
                           Target = indelphi_file$genomeSeq,
                           "PAM Index" = indelphi_file$cutpos+ 3)
ToNX::write_tb(forcast_file, file="~/118_for_forecast.txt", col.names = F)

###后面是读取预测的结果
setwd('~/indelphi_res/')
indelphi_118_res <- list()
for(file in indelphi_file$id2){
  indelphi_118_res[[file]] <- read.table(paste0(file, ".csv"), sep = ",", header = T)
}
indelphi_118_res <- lapply(indelphi_118_res, function(x){
  x$Length[x$Category == "del"] <- -x$Length[x$Category == "del"]
  res <- lapply(split(x$Predicted.frequency, x$Length), sum)
  res <- data.frame(indel_size = names(res), fq = unlist(res))
  res$Pct <- res$fq
  res <- res[gtools::mixedorder(res$indel_size),]
  res
})
save(indelphi_118_res, file="~/data/project/ear_project/gene_therapy_ll/Previews/118_indelphi.rda")

setwd('~/forecast_result/')
forecast_118_res <- list()
for(file in forcast_file$ID){
  forecast_118_res[[file]] <- read.table(paste0(file, "_predictedindelsummary.txt"), 
                                         sep = "\t", header = T)
}
forecast_118_res <- lapply(forecast_118_res, function(x){
  edit <- unlist(lapply(x[,1], function(y){
    unlist(strsplit(y, "[_]"))[1]
  }))
  edit <- data.frame(type = str_sub(edit, 1, 1), len = str_sub(edit, 2), counts = x[,3])
  edit$len <- as.numeric(edit$len)
  edit$len[edit$type == "D"] <- -edit$len[edit$type == "D"]
  res <- lapply(split(edit$counts, edit$len), sum)
  res <- data.frame(indel_size = names(res), fq = unlist(res))
  res$Pct <- res$fq / sum(res$fq) * 100
  res <- res[gtools::mixedorder(res$indel_size),]
  res
})
save(forecast_118_res, file="~/data/project/ear_project/gene_therapy_ll/Previews/118_forecast.rda")
