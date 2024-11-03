setwd('~/indelphi_res2/')
indelphi_132_res <- list()
indelphi_file <- read.csv("~/for_indelphi.csv")
for(file in indelphi_file$id2){
  indelphi_132_res[[file]] <- read.table(paste0(file, ".csv"), sep = ",", header = T)
}
indelphi_132_res <- lapply(indelphi_132_res, function(x){
  x$Length[x$Category == "del"] <- -x$Length[x$Category == "del"]
  res <- lapply(split(x$Predicted.frequency, x$Length), sum)
  res <- data.frame(indel_size = names(res), fq = unlist(res))
  res$Pct <- res$fq
  res <- res[gtools::mixedorder(res$indel_size),]
  res
})
save(indelphi_132_res, file="~/data/project/ear_project/gene_therapy_ll/Result/indelphi_132_result.rda")


setwd('~/forecast132_result/')
forecast_132_res <- list()
forcast_file <- read.table("~/132_for_forecast.txt")
for(file in forcast_file$V1){
  forecast_132_res[[file]] <- read.table(paste0(file, "_predictedindelsummary.txt"), 
                                         sep = "\t", header = T)
}
forecast_132_res <- lapply(forecast_132_res, function(x){
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
save(forecast_132_res, file="~/data/project/ear_project/gene_therapy_ll/Result/132_forecast.rda")
