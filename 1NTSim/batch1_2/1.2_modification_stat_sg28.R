setwd("~/data/project/ear_project/gene_therapy_ll/batch1")
samples <- c("./Sg1-28-1-LFM17947", "./Sg1-28-2-LFM17948", "./Sg1-28-3-LFM17949")
total_sample_stat <- list()
for(sample in samples){
  setwd(sample)
  total_sample_stat[[sample]] <- list()
  sgs <- list.files(pattern = "^sg[0-9]+_split")
  for(sg in sgs){
    setwd(sg)
    setwd("CRISPResso_on_nhej")
    mapping_rate <- read.table("CRISPResso_quantification_of_editing_frequency.txt", 
                               sep="\t", header = T)
    total_sample_stat[[sample]][[sg]] <- c(unlist(mapping_rate[1,c(7:12)]))
    setwd("../..")
  }
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples, sg)
total_sample_stat <- lapply(total_sample_stat, function(x){
  x <- data.frame(do.call(rbind, x))
  x$sg <- rownames(x)
  colnames(x) <- c("Unmodified", "Modified", "Discarded", "Insertions","Deletions", "Substitutions", "ID" )
  x
})
total_sample_stat <- lapply(total_sample_stat, function(x){
  x <- x[,c(-3)]
})

#output Excel
output_res<- data.frame(do.call(rbind, total_sample_stat))
output_res$Sample <- rownames(output_res)
output_res$Sample <- str_remove(str_remove(output_res$Sample, "./"),"[.].*")
output_res <- output_res[gtools::mixedorder(output_res$ID),]
output_res$Sample <- str_remove(output_res$Sample, "CRISPRessoPooled_on_")
output_res <- output_res[,c(6,7,1:5)]

total_sample_stat_pct <- lapply(total_sample_stat, function(x){
  total_counts <- unlist(rowSums(x[,-ncol(x)]))
  for(i in 1 : 5){
    x[,i] <- x[, i] / total_counts
  }
  x
})
output_res2<- data.frame(do.call(rbind, total_sample_stat_pct))
output_res2$Sample <- rownames(output_res2)
output_res2$Sample <- str_remove(str_remove(output_res2$Sample, "./"),"[.].*")
output_res2 <- output_res2[gtools::mixedorder(output_res2$ID),]
output_res2$Sample <- str_remove(output_res2$Sample, "CRISPRessoPooled_on_")
output_res2 <- output_res2[,c(6,7,1:5)]
colnames(output_res2) <- paste0(colnames(output_res2), "(Pct)")
output_res <- data.frame(cbind(output_res, output_res2[,c(-1,-2)]))
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/First1NT/2024_1_12_integrated_design_result.xlsx")
sgRNA$id2 <- paste("Sg", sgRNA$new_id, sep = "-")
sgRNA$id2 <- str_replace_all(sgRNA$id2, "-", "_")
output_res$ID <- str_remove(output_res$ID, "_split")
for(i in 1 : length(output_res$ID)){
  n <- as.numeric(str_remove(output_res$ID[i], "sg"))
  output_res$ID[i] <- sgRNA$id2[n]
}
openxlsx::write.xlsx(list(Res = output_res), "Result/sg28_edit_stat1.xlsx",
                     colNames=T, rowNames=F)
