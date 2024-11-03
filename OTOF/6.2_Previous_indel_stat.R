print(load("~/data/project/ear_project/gene_therapy_ll/Previews/Result/Rep123_indel_stat.rda"))

indel_pct <- list(Rep1 = old_rep1, Rep2 = old_rep2, Rep3 = old_rep3)

indel_pct <- lapply(indel_pct, function(x){
  x <- x[names(x ) %in% names(old_tissue)]
  res <- data.frame(name = names(old_tissue), Pct = NA)
  rownames(res) <- res$name
  tmp <- unlist(lapply(x, function(x){
    sum(x[x[,1] != 0,2]) / sum(x[,2]) 
  }))
  res[names(x), "Pct"] <- tmp
  res
})
indel_pct_result <- data.frame(name = names(old_tissue), 
                               Rep1 = indel_pct$Rep1$Pct * 100,
                               Rep2 = indel_pct$Rep2$Pct * 100, 
                               Rep3 = indel_pct$Rep3$Pct * 100)
openxlsx::write.xlsx(indel_pct_result, file="~/Nutstore Files/Tobin/Previous/indel_pct_in_tissue.xlsx", rowNames=F, colNames=T)
