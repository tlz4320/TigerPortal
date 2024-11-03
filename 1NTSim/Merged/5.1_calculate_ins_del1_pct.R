load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_ins_table_data.rda")
      
ins1_pct <- lapply(select_edit_pattern, function(x){
  lapply(x, function(y){
    sum(y$Pct)
  })
  })
load("~/data/project/ear_project/gene_therapy_ll/Result/sgCmp132.rda")
ins1_pct <- ins1_pct[sgCmp132$id2]
ins1_pct <- lapply(ins1_pct, function(x){
  name <- names(x)
  name <- lapply(name, function(name){
    if(name == "Sg1-28-1-LFM17947"){
      name <- "A"
    }
    if(name == "Sg1-28-2-LFM17948" ){
      name <- "B"
    }
    if(name == "Sg1-28-3-LFM17949"  ){
      name <- "C"
    }
    name <- str_remove(name, "CRISPRessoPooled_on_")
    name <- str_remove(name, "[0-9]+")
    name
  })
  res <- data.frame(A = NA, B = NA, C = NA)
  name <- unlist(name)
  for(i in 1 : length(name)){
    n <- name[i]
    res[1, n] <- x[[i]]
  }
  res
})
for(name in names(ins1_pct)){
  ins1_pct[[name]]$sg <- name
}
ins1_pct <- data.frame(do.call(rbind, ins1_pct))
ins1_pct$Mean <- rowMeans(ins1_pct[,-4], na.rm = T)


openxlsx::write.xlsx(ins1_pct,"~/data/project/ear_project/gene_therapy_ll/Result/132_ins1_pct_in_indel.xlsx", colNames=T, rowNames=F)




load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_del_table_data.rda")

del1_pct <- lapply(select_edit_pattern, function(x){
  lapply(x, function(y){
    sum(y$Pct)
  })
})
load("~/data/project/ear_project/gene_therapy_ll/Result/sgCmp132.rda")
del1_pct <- del1_pct[sgCmp132$id2]
del1_pct <- lapply(del1_pct, function(x){
  name <- names(x)
  name <- lapply(name, function(name){
    if(name == "Sg1-28-1-LFM17947"){
      name <- "A"
    }
    if(name == "Sg1-28-2-LFM17948" ){
      name <- "B"
    }
    if(name == "Sg1-28-3-LFM17949"  ){
      name <- "C"
    }
    name <- str_remove(name, "CRISPRessoPooled_on_")
    name <- str_remove(name, "[0-9]+")
    name
  })
  res <- data.frame(A = NA, B = NA, C = NA)
  name <- unlist(name)
  for(i in 1 : length(name)){
    n <- name[i]
    res[1, n] <- x[[i]]
  }
  res
})
for(name in names(del1_pct)){
  del1_pct[[name]]$sg <- name
}
del1_pct <- data.frame(do.call(rbind, del1_pct))
del1_pct$Mean <- rowMeans(del1_pct[,-4], na.rm = T)


openxlsx::write.xlsx(del1_pct,"~/data/project/ear_project/gene_therapy_ll/Result/132_del1_pct_in_indel.xlsx", colNames=T, rowNames=F)
