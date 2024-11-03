mutation_info <- lapply(split(mutation_info, mutation_info$V6), function(x){
  len <- str_length(x[,27])
  x[which.max(len),]
})
mutation_info <- data.frame(do.call(rbind, mutation_info))

total_design <- mclapply(1 : nrow(mutation_info), function(index){
    inputseq <- mutation_info$V7[index]
    mutid <- mutation_info$V6[index]
    checked <- formatCheck(inputseq)
    tmp_res <- list()
    if(class(checked) == "logical"){
      if(!checked){
        shinyalert("Not a valid input", "Not a valid input", type="error")
        return("")
      }
    }
    indel_pos <- checked$pos
    seq <- checked$seq
    isIndel <- checked$isIndel
    indelNT <- checked$indelNT
    possiblePAM <- findPAM(seq, pam = "NGG")
    if(T){
      reverseSeq <- as.character(reverseComplement(DNAString(seq)))
      possiblePAM2 <- findPAM(reverseSeq, pam = "NGG")
    } else {
      #只是为了方便下面判断而已
      possiblePAM2 <- possiblePAM
    }
    if(nrow(possiblePAM) == 0 & nrow(possiblePAM2) == 0){
      return("")
    }
    if(nrow(possiblePAM) != 0){
      possiblePAM$dis <- abs(possiblePAM$start - 4 - indel_pos)
      coloredSeq_f <- colorSeq(seq, possiblePAM, indel_pos, isIndel,indelNT)
      coloredSeq_f$strand <- "+"
    } else {
      coloredSeq_f <- data.frame()
    }
    if(T & nrow(possiblePAM2) != 0){
      possiblePAM2$dis <- abs(possiblePAM2$start - 4 - (str_length(seq) - indel_pos + 1))
      coloredSeq_r <- colorSeq(reverseSeq, possiblePAM2, str_length(reverseSeq) - indel_pos + 1, isIndel,indelNT)
      coloredSeq_r$strand <- "-"
      tmp_res[["reverseSeq"]] <- reverseSeq
    } else {
      coloredSeq_r <- data.frame()
    }
    
    if(nrow(coloredSeq_r) != 0 & nrow(coloredSeq_f) != 0){
      coloredSeq <- data.frame(rbind(coloredSeq_f, coloredSeq_r))
    } else if(nrow(coloredSeq_r) == 0){
      coloredSeq <- coloredSeq_f
    } else {
      coloredSeq <- coloredSeq_r
    }
    if(isIndel != 0){
      coloredSeq <- coloredSeq[order(coloredSeq$dis),]
      coloredSeq <- coloredSeq[coloredSeq$dis < 6,]
    }
    
    if(nrow(possiblePAM) == 0 & nrow(possiblePAM2) == 0){
      return("")
    }
    if(nrow(coloredSeq) == 0){
      return("")
    }
    tmp_res[["infos"]] <- coloredSeq
    tmp_res[["seq"]] <- seq
    tmp_res[["infos2"]] <- checked
    tmp_res[["mut"]] <- mutid
    tmp_res
}, mc.cores = 10)
should_keep <- unlist(lapply(total_design, function(x){
  class(x) != "character"
}))
total_design <- total_design[should_keep]
# save(total_design, file="~/data/project/ear_project/gene_therapy_ll/tmp_total_design.rda")
total_for_indelphi <- mclapply(total_design, function(tmp){
  sgInfos <- tmp$infos
  mutid <- tmp$mut
  
  sel <- rep(T, nrow(sgInfos))
  sgInfos$id <- paste0("id_", 1 : nrow(sgInfos))
  sgInfos_sel <- sgInfos[sel,]
  for_indelphi <- lapply(1 : nrow(sgInfos_sel), function(i){
    
    strand <- sgInfos_sel$strand[i]
    start <- sgInfos_sel$start[i]
    if(strand == "+"){
      seq <- tmp$seq
    } else {
      seq <- tmp$reverseSeq
    }
    data.frame(id = paste0(mutid, "-", sgInfos_sel$id[i]), 
               seq = seq, pos = start - 4)
    
  })
  data.frame(do.call(rbind, for_indelphi))
}, mc.cores = 10)
total_for_indelphi <- data.frame(do.call(rbind, total_for_indelphi))
steps <- seq(0, nrow(total_for_indelphi), floor(nrow(total_for_indelphi) / 20))
steps[length(steps)] <- nrow(total_for_indelphi)
for(i in 2 : length(steps)){
  write.table(total_for_indelphi[(steps[i - 1] + 1) : steps[i],], 
              paste0("~/data/all_mutation_predict/indelphi_table", i, ".csv"), col.names = T, row.names = F,quote = T,sep = ",")
  
}


sgs_list <- list()
infos_list <- list()
region_list <- list()
for(i in 1 : length(total_design)){
  tmp <- total_design[[i]]
  sgs_list[[tmp$mut]] <- tmp$infos
  infos_list[[tmp$mut]] <- tmp$infos2
  region_list[[tmp$mut]] <- unlist(mutation_info[mutation_info[,6] == tmp$mut,c(8,9)])
}





