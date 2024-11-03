setwd("~/data/all_mutation_predict/")
csv_file <- list.files(pattern = "csv")
total_result <- list()
for(file in csv_file){
  indelphi_table <- read.csv(file)
  setwd(str_remove(file, ".csv"))
  indelphi_table$mut <- str_remove(indelphi_table$id, "[-].*")
  mutids <- unique(indelphi_table$mut)
  tmp_total_result <- mclapply(mutids, function(mutid){
    sgs <- sgs_list[[mutid]]
    info <- infos_list[[mutid]]
    cds_region <- region_list[[mutid]]
    ids <- indelphi_table$id[indelphi_table$mut == mutid]
    sgs$id <- paste0(mutid, "-id_", 1 : nrow(sgs))
    sgs <- sgs[sgs$id %in% ids,]
    indelphi_res <- list()
    for(file in ids){
      indelphi_res[[file]] <- read.table(paste0(file, ".csv"), sep = ",", header = T)
    }
    result_list <- list()
    wt <- info$WT
    rev_wt <- as.character(reverseComplement(DNAString(wt)))
    isIndel <- info$isIndel
    indelLen <- str_length(info$indelNT)
    ### Inframe and 0AA stat
    inframe <- lapply(ids, function(id){
      tables <- indelphi_res[[id]]
      indel <- ifelse(isIndel == 1, indelLen, -indelLen)
      tables$indel <- ifelse(tables$Category == "del", -tables$Length, tables$Length)
      tables$inframe <- tables$indel + indel
      res <- tables[tables$inframe %% 3 == 0, ]
      res <- res[res$Predicted.frequency != 0,]
      res$aa <- res$inframe / 3
      res$aa[abs(res$aa) > 2] <- "Other"
      res <- lapply(split(res$Predicted.frequency, res$aa), sum)
      res <- data.frame(aa = names(res), Pct = unlist(res), id = id)
      res$Pct2 <- res$Pct / sum(res$Pct) * 100
      res
    })
    inframe <- data.frame(do.call(rbind, inframe))
    result_list[["inframe"]] <- inframe
    sgs$inframe <- unlist(lapply(sgs$id, function(x){
      sum(inframe$Pct[inframe$id == x])
    }))
    aa0 <- unlist(lapply(ids, function(id){
      tables <- indelphi_res[[id]]
      indel <- ifelse(isIndel == 1, indelLen, -indelLen)
      tables$indel <- ifelse(tables$Category == "del", -tables$Length, tables$Length)
      tables$inframe <- tables$indel + indel
      res <- tables[tables$inframe == 0, ]
      res <- res[res$Predicted.frequency != 0,]
      sum(res$Predicted.frequency)
    }))
    sgs$aa0 <- aa0
    wts <- unlist(strsplit(str_to_upper(wt), "*"))
    rev_wts <- unlist(strsplit(str_to_upper(rev_wt), "*"))
    ##calculate NT diff of 0AA
    diff <- unlist(lapply(ids, function(id){
      tables <- indelphi_res[[id]]
      strand <- sgs$strand[sgs$id ==id]
      indel <- ifelse(isIndel == 1, indelLen, -indelLen)
      tables$indel <- ifelse(tables$Category == "del", -tables$Length, tables$Length)
      tables$inframe <- tables$indel + indel
      res <- tables[tables$inframe == 0, ]
      res <- res[res$Predicted.frequency != 0,]
      if(nrow(res) == 0){
        return(-1)
      }
      res$diff <- -1
      if(strand == "+"){
        wt_seqs = wts
      } else {
        wt_seqs = rev_wts
      }
      for(i in 1 : nrow(res)){
        res$diff[i] <- 
          sum(wt_seqs != unlist(strsplit(str_to_upper(res$Genotype[i]), "*")))
      }
      sum(res$diff * res$Predicted.frequency / sum(res$Predicted.frequency))
    }))
    sgs$diff <- diff
    sgs$aa0_inframe <- sgs$aa0 / sgs$inframe
    wtLen <- length(wts)
    sgs$cut_mut_dis <- ifelse(sgs$strand == "+", info$pos - sgs$start + 4, 
                              wtLen - info$pos - sgs$start + 5)
    proteins <- seqinr::translate(seqinr::s2c(str_sub(wt, cds_region[1] + 1, -cds_region[2] - 1)), ambiguous = T)
    wt <- lapply(sgs$id[sgs$aa0 > 0], function(id){
      tables <- indelphi_res[[id]]
      indel <- ifelse(isIndel == 1, indelLen, -indelLen)
      tables$indel <- ifelse(tables$Category == "del", -tables$Length, tables$Length)
      tables$inframe <- tables$indel + indel
      res <- tables[tables$inframe == 0, ]
      res <- res[res$Predicted.frequency != 0,]
      res$protein = ""
      res$proteindiff <- 0
      strand <- sgs$strand[sgs$id ==id]
      for(i in 1 : nrow(res)){
        genotype <- res$Genotype[i]
        if(strand == "-"){
          genotype <- as.character(reverseComplement(DNAString(genotype)))
        }
        genotype <- str_sub(genotype, cds_region[1] + 1, -cds_region[2] - 1)
        tmpp <- seqinr::translate(seqinr::s2c(genotype), ambiguous = T)
        res$protein[i] <- paste0(tmpp, collapse = "")
        res$proteindiff[i] <- sum(proteins != tmpp)
      }
      res$NTdiff <- -1
      
      if(strand == "+"){
        wt_seqs = wts
      } else {
        wt_seqs = rev_wts
      }
      for(i in 1 : nrow(res)){
        res$NTdiff[i] <- 
          sum(wt_seqs != unlist(strsplit(str_to_upper(res$Genotype[i]), "*")))
      }
      res$id <- id
      res
    })
    aa0_result <- data.frame(do.call(rbind, wt))
    result_list[["0AA"]] <- aa0_result
    # protein <- paste0(proteins, collapse = "")
    # inframe_result <- lapply(sgs$id, function(id){
    #   tables <- indelphi_res[[id]]
    #   indel <- ifelse(isIndel == 1, indelLen, -indelLen)
    #   tables$indel <- ifelse(tables$Category == "del", -tables$Length, tables$Length)
    #   tables$inframe <- tables$indel + indel
    #   res <- tables[tables$inframe %% 3 == 0 & abs(tables$inframe / 3) <= 1, ]
    #   res <- res[res$Predicted.frequency != 0,]
    #   res$protein = ""
    #   res$proteindiff <- 0
    #   res$aa <- res$inframe / 3
    #   strand <- sgs$strand[sgs$id ==id]
    #   for(i in 1 : nrow(res)){
    #     genotype <- res$Genotype[i]
    #     if(strand == "-"){
    #       genotype <- as.character(reverseComplement(DNAString(genotype)))
    #     }
    #     genotype <- str_sub(genotype, cds_region[1] + 1, -cds_region[2] - 1)
    #     tmpp <- seqinr::translate(seqinr::s2c(genotype), ambiguous = T)
    #     tmpp <- paste0(tmpp, collapse = "")
    #     res$protein[i] <- tmpp
    #     comp <- Biostrings::pairwiseAlignment(protein, tmpp)
    #     pt <- Biostrings::alignedPattern(comp)
    #     sj <- Biostrings::alignedSubject(comp)
    #     pt <- unlist(strsplit(as.character(pt), "*"))
    #     sj <- unlist(strsplit(as.character(sj), "*"))
    #     res$proteindiff[i] <- sum(pt != sj) - abs(res$aa[i])
    #   }
    #   res$id <- id
    #   res
    # })
    # inframe_result <- data.frame(do.call(rbind, inframe_result))
    # result_list[["1AA"]] <- inframe_result
    result_list[["id"]] <- mutid
    result_list
  }, mc.cores = 20)
  names(tmp_total_result) <- mutids
  for(name in names(tmp_total_result)){
    if(name %in% names(total_result)){
      tmp <- total_result[[name]]
      tmp2 <- tmp_total_result[[name]]
      tmp$inframe <- data.frame(rbind(tmp$inframe, tmp2$inframe))
      tmp[["0AA"]] <- data.frame(rbind(tmp[["0AA"]], tmp2[['0AA']]))
      total_result[[name]] <- tmp
    } else {
      total_result[[name]] <- tmp_total_result[[name]]
    }

  }
  print(paste0(file, " Processed"))
  setwd("..")
}
# save(total_result, file="~/data/all_mutation_predict/total_result.rda")
best_wt_percent <- mclapply(total_result, function(res){
  assign("tmp",res, envir = .GlobalEnv)
  aa0_res <- res[["0AA"]]
  aa0_pct <- res$inframe[res$inframe$aa == 0,]
  top_aa0_res <- lapply(split(aa0_res, aa0_res$id), function(tmp){
    pct <- aa0_pct$Pct2[aa0_pct$id == unique(tmp$id)]
    tmp$Pct <- tmp$Predicted.frequency / sum(tmp$Predicted.frequency) * pct
    tmp2 <- tmp[tmp$proteindiff == 0,]
    if(nrow(tmp2) == 0){
      return(data.frame(pct = 0, pct2 = 0, id = unique(tmp$id)))
    }
    data.frame(pct = sum(tmp2$Predicted.frequency), pct2 = sum(tmp2$Pct), id = unique(tmp$id))
  })
  top_aa0_res <- data.frame(do.call(rbind, top_aa0_res))
  if(all(top_aa0_res$pct2 == 0)){
    return("No")
  }
  return(top_aa0_res[which.max(top_aa0_res$pct2),])
}, mc.cores = 10)
best_wt_percent_sel <- unlist(lapply(best_wt_percent, function(x){
  class(x) != "character"
}))
best_wt_percent_sel <- best_wt_percent[best_wt_percent_sel]
best_wt_percent_sel <- data.frame(do.call(rbind, best_wt_percent_sel))
best_wt_percent_sel$mut <- str_remove(best_wt_percent_sel$id, "[-].*")
mutation_info_sel <- mutation_info[,c(1:6)]
best_wt_percent_sel <- merge(mutation_info_sel, best_wt_percent_sel, by.x="V6", by.y='mut')
best_wt_percent_sel <- best_wt_percent_sel[order(best_wt_percent_sel$pct2, decreasing = T),]
best_wt_percent_sel <- best_wt_percent_sel[!duplicated(best_wt_percent_sel$id),]


clinvar <- read.table("~/data/all_mutation_predict/clinvar_filtered.txt", sep="\t", header = F)
clinvar <- clinvar[clinvar$V4 == "-" | clinvar$V5 == "-",]
clinvar <- clinvar[str_length(clinvar$V4) == 1 & str_length(clinvar$V5) == 1,]
clinvar_sel <- clinvar[grep("pathogenic", str_to_lower(clinvar$V10)), ]
mutation_id <- paste(mutation_info$V1, mutation_info$V2, mutation_info$V4, mutation_info$V5, sep="-")
clinvar_sel$id <- paste(paste0("chr", clinvar_sel[,1]), clinvar_sel[,2], clinvar_sel[,4],clinvar_sel[,5], sep="-")
clinvar_sel <- clinvar_sel[clinvar_sel$id %in% mutation_id,]
ToNX::write_tb(clinvar_sel[,c(1:6)], file="~/data/all_mutation_predict/clinvar_sel_info.vcf")
ToNX::write_tb(clinvar_sel$id, file="~/data/all_mutation_predict/clivar_sel_id.txt")
clinvar_anno <- read.table("~/data/all_mutation_predict/annotation_snp.hg38_multianno.txt", sep="\t", header = T)
clinvar_anno$gnomad312_AF <- as.numeric(clinvar_anno$gnomad312_AF)
clinvar_anno$gnomad312_AF[is.na(clinvar_anno$gnomad312_AF)] <- 0

clinvar_anno_info <- mutation_info[mutation_info$V6 %in% clinvar_anno$avsnp150,]
#去掉接近3' 5'的例子 主要是rs764994176这种

clinvar_anno_info <- clinvar_anno_info[clinvar_anno_info[,8] == 0 & 
                                         clinvar_anno_info[,9] == 0,]
clinvar_anno <- clinvar_anno[clinvar_anno$avsnp150 %in% clinvar_anno_info$V6,]
clinvar_anno <- clinvar_anno[order(clinvar_anno$gnomad312_AF, decreasing = T),]










