###前面是把所有突变都跑了一下 现在取出了top100的结果
###需要对top100进行更加精细的运行  所以接下来是重新跑一遍
library(parallel)
top100_design <- list() 


for(id in clinvar_anno$avsnp150){
  index <- which(mutation_info$V6 == id)
  inputseq <- mutation_info$V7[index]
  mutid <- mutation_info$V6[index]
  checked <- formatCheck(inputseq)
  tmp_res <- list()
  if(class(checked) == "logical"){
    if(!checked){
      shinyalert("Not a valid input", "Not a valid input", type="error")
      next
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
    next
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
    coloredSeq <- coloredSeq[coloredSeq$dis < 20,]
  }
  
  if(nrow(possiblePAM) == 0 & nrow(possiblePAM2) == 0){
    next
  }
  if(nrow(coloredSeq) == 0){
    next
  }
  tmp_res[["infos"]] <- coloredSeq
  tmp_res[["seq"]] <- seq
  tmp_res[["infos2"]] <- checked
  tmp_res[["mut"]] <- mutid
  top100_design[[length(top100_design) + 1]] <- tmp_res
  if(length(top100_design) == 100){
    break
  }
} 
top100_for_indelphi <- mclapply(top100_design, function(tmp){
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
top100_for_indelphi <- data.frame(do.call(rbind, top100_for_indelphi))
# tmp <- read.table("~/data/all_mutation_predict/top100_indelphi.csv", sep=",", header = T)
# table(top100_for_indelphi$seq == tmp$seq)
# write.table(top100_for_indelphi,"~/data/all_mutation_predict/top100_indelphi.csv",
#             col.names = T, row.names = F,quote = T,sep = ",")
  
setwd("~/data/all_mutation_predict/top100_result/")
csv_file <- list.files(pattern = "csv")
top100_result <- list()

sgs_list <- list()
infos_list <- list()
region_list <- list()
for(i in 1 : length(top100_design)){
  tmp <- top100_design[[i]]
  sgs_list[[tmp$mut]] <- tmp$infos
  infos_list[[tmp$mut]] <- tmp$infos2
  region_list[[tmp$mut]] <- unlist(mutation_info[mutation_info[,6] == tmp$mut,c(8,9)])
}


top100_for_indelphi$mut <- str_remove(top100_for_indelphi$id, "[-].*")
mutids <- unique(top100_for_indelphi$mut)
tmp_total_result <- mclapply(mutids, function(mutid){
  sgs <- sgs_list[[mutid]]
  info <- infos_list[[mutid]]
  cds_region <- region_list[[mutid]]
  ids <- top100_for_indelphi$id[top100_for_indelphi$mut == mutid]
  sgs$id <- paste0(mutid, "-id_", 1 : nrow(sgs))
  sgs <- sgs[sgs$id %in% ids,]
  indelphi_res <- list()
  for(file in ids){
    indelphi_edit_table <- read.table(paste0(file, ".csv"), sep = ",", header = T)
    indelphi_res[[file]] <- fix_indel1_func(indelphi_edit_table, sgs$gRNA[sgs$id == file])
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
  protein <- paste0(proteins, collapse = "")
  inframe_result <- lapply(sgs$id, function(id){
    tables <- indelphi_res[[id]]
    indel <- ifelse(isIndel == 1, indelLen, -indelLen)
    tables$indel <- ifelse(tables$Category == "del", -tables$Length, tables$Length)
    tables$inframe <- tables$indel + indel
    res <- tables[tables$inframe %% 3 == 0 & abs(tables$inframe / 3) <= 1, ]
    res <- res[res$Predicted.frequency != 0,]
    res$protein = ""
    res$proteindiff <- 0
    res$aa <- res$inframe / 3
    strand <- sgs$strand[sgs$id ==id]
    for(i in 1 : nrow(res)){
      genotype <- res$Genotype[i]
      if(strand == "-"){
        genotype <- as.character(reverseComplement(DNAString(genotype)))
      }
      genotype <- str_sub(genotype, cds_region[1] + 1, -cds_region[2] - 1)
      tmpp <- seqinr::translate(seqinr::s2c(genotype), ambiguous = T)
      tmpp <- paste0(tmpp, collapse = "")
      res$protein[i] <- tmpp
      comp <- Biostrings::pairwiseAlignment(protein, tmpp)
      pt <- Biostrings::alignedPattern(comp)
      sj <- Biostrings::alignedSubject(comp)
      pt <- unlist(strsplit(as.character(pt), "*"))
      sj <- unlist(strsplit(as.character(sj), "*"))
      res$proteindiff[i] <- sum(pt != sj) - abs(res$aa[i])
    }
    res$id <- id
    res
  })
  inframe_result <- data.frame(do.call(rbind, inframe_result))
  result_list[["1AA"]] <- inframe_result
  result_list[["id"]] <- mutid
  result_list[["sgs"]] <- sgs
  result_list[["indelphi"]] <- indelphi_res
  result_list
}, mc.cores = 20)
names(tmp_total_result) <- mutids
for(name in names(tmp_total_result)){
    top100_result[[name]] <- tmp_total_result[[name]]
}
rm(tmp_total_result)

top100_sgs <- lapply(top100_result, function(x){
  x[['sgs']]
})
top100_sgs <- data.frame(do.call(rbind, top100_sgs))



p <- ggplot(top100_sgs) + geom_point(aes(x = cut_mut_dis, y =diff, color = inframe, 
                                  fill = aa0_inframe), 
                              position = position_jitter(seed = 1), 
                              shape = 21, size = 3, stroke = 1) + 
  scale_color_gradientn(colours = c("blue",  "red"),limits = c(0,100),
                        name = "Inframe%") + theme_bw() + 
  ggnewscale::new_scale_colour() + 
  # geom_text_repel(aes(x = cut_mut_dis, y = diff, label = id), 
  #                 position = position_jitter(seed = 1), 
  #                 min.segment.length = unit(0, 'lines')) + 
  xlab("Distance cleavage-Mut (NT)")+
  ylab("Number Change (NT)") + 
  scale_fill_gradientn(colours = c("#FFEBCD", "#F09B28"), name = "0AA/Inframe") + 
  scale_x_continuous(breaks = (min(c(top100_sgs$cut_mut_dis, -1)) : max(c(top100_sgs$cut_mut_dis, 1))), 
                     labels = c((min(c(top100_sgs$cut_mut_dis, -1)) - 1) : -1, 1 : max(c(top100_sgs$cut_mut_dis, 1))),
                     limits = c(min(c(top100_sgs$cut_mut_dis, -1)) - 1 , max(c(top100_sgs$cut_mut_dis, 1)) + 1)) + 
  scale_y_continuous(breaks = c(0 : ceiling(max(top100_sgs$diff))), limits = c(0, ceiling(max(top100_sgs$diff)))) + 
  scale_color_manual(values = c("black", "red")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf("~/Nutstore Files/Tobin/Predict/volcano_plot_of_distance_ntdiff_fixed.pdf", width = 9, height = 6)
print(p)
dev.off()


aa0_result <- lapply(top100_result, function(x){
  x[['0AA']]
})
aa0_result <- data.frame(do.call(rbind, aa0_result))
inframe_result <- lapply(top100_result, function(x){
  x[["1AA"]]
})
inframe_result <- data.frame(do.call(rbind, inframe_result))

wt_0aa <- lapply(split(aa0_result,aa0_result$id), function(x){
  x$Pct <- x$Predicted.frequency / 
    sum(x$Predicted.frequency)
  sum(x$Pct[x$proteindiff == 0])
})
wt_0aa <- data.frame(id = names(wt_0aa), Pct = unlist(wt_0aa))

id_change <- lapply(split(inframe_result, inframe_result$id), function(x){
  sum(x$proteindiff * x$Predicted.frequency / sum(x$Predicted.frequency))
})
id_change <- data.frame(id = names(id_change), ids = unlist(id_change))
plot_data <- merge(id_change, wt_0aa, by="id")
plot_data <- merge(plot_data, top100_sgs, by="id")
plot_data <- plot_data[plot_data$Pct != 0,]
p <- ggplot(plot_data, aes(x = Pct, y = ids)) + 
  geom_point(aes(color = inframe, fill = aa0_inframe), shape = 21, size = 3, strock = 1) + 
  scale_color_gradientn(colours = c("blue",  "red"),limits = c(0,100),
                        name = "Inframe%") + theme_bw() + 
  ggnewscale::new_scale_colour() + 
  scale_fill_gradientn(colours = c("#FFEBCD", "#F09B28"), name = "0AA/Inframe") + 
  xlab("WT% In Ins1") + ylab("Total Identity Change") + 
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1), 
                                                           labels = c(0, 0.25, 0.5, 0.75, "1")) + 
  theme_bw()
pdf("~/Nutstore Files/Tobin/Predict/WT_of_AA_IDchange_fixed.pdf", width = 7.5, height = 6)
print(p)
dev.off()





arround_len <- 6


top100_sgs <- top100_sgs[order(top100_sgs$aa0_inframe, decreasing = T),]

top10_wt_0aa <- aa0_result
top10_wt_0aa <- top10_wt_0aa[top10_wt_0aa$inframe == 0& top10_wt_0aa$proteindiff == 0,]
top10_wt_0aa <- lapply(split(top10_wt_0aa$Predicted.frequency, top10_wt_0aa$id), sum)
top10_wt_0aa <- data.frame(pct = unlist(top10_wt_0aa), id = names(top10_wt_0aa))
top10_wt_0aa$mutid <- str_remove(top10_wt_0aa$id, "[-].*")
top10_wt_0aa <- top10_wt_0aa[order(top10_wt_0aa$pct, decreasing = T),]

top10_wt_0aa <- top10_wt_0aa[!duplicated(top10_wt_0aa$mutid),]
top10_wt_0aa <- top10_wt_0aa[1:10,]
top10_wt_sgs <- top100_sgs[top100_sgs$id %in% top10_wt_0aa$id,]
seq_mat <- matrix('', ncol = arround_len * 2, nrow = nrow(top10_wt_sgs))
colnames(seq_mat) <- as.character(setdiff(c(-arround_len:arround_len), 0))
rownames(seq_mat) <- top10_wt_sgs$id
top10_wt_sgs$mut <- str_remove(top10_wt_sgs$id, "[-].*")
top10_wt_sgs <- top10_wt_sgs[order(top10_wt_sgs$inframe, decreasing = T),]
top10_del1_mut <- top10_wt_sgs$id
cds_pos <- list()
pam_mat <- matrix(0, ncol = arround_len * 2, nrow = nrow(top10_wt_sgs))
colnames(seq_mat) <- as.character(setdiff(c(-arround_len:arround_len), 0))
rownames(seq_mat) <- top10_wt_sgs$id
for(i in 1 : nrow(top10_wt_sgs)){
  strand <- top10_wt_sgs$strand[i]
  mutid <- top10_wt_sgs$mut[i]
  infos <- infos_list[[mutid]]
  wtseq <- infos[["WT"]]
  if(strand == "-"){
    wtseq <- as.character(reverseComplement(DNAString(wtseq)))
  }
  
  plotseq <- str_locate(wtseq, top10_wt_sgs$right[i])[[1]]
  plotseq <- plotseq - str_length(top10_wt_sgs$PAM[i])
  cutpos <- plotseq - ifelse(top10_wt_sgs$cut_mut_dis[i] > 0, 5, 4)
  wtLen <- str_length(wtseq)
  cdspos <- 1
  while(cdspos < cutpos){
    cdspos <- cdspos + 3
    if(cdspos + 3 > cutpos){
      cds_pos[[mutid]] <- cdspos - cutpos
      break
    }
  }
  
  
  
  plotseq <- str_sub(wtseq, plotseq - 3 - arround_len,plotseq + 2)
  # left <- paste0(top10_wt_sgs$left[i], top10_wt_sgs$gRNA[i])
  # left <- strsplit(str_sub(left, -3 - arround_len , -4), "*")
  # right <- paste0(top10_wt_sgs$gRNA[i],top10_wt_sgs$PAM[i], top10_wt_sgs$right[i])
  # right <- strsplit(str_sub(right, str_length(top10_wt_sgs$gRNA)[i] - 2, str_length(top10_wt_sgs$gRNA)[i] + 3), "*")
  # if(strand == "-"){
  #   rev <- as.character(reverseComplement(DNAString(paste0(unlist(c(left, right)), collapse = ""))))
  #   seq_mat[i,] <- unlist(strsplit(rev, "*"))
  # } else {
  # mutseq <- unlist(c(left, right))
  
  seq_mat[i,] <- unlist(strsplit(plotseq, "*"))
  # }
}

indelphi_res <- list()
for(name in names(top100_result)){
  tmp <- top100_result[[name]][["indelphi"]]
  for(name2 in names(tmp)){
    indelphi_res[[name2]] <- tmp[[name2]]
  }
}



indel_mat <- matrix(0, nrow = nrow(top10_wt_sgs), ncol = 4)
colnames(indel_mat) <- c( "A", "T", "C", "G")
rownames(indel_mat) <- top10_wt_sgs$id


for(id in top10_wt_sgs$id){
  tables <- indelphi_res[[id]]
  gRNA <- top10_wt_sgs$gRNA[top10_wt_sgs$id == id]
  isIndel <- infos_list[[str_remove(id, "[-].*")]]$isIndel
  #If insert 1 mut cal del
  if(isIndel == 1){
    tables <- tables[tables$Length == 1 & tables$Category == "del",]
    tables$Pct <- tables$Predicted.frequency / sum(tables$Predicted.frequency) * 100
    leftnt <- str_sub(gRNA, -4, -4)
    rightnt <- str_sub(gRNA, -3, -3)
    if(0 %in% tables$Genotype.position){
      indel_mat[id, leftnt] <- tables$Pct[tables$Genotype.position == 0]
    }
    if(1 %in% tables$Genotype.position){
      indel_mat[id, rightnt] <- tables$Pct[tables$Genotype.position == 1]
    }
  } else {
    #for del1 mutation cal ins1
    tables <- tables[tables$Length == 1 & tables$Category == "ins",]
    tables$Pct <- tables$Predicted.frequency / sum(tables$Predicted.frequency) * 100
    for(i in 1 : nrow(tables)){
      indel_mat[id, tables$Inserted.Bases[i]] <- tables$Pct[i]
    }
  }
}
top10_inframe <- lapply(top100_result, function(x){
  x$inframe
})
top10_inframe <- data.frame(do.call(rbind, top10_inframe))
top10_inframe <- top10_inframe[top10_inframe$id %in% top10_wt_sgs$id,]
inframe_spec <- reshape2::dcast(top10_inframe, id ~ aa, value.var = "Pct2")
rownames(inframe_spec) <- inframe_spec[,1]
inframe_spec <- inframe_spec[,-1]
inframe_spec[is.na(inframe_spec)] <- 0

inframe_spec <- inframe_spec[, c(as.character(sort(as.numeric(colnames(inframe_spec)[colnames(inframe_spec) != "Other"]))), "Other")]


cols2 <- c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")
names(cols2) <- c("-2", "-1", "0", "1", "2", "Others")
cols = c( "#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(cols) <- c( "A", "T", "C", "G")

top10_aa0 <- aa0_result
top10_aa0 <- top10_aa0[top10_aa0$id %in% top10_wt_sgs$id,]
wt_0aa <- lapply(split(top10_aa0,top10_aa0$id), function(x){
  x$Pct <- x$Predicted.frequency / 
    sum(x$Predicted.frequency) * 100
  sum(x$Pct[x$proteindiff == 0])
})
wt_0aa <- data.frame(id = names(wt_0aa), Pct = unlist(wt_0aa))
rownames(wt_0aa) <- wt_0aa$id
wt_0aa <- wt_0aa[top10_wt_sgs$id,]
right_anno <- rowAnnotation(
                        "Inframe/\nIndel%" = anno_barplot(x = top10_wt_sgs$inframe, 
                                    gp = gpar(lwd = 0.1, fill = "#BBBBBB"),
                                    width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 100), 
                                    axis_param = c("side"="top")), 
                        emp = anno_empty(width = unit(2, "mm"), border = F),
                        "Inframe\nSpec" = anno_barplot2(x = inframe_spec, width = unit(2, "cm"), 
                                                        gp = gpar(lwd = 0.1, fill = c(cols2)),
                                                        bar_width = 0.5, ylim = ceiling(c(0, 100)), 
                                                        axis_param = c("side"="top")), 
                        emp2 = anno_empty(width = unit(2, "mm"), border = F),
                        InsertNT = anno_barplot2(indel_mat, 
                                                     gp = gpar(fill = cols, lwd = 0.1), 
                                                     bar_width = 0.5, 
                                                     width = unit(2, "cm"), 
                                                     axis_param = c("side"="top")), 


                            "WT/\n0AA%" = anno_barplot(x = wt_0aa$Pct, 
                                                       gp = gpar(lwd = 0.1, fill = "#BBBBBB"),
                                                       width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 100), 
                                                       axis_param = c("side"="top")))

plot_mat <- matrix(runif(arround_len * 2 * nrow(top10_wt_sgs)), nrow = nrow(top10_wt_sgs), ncol = arround_len * 2)
rownames(plot_mat) <- top10_wt_sgs$id

ht <- Heatmap(plot_mat, 
              name = "EditDis",
              column_split = c(rep(" ", arround_len), rep("   ", arround_len)),
              column_title=NULL,
              row_title = NULL,
              show_column_names = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 24),
              column_names_side = "top",
              right_annotation = right_anno,
              row_names_side = "left",
              show_row_names = T,
              cell_fun = function(j, i, x, y, w, h, fill) {
                nt = pindex(seq_mat, i, j)
                strand <- top10_wt_sgs$strand[i]
                cut_mut <- top10_wt_sgs$cut_mut_dis[i]
                
                col2 = "white"
                grid.rect(x, y, w, h, gp = gpar(fill = col2, col = col2))
                x_from <- x - unit(as.numeric(w) / 2, "npc")
                x_to <- x + unit(as.numeric(w) / 2, "npc")
                y_from <- y - unit(as.numeric(h) / 2, "npc")
                y_to <- y + unit(as.numeric(h) / 2, "npc")
                grid.polyline(x = c(x_from, x_from, x_to , x_to), 
                              y = c(y_from,y_to ,y_from, y_to), 
                              gp = gpar(col = "#D0D0D0"), 
                              id = rep(1 : (length(x) * 2), 2))
                
                mutpos <- top10_wt_sgs$cut_mut_dis[i] + 5
                isIndel <- infos_list[[str_remove(top10_wt_sgs$id[i], "[-].*")]]$isIndel
                if(j == mutpos){
                  if(isIndel == 1){
                    col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                    names(col) <- c("A", 'T','C','G')
                    grid.text(nt, x, y,
                              gp = gpar(fontsize = 20, col = col))
                  } else if(isIndel == 2){
                    col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                    names(col) <- c("A", 'T','C','G')
                    grid.rect(x = x, y= y, width = as.numeric(w)* 0.88, height = as.numeric(h)* 0.95, 
                              gp = gpar(col = "red"))
                    grid.lines(x = c(x - unit(as.numeric(w) * 0.44, "npc"), 
                                     x + unit(as.numeric(w) * 0.44, "npc")), 
                               y = c(y - unit(as.numeric(h) * 0.475, "npc"),
                                     y +  unit(as.numeric(h) * 0.475, "npc")), 
                               gp = gpar(col = "red"))
                    grid.text(nt, x, y,
                              gp = gpar(fontsize = 20, col = col))
                  } else {
                    grid.text(nt, x, y,
                              gp = gpar(fontsize = 18))
                  }

                } else if(j == 9 & cut_mut >= 5){
                  grid.text(nt, x, y,
                            gp = gpar(fontsize = 18, col = "#FFB829"))
                } else if(j > 9){
                  grid.text(nt, x, y,
                            gp = gpar(fontsize = 18, col = "#FFB829"))
                } else{
                  grid.text(nt, x, y,
                            gp = gpar(fontsize = 18))
                }
                if(j == (7 - sum(cut_mut > 0))){
                  if(strand == "+"){
                    grid.polygon(x = c(x_from - unit(as.numeric(w) / 8, "npc"), 
                                       x_from, 
                                       x_from + unit(as.numeric(w) / 8, "npc")),
                                 y = c(y + unit(as.numeric(h) / 3, "npc"), 
                                       y - unit(as.numeric(h) / 3, "npc"), 
                                       y + unit(as.numeric(h) / 3, "npc")), 
                                 gp = gpar(col = "red", fill = "red"))
                    
                  } else{
                    grid.polygon(x = c(x_from - unit(as.numeric(w) / 8, "npc"), 
                                       x_from, 
                                       x_from + unit(as.numeric(w) / 8, "npc")),
                                 y = c(y - unit(as.numeric(h) / 3, "npc"), 
                                       y + unit(as.numeric(h) / 3, "npc"), 
                                       y - unit(as.numeric(h) / 3, "npc")), 
                                 gp = gpar(col = "red", fill = "red"))
                  }
                }
                
                
              },
              col=c("white", "white"),show_heatmap_legend = F,
              border = F)
lgd_title <- "InsertNT"
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = lgd_title, type = "grid", pch = 16, 
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF"))),
  Legend(labels = c("-2", "-1", "0", "1", "2", "Others"), title = "Inframe Spec", type = "grid", pch = 16, 
         legend_gp = gpar(fill = c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")))
)


pdf("~/Nutstore Files/Tobin/Predict/Del1_top10_edit_LogoPlot_fixed.pdf", width = 11, height = 9 / 13 * 6)
draw(ht, annotation_legend_list = lgd_list)
tmp_top10_wt_sgs <- top10_wt_sgs[nrow(top10_wt_sgs) : 1,]
decorate_heatmap_body("EditDis", {
  cell_width <- 1 / 6
  cell_height <- 1 / nrow(tmp_top10_wt_sgs)
  
  for(i in 1 : nrow(tmp_top10_wt_sgs)){
    mutid <- tmp_top10_wt_sgs$mut[i]
    cdspos <- cds_pos[[mutid]]
    strand <- tmp_top10_wt_sgs$strand[i]
    cutpos <- ifelse(tmp_top10_wt_sgs$cut_mut_dis[[i]] <= 0, 6, 5)
    grid.segments(cell_width * (cutpos + cdspos - 1) - cell_width / 50 , cell_height * (i - 1) + cell_height / 20, 
                  cell_width * (cutpos + 2 + cdspos) + cell_width / 50, cell_height * (i - 1) + cell_height / 20,
                  gp = gpar(col = "green", lwd = 2))
  }
  

})

dev.off()




####Ins mut

top100_mut <- mutation_info[mutation_info$V6 %in% top100_for_indelphi$mut,]
top100_mut_ins <- top100_mut[top100_mut$V4 == "-",]
arround_len <- 6


top10_wt_0aa <- aa0_result
top10_wt_0aa <- top10_wt_0aa[top10_wt_0aa$inframe == 0& top10_wt_0aa$proteindiff == 0,]
top10_wt_0aa <- lapply(split(top10_wt_0aa$Predicted.frequency, top10_wt_0aa$id), sum)
top10_wt_0aa <- data.frame(pct = unlist(top10_wt_0aa), id = names(top10_wt_0aa))
top10_wt_0aa$mutid <- str_remove(top10_wt_0aa$id, "[-].*")
top10_wt_0aa <- top10_wt_0aa[order(top10_wt_0aa$pct, decreasing = T),]

top10_wt_0aa <- top10_wt_0aa[!duplicated(top10_wt_0aa$mutid),]
top10_wt_0aa <- top10_wt_0aa[top10_wt_0aa$mutid %in% top100_mut_ins$V6,]
top100_sgs <- top100_sgs[abs(top100_sgs$cut_mut_dis) < 6,]
top10_wt_0aa <- top10_wt_0aa[top10_wt_0aa$id %in% top100_sgs$id,]
top10_wt_0aa <- top10_wt_0aa[1:10,]
top10_wt_sgs <- top100_sgs[top100_sgs$id %in% top10_wt_0aa$id,]
seq_mat <- matrix('', ncol = arround_len * 2, nrow = nrow(top10_wt_sgs))
colnames(seq_mat) <- as.character(setdiff(c(-arround_len:arround_len), 0))
rownames(seq_mat) <- top10_wt_sgs$id
top10_wt_sgs$mut <- str_remove(top10_wt_sgs$id, "[-].*")
top10_wt_sgs <- top10_wt_sgs[order(top10_wt_sgs$inframe, decreasing = T),]
top10_ins1_mut <- top10_wt_sgs$id
cds_pos <- list()
pam_mat <- matrix(0, ncol = arround_len * 2, nrow = nrow(top10_wt_sgs))
colnames(seq_mat) <- as.character(setdiff(c(-arround_len:arround_len), 0))
rownames(seq_mat) <- top10_wt_sgs$id
for(i in 1 : nrow(top10_wt_sgs)){
  strand <- top10_wt_sgs$strand[i]
  mutid <- top10_wt_sgs$mut[i]
  mutpos <- top10_wt_sgs$cut_mut_dis[i]
  infos <- infos_list[[mutid]]
  wtseq <- infos[["seq"]]
  if(strand == "-"){
    wtseq <- as.character(reverseComplement(DNAString(wtseq)))
  }
  
  plotseq <- str_locate(wtseq, top10_wt_sgs$right[i])[[1]]
  plotseq <- plotseq - str_length(top10_wt_sgs$PAM[i])
  cutpos <- plotseq - 4
  mutpos <- mutpos + cutpos
  cdspos <- 1
  while(cdspos < cutpos){

    cdspos <- cdspos + 3
    if(cdspos + 3 > cutpos){
      cds_pos[[mutid]] <- cdspos - cutpos
      break
    }
    if(mutpos %in% c(cdspos, cdspos + 1, cdspos + 2)){
      cdspos <- cdspos + 1
    }
  }
  
  
  
  left <- paste0(top10_wt_sgs$left[i], top10_wt_sgs$gRNA[i])
  left <- strsplit(str_sub(left, -3 - arround_len , -4), "*")
  right <- paste0(top10_wt_sgs$gRNA[i],top10_wt_sgs$PAM[i], top10_wt_sgs$right[i])
  right <- strsplit(str_sub(right, str_length(top10_wt_sgs$gRNA)[i] - 2, str_length(top10_wt_sgs$gRNA)[i] + 3), "*")
  # if(strand == "-"){
  #   rev <- as.character(reverseComplement(DNAString(paste0(unlist(c(left, right)), collapse = ""))))
  #   seq_mat[i,] <- unlist(strsplit(rev, "*"))
  # } else {
  mutseq <- unlist(c(left, right))
  # }
  seq_mat[i,] <- mutseq
}

indelphi_res <- list()
for(name in names(top100_result)){
  tmp <- top100_result[[name]][["indelphi"]]
  for(name2 in names(tmp)){
    indelphi_res[[name2]] <- tmp[[name2]]
  }
}



indel_mat <- matrix(0, nrow = nrow(top10_wt_sgs), ncol = 4)
colnames(indel_mat) <- c( "A", "T", "C", "G")
rownames(indel_mat) <- top10_wt_sgs$id


for(id in top10_wt_sgs$id){
  tables <- indelphi_res[[id]]
  gRNA <- top10_wt_sgs$gRNA[top10_wt_sgs$id == id]
  isIndel <- infos_list[[str_remove(id, "[-].*")]]$isIndel
  #If insert 1 mut cal del
  if(isIndel == 1){
    tables <- tables[tables$Length == 1 & tables$Category == "del",]
    tables$Pct <- tables$Predicted.frequency / sum(tables$Predicted.frequency) * 100
    leftnt <- str_sub(gRNA, -4, -4)
    rightnt <- str_sub(gRNA, -3, -3)
    if(0 %in% tables$Genotype.position){
      indel_mat[id, leftnt] <- tables$Pct[tables$Genotype.position == 0]
    }
    if(1 %in% tables$Genotype.position){
      indel_mat[id, rightnt] <- tables$Pct[tables$Genotype.position == 1]
    }
  } else {
    #for del1 mutation cal ins1
    tables <- tables[tables$Length == 1 & tables$Category == "ins",]
    tables$Pct <- tables$Predicted.frequency / sum(tables$Predicted.frequency) * 100
    for(i in 1 : nrow(tables)){
      indel_mat[id, tables$Inserted.Bases[i]] <- tables$Pct[i]
    }
  }
}
top10_inframe <- lapply(top100_result, function(x){
  x$inframe
})
top10_inframe <- data.frame(do.call(rbind, top10_inframe))
top10_inframe <- top10_inframe[top10_inframe$id %in% top10_wt_sgs$id,]
inframe_spec <- reshape2::dcast(top10_inframe, id ~ aa, value.var = "Pct2")
rownames(inframe_spec) <- inframe_spec[,1]
inframe_spec <- inframe_spec[,-1]
inframe_spec[is.na(inframe_spec)] <- 0

inframe_spec <- inframe_spec[, c(as.character(sort(as.numeric(colnames(inframe_spec)[colnames(inframe_spec) != "Other"]))), "Other")]


cols2 <- c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")
names(cols2) <- c("-2", "-1", "0", "1", "2", "Others")
cols = c( "#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(cols) <- c( "A", "T", "C", "G")

top10_aa0 <- aa0_result
top10_aa0 <- top10_aa0[top10_aa0$id %in% top10_wt_sgs$id,]
wt_0aa <- lapply(split(top10_aa0,top10_aa0$id), function(x){
  x$Pct <- x$Predicted.frequency / 
    sum(x$Predicted.frequency) * 100
  sum(x$Pct[x$proteindiff == 0])
})
wt_0aa <- data.frame(id = names(wt_0aa), Pct = unlist(wt_0aa))
rownames(wt_0aa) <- wt_0aa$id
wt_0aa <- wt_0aa[top10_wt_sgs$id,]
right_anno <- rowAnnotation(
  "Inframe/\nIndel%" = anno_barplot(x = top10_wt_sgs$inframe, 
                                    gp = gpar(lwd = 0.1, fill = "#BBBBBB"),
                                    width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 100), 
                                    axis_param = c("side"="top")), 
  emp = anno_empty(width = unit(2, "mm"), border = F),
  "Inframe\nSpec" = anno_barplot2(x = inframe_spec, width = unit(2, "cm"), 
                                  gp = gpar(lwd = 0.1, fill = c(cols2)),
                                  bar_width = 0.5, ylim = ceiling(c(0, 100)), 
                                  axis_param = c("side"="top")), 
  emp2 = anno_empty(width = unit(2, "mm"), border = F),
  DeleteNT = anno_barplot2(indel_mat, 
                           gp = gpar(fill = cols, lwd = 0.1), 
                           bar_width = 0.5, 
                           width = unit(2, "cm"), 
                           axis_param = c("side"="top")), 


  "WT/\n0AA%" = anno_barplot(x = wt_0aa$Pct, 
                             gp = gpar(lwd = 0.1, fill = "#BBBBBB"),
                             width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 100), 
                             axis_param = c("side"="top")))

plot_mat <- matrix(runif(arround_len * 2 * nrow(top10_wt_sgs)), nrow = nrow(top10_wt_sgs), ncol = arround_len * 2)
rownames(plot_mat) <- top10_wt_sgs$id

ht <- Heatmap(plot_mat, 
              name = "EditDis",
              column_title=NULL,
              row_title = NULL,
              show_column_names = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 24),
              column_names_side = "top",
              right_annotation = right_anno,
              row_names_side = "left",
              show_row_names = T,
              cell_fun = function(j, i, x, y, w, h, fill) {
                nt = pindex(seq_mat, i, j)
                strand <- top10_wt_sgs$strand[i]
                cut_mut <- top10_wt_sgs$cut_mut_dis[i]
                
                col2 = "white"
                grid.rect(x, y, w, h, gp = gpar(fill = col2, col = col2))
                x_from <- x - unit(as.numeric(w) / 2, "npc")
                x_to <- x + unit(as.numeric(w) / 2, "npc")
                y_from <- y - unit(as.numeric(h) / 2, "npc")
                y_to <- y + unit(as.numeric(h) / 2, "npc")
                grid.polyline(x = c(x_from, x_from, x_to , x_to), 
                              y = c(y_from,y_to ,y_from, y_to), 
                              gp = gpar(col = "#D0D0D0"), 
                              id = rep(1 : (length(x) * 2), 2))
                
                mutpos <- top10_wt_sgs$cut_mut_dis[i] + 6
                if(strand == "-"){

                    mutpos <- mutpos + 1

                }
                isIndel <- infos_list[[str_remove(top10_wt_sgs$id[i], "[-].*")]]$isIndel
                if(j == mutpos){
                  if(isIndel == 1){
                    col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                    names(col) <- c("A", 'T','C','G')
                    grid.text(nt, x, y,
                              gp = gpar(fontsize = 20, col = col))
                  } else if(isIndel == 2){
                    col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                    names(col) <- c("A", 'T','C','G')
                    grid.rect(x = x, y= y, width = as.numeric(w)* 0.88, height = as.numeric(h)* 0.95, 
                              gp = gpar(col = "red"))
                    grid.lines(x = c(x - unit(as.numeric(w) * 0.44, "npc"), 
                                     x + unit(as.numeric(w) * 0.44, "npc")), 
                               y = c(y - unit(as.numeric(h) * 0.475, "npc"),
                                     y +  unit(as.numeric(h) * 0.475, "npc")), 
                               gp = gpar(col = "red"))
                    grid.text(nt, x, y,
                              gp = gpar(fontsize = 20, col = col))
                  } else {
                    grid.text(nt, x, y,
                              gp = gpar(fontsize = 18))
                  }
                  
                } else if(j > 9){
                  grid.text(nt, x, y,
                            gp = gpar(fontsize = 18, col = "#FFB829"))
                } else{
                  grid.text(nt, x, y,
                            gp = gpar(fontsize = 18))
                }
                if(j == 7){
                  if(strand == "+"){
                    grid.polygon(x = c(x_from - unit(as.numeric(w) / 8, "npc"), 
                                       x_from, 
                                       x_from + unit(as.numeric(w) / 8, "npc")),
                                 y = c(y + unit(as.numeric(h) / 3, "npc"), 
                                       y - unit(as.numeric(h) / 3, "npc"), 
                                       y + unit(as.numeric(h) / 3, "npc")), 
                                 gp = gpar(col = "red", fill = "red"))
                    
                  } else{
                    grid.polygon(x = c(x_from - unit(as.numeric(w) / 8, "npc"), 
                                       x_from, 
                                       x_from + unit(as.numeric(w) / 8, "npc")),
                                 y = c(y - unit(as.numeric(h) / 3, "npc"), 
                                       y + unit(as.numeric(h) / 3, "npc"), 
                                       y - unit(as.numeric(h) / 3, "npc")), 
                                 gp = gpar(col = "red", fill = "red"))
                  }
                }
                
                
              },
              col=c("white", "white"),show_heatmap_legend = F,
              border = F)
lgd_title <- ifelse(isIndel == 2, "InsertNT","DeleteNT")
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = lgd_title, type = "grid", pch = 16, 
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF"))),
  Legend(labels = c("-2", "-1", "0", "1", "2", "Others"), title = "Inframe Spec", type = "grid", pch = 16, 
         legend_gp = gpar(fill = c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")))
)


pdf("~/Nutstore Files/Tobin/Predict/Ins1_top10_edit_LogoPlot_fixed.pdf", width = 11, height = 9 / 13 * 6)
draw(ht, annotation_legend_list = lgd_list)
tmp_top10_wt_sgs <- top10_wt_sgs[nrow(top10_wt_sgs) : 1,]

decorate_heatmap_body("EditDis", {
  cell_width <- 1 / 12
  cell_height <- 1 / nrow(tmp_top10_wt_sgs)
  
  for(i in 1 : nrow(tmp_top10_wt_sgs)){
    mutid <- tmp_top10_wt_sgs$mut[i]
    cdspos <- cds_pos[[mutid]]
    strand <- tmp_top10_wt_sgs$strand[i]
    mutpos <- tmp_top10_wt_sgs$cut_mut_dis[i]
    if(strand == "-"){
      mutpos <- mutpos + 1
    }

    mutpos <- mutpos + 6
    add_color <- c(6 : 8) + cdspos
    apdl <- apdr <- 0
    if(mutpos %in% add_color){
      if(mutpos == add_color[1]){
        apdl <- 1
      } 
        apdr <- 1
      
    }
    grid.segments(cell_width * (cdspos + 5 + apdl) - cell_width / 50 , cell_height * (i - 1) + cell_height / 20, 
                  cell_width * (8 + cdspos + apdr) + cell_width / 50, cell_height * (i - 1) + cell_height / 20,
                  gp = gpar(col = "green", lwd = 2))
  }
  
  
})


dev.off()


mutgene <- mutation_info[,c(6, 22)]
mutgene <- mutgene[!duplicated(mutgene$V6),]
plot_data$shape <- "Other"
plot_data$shape[plot_data$id %in% top10_del1_mut] <- "Top10-Del1"
plot_data$shape[plot_data$id %in% top10_ins1_mut] <- "Top10-Ins1"

plot_data2 <- plot_data
plot_data2$mutid <- str_remove(plot_data2$id, "[-].*") 
plot_data2 <- merge(mutgene, plot_data2, by.x="V6", by.y="mutid")

plot_data2$addName <- paste(plot_data2[,2], plot_data2[,1], sep="-")
plot_data2$addName[plot_data2$shape == "Other"] <- ""
plot_data3 <- plot_data2[plot_data2$addName != "" & plot_data2$inframe > 50,]
p <- ggplot(plot_data, aes(x = Pct, y = inframe)) + 
  geom_point(aes(fill = aa0_inframe, 
                 shape= shape), size = 3, strock = 1) + 
  geom_text_repel(data = plot_data3, 
                  aes(label = addName), min.segment.length = unit(0, 'lines'),
                  max.overlaps = 10000) + 
  scale_color_gradientn(colours = c("blue",  "red"),limits = c(0,100),
                        name = "Inframe%") + theme_bw() + 
  scale_shape_manual(values = c(21, 23, 24)) + 
  ggnewscale::new_scale_colour() + 
  scale_fill_gradientn(colours = c("#FFEBCD", "#F09B28"), name = "0AA/Inframe") + 
  xlab("WT% In Indel1") + ylab("Inframe%") + 
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,100)) + 
  theme_bw()
pdf("~/Nutstore Files/Tobin/Predict/WT_of_AA_IDchange_v2_fixed.pdf", width = 7.5, height = 6)
print(p)
dev.off()




