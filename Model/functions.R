checkRS <- function(rs){
  if(str_length(rs) < 3){
    return(F)
  }
  res <- mutation_info[mutation_info$V6 == rs,]
  if(nrow(res) == 0){
    return(F)
  }
  return(res)
}

#c.1236del C
checkMut <- function(id){
  id <- unlist(strsplit(str_to_lower(id), "*"))
  if(id[1] != 'c' || id[2] != '.'){
    return(F)
  }
  step <- 3
  pos <- 0
  while(step < length(id)){
    if(str_detect(id[step], "[0-9]")){
      pos <- pos * 10 + as.numeric(id[step])
    } else {
      break
    }
    step <- step + 1
  }
  if(step == length(id) || pos == 0 || length(id) > (step + 3)){
    return(F)
  }
  muttype <- "ins"
  if(id[step] == 'd'){
    muttype <- "del"
  }
  indelNT <- id[length(id)]
  return(list(pos = pos, type= muttype, nt = indelNT))
}
generateCDS <- function(gene, mutinfo){
  gene_sel <- largest_cds[largest_cds$V13 == gene,]
  if(nrow(gene_sel) == 0){
    return(F)
  }
  cds <- gene_sel$Seq
  if(mutinfo$type == "del"){
    res <- paste0(str_sub(cds, 1, mutinfo$pos - 1), 
                  "{", 
                  str_sub(cds, mutinfo$pos, mutinfo$pos),
                  "}", 
                  str_sub(cds, mutinfo$pos + 1))
  } else {
    res <- paste0(str_sub(cds, 1, mutinfo$pos), 
                  "[", 
                  mutinfo$nt,
                  "]", 
                  str_sub(cds, mutinfo$pos + 1))
  }
  ###可能部分突变发生的位置比较靠两端 就需要续上两边的UTR做一些更全面的计算
  #假如连UTR都没有了  就暂时先不考虑了
  leftLen <- str_length(gene_sel$Left)
  leftApp <- ""
  if(mutinfo$pos < 50 & leftLen > 0){
    len <- as.integer(min(50 - mutinfo$pos, leftLen) / 3)
    if(len > 0){
      leftApp <- str_sub(gene_sel$Left, -len * 3,-1)
    }
  }
  cdsLen <- str_length(cds)
  rightLen <- str_length(gene_sel$Right)
  rightApp <- ""
  if((cdsLen - mutinfo$pos) < 50 & rightLen > 0){
    len <- as.integer(min(50 - (cdsLen - mutinfo$pos), rightLen) / 3)
    if(len > 0){
      leftApp <- str_sub(gene_sel$Right, 1, len * 3)
    }
  }
  res <- paste0(leftApp, res, rightApp)
  if(mutinfo$type == "del"){
    start <- max(as.integer((mutinfo$pos + str_length(leftApp) - 50) / 3) * 3, 0) + 1
    end <- min(as.integer((mutinfo$pos + 50) / 3) * 3 + 2, str_length(res))
  } else {
    start <- max(as.integer((mutinfo$pos + str_length(leftApp) - 50) / 3) * 3, 0) + 1
    end <- min(as.integer((mutinfo$pos + 50) / 3) * 3, str_length(res))
  }
  return(list(seq = str_sub(res, start, end), pos = str_length(leftApp)))
}



#[insert] {delete}
formatCheck <- function(inputseq){
  if(!str_detect(inputseq, "[A-Z|\\[|\\]\\{|\\}]")){
    return(F)
  }
  seq <- unlist(strsplit(inputseq, "*"))
  indel_nt <- c()
  indel_pos <- 0
  insert <- F
  delete <- F
  #1 INSERT 2 DELETE
  isIndel <- 0
  pos <- 0
  for(ch in seq){
    pos <- pos + 1
    if(ch == '['){
      indel_pos <- pos
      if(insert | delete | length(indel_nt) != 0){
        return(F)
      }
      isIndel <- 1
      insert <- T
    } else if(ch == '{'){
      indel_pos <- pos
      if(insert | delete | length(indel_nt) != 0){
        return(F)
      }
      isIndel <- 2
      delete <- T
    } else if(ch == ']'){
      if(!insert | length(indel_nt) == 0){
        return(F)
      }
      insert <- F
    } else if(ch == '}'){
      if(!delete | length(indel_nt) == 0){
        return(F)
      }
      delete <- F
    } else {
      if(insert | delete){
        indel_nt <- c(indel_nt, ch)
      }
    }
  }
  if(delete |insert){
    return(F)
  }
  seq <- str_remove_all(str_remove_all(inputseq, "\\[|\\]"), "\\{.*\\}")
  wt <- str_remove_all(str_remove_all(inputseq, "\\[.*\\]"), "\\{|\\}")
  checkinfo <- data.frame(pos = indel_pos, seq = seq, isIndel = isIndel, 
                          input = inputseq,
                          WT = wt,
                          indelNT = paste0(indel_nt,collapse = ""))
  return(checkinfo)
}


findPAM <- function(seq, pam="NGG"){
  pam <- Biostrings::DNAString(pam)
  matched <- Biostrings::matchPattern(pattern= pam,
                                       subject = DNAString(seq),
                                       fixed = "subject",
                                       max.mismatch = 0)
  pam_table <- data.frame(matched@ranges)
  pam_table[pam_table$start > 21,]
}
colorRegion <- function(seq, color, start, end){
  seq <- unlist(strsplit(seq, "*"))
  result <- c()
  step <- 0
  i <- 1
  while(T){
    while(seq[i] == "<"){
      while(seq[i] != ">"){
        i <- i + 1
      }
      i <- i + 1
    }
    step <- step + 1
    if(step == start){
      if(start != 1){
        result <- c(seq[1 : (i - 1)])
      }
      result <- c(result, paste0("<span style = 'color:", color, ";'>"))
    }
    if(step == end){
      result <- c(result, seq[(i + start - end) : i], "</span>", seq[(i + 1) : length(seq)])
      break
    }
    i <- i + 1
  }
  paste(result, collapse = "")
}

colorSeq <- function(seq, pams, pos, isIndel, indelNT){
  res <- lapply(1 : nrow(pams), function(i){
    start <- pams$start[i]
    end <- pams$end[i]
    
    left <- str_sub(seq, 1, start - 21)
    gRNA <- str_sub(seq, start - 20, start - 1)
    PAM <- str_sub(seq, start, end)
    right <- str_sub(seq, end + 1)
    colored <- paste0(left, "<span style = 'color:#8EC22F;'>", gRNA, "</span>", 
                        "<span style = 'color:#FFB829;'>", PAM, "</span>", right)
    if(isIndel == 1){
      colored <- colorRegion(colored, color = "red", start = pos, end = pos + str_length(indelNT) - 1)
    }
    data.frame(start = start, end = end, dis = pams$dis[i],
               left = left, gRNA = gRNA, PAM = PAM, right = right, colored = colored)
  })
  res <- data.frame(do.call(rbind, res))
  res
} 

callIndelPhi <- function(seq, pos, celltype="U2OS", desdir = ""){
  if(desdir == ""){
    tmpdir <- tempdir()
  } else {
    tmpdir = desdir
  }
  setwd("/home/bioinfo/code/inDelphi-model/")
  system(paste("/home/bioinfo/micromamba/envs/bio/bin/python script4R.py", 
                celltype, seq, pos, tmpdir, sep = " "))
  desdir
}

readPredict <- function(td,sgs, info, cds_region){
  ids <- sgs$id
  indelphi_res <- list()
  for(file in ids){
    indelphi_edit_table <- read.table(paste0(td, "/", file, ".csv"), sep = ",", header = T)
    indelphi_res[[file]] <- fix_indel1_func(indelphi_edit_table, sgs$gRNA[sgs$id == file])
  }
  
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
  openxlsx::write.xlsx(inframe, file="Inframe_ratio.xlsx", rowNames=F, colNames=T)
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
  aa0_result$proteindiff[grep("[*]", aa0_result$protein)] <- -1
  
  openxlsx::write.xlsx(aa0_result, file=paste0(td, "/0AA_result.xlsx"), rowNames=F, colNames=T)
  protein <- paste0(proteins, collapse = "")
  inframe_result <- lapply(sgs$id, function(id){
    tables <- indelphi_res[[id]]
    indel <- ifelse(isIndel == 1, indelLen, -indelLen)
    tables$indel <- ifelse(tables$Category == "del", -tables$Length, tables$Length)
    tables$inframe <- tables$indel + indel
    #res <- tables[tables$inframe %% 3 == 0 & abs(tables$inframe / 3) <= 1, ]
    res <- tables[tables$inframe %% 3 == 0 & abs(tables$inframe / 3) == 0, ]
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
      # tmpp <- paste0(tmpp, collapse = "")
      res$protein[i] <- paste0(tmpp, collapse = "")
      # comp <- Biostrings::pairwiseAlignment(protein, tmpp)
      # pt <- Biostrings::alignedPattern(comp)
      # sj <- Biostrings::alignedSubject(comp)
      # pt <- unlist(strsplit(as.character(pt), "*"))
      # sj <- unlist(strsplit(as.character(sj), "*"))
      res$proteindiff[i] <- sum(tmpp != proteins)
    }
    res$id <- id
    res
  })
  inframe_result <- data.frame(do.call(rbind, inframe_result))
  inframe_result$proteindiff[grep("[*]", inframe_result$protein)] <- -1
  openxlsx::write.xlsx(inframe_result, paste0(td, "/1AA_result.xlsx"), rowNames=F, colNames=T)
  tmp_sgs <- sgs
  tmp_sgs$cut_mut_dis[tmp_sgs$cut_mut_dis <= 0] <- tmp_sgs$cut_mut_dis[tmp_sgs$cut_mut_dis <= 0] - 1
  openxlsx::write.xlsx(sgs, paste0(td, "/sgRNA_information.xlsx"), rowNames=F, colNames=T)
  drawVolcano(sgs, td)
  if(isIndel == 2){
    drawLogo_v2_del(sgs, info, indelphi_res, inframe, aa0_result, td)
  } else {
    drawLogo_v2_ins(sgs, info, indelphi_res, inframe, aa0_result, td)
  }
  drawProteinStat(sgs, aa0_result, inframe_result, td)
}
drawProteinStat <- function(sgs, aa0_result, inframe_result, td){
  wt_0aa <- lapply(split(aa0_result,aa0_result$id), function(x){
    x$Pct <- x$Predicted.frequency / 
      sum(x$Predicted.frequency) * 100
    sum(x$Pct[x$proteindiff == 0])
  })
  wt_0aa <- data.frame(id = names(wt_0aa), Pct = unlist(wt_0aa))
  inframe_result <- inframe_result[inframe_result$proteindiff != -1,]
  id_change <- lapply(split(inframe_result, inframe_result$id), function(x){
    sum(x$proteindiff * x$Predicted.frequency / sum(x$Predicted.frequency))
  })
  id_change <- data.frame(id = names(id_change), ids = unlist(id_change))
  plot_data <- merge(id_change, wt_0aa, by="id")
  plot_data <- merge(plot_data, sgs, by="id")
  #plot_data <- plot_data[plot_data$ids < 1,]
  if(nrow(plot_data) > 0){
    p <- ggplot(plot_data, aes(x = Pct, y = inframe)) + 
      geom_point(aes(fill = aa0_inframe), shape = 21,size = 3) + 
      geom_text_repel(aes(x = Pct, y = inframe,label = id), 
                      min.segment.length = unit(0, 'lines'), 
                      max.overlaps = 10000) + 
      xlab("WT% In 0AA") + ylab("Inframe/Indel%") + 
      scale_x_continuous(limits = c(-1,101)) + scale_y_continuous(limits = c(-1, 101)) + 
      scale_fill_gradientn(colours = c("#FFEBCD", "#F09B28"), name = "0AA/Inframe") + 
      theme_bw()
    pdf(paste0(td,"/WT_of_AA_IDchange.pdf"), width = 7.5, height = 6)
    print(p)
    dev.off()
    png(paste0(td,"/WT_of_AA_IDchange.png"), width = 750, height = 600)
    print(p)
    dev.off()
  }
}

drawVolcano <- function(sgs, td){
  pdf(paste0(td,"/volcano_plot_of_distance_ntdiff.pdf"), width = 9, height = 6)
  assign("tmp", sgs, envir= .GlobalEnv)
  
  p <- ggplot(sgs) + geom_point(aes(x = cut_mut_dis, y =diff, color = inframe, 
                               fill = aa0_inframe), 
                           position = position_jitter(seed = 1), 
                           shape = 21, size = 3, stroke = 1) + 
    scale_color_gradientn(colours = c("blue",  "red"),limits = c(0,101),
                          name = "Inframe%") + theme_bw() + 
    ggnewscale::new_scale_colour() + 
    geom_text_repel(aes(x = cut_mut_dis, y = diff, label = id), 
                    position = position_jitter(seed = 1), 
                    max.overlaps = 10000,
                    min.segment.length = unit(0, 'lines')) + 
    xlab("Distance cleavage-Mut (NT)")+
    ylab("Number Change (NT)") + 
    scale_fill_gradientn(colours = c("#FFEBCD", "#F09B28"), name = "0AA/Inframe") + 
    scale_x_continuous(breaks = (min(c(sgs$cut_mut_dis, -1)) : max(c(sgs$cut_mut_dis, 1))), 
                       labels = c((min(c(sgs$cut_mut_dis, -1)) - 1) : -1, 1 : max(c(sgs$cut_mut_dis, 1))),
                       limits = c(min(c(sgs$cut_mut_dis, -1)) - 1 , max(c(sgs$cut_mut_dis, 1)) + 1)) + 
    scale_y_continuous(breaks = c(0 : ceiling(max(sgs$diff))), limits = c(-0.2, ceiling(max(sgs$diff)) + 0.2)) + 
    scale_color_manual(values = c("black", "red")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(p)
  dev.off()
  png(paste0(td, "/volcano_plot_of_distance_ntdiff.png"), width = 900, height = 600)
  print(p)
  dev.off()
}
drawLogo_v2_del <- function(sgs, info, res, inframe, aa0_result, td){
  
  arround_len <- 6
  seq_mat <- matrix('', ncol = arround_len * 2, nrow = nrow(sgs))
  colnames(seq_mat) <- as.character(setdiff(c(-arround_len:arround_len), 0))
  rownames(seq_mat) <- sgs$id
  isIndel <- info$isIndel
  cds_pos <- list()
  pam_mat <- matrix(0, ncol = arround_len * 2, nrow = nrow(sgs))
  for(i in 1 : nrow(sgs)){
    strand <- sgs$strand[i]
    sgsid <- sgs$id[i]
    wtseq <- info[["WT"]]
    if(strand == "-"){
      wtseq <- as.character(reverseComplement(DNAString(wtseq)))
    }
    
    plotseq <- str_locate(wtseq, sgs$right[i])[[1]]
    if(is.na(plotseq)){
      plotseq <- str_locate(wtseq, sgs$left[i])[[2]]
      cutpos <- plotseq + 17

      plotseq <- plotseq + 21
    } else{
      plotseq <- plotseq - str_length(sgs$PAM[i])
      cutpos <- plotseq - ifelse(sgs$cut_mut_dis[i] > 0, 5, 4)
    }

    wtLen <- str_length(wtseq)
    cdspos <- 1
    while(cdspos < cutpos){
      cdspos <- cdspos + 3
      if(cdspos + 3 > cutpos){
        cds_pos[[sgsid]] <- cdspos - cutpos
        break
      }
    }
    plotseq <- str_sub(wtseq, plotseq - 3 - arround_len,plotseq + 2)
    seq_mat[i,] <- unlist(strsplit(plotseq, "*"))
  }
  
  
  indel_mat <- matrix(0, nrow = nrow(sgs), ncol = 4)
  colnames(indel_mat) <- c( "A", "T", "C", "G")
  rownames(indel_mat) <- sgs$id
  
  
  for(id in sgs$id){
    tables <- res[[id]]
    gRNA <- sgs$gRNA[sgs$id == id]
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

  
  inframe_spec <- reshape2::dcast(inframe, id ~ aa, value.var = "Pct2")
  rownames(inframe_spec) <- inframe_spec[,1]
  inframe_spec <- inframe_spec[,-1]
  inframe_spec[is.na(inframe_spec)] <- 0
  
  inframe_spec <- inframe_spec[, c(as.character(sort(as.numeric(colnames(inframe_spec)[colnames(inframe_spec) != "Other"]))), "Other")]
  
  
  cols2 <- c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")
  names(cols2) <- c("-2", "-1", "0", "1", "2", "Others")
  cols = c( "#2F89FC","#30E3CA","#66CD00", "#98ABEF")
  names(cols) <- c( "A", "T", "C", "G")
  
  wt_0aa <- lapply(split(aa0_result,aa0_result$id), function(x){
    x$Pct <- x$Predicted.frequency / 
      sum(x$Predicted.frequency) * 100
    sum(x$Pct[x$proteindiff == 0])
  })
  wt_0aa <- data.frame(id = names(wt_0aa), Pct = unlist(wt_0aa))
  rownames(wt_0aa) <- wt_0aa$id
  wt_0aa <- wt_0aa[sgs$id,]
  right_anno <- rowAnnotation(
    "Inframe/\nIndel%" = anno_barplot(x = sgs$inframe, 
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
  
  plot_mat <- matrix(runif(arround_len * 2 * nrow(sgs)), nrow = nrow(sgs), ncol = arround_len * 2)
  rownames(plot_mat) <- sgs$id
  
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
                  strand <- sgs$strand[i]
                  cut_mut <- sgs$cut_mut_dis[i]
                  
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
                  
                  mutpos <- sgs$cut_mut_dis[i] + 5
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
                    
                  } else if(j == 9 & cut_mut >= 5 & cut_mut <= 7){
                    grid.text(nt, x, y,
                              gp = gpar(fontsize = 18, col = "#FFB829"))
                  } else if(j > 9){
                    grid.text(nt, x, y,
                              gp = gpar(fontsize = 18, col = "#FFB829"))
                  } else{
                    grid.text(nt, x, y,
                              gp = gpar(fontsize = 18))
                  }
                  if(j == (7 - sum(cut_mut > 0 & cut_mut < 7))){
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
  
  # draw(ht, annotation_legend_list = lgd_list)
  tmp_sgs <- sgs[nrow(sgs) : 1,]
  pdf(paste0(td,"/Indel_edit_LogoPlot.pdf"), width = 11, height = 8 / 13 * nrow(sgs))
  draw(ht, annotation_legend_list = lgd_list)
  decorate_heatmap_body("EditDis", {
    cell_width <- 1 / 6
    cell_height <- 1 / nrow(tmp_sgs)
    
    for(i in 1 : nrow(tmp_sgs)){
      mutid <- tmp_sgs$id[i]
      cdspos <- cds_pos[[mutid]]
      strand <- tmp_sgs$strand[i]
      cutpos <- ifelse(tmp_sgs$cut_mut_dis[[i]] > 0 & tmp_sgs$cut_mut_dis[[i]] < 7, 5, 6)
      grid.segments(cell_width * (cutpos + cdspos - 1) - cell_width / 50 , cell_height * (i - 1) + cell_height / 20, 
                    cell_width * (cutpos + 2 + cdspos) + cell_width / 50, cell_height * (i - 1) + cell_height / 20,
                    gp = gpar(col = "green", lwd = 2))
    }
    
    
  })
  dev.off()
  png(paste0(td,"/Indel_edit_LogoPlot.png"), width = 1100, height = 8 / 13 * nrow(sgs) * 100)
  draw(ht, annotation_legend_list = lgd_list)
  decorate_heatmap_body("EditDis", {
    cell_width <- 1 / 6
    cell_height <- 1 / nrow(tmp_sgs)
    
    for(i in 1 : nrow(tmp_sgs)){
      mutid <- tmp_sgs$id[i]
      cdspos <- cds_pos[[mutid]]
      strand <- tmp_sgs$strand[i]
      cutpos <- ifelse(tmp_sgs$cut_mut_dis[[i]] > 0 & tmp_sgs$cut_mut_dis[[i]] < 7, 5, 6)
      grid.segments(cell_width * (cutpos + cdspos - 1) - cell_width / 50 , cell_height * (i - 1) + cell_height / 20, 
                    cell_width * (cutpos + 2 + cdspos) + cell_width / 50, cell_height * (i - 1) + cell_height / 20,
                    gp = gpar(col = "green", lwd = 2))
    }
    
    
  })
  dev.off()
}


drawLogo_v2_ins <- function(sgs, info, res, inframe, aa0_result, td){
  
  arround_len <- 6
  seq_mat <- matrix('', ncol = arround_len * 2, nrow = nrow(sgs))
  colnames(seq_mat) <- as.character(setdiff(c(-arround_len:arround_len), 0))
  rownames(seq_mat) <- sgs$id
  isIndel <- info$isIndel
  cds_pos <- list()
  pam_mat <- matrix(0, ncol = arround_len * 2, nrow = nrow(sgs))
  for(i in 1 : nrow(sgs)){
    strand <- sgs$strand[i]
    mutid <- sgs$id[i]
    wtseq <- info[["seq"]]
    mutpos <- sgs$cut_mut_dis[i]
    if(strand == "-"){
      wtseq <- as.character(reverseComplement(DNAString(wtseq)))
    }
    
    plotseq <- str_locate(wtseq, sgs$right[i])[[1]]
    if(is.na(plotseq)){
      plotseq <- str_locate(wtseq, sgs$left[i])[[2]]
      cutpos <- plotseq + 17
      plotseq <- plotseq + 21
    } else{
      plotseq <- plotseq - str_length(sgs$PAM[i])
      cutpos <- plotseq - 4
    }
    #ins1 yima
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
    
    
    
    left <- paste0(sgs$left[i], sgs$gRNA[i])
    left <- strsplit(str_sub(left, -3 - arround_len , -4), "*")
    right <- paste0(sgs$gRNA[i],sgs$PAM[i], sgs$right[i])
    right <- strsplit(str_sub(right, str_length(sgs$gRNA)[i] - 2, str_length(sgs$gRNA)[i] + 3), "*")
    mutseq <- unlist(c(left, right))
    seq_mat[i,] <- mutseq
  }
  
  
  indel_mat <- matrix(0, nrow = nrow(sgs), ncol = 4)
  colnames(indel_mat) <- c( "A", "T", "C", "G")
  rownames(indel_mat) <- sgs$id
  
  
  for(id in sgs$id){
    tables <- res[[id]]
    gRNA <- sgs$gRNA[sgs$id == id]
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
  
  
  inframe_spec <- reshape2::dcast(inframe, id ~ aa, value.var = "Pct2")
  rownames(inframe_spec) <- inframe_spec[,1]
  inframe_spec <- inframe_spec[,-1]
  inframe_spec[is.na(inframe_spec)] <- 0
  
  inframe_spec <- inframe_spec[, c(as.character(sort(as.numeric(colnames(inframe_spec)[colnames(inframe_spec) != "Other"]))), "Other")]
  
  
  cols2 <- c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")
  names(cols2) <- c("-2", "-1", "0", "1", "2", "Others")
  cols = c( "#2F89FC","#30E3CA","#66CD00", "#98ABEF")
  names(cols) <- c( "A", "T", "C", "G")
  
  wt_0aa <- lapply(split(aa0_result,aa0_result$id), function(x){
    x$Pct <- x$Predicted.frequency / 
      sum(x$Predicted.frequency) * 100
    sum(x$Pct[x$proteindiff == 0])
  })
  wt_0aa <- data.frame(id = names(wt_0aa), Pct = unlist(wt_0aa))
  rownames(wt_0aa) <- wt_0aa$id
  wt_0aa <- wt_0aa[sgs$id,]
  right_anno <- rowAnnotation(
    "Inframe/\nIndel%" = anno_barplot(x = sgs$inframe, 
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
  
  plot_mat <- matrix(runif(arround_len * 2 * nrow(sgs)), nrow = nrow(sgs), ncol = arround_len * 2)
  rownames(plot_mat) <- sgs$id
  
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
                  strand <- sgs$strand[i]
                  cut_mut <- sgs$cut_mut_dis[i]
                  
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
                  
                  mutpos <- sgs$cut_mut_dis[i] + 6
                  if(strand == "-"){
                    
                    mutpos <- mutpos + 1
                    
                  }
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
  lgd_title <- "DeleteNT"
  lgd_list = list(
    Legend(labels = c("A", 'T','C','G'), title = lgd_title, type = "grid", pch = 16, 
           legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF"))),
    Legend(labels = c("-2", "-1", "0", "1", "2", "Others"), title = "Inframe Spec", type = "grid", pch = 16, 
           legend_gp = gpar(fill = c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")))
  )
  
  
  # draw(ht, annotation_legend_list = lgd_list)
  tmp_sgs <- sgs[nrow(sgs) : 1,]
  pdf(paste0(td,"/Indel_edit_LogoPlot.pdf"), width = 11, height = 8 / 13 * nrow(sgs))
  draw(ht, annotation_legend_list = lgd_list)
  decorate_heatmap_body("EditDis", {
    cell_width <- 1 / 12
    cell_height <- 1 / nrow(tmp_sgs)
    
    for(i in 1 : nrow(tmp_sgs)){
      mutid <- tmp_sgs$id[i]
      cdspos <- cds_pos[[mutid]]
      strand <- tmp_sgs$strand[i]
      mutpos <- tmp_sgs$cut_mut_dis[i]
      if(strand == "-")
        mutpos <- mutpos + 1
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
  png(paste0(td,"/Indel_edit_LogoPlot.png"), width = 1100, height = 8 / 13 * nrow(sgs) * 100)
  draw(ht, annotation_legend_list = lgd_list)
  decorate_heatmap_body("EditDis", {
    cell_width <- 1 / 12
    cell_height <- 1 / nrow(tmp_sgs)
    
    for(i in 1 : nrow(tmp_sgs)){
      mutid <- tmp_sgs$id[i]
      cdspos <- cds_pos[[mutid]]
      strand <- tmp_sgs$strand[i]
      mutpos <- tmp_sgs$cut_mut_dis[i]
      if(strand == "-")
        mutpos <- mutpos + 1
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
}


drawLogo <- function(sgs, info, res, inframe, aa0_result, td){
  
  arround_len <- 6
  seq_mat <- matrix('', ncol = arround_len * 2, nrow = nrow(sgs))
  colnames(seq_mat) <- as.character(setdiff(c(-arround_len:arround_len), 0))
  rownames(seq_mat) <- sgs$id
  isIndel <- info$isIndel
  for(i in 1 : nrow(sgs)){
    strand <- sgs$strand[i]
    left <- paste0(sgs$left[i], sgs$gRNA[i])
    left <- strsplit(str_sub(left, -3 - arround_len , -4), "*")
    right <- paste0(sgs$gRNA[i],sgs$PAM[i], sgs$right[i])
    right <- strsplit(str_sub(right, str_length(sgs$gRNA)[i] - 2, str_length(sgs$gRNA)[i] + 3), "*")
    # if(strand == "-"){
    #   rev <- as.character(reverseComplement(DNAString(paste0(unlist(c(left, right)), collapse = ""))))
    #   seq_mat[i,] <- unlist(strsplit(rev, "*"))
    # } else {
      seq_mat[i,] <- unlist(c(left, right))
    # }
  }
  
  indel_mat <- matrix(0, nrow = nrow(sgs), ncol = 4)
  colnames(indel_mat) <- c( "A", "T", "C", "G")
  rownames(indel_mat) <- sgs$id
  for(id in sgs$id){
    tables <- res[[id]]
    gRNA <- sgs$gRNA[sgs$id == id]
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
  
  inframe_spec <- reshape2::dcast(inframe, id ~ aa, value.var = "Pct2")
  rownames(inframe_spec) <- inframe_spec[,1]
  inframe_spec <- inframe_spec[,-1]
  inframe_spec[is.na(inframe_spec)] <- 0
  
  inframe_spec <- inframe_spec[, c(as.character(sort(as.numeric(colnames(inframe_spec)[colnames(inframe_spec) != "Other"]))), "Other")]
  
  
  cols2 <- c("#21913C","#F0BE5A", "#EE2626", "#F0AAAA", "#356FBA", "gray")
  names(cols2) <- c("-2", "-1", "0", "1", "2", "Others")
  cols = c( "#2F89FC","#30E3CA","#66CD00", "#98ABEF")
  names(cols) <- c( "A", "T", "C", "G")
  
  
  wt_0aa <- lapply(split(aa0_result,aa0_result$id), function(x){
    x$Pct <- x$Predicted.frequency / 
      sum(x$Predicted.frequency) * 100
    sum(x$Pct[x$proteindiff == 0])
  })
  wt_0aa <- data.frame(id = names(wt_0aa), Pct = unlist(wt_0aa))
  rownames(wt_0aa) <- wt_0aa$id
  wt_0aa <- wt_0aa[sgs$id,]
  right_anno <- rowAnnotation(InsertNT = anno_barplot2(indel_mat, 
                                                       gp = gpar(fill = cols, lwd = 0.1), 
                                                       bar_width = 0.5, 
                                                       width = unit(2, "cm")), 
                              emp = anno_empty(width = unit(2, "mm"), border = F),
                              "Inframe\nSpec" = anno_barplot2(x = inframe_spec, width = unit(2, "cm"), 
                                                             gp = gpar(lwd = 0.1, fill = c(cols2)),
                                                             bar_width = 0.5, ylim = ceiling(c(0, 100))), 
                              emp2 = anno_empty(width = unit(2, "mm"), border = F),
                              "Inframe/\nIndel%" = anno_barplot(x = sgs$inframe, 
                                                                 gp = gpar(lwd = 0.1, fill = "#BBBBBB"),
                                                                 width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 100)), 
                              "WT/\n0AA%" = anno_barplot(x = wt_0aa$Pct, 
                                                          gp = gpar(lwd = 0.1, fill = "#BBBBBB"),
                                                          width = unit(2, "cm"), bar_width = 0.5, ylim = c(0, 100)))
  
  plot_mat <- matrix(runif(arround_len * 2 * nrow(sgs)), nrow = nrow(sgs), ncol = arround_len * 2)
  rownames(plot_mat) <- sgs$id
  
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
                  strand <- sgs$strand[i]
                  
                  
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

                    mutpos <- sgs$cut_mut_dis[i] + 6
                    
                    if(j == mutpos && isIndel == 1){
                      col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
                      names(col) <- c("A", 'T','C','G')
                      grid.text(nt, x, y,
                                gp = gpar(fontsize = 20, col = col))
                    }else{
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
  pdf(paste0(td,"/Indel_edit_LogoPlot.pdf"), width = 11, height = 6 / 13 * nrow(sgs))
  draw(ht, annotation_legend_list = lgd_list)
  dev.off()
  png(paste0(td,"/Indel_edit_LogoPlot.png"), width = 1100, height = 6 / 13 * nrow(sgs) * 100)
  draw(ht, annotation_legend_list = lgd_list)
  dev.off()
  
}
renderResult <- function(output, td){
  output$resultList <- renderUI({
    tagList(
      fluidPage(fluidRow(column(2),
                         column(1, "Volcano Plot:"),
                         column(4, img(src = paste0("tmp/",basename(td), "/volcano_plot_of_distance_ntdiff.png"))),
                         height = "30%", align = "center"
      ), 
      fluidRow(column(2),
                         column(1, "Protein Plot:"),
                         column(4, img(src = paste0("tmp/",basename(td), "/WT_of_AA_IDchange.png"))),
                         height = "30%", align = "center"
      ),
      fluidRow(column(2),column(1, "Logo Plot:"),
               column(4, img(src = paste0("tmp/",basename(td), "/Indel_edit_LogoPlot.png"))),
                                  height = "40%", align = "center"
      )
    ))
  })
}



