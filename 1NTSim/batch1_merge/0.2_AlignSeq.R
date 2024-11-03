align_edit_table <- function(edit_table, wtSeq, len = 12){
  edit_table$AlignedX = ""
  edit_table$AlignedY = ""
  align_res <- lapply(1 : nrow(edit_table), function(i, wtSeq){
    x <- edit_table[i, 1]
    y <- edit_table[i, 2]
    x <- unlist(strsplit(x, "*"))
    y <- unlist(strsplit(y, "*"))
    i <- 20
    j <- 21
    acc <- 0
    wt_i <- 20
    wt_j <- 21
    while(i > 0){
      if(y[i] != '-' & x[i] != '-'){
        acc <- acc + 1
        x[i] <- y[i] <- wtSeq[wt_i]
      }
      if(acc > len){
        break
      }

      if(y[i] != '-'){
        wt_i <- wt_i - 1
      }
      i <- i - 1
    }
    if(i != 0){
      while(i > 0){
        x[i] <- wtSeq[wt_i]
        i <- i - 1
        wt_i <- wt_i - 1
      }
    }
    acc <- 0
    while(j < 41){
      if(y[j] != '-' & x[j] != '-'){
        acc <- acc + 1
        x[j] <- y[j] <- wtSeq[wt_j]
      }
      if(acc > len){
        break
      }
      if(y[j] != '-'){
        wt_j <- wt_j + 1
      }
      j <- j + 1
    }
    if(j != 41){
      while(j < 41){
        x[j] <- wtSeq[wt_j]
        j <- j + 1
        wt_j <- wt_j + 1
      }
    }

    if(wt_i != 0){
      x <- unlist(c(wtSeq[1 : wt_i], x))
      y <- unlist(c(wtSeq[1 : wt_i], y))
    }
    x_align <- paste0(x[1:40], collapse = "")
    y_align <- paste0(y[1:40], collapse = "")

    data.frame(x = x_align, y = y_align)
  }, wtSeq)
  align_res <- data.frame(do.call(rbind, align_res))
  edit_table$AlignedX <- align_res$x
  edit_table$AlignedY <- align_res$y
  edit_table
}


getStat <- function(edit_table, insert){
  if(insert){
    res <- lapply(1 : nrow(edit_table), function(i){
      x <- edit_table[i, 11]
      
      x <- unlist(strsplit(x, "*"))
      if(!"-" %in% x){
        return(data.frame(pos = -1, 
                          Pct = edit_table$Pct[i]))
      }
      data.frame(pos = which(x == '-'), 
                 Pct = edit_table$Pct[i])
    })
  }else{
    res <- lapply(1 : nrow(edit_table), function(i){
      x <- edit_table[i, 10]
      x <- unlist(strsplit(x, "*"))
      if(!"-" %in% x){
        return(data.frame(pos = -1, 
                          Pct = edit_table$Pct[i]))
      }
      data.frame(pos = which(x == '-'), 
                 Pct = edit_table$Pct[i])
    })
  }
  res <- data.frame(do.call(rbind, res))
  res <- res[res$pos != -1,]
  res <- lapply(split(res$Pct, res$pos), sum)
  res <- data.frame(pos = names(res), Pct = unlist(res))
  if(sum(c(1 : 40) %in% res$pos) != 40){
    res_append <- data.frame(pos = setdiff(as.character(1:40), res$pos), 
                             Pct = 0)
    res <- rbind(res, res_append)
  }
  res <- res[gtools::mixedorder(res$pos),]
  res$GroupPct <- res$Pct / sum(edit_table$Pct)
  res$Pct <- res$Pct * 100
  res$GroupPct <- res$GroupPct * 100
  res
}


getWtSeq <- function(edit_table){
  wtSeq <- unlist(lapply(split(edit_table[,7], edit_table[,2]), sum))
  wtSeq <- names(wtSeq)[which.max(wtSeq)]
  unlist(strsplit(wtSeq, "*"))
}

mkColorfulSeq <- function(seq1, seq2, cmp1, cmp2){
  cmp1 <- unlist(strsplit(cmp1, "*"))
  cmp2 <- unlist(strsplit(cmp2, "*"))
  if(cmp1[1] != '-' & cmp2[1] != '-'){
    pos <- min(c(which(cmp1 == '-'), which(cmp2 == '-'))) - 3
  } 
  else{
    pos <- max(c(which(cmp1 == '-'), which(cmp2 == '-'))) - 4
  }
  for(i in length(seq1) : 1){
    if(i > length(seq1) - 3){
      seq1[i] <- paste0("<span style = 'color:#B40000;'>", seq1[i], "</span>")
      seq2[i] <- paste0("<span style = 'color:#B40000;'>", seq2[i], "</span>")
    }
    if(i == pos){
      seq1[i] <- paste0("<span style = 'color:#0072B2;'>", seq1[i], "</span>")
    }
      #0072B2
  }
  return(list(seq1 = seq1, seq2 = seq2))
}




getStat2 <- function(edit_table, insert){
  if(insert){
    res <- lapply(1 : nrow(edit_table), function(i){
      x <- edit_table[i, "Reference_Sequence"]
      y <- edit_table[i, "Aligned_Sequence"]
      x <- unlist(strsplit(x, "*"))
      y <- unlist(strsplit(y, "*"))
      
      if(!"-" %in% x){
        return(data.frame(pos = -1, 
                          NT= "",
                          Pct = edit_table$Pct[i]))
      }
      pos <- which(x == '-')
      pos <- pos[which.min(abs(pos - 20))]
      data.frame(pos = pos, 
                 counts = edit_table$X.Reads[i],
                 NT = y[pos],
                 Pct = edit_table$Pct[i])
    })
  }else{
    res <- lapply(1 : nrow(edit_table), function(i){
      x <- edit_table[i, "Aligned_Sequence"]
      y <- edit_table[i, "Reference_Sequence"]
      x <- unlist(strsplit(x, "*"))
      y <- unlist(strsplit(y, "*"))
      if(!"-" %in% x){
        return(data.frame(pos = -1, 
                          NT= "",
                          Pct = edit_table$Pct[i]))
      }
      
      pos <- which(x == '-')
      pos <- pos[which.min(abs(pos - 20))]
      data.frame(pos = pos, 
                 counts = edit_table$X.Reads[i],
                 NT = y[pos],
                 Pct = edit_table$Pct[i])
    })
  }
  res <- data.frame(do.call(rbind, res))
  res <- res[res$pos != -1,]
  
  counts <- lapply(split(res$counts, paste(res$pos, res$NT, sep = "-")), sum)
  
  res <- lapply(split(res$Pct, paste(res$pos, res$NT, sep = "-")), sum)
  pos <- unlist(lapply(names(res), function(x){
    unlist(strsplit(x, "[-]"))[1]
  }))
  NT <- unlist(lapply(names(res), function(x){
    unlist(strsplit(x, "[-]"))[2]
  }))
  res <- data.frame(pos = as.numeric(pos), NT = NT, Pct = unlist(res), counts = unlist(counts))
  res <- res[gtools::mixedorder(res$pos),]
  res$GroupPct <- res$Pct / sum(edit_table$Pct)
  res$Pct <- res$Pct
  res$GroupPct <- res$GroupPct * 100
  res
}
