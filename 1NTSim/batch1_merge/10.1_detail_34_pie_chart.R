#接下来选择出切割位点左边碱基相同但是右边碱基不相同的例子
bulge_pos$sameL <- unlist(lapply(1:nrow(bulge_pos), function(i){
  seq1 <- unlist(strsplit(bulge_pos[i,1], "*"))
  seq1 <- seq1[seq1 != "-"]
  seq2 <- unlist(strsplit(bulge_pos[i,2], "*"))
  seq2 <- seq2[seq2 != "-"]
  
  return(seq1[17] == seq2[17])
}))

bulge_pos$diffR <- unlist(lapply(1:nrow(bulge_pos), function(i){
  seq1 <- unlist(strsplit(bulge_pos[i,1], "*"))
  seq1 <- seq1[seq1 != "-"]
  seq2 <- unlist(strsplit(bulge_pos[i,2], "*"))
  seq2 <- seq2[seq2 != "-"]
  
  return(seq1[18] != seq2[18])
}))

table(bulge_pos$sameL & bulge_pos$diffR)

detail_34_bulge_pos <- bulge_pos[bulge_pos$sameL & bulge_pos$diffR,]
detail_34_bulge_pos$del_id <-  unlist(lapply(detail_34_bulge_pos$V4, function(x){
  tmp <- sgCmp[sgCmp$id == x,]
  #把Ins1放在前面
  tmp <- tmp[order(tmp$isInsert, decreasing = T),]
  tmp$id2[2]
}))
detail_34_bulge_pos$ins_id <- unlist(lapply(detail_34_bulge_pos$V4, function(x){
  tmp <- sgCmp[sgCmp$id == x,]
  #把Ins1放在前面
  tmp <- tmp[order(tmp$isInsert, decreasing = T),]
  tmp$id2[1]
}))
detail_34_bulge_pos <- detail_34_bulge_pos[detail_34_bulge_pos$del_id %in% names(nt_counts_mats) & 
                                             detail_34_bulge_pos$ins_id %in% names(nt_counts_mats), ]
del_one_stat <- list()

for(id in detail_34_bulge_pos$del_id){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat[id,])
  for(i in 12){
    #还是过滤1一下的结果
    if(sum(tmp_pct_mat[,i]) > 1){
      leftnt <- tmp_seq_mat[i - 1]
      rightnt <- tmp_seq_mat[i + 1]
      if(!leftnt %in% names(del_one_stat)){
        del_one_stat[[leftnt]] <- list()
      }
      tmp_besides_nt_stat <- del_one_stat[[leftnt]]
      if(!rightnt %in% names(tmp_besides_nt_stat)){
        tmp_nt_stat <- c(0,0,0,0)
        names(tmp_nt_stat) <- c("A", "T", "C", "G")
        tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      }
      tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
      
      for(nt in rownames(tmp_counts_mat)){
        tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_counts_mat[nt, i]
      }
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      del_one_stat[[leftnt]] <- tmp_besides_nt_stat
    }
  }
}

ins_one_stat <- list()

for(id in detail_34_bulge_pos$ins_id){
  tmp_counts_mat <- nt_counts_mats[[id]]
  tmp_pct_mat <- nt_mats[[id]]
  tmp_seq_mat <- unlist(seq_mat[id,])
  del_id <- detail_34_bulge_pos$del_id[detail_34_bulge_pos$ins_id == id]
  del_right <- unlist(seq_mat[del_id,13])
  for(i in 12){
    #还是过滤1一下的结果
    if(sum(tmp_counts_mat[,i]) > 500){
      leftnt <- tmp_seq_mat[i - 1]
      rightnt <- paste0(del_right,tmp_seq_mat[i + 1])
      if(!leftnt %in% names(ins_one_stat)){
        ins_one_stat[[leftnt]] <- list()
      }
      tmp_besides_nt_stat <- ins_one_stat[[leftnt]]
      if(!rightnt %in% names(tmp_besides_nt_stat)){
        tmp_nt_stat <- c(0,0,0,0)
        names(tmp_nt_stat) <- c("A", "T", "C", "G")
        tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      }
      tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
      
      for(nt in rownames(tmp_counts_mat)){
        tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_counts_mat[nt, i]
      }
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
      ins_one_stat[[leftnt]] <- tmp_besides_nt_stat
    }
  }
}


fake_plot_mat <- matrix(runif(8 * 16), nrow = 8, ncol = 16)

nt_plot_mat <- matrix("", nrow = 8, ncol = 16)
for(i in seq(1, 16, 4)){
  nt_plot_mat[1, (i : (i + 3))] <- c("A", "G", "C", "T")
  nt_plot_mat[3, (i : (i + 3))] <- c("G", "A", "C", "T")
  nt_plot_mat[5, (i : (i + 3))] <- c("C", "A", "G", "T")
  nt_plot_mat[7, (i : (i + 3))] <- c("T", "A", "C", "G")
}

top_anno <- HeatmapAnnotation(
  "nt" =  anno_empty(border = FALSE, height = unit(10, "mm")),
  "lines" = anno_empty(border = FALSE, height = unit(1, "mm")), 
  border = F
)
pdf("~/data/project/ear_project/gene_therapy_ll/Result/left_right_detail_edit_result.pdf", width = 7, height = 4)
Heatmap(fake_plot_mat, 
        name = "STAT",
        column_split = rep(c(" ", "  ", "   ", "    "), c(4,4,4,4)),
        column_title=NULL,
        row_title = NULL,
        row_split = NULL,
        show_column_names = F, 
        cluster_rows = F, 
        cluster_columns = F, 
        column_names_rot = 0,
        column_names_centered = T,
        column_names_gp = gpar(fontsize = 24),
        column_names_side = "top", 
        top_annotation = top_anno,
        left_annotation = NULL,
        right_annotation = NULL,
        show_row_names = T,
        row_names_side = "left",
        cell_fun = function(j, i, x, y, w, h, fill) {
          order_seq <- c("A", "T", "C", "G")
          nt_left <- order_seq[ceiling(j / 4)]
          nt = pindex(nt_plot_mat, i, j)
          col = c(2 : 5)
          names(col) <- c("A", 'T','C','G')
          sg <- rownames(seq_mat)[i]
          if(i %% 2 == 1){
            grid.text(nt, x, y, 
                      gp = gpar(fontsize = 20, col = col[nt]))
          }
          nt_right = pindex(nt_plot_mat, i - 1, j)
          if(nt == ""){
            if(j %% 4 == 1){
              sel_stat = del_one_stat
            } else {
              sel_stat = ins_one_stat
              nt_right <- paste0(pindex(nt_plot_mat, i - 1, ceiling(j / 4) * 4 - 3),nt_right)
            }
            if(nt_left %in% names(sel_stat)){
              if(nt_right %in% names(sel_stat[[nt_left]])){
                pct <- sel_stat[[nt_left]][[nt_right]]
                labels <- names(pct)
                pct <- unlist(pct)
                col = c(2 : 5)
                names(col) <- c("A", 'T','C','G')
                pct <- c(0, cumsum(pct)/sum(pct))
                dx <- diff(pct)
                nx <- length(dx)
                twopi <- 2 * pi
                init.angle <- 0
                t2xy <- function(t) {
                  t2p <- twopi * t + init.angle * pi/180
                  list(x = 0.8 * cos(t2p), y = radius * sin(t2p))
                }
                maxx <- max(unlist(lapply(1 : nx, function(i){
                  n <- max(10, floor(200 * dx[i]))
                  P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                  max(abs(P$x))
                })))
                maxy <- max(unlist(lapply(1 : nx, function(i){
                  n <- max(10, floor(200 * dx[i]))
                  P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                  max(abs(P$y))
                })))
                facy <- as.numeric(h) / maxy * 0.3
                facx <- as.numeric(w) / maxx * 0.3
                for (i in 1L:nx) {
                  n <- max(10, floor(200 * dx[i]))
                  
                  P <- t2xy(seq.int(pct[i], pct[i + 1], length.out = n))
                  
                  grid.polygon(unit(as.numeric(x) + c(P$x * facx, 0), units = "npc"), 
                               unit(as.numeric(y) + c(P$y * facy, 0), units = "npc"), 
                               gp = gpar(fill = col[labels[i]]))
              }
            }
              
            }
          }
          
        },layer_fun = function(j, i, x, y, w, h, fill) {
          grid.lines(x = unit(c(0.25, 0.25),units = "npc"))
        }, col=c("white", "white"),show_heatmap_legend = F,
        
        border = T)



for(i in 1 : 4){
  decorate_annotation("nt", slice = i, {
    tg <- richtext_grob(gt_render(c("A", "T", "C", "G")[i]), 
                        rot = 0, 
                        x = unit(0.5, "npc"),
                        y=unit(0.5, "npc"), hjust = 0.5, gp = gpar(fontsize = 20))
    grid.draw(tg)
    invisible(tg)
  })
  
}
dev.off()






