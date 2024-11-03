###在4.1之后seqs_table

old_wt_seq <- lapply(total_old_edit_table_rev, function(tmp){
  tmp <- tmp[[1]]
  tmp <- tmp[tmp$Unedited == 'True',]
  getWtSeq(tmp)
})
noAA_change_seq <- list()

for(i in 1 : nrow(total_used_seq_tissue)){
  id <- total_used_seq_tissue$result[i]
  edit_tables <- total_old_edit_table_rev[[id]]
  edit_table <- lapply(edit_tables, function(edit_table){
    edit_table <- edit_table[edit_table$n_inserted + edit_table$n_deleted != 0,]
    edit_table$Pct <- edit_table[,7] / sum(edit_table[,7]) * 100
    edit_table <- lapply(split(edit_table, edit_table$Aligned_Sequence), 
                         function(x){
                           if(nrow(x) == 1){
                             return(x)
                           }
                           x[,9] <- sum(x[,9])
                           return(x[1,])
                         })
    edit_table <- data.frame(do.call(rbind, edit_table))
    edit_table <- edit_table[order(edit_table$Pct, decreasing = T),]
    ##选出只有Ins1的结果
    edit_table <- edit_table[edit_table$n_deleted == 0 & 
                               edit_table$n_inserted == 1& 
                               edit_table$n_mutated <= 1,]
    edit_table
  })
  noAA_change_seq[[id]] <- edit_table
  
}





seqs_table <- list()
for(id in names(noAA_change_seq)){
  edit_table <- noAA_change_seq[[id]]
  
  edit_table <- data.frame(do.call(rbind, edit_table))
  if(nrow(edit_table) == 0)
    next
  if(sum(edit_table[,8]) < 0.1)
    next
  edit_table <- data.frame(do.call(rbind, lapply(split(edit_table, edit_table$Aligned_Sequence), function(x){
    x[,7] <- mean(x[,7])
    x[,8] <- mean(x[,8])
    x[,9] <- mean(x[,9])
    x
    x[1,]
  })))
  edit_table$id <- id
  seqs_table[[id]] <- edit_table[,c(1,2, 7,8, 9, 10)]
}
seqs_table <- data.frame(do.call(rbind, seqs_table))

old_ins1_pct_table <- lapply(unique(seqs_table$id), function(id){
  ins1_table <- seqs_table[seqs_table$id == id,]
  wt_seq <- old_wt_seq[[id]]
  tmp <- getStat2(ins1_table, T)
  tmp$id <- id
  tmp
})
names(old_ins1_pct_table) <- unique(seqs_table$id)

besides_nt_stat <- list()
####这个是由于这个里面有不少是负链上面的 需要反向互补  不然左右碱基找的是错的
revComp <- function(x){
  as.character(reverseComplement(DNAString(x)))
}
for(id in names(old_ins1_pct_table)){
  tmp_pct_table <- old_ins1_pct_table[[id]]
  tmp_pct_table <- tmp_pct_table[tmp_pct_table$NT != 'N',]
  tmp_pct_table_sel <- tmp_pct_table[tmp_pct_table$pos %in% c(20, 21),]
  strand <- total_used_seq$strand[total_used_seq$result == id]
  wt_seq <- old_wt_seq[[id]]
  if(sum(tmp_pct_table_sel$Pct) > 1){
    if(strand == "+"){
      leftnt <- wt_seq[20]
      rightnt <- wt_seq[21]
    }else{
      leftnt <- revComp(wt_seq[21])
      rightnt <- revComp(wt_seq[20])
    }
    if(leftnt == "G" & rightnt == "C"){
      print(id)
    }
    if(!leftnt %in% names(besides_nt_stat)){
      besides_nt_stat[[leftnt]] <- list()
    }
    tmp_besides_nt_stat <- besides_nt_stat[[leftnt]]
    if(!rightnt %in% names(tmp_besides_nt_stat)){
      tmp_nt_stat <- c(0,0,0,0)
      names(tmp_nt_stat) <- c("A", "T", "C", "G")
      tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
    }
    tmp_nt_stat <- tmp_besides_nt_stat[[rightnt]]
    
    for(index in 1 : nrow(tmp_pct_table_sel)){
      nt <- tmp_pct_table_sel$NT[index]
      if(strand != "+"){
        nt <- revComp(nt)
      }
      
      tmp_nt_stat[nt] <- tmp_nt_stat[nt] + tmp_pct_table_sel$Pct[index]
    }
    tmp_besides_nt_stat[[rightnt]] <- tmp_nt_stat
    besides_nt_stat[[leftnt]] <- tmp_besides_nt_stat
  }
}
library(ComplexHeatmap)
library(gridtext)
fake_plot_mat <- matrix(runif(4*4), nrow = 4, ncol = 4)
rownames(fake_plot_mat) <- c("A", "T", "C", "G")
colnames(fake_plot_mat) <- c("A", "T", "C", "G")
pdf("~/data/project/ear_project/gene_therapy_ll/Previews/Result/old_del_mut_left_right_nt_pie_chart.pdf", width = 6, height = 5)
ht <- Heatmap(fake_plot_mat, 
              name = "Name",
              column_split = NULL,
              column_title=NULL,
              row_title = NULL,
              row_split = NULL,
              show_column_names = T, 
              cluster_rows = F, 
              cluster_columns = F, 
              row_names_gp = gpar(fontsize = 24),
              column_names_rot = 0,
              column_names_centered = T,
              column_names_gp = gpar(fontsize = 24),
              column_names_side = "top", top_annotation = NULL,
              left_annotation = NULL,
              right_annotation = NULL,
              show_row_names = T,
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                leftnt <- rownames(fake_plot_mat)[i]
                rightnt <- colnames(fake_plot_mat)[j]
                if(leftnt %in% names(besides_nt_stat) & rightnt %in% names(besides_nt_stat[[leftnt]])){
                  pct <- besides_nt_stat[[leftnt]][[rightnt]]
                  
                  labels <- names(pct)
                  pct <- unlist(pct)
                  col = c(2 : 5)
                  names(col) <- c("A", 'T','C','G')
                  grid.pie(pct,x = x,y = y, w= w, h=h, colors = col)
                }

                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16, 
         legend_gp = gpar(fill = 2:5))
)
draw(ht, annotation_legend_list = lgd_list)
dev.off()
