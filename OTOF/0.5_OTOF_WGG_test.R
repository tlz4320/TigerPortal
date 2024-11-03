setwd("~/data/project/ear_project/gene_therapy_wgg/batch6/fastq/")

indel_result <- read.table("withNGG_res/CRISPResso_on_withNGG/Alleles_frequency_table_around_sgRNA_agggcgtgccccccgaacg.txt", 
                         sep="\t", header = T, comment.char = "")
indel <- indel_result$n_inserted + indel_result$n_deleted
indel_result <- indel_result[indel != 0 & indel < 15,]
indel_result$Pct <- indel_result[,8] / sum(indel_result[,8]) * 100
indel_result <- indel_result[indel_result$Pct > 0.1,]
indel_matrix <- matrix(0, nrow = nrow(indel_result), ncol = 40)

for(i in 1 : nrow(indel_result)){
  if(indel_result$n_inserted[i] != 0){
    pos <- unlist(strsplit(indel_result$Reference_Sequence[i], "*"))
    pos <- which(pos == "-")
    indel_matrix[i, pos] <- 1
  }
  if(indel_result$n_deleted[i] != 0){
    pos <- unlist(strsplit(indel_result$Aligned_Sequence[i], "*"))
    pos <- which(pos == "-")
    indel_matrix[i, pos] <- -1
  }
}

seq_matrix <- matrix("", nrow = nrow(indel_result), ncol = 40)
for(i in 1 : nrow(indel_result)){
  seq <- unlist(strsplit(indel_result$Aligned_Sequence[i], "*"))
  seq_matrix[i,] <- seq
}
fake_mat <- matrix(runif(nrow(seq_matrix) * ncol(seq_matrix)), nrow = nrow(seq_matrix))


right_anno <- rowAnnotation(pct2 = anno_barplot(indel_result$Pct, border = T, 
                                                gp = gpar(fill = "#CA2D26"), 
                                                axis_param = c("side"="top"), bar_width = 0.5,
                                                ylim = c(0,100), width = unit(2, "cm")), 
                            show_annotation_name = F)
# tmp <- ins_result[,c("Pct", "se")]
# openxlsx::write.xlsx(tmp,
#                      file= "~/Nutstore Files/Tobin/OTOF_sg1_wt_ins1_row_stat.xlsx", 
#                      rowNames=F)

ht <- Heatmap(fake_mat, rect_gp = gpar(col = "gray"), col = c("white", "white"), 
              cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, 
              right_annotation = right_anno,
              cell_fun = function(j, i, x, y, w, h, fill) {
                nt = pindex(seq_matrix, i, j)
                ins <- pindex(indel_matrix, i, j)
                if(ins > 0){
                  col <- "#CA2D26"
                  grid.rect(x , y, w, h, gp = gpar(col = col, fill = "white"))
                }
                if(ins < 0){
                  col <- "#5A962F"
                  grid.rect(x , y, w, h, gp = gpar(col = col, fill = "white"))
                }
                if(j %in% c(25 : 27)){
                  grid.text(nt, x, y,  gp = gpar(fontsize = 10, col = "#CA2D26"))
                }
                else{
                  grid.text(nt, x, y,  gp = gpar(fontsize = 10, col = "black"))
                }
                
              }
)
pdf("~/Nutstore Files/Tobin/Merged1NT/OTOF_withNGG.pdf", 
    width = 8, height = 2)
draw(ht)
dev.off()

setwd("~/data/project/ear_project/gene_therapy_wgg/batch6/fastq/")

indel_result <- read.table("noNGG_res/CRISPResso_on_noNGG/Alleles_frequency_table_around_sgRNA_agggcgtgccccccgaacgg.txt", 
                           sep="\t", header = T, comment.char = "")
indel <- indel_result$n_inserted + indel_result$n_deleted
indel_result <- indel_result[indel != 0 & indel < 15,]
indel_result$Pct <- indel_result[,8] / sum(indel_result[,8]) * 100
indel_result <- indel_result[indel_result$Pct > 0.1,]
indel_matrix <- matrix(0, nrow = nrow(indel_result), ncol = 40)

for(i in 1 : nrow(indel_result)){
  if(indel_result$n_inserted[i] != 0){
    pos <- unlist(strsplit(indel_result$Reference_Sequence[i], "*"))
    pos <- which(pos == "-")
    indel_matrix[i, pos] <- 1
  }
  if(indel_result$n_deleted[i] != 0){
    pos <- unlist(strsplit(indel_result$Aligned_Sequence[i], "*"))
    pos <- which(pos == "-")
    indel_matrix[i, pos] <- -1
  }
}

seq_matrix <- matrix("", nrow = nrow(indel_result), ncol = 40)
for(i in 1 : nrow(indel_result)){
  seq <- unlist(strsplit(indel_result$Aligned_Sequence[i], "*"))
  seq_matrix[i,] <- seq
}
fake_mat <- matrix(runif(nrow(seq_matrix) * ncol(seq_matrix)), nrow = nrow(seq_matrix))


right_anno <- rowAnnotation(pct2 = anno_barplot(indel_result$Pct, border = T, 
                                                gp = gpar(fill = "#CA2D26"), 
                                                axis_param = c("side"="top"), bar_width = 0.5,
                                                ylim = c(0,100), width = unit(2, "cm")), 
                            show_annotation_name = F)
# tmp <- ins_result[,c("Pct", "se")]
# openxlsx::write.xlsx(tmp,
#                      file= "~/Nutstore Files/Tobin/OTOF_sg1_wt_ins1_row_stat.xlsx", 
#                      rowNames=F)

ht <- Heatmap(fake_mat, rect_gp = gpar(col = "gray"), col = c("white", "white"), 
              cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, 
              right_annotation = right_anno,
              cell_fun = function(j, i, x, y, w, h, fill) {
                nt = pindex(seq_matrix, i, j)
                ins <- pindex(indel_matrix, i, j)
                if(ins > 0){
                  col <- "#CA2D26"
                  grid.rect(x , y, w, h, gp = gpar(col = col, fill = "white"))
                }
                if(ins < 0){
                  col <- "#5A962F"
                  grid.rect(x , y, w, h, gp = gpar(col = col, fill = "white"))
                }
                if(j %in% c(25 : 27)){
                  grid.text(nt, x, y,  gp = gpar(fontsize = 10, col = "#CA2D26"))
                }
                else{
                  grid.text(nt, x, y,  gp = gpar(fontsize = 10, col = "black"))
                }
                
              }
)
pdf("~/Nutstore Files/Tobin/Merged1NT/OTOF_noNGG.pdf", 
    width = 8, height = 1.6)
draw(ht)
dev.off()



tmp <- read.table("noNGG_res/CRISPResso_on_noNGG/Indel_histogram.txt", header = T)
tmp <- tmp[tmp$indel_size != 0,]
sum(tmp$fq[abs(tmp$indel_size) == 1]) / sum(tmp$fq) * 100
tmp <- read.table("withNGG_res/CRISPResso_on_withNGG//Indel_histogram.txt", header = T)
tmp <- tmp[tmp$indel_size != 0,]
sum(tmp$fq[abs(tmp$indel_size) == 1]) / sum(tmp$fq) * 100

