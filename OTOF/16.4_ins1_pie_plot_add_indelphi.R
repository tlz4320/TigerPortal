tissue_data <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/tissue_ins1_left_right_nt_pie_data.xlsx")
cell_data <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Previous/cell_ins1_left_right_nt_pie_data.xlsx")
hek293_data <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Merged1NT/merged_ins1_left_right_nt_pie_data.xlsx")
print(load("~/Nutstore Files/Tobin/Predict/tissue_u2os_stat_data.rda"))


fake_plot_mat <- matrix(runif(12*4), nrow = 4, ncol = 4 * 4)
colnames(fake_plot_mat) <- c("A_tissue","A_cell","A_cell2", "A_merge",
                             "T_tissue","T_cell", "T_cell2", "T_merge",
                             "C_tissue","C_cell", "C_cell2", "C_merge",
                             "G_tissue", "G_cell", "G_cell2", "G_merge")
rownames(fake_plot_mat) <- c("A", "T", "C", "G")


col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")
names(col) <-  c("A", 'T','C','G')
ht <- Heatmap(fake_plot_mat,
              name = "Name",
              column_split = c(rep(c("", " ", "  ", "   "), rep(4, 4))),
              column_title=NULL,
              row_title = NULL,
              row_split = NULL,
              show_column_names = F, 
              cluster_rows = F, 
              cluster_columns = F, 
              top_annotation = HeatmapAnnotation("colname" = anno_empty(height = unit(2, "cm"), border = F)),
              left_annotation = rowAnnotation("rowname" = anno_empty(height = unit(2, "cm"), border = F)),
              right_annotation = NULL,
              show_row_names = F,
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                
                leftnt <- rownames(fake_plot_mat)[i]
                rightnt <- colnames(fake_plot_mat)[j]
                rightnt <- unlist(strsplit(rightnt, "[_]"))
                if(rightnt[2] == "tissue"){
                  pct <- unlist(tissue_data[tissue_data$leftNT == leftnt & 
                                              tissue_data$rightNT == rightnt[1], c(3:6)])
                } else if(rightnt[2] == "cell"){
                  pct <- unlist(cell_data[cell_data$leftNT == leftnt & 
                                            cell_data$rightNT == rightnt[1], c(3:6)])
                } else if(rightnt[2] == "cell2"){
                  pct <- unlist(hek293_data[hek293_data$leftNT == leftnt & 
                                              hek293_data$rightNT == rightnt[1], c(3:6)])
                } else {
                  pct <- ins1_besides_1nt_stat[[leftnt]][[rightnt[1]]]
                }
                pct <- pct[c(which(names(pct) == leftnt), which(names(pct) != leftnt))]
                if(length(pct) != 0){
                  ToNX::grid.pie(pct,x = x,y = y, w= w, h=h, 
                                 colors = col,size.scales = 0.3,
                                 init.angle = 90)
                }
                
                
              }, col = c("white", "white"), show_heatmap_legend = F,
              
              border = F)
lgd_list = list(
  Legend(labels = c("A", 'T','C','G'), title = "InsertNT", type = "grid", pch = 16,
         legend_gp = gpar(fill = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")))
)
pdf("~/Nutstore Files/Tobin/Previous/tissue_cell_293T_ins1_left_right_nt_addindelphi.pdf", 
    width = 21.5, height = 6)
draw(ht, annotation_legend_list = lgd_list)
library(gridtext)
for(i in 1 : 4){
  decorate_annotation("colname", slice = i, {
    tg <- richtext_grob(gt_render(c("A", "T", "C", "G")[i]), 
                        rot = 0, 
                        x = unit(0.5, "npc"),
                        y=unit(0.5, "npc"), hjust = 0.5, 
                        gp = gpar(fontsize = 50, 
                                  col = c("#2F89FC","#30E3CA","#66CD00", "#98ABEF")[i]))
    grid.draw(tg)
    invisible(tg)
  })
}
decorate_annotation("rowname", slice = 1, {
  tg <- richtext_grob(gt_render(rev(c("A", "T", "C", "G"))), 
                      rot = 0, 
                      y = unit(1 : 4 / 4 - 0.5 / 4, "npc"),
                      x=unit(0, "npc"), hjust = 0, 
                      gp = gpar(fontsize = 50, 
                                col = rev(c("#2F89FC","#30E3CA","#66CD00", "#98ABEF"))))
  grid.draw(tg)
  invisible(tg)
})

dev.off()

