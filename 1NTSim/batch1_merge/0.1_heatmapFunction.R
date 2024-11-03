plotStat_merge <- function(stat, region = -9 : 6, 
                     noAA, plotOut = F, split_col, title,
                     convertGene = function(x){x}){
  library(ComplexHeatmap)
  sample_name <- names(stat)
  region_keep <- (min(region) + 1) : (max(region) - 1)
  total_stat <- lapply(stat, function(y){
    long_del <- sum(y[y[,1] <= min(region), 2])
    long_ins <- sum(y[y[,1] >= max(region), 2])
    
    long_res <- data.frame(indel_size = c(paste0(">=", min(region)), paste0(">=", max(region))), 
                           fq = c(long_del, long_ins))
    tmp_region <- region_keep[!region_keep %in% y[,1]]
    y <- y[y[,1] %in% region_keep,]
    if(length(tmp_region) != 0){
      y <- data.frame(rbind(y, data.frame(indel_size = tmp_region, fq = 0)))
      y <- y[order(y[,1], decreasing = F),]
    }
    if(plotOut){
      y <- data.frame(rbind(y, long_res))
    }
    colnames(y) <- c("Type","Freq")
    y
  })
  
  total_stat <- data.frame(do.call(cbind, total_stat))
  total_stat <- total_stat[,c(1, seq(2, ncol(total_stat), 2))]
  rownames(total_stat) <- total_stat[,1]
  total_stat <- total_stat[,-1]
  
  for(i in 1 : ncol(total_stat)){
    total_stat[,i] <- as.numeric(total_stat[,i])
  }
  total_stat_pct <- total_stat
  total_stat_pct["0",] <- 0
  for(i in 1 : ncol(total_stat_pct)){
    total_stat_pct[,i] <- total_stat_pct[,i] / sum(total_stat_pct[,i])
  }
  if(plotOut){
    total_stat_insert <- total_stat_pct[c(paste0(">=", max(region)),max(region_keep) : 1),]
    total_stat_del <- total_stat_pct[c(-1 : min(region_keep), paste0(">=", min(region))),]
  } else{
    total_stat_insert <- total_stat_pct[as.character(c(max(region_keep) : 1)),]
    total_stat_del <- total_stat_pct[as.character(c(-1 : min(region_keep))),]
  }
  total_stat_insert <- as.matrix(total_stat_insert)
  total_stat_del <- as.matrix(total_stat_del)
  total_plot <- rbind(total_stat_insert, total_stat_del)
  colnames(total_plot) <- sample_name
  library(circlize)
  plot_data2 <- colSums(total_stat)
  plot_data2 <- data.frame(sample = names(plot_data2), counts = plot_data2)
  rownames(plot_data2) <- plot_data2$sample
  plot_data2 <- plot_data2[colnames(total_plot),]
  heatanno <- HeatmapAnnotation(
    "line" =  anno_empty(border = FALSE, height = unit(2, "mm")),
    "Log Ave\nCounts" = anno_barplot(x = log1p(plot_data2$counts), border = F), 
    border = F
  )
  heatanno2 <- HeatmapAnnotation(
    "Percent" = anno_boxplot(x = noAA, border = F), border = F
  )
  
  split_col <- as.numeric(factor(split_col, levels = unique(split_col)))
  ht <- Heatmap(total_plot, 
          name = title,
          column_split = split_col,
          column_title=NULL,
          show_column_names = T, 
          cluster_rows = F, cluster_columns = F, 
          column_names_side = "top",
          top_annotation = heatanno,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x = x, y = y, width = width, height = height, 
                      gp = gpar(fill = "transparent", col = "white"))
          },
          bottom_annotation = heatanno2,
          col=c("white", "blue"),
          show_heatmap_legend = T,
          border = F)
  draw(ht)
  lens <- unique(split_col)
    
  for(i in 1:length(unique(split_col))) {
    len <- sum(split_col == lens[i])
    decorate_annotation("line", slice = i, {
      grid.lines(x = unit(c(1 / len / 2, 1 / len / 2), "npc"), 
                 y=unit(c(0, 1), "npc"))
      grid.lines(x = unit(c(1 / len / 2, 1 - 1 / len / 2), "npc"), 
                 y=unit(c(1, 1), "npc"))
      grid.lines(x = unit(c(1 - 1 / len / 2, 1 - 1 / len / 2), "npc"), 
                 y=unit(c(0, 1), "npc"))
    })
  }
  
  # heatmap_plot <- ggheatmap(total_plot, color = colorRampPalette(c( "white", "blue"))(100), 
  #           legendName = "% in Indel Reads", levels_rows = rev(rownames(total_plot)), 
  #           levels_cols = colnames(total_plot), border = "white", text_show_cols = F)
  # noAA$ID <- factor(noAA$ID, levels = colnames(total_plot))
  # noAA_plot <- ggplot(noAA, aes(x=ID, y=Percent)) + 
  #   geom_bar(stat="identity") +
  #   geom_errorbar(aes(ymin=Percent-se, ymax=Percent+se))  + theme_classic2() +  
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
  #         axis.title.x = element_blank())
  # aplot::insert_bottom(aplot::insert_bottom(top_anno, heatmap_plot, height = 5), noAA_plot)
}

plotStat_merge_v2 <- function(stat, region = -9 : 6, 
                           noAA, plotOut = F, title,
                           convertGene = function(x){x}){
  library(ComplexHeatmap)
  sample_name <- names(stat)
  region_keep <- (min(region) + 1) : (max(region) - 1)
  total_stat <- lapply(stat, function(y){
    long_del <- sum(y[y[,1] <= min(region), 2])
    long_ins <- sum(y[y[,1] >= max(region), 2])
    
    long_res <- data.frame(indel_size = c(paste0(">=", min(region)), paste0(">=", max(region))), 
                           fq = c(long_del, long_ins))
    tmp_region <- region_keep[!region_keep %in% y[,1]]
    y <- y[y[,1] %in% region_keep,]
    if(length(tmp_region) != 0){
      y <- data.frame(rbind(y, data.frame(indel_size = tmp_region, fq = 0)))
      y <- y[order(y[,1], decreasing = F),]
    }
    if(plotOut){
      y <- data.frame(rbind(y, long_res))
    }
    colnames(y) <- c("Type","Freq")
    y
  })
  
  total_stat <- data.frame(do.call(cbind, total_stat))
  total_stat <- total_stat[,c(1, seq(2, ncol(total_stat), 2))]
  rownames(total_stat) <- total_stat[,1]
  total_stat <- total_stat[,-1]
  
  for(i in 1 : ncol(total_stat)){
    total_stat[,i] <- as.numeric(total_stat[,i])
  }
  total_stat_pct <- total_stat
  total_stat_pct["0",] <- 0
  for(i in 1 : ncol(total_stat_pct)){
    total_stat_pct[,i] <- total_stat_pct[,i] / sum(total_stat_pct[,i])
  }
  if(plotOut){
    total_stat_insert <- total_stat_pct[c(paste0(">=", max(region)),max(region_keep) : 1),]
    total_stat_del <- total_stat_pct[c(-1 : min(region_keep), paste0(">=", min(region))),]
  } else{
    total_stat_insert <- total_stat_pct[as.character(c(max(region_keep) : 1)),]
    total_stat_del <- total_stat_pct[as.character(c(-1 : min(region_keep))),]
  }
  total_stat_insert <- as.matrix(total_stat_insert)
  total_stat_del <- as.matrix(total_stat_del)
  total_plot <- rbind(total_stat_insert, total_stat_del)
  colnames(total_plot) <- sample_name
  library(circlize)
  plot_data2 <- colSums(total_stat)
  plot_data2 <- data.frame(sample = names(plot_data2), counts = plot_data2)
  rownames(plot_data2) <- plot_data2$sample
  plot_data2 <- plot_data2[colnames(total_plot),]
  heatanno <- HeatmapAnnotation(
    "line" =  anno_empty(border = FALSE, height = unit(2, "mm")),
    "Log Ave\nCounts" = anno_barplot(x = log1p(plot_data2$counts), border = F), 
    border = F
  )
  heatanno2 <- HeatmapAnnotation(
    "Percent" = anno_boxplot(x = noAA, border = F), border = F
  )
  
  ht <- Heatmap(total_plot, 
                name = title,
                show_column_names = T, 
                cluster_rows = F, cluster_columns = F, 
                column_names_side = "top",
                top_annotation = heatanno,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.rect(x = x, y = y, width = width, height = height, 
                            gp = gpar(fill = "transparent", col = "white"))
                },
                bottom_annotation = heatanno2,
                col=c("white", "blue"),
                show_heatmap_legend = T,
                border = F)
  
  draw(ht)
  decorate_annotation("Percent",{
    grid.function(function(x) list(x=x, y=0.2),
                  range="x", gp=gpar(lty=2) )
  })
  # heatmap_plot <- ggheatmap(total_plot, color = colorRampPalette(c( "white", "blue"))(100), 
  #           legendName = "% in Indel Reads", levels_rows = rev(rownames(total_plot)), 
  #           levels_cols = colnames(total_plot), border = "white", text_show_cols = F)
  # noAA$ID <- factor(noAA$ID, levels = colnames(total_plot))
  # noAA_plot <- ggplot(noAA, aes(x=ID, y=Percent)) + 
  #   geom_bar(stat="identity") +
  #   geom_errorbar(aes(ymin=Percent-se, ymax=Percent+se))  + theme_classic2() +  
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
  #         axis.title.x = element_blank())
  # aplot::insert_bottom(aplot::insert_bottom(top_anno, heatmap_plot, height = 5), noAA_plot)
}



plotStat_merge_v3 <- function(stat, region = -9 : 6, 
                              noAA, plotOut = F, title,
                              convertGene = function(x){x}){
  library(ComplexHeatmap)
  sample_name <- names(stat)
  region_keep <- (min(region) + 1) : (max(region) - 1)
  total_stat <- lapply(stat, function(y){
    long_del <- sum(y[y[,1] <= min(region), 2])
    long_ins <- sum(y[y[,1] >= max(region), 2])
    
    long_res <- data.frame(indel_size = c(paste0(">=", min(region)), paste0(">=", max(region))), 
                           fq = c(long_del, long_ins))
    tmp_region <- region_keep[!region_keep %in% y[,1]]
    y <- y[y[,1] %in% region_keep,]
    if(length(tmp_region) != 0){
      y <- data.frame(rbind(y, data.frame(indel_size = tmp_region, fq = 0)))
      y <- y[order(y[,1], decreasing = F),]
    }
    y <- data.frame(rbind(y, long_res))
    colnames(y) <- c("Type","Freq")
    y
  })
  
  total_stat <- data.frame(do.call(cbind, total_stat))
  total_stat <- total_stat[,c(1, seq(2, ncol(total_stat), 2))]
  rownames(total_stat) <- total_stat[,1]
  total_stat <- total_stat[,-1]
  
  for(i in 1 : ncol(total_stat)){
    total_stat[,i] <- as.numeric(total_stat[,i])
  }
  total_stat_pct <- total_stat
  total_stat_pct["0",] <- 0
  for(i in 1 : ncol(total_stat_pct)){
    total_stat_pct[,i] <- total_stat_pct[,i] / sum(total_stat_pct[,i])
  }

  if(plotOut){
    total_stat_insert <- total_stat_pct[c(paste0(">=", max(region)),max(region_keep) : 1),]
    total_stat_del <- total_stat_pct[c(-1 : min(region_keep), paste0(">=", min(region))),]
  } else{
    total_stat_insert <- total_stat_pct[as.character(c(max(region_keep) : 1)),]
    total_stat_del <- total_stat_pct[as.character(c(-1 : min(region_keep))),]
  }
  total_stat_insert <- as.matrix(total_stat_insert)
  total_stat_del <- as.matrix(total_stat_del)
  total_plot <- rbind(total_stat_insert, total_stat_del)
  colnames(total_plot) <- sample_name
  library(circlize)
  plot_data2 <- colSums(total_stat)
  plot_data2 <- data.frame(sample = names(plot_data2), counts = plot_data2)
  rownames(plot_data2) <- plot_data2$sample
  plot_data2 <- plot_data2[colnames(total_plot),]
  heatanno <- HeatmapAnnotation(
    "line" =  anno_empty(border = FALSE, height = unit(2, "mm")),
    "Log Ave\nCounts" = anno_barplot(x = log1p(plot_data2$counts), border = F), 
    border = F
  )
  se <- unlist(apply(noAA, 2, function(x){
    x <- unlist(x)
    x <- x[!is.na(x)]
    plotrix::std.error(x)
  }))
  noAA <- unlist(apply(noAA, 2, function(x){
    mean(x, na.rm = T)
  }))
  heatanno2 <- HeatmapAnnotation(
    "Percent" = anno_barplot(x = noAA, border = F, ylim = c(0,1)), border = F
  )
  insert_color <- colorRamp2(c(0, 1), c("white", "#FF00FF"))
  del_color = colorRamp2(c(0, 1), c("white", "#0000ff"))
  lgd_insert = Legend(col_fun = insert_color, title = "insert")
  lgd_del = Legend(col_fun = del_color, title = "delete")
  pd = packLegend(lgd_insert, lgd_del)
  ht <- Heatmap(total_plot, 
                row_split = c(rep("", nrow(total_stat_insert)), 
                              rep(" ", nrow(total_stat_del))),
                name = title,
                show_column_names = T, 
                cluster_rows = F, cluster_columns = F, 
                column_names_side = "top",
                top_annotation = heatanno,
                layer_fun = function(j, i, x, y, w, h, fill) {
                  
                  # transform the matrix into a long vector
                  v = pindex(total_plot, i, j) 
                  # `j` here is also a vector with the same length of `v`
                  col = ifelse(i > nrow(total_stat_insert) , del_color(v), insert_color(v))
                  col2 = "white"
                  grid.rect(x, y, w, h, gp = gpar(fill = col, col = col2))
                },
                bottom_annotation = heatanno2,
                col=c("white", "blue"),
                show_heatmap_legend = F,
                border = F)
  
  ht <- draw(ht, annotation_legend_list = pd, 
             column_title=title,
             column_title_gp=grid::gpar(fontsize=16))
  co = column_order(ht)
  decorate_annotation("Percent", {
    od = co
    #pushViewport(viewport(xscale = c(0.5, length(od) + 0.5), yscale = c(0, 1)))
    grid.segments(seq_along(od), noAA[od] - se[od], seq_along(od), noAA[od] + se[od], default.units = "native")
    grid.segments(seq_along(od) - 0.1, noAA[od] - se[od], seq_along(od) + 0.1, noAA[od] - se[od], default.units = "native")
    grid.segments(seq_along(od) - 0.1, noAA[od] + se[od], seq_along(od) + 0.1, noAA[od] + se[od], default.units = "native")
    grid.function(function(x) {list(x=x, y=0.2)},
                  range="x", gp=gpar(lty=2) )
    #grid.yaxis()
    #popViewport()
  })
}
