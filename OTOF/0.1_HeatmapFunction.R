plotStat <- function(stat, region = -18 : 18, showname = F, title = "",
                     shownumber = F, split_col = NULL, restore = NULL, 
                     colNames = NULL,
                     maxRestore = 1, convertGene = function(x){x}){
  library(ComplexHeatmap)
  library(patchwork)
  library(gridtext)
  if(!is.null(split_col)){
    showname <- F
  }
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
  total_stat_insert <- total_stat_pct[c(paste0(">=", max(region)),max(region_keep) : 1),]
  total_stat_del <- total_stat_pct[c(-1 : min(region_keep), paste0(">=", min(region))),]
  
  total_stat_insert <- as.matrix(total_stat_insert)
  total_stat_del <- as.matrix(total_stat_del)
  total_plot <- rbind(total_stat_insert, total_stat_del)
  colnames(total_plot) <- sample_name
  if(is.null(restore)){
    rowsplit <- c(rep(" ", nrow(total_stat_insert)), 
                  rep("  ",nrow(total_stat_del)))
  }else{
    rowsplit <- c("", 
                  rep(" ", nrow(total_stat_insert)), 
                  rep("  ",nrow(total_stat_del)))
  }
  library(circlize)
  insert_color <- colorRamp2(c(0, 1), c("white", "#285A32"))
  del_color = colorRamp2(c(0, 1), c("white", "#1432EB"))
  restore_color <- colorRamp2(c(0, maxRestore), c("white", "#E13C32"))
  lgd_restore = Legend(col_fun = restore_color, title = "restore")
  lgd_insert = Legend(col_fun = insert_color, title = "insert")
  lgd_del = Legend(col_fun = del_color, title = "delete")
  pd = packLegend(lgd_insert, lgd_del, lgd_restore)
  unique_split <- unique(split_col)
  if(is.null(split_col)){
    heatanno <- NULL
  }else{
    heatanno <- 
      HeatmapAnnotation(group = anno_empty(border = F, height = unit(30, "mm")),
                        "line" = 
                          anno_empty(border = FALSE, height = unit(2, "mm"))
      )
  }
  
  if(!is.null(split_col)){
    split_col <- factor(split_col, levels = unique(split_col))
  }
  
  total_plot <- rbind(restore, total_plot)
  ht <- Heatmap(total_plot, 
                name = title,
                column_split = split_col,
                column_title=NULL,
                row_split = rowsplit,
                show_column_names = showname, 
                cluster_rows = F, cluster_columns = F, 
                column_title_rot = 45,
                column_names_side = "top",
                top_annotation = heatanno,
                layer_fun = function(j, i, x, y, w, h, fill) {
                  
                  # transform the matrix into a long vector
                  v = pindex(total_plot, i, j) 
                  
                  # `j` here is also a vector with the same length of `v`
                  col = ifelse(rowsplit[i] == " " , insert_color(v),
                               ifelse(rowsplit[i] == "", restore_color(v), del_color(v)))
                  col2 = "white"
                  if(!is.null(split_col)){
                    col2 = col
                  }
                  grid.rect(x, y, w, h, gp = gpar(fill = col, col = col2))
                  x_from <- x - unit(as.numeric(w) / 2, "npc")
                  x_to <- x + unit(as.numeric(w) / 2, "npc")
                  y_from <- y - unit(as.numeric(h) / 2, "npc")
                  y_to <- y + unit(as.numeric(h) / 2, "npc")
                  grid.polyline(x = c(x_from, x_to, x_from, x_to), 
                             y = c(y_from,y_from,y_to , y_to), 
                             gp = gpar(col = "#D0D0D0"), 
                             id = rep(1 : (length(x) * 2), 2))
                  if(shownumber)
                  {
                    grid.text(sprintf("%.2f", v), x, y, 
                              gp = gpar(fontsize = 3))
                  }
                  
                }, col=c("blue", "white", "red"),show_heatmap_legend = F,
                border = F)
  
  draw(ht, annotation_legend_list = pd)
  if(!is.null(split_col)){
    for(i in 1 : length(unique(split_col))){
      for(j in 1 : ifelse(is.null(restore), 2, 3)){
        if(j != 1){
          decorate_heatmap_body(heatmap = title, row_slice = j, column_slice = i,
                                {
                                  grid.lines(x = unit(c(0, 1), "npc"),
                                             y=unit(c(1, 1), "npc"), 
                                             gp = gpar(lwd = 2))
                                })
        }
        decorate_heatmap_body(heatmap = title, row_slice = j, column_slice = i,
                              {
                                grid.lines(x = unit(c(0, 0), "npc"),
                                           y=unit(c(0, 1), "npc"), 
                                           gp = gpar(lwd = 1))
                              })
        decorate_heatmap_body(heatmap = title, row_slice = j, column_slice = i,
                              {
                                grid.lines(x = unit(c(1, 1), "npc"),
                                           y=unit(c(0, 1), "npc"), 
                                           gp = gpar(lwd = 1))
                              })
        # }
      }
      
    }
    lens <- unique(split_col)
    
    show_name <- lapply(unique_split, convertGene)
    
    
    for(i in 1:length(unique(split_col))) {
      len <- sum(split_col == lens[i])
      decorate_annotation("group", slice = i, {
        tg <- richtext_grob(gt_render(show_name[i]), 
                            rot = 45, 
                            x = unit(0.5, "npc"),
                            y=unit(0, "npc"), hjust = 0)
        grid.draw(tg)
        invisible(tg)
      })
      if(len == 1){
        decorate_annotation("line", slice = i, {
          grid.lines(x = unit(c(1 / len / 2, 1 / len / 2), "npc"), 
                     y=unit(c(0, 1), "npc"))
        })
      }
      else{
        decorate_annotation("line", slice = i, {
          grid.lines(x = unit(c(1 / len / 2, 1 / len / 2), "npc"), 
                     y=unit(c(0, 1), "npc"))
          grid.lines(x = unit(c(1 / len / 2, 1 - 1 / len / 2), "npc"), 
                     y=unit(c(1, 1), "npc"))
          grid.lines(x = unit(c(1 - 1 / len / 2, 1 - 1 / len / 2), "npc"), 
                     y=unit(c(0, 1), "npc"))
        })
      }
      
    }
  }
}




plotStat_AA <- function(stat, region = 5, showname = F, title = "",
                        shownumber = F, split_col = NULL, 
                        convertGene = function(x){x}){
  library(ComplexHeatmap)
  library(patchwork)
  library(gridtext)
  if(!is.null(split_col)){
    showname <- F
  }
  sample_name <- names(stat)
  region_keep <- 0 : c(region - 1)
  total_stat <- lapply(stat, function(y){
    long_change <- sum(y[y[,1] >= region, 2])
    long_res <- data.frame(aa = paste0(">=", region), 
                           Counts = long_change)
    tmp_region <- region_keep[!region_keep %in% y[,1]]
    y <- y[y[,1] %in% region_keep,]
    if(length(tmp_region) != 0){
      y <- data.frame(rbind(y, data.frame(aa = tmp_region, Counts = 0)))
      y <- y[order(y[,1], decreasing = F),]
    }
    y <- data.frame(rbind(y, long_res))
    colnames(y) <- c("AA","Counts")
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
  for(i in 1 : ncol(total_stat_pct)){
    sum_value <- max(sum(total_stat_pct[,i]), 1)
    total_stat_pct[,i] <- total_stat_pct[,i] / sum_value
  }
  
  total_plot <- total_stat_pct
  colnames(total_stat) <- sample_name
  colnames(total_plot) <- sample_name
  library(circlize)
  aa_color <- colorRamp2(c(0, max(unlist(total_plot))), c("white", "#E13C32"))
  lgd_aa = Legend(col_fun = aa_color, title = "AA change")
  pd = packLegend(lgd_aa)
  unique_split <- unique(split_col)
  if(is.null(split_col)){
    heatanno <- NULL
  }else{
    heatanno <- 
      HeatmapAnnotation(group = anno_empty(border = F, height = unit(30, "mm")),
                        "line" = 
                          anno_empty(border = FALSE, height = unit(2, "mm"))
      )
  }
  
  if(!is.null(split_col)){
    split_col <- factor(split_col, levels = unique(split_col))
  }
  total_plot <- as.matrix(total_plot)
  ht <- Heatmap(total_plot, 
                name = title,
                column_split = split_col,
                column_title=NULL,
                show_column_names = showname, 
                cluster_rows = F, cluster_columns = F, 
                column_title_rot = 45,
                column_names_side = "top",
                top_annotation = heatanno,
                layer_fun = function(j, i, x, y, w, h, fill) {
                  v = pindex(total_plot, i, j) 
                  
                  # `j` here is also a vector with the same length of `v`
                  col = aa_color(v)
                  col2 = "white"
                  if(!is.null(split_col)){
                    col2 = col
                  }
                  grid.rect(x, y, w, h, gp = gpar(fill = col, col = col2))
                  if(shownumber)
                  {
                    grid.text(sprintf("%.2f", v), x, y, 
                              gp = gpar(fontsize = 3))
                  }
                },, col=c("blue", "white", "red"),show_heatmap_legend = F,
                border = F)
  
  draw(ht, annotation_legend_list = pd)
  if(!is.null(split_col)){
    for(i in 1 : length(unique(split_col))){
      
      decorate_heatmap_body(heatmap = title, row_slice = 1, column_slice = i,
                            {
                              grid.lines(x = unit(c(0, 0), "npc"),
                                         y=unit(c(0, 1), "npc"))
                            })
      decorate_heatmap_body(heatmap = title, row_slice = 1, column_slice = i,
                            {
                              grid.lines(x = unit(c(1, 1), "npc"),
                                         y=unit(c(0, 1), "npc"))
                            })
    }
    lens <- unique(split_col)
    
    show_name <- lapply(unique_split, convertGene)
    
    
    for(i in 1:length(unique(split_col))) {
      len <- sum(split_col == lens[i])
      decorate_annotation("group", slice = i, {
        tg <- richtext_grob(gt_render(show_name[i]), 
                            rot = 45, 
                            x = unit(0.5, "npc"),
                            y=unit(0, "npc"), hjust = 0)
        grid.draw(tg)
        invisible(tg)
      })
      if(len == 1){
        decorate_annotation("line", slice = i, {
          grid.lines(x = unit(c(1 / len / 2, 1 / len / 2), "npc"), 
                     y=unit(c(0, 1), "npc"))
        })
      }
      else{
        decorate_annotation("line", slice = i, {
          grid.lines(x = unit(c(1 / len / 2, 1 / len / 2), "npc"), 
                     y=unit(c(0, 1), "npc"))
          grid.lines(x = unit(c(1 / len / 2, 1 - 1 / len / 2), "npc"), 
                     y=unit(c(1, 1), "npc"))
          grid.lines(x = unit(c(1 - 1 / len / 2, 1 - 1 / len / 2), "npc"), 
                     y=unit(c(0, 1), "npc"))
        })
      }
      
    }
  }
  return(list(Counts = total_stat, Pct = total_plot))
}

