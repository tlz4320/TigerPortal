setwd('~/indelphi_res/')
indelphi_118_res <- list()
indelphi_file <- read.csv("~/118_for_indelphi.csv")
for(file in indelphi_file$id2){
  indelphi_118_res[[file]] <- read.table(paste0(file, ".csv"), sep = ",", header = T)
}
indelphi_118_res <- lapply(indelphi_118_res, function(x){
  x$Length[x$Category == "del"] <- -x$Length[x$Category == "del"]
  res <- lapply(split(x$Predicted.frequency, x$Length), sum)
  res <- data.frame(indel_size = names(res), fq = unlist(res))
  res$Pct <- res$fq
  res <- res[gtools::mixedorder(res$indel_size),]
  res
})

setwd('~/forecast_result/')
forecast_118_res <- list()
forcast_file <- read.table("~/118_for_forecast.txt")
for(file in forcast_file$V1){
  forecast_118_res[[file]] <- read.table(paste0(file, "_predictedindelsummary.txt"), 
                                         sep = "\t", header = T)
}
forecast_118_res <- lapply(forecast_118_res, function(x){
  edit <- unlist(lapply(x[,1], function(y){
    unlist(strsplit(y, "[_]"))[1]
  }))
  edit <- data.frame(type = str_sub(edit, 1, 1), len = str_sub(edit, 2), counts = x[,3])
  edit$len <- as.numeric(edit$len)
  edit$len[edit$type == "D"] <- -edit$len[edit$type == "D"]
  res <- lapply(split(edit$counts, edit$len), sum)
  res <- data.frame(indel_size = names(res), fq = unlist(res))
  res$Pct <- res$fq / sum(res$fq) * 100
  res <- res[gtools::mixedorder(res$indel_size),]
  res
})

load("~/data/project/ear_project/gene_therapy_ll/Previews/total_used_seq_split.rda")
total_ids <- str_remove(total_used_seq_tissue$result, "[-].*")
total_used_seq_tissue$id <- total_ids
inframe_result_indelphi <- lapply(total_ids, function(id){
  tmp <- indelphi_118_res[[id]]
  info <- total_used_seq_tissue[total_used_seq_tissue$id == id,]
  total_reads <- sum(tmp$fq)
  
  if(info$mutType == "d"){
    tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 1, ]
    tmp$inframe <- as.numeric(tmp$indel_size) - 1
  } else {
    tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 2, ]
    tmp$inframe <- as.numeric(tmp$indel_size) + 1
  }
  
  
  tmp <- lapply(split(tmp$fq, tmp$inframe), sum)
  tmp <- data.frame(indel = as.numeric(names(tmp)), fq = unlist(tmp))
  tmp$pct <- tmp$fq / total_reads * 100
  res <- tmp
  res <- res[res$pct != 0,]
  res$aa <- as.numeric(res$indel) / 3
  res
})
names(inframe_result_indelphi) <- total_ids
inframe_result_indelphi <- lapply(inframe_result_indelphi, function(x){
  if(nrow(x) == 0)
    return(x)
  x <- x[order(abs(x$aa)),]
  x2 <- x
  x2$aa[abs(x2$aa) > 5] <- "Others"
  x2$fq[abs(x$aa) > 5] <- sum(x2$fq[abs(x$aa) > 5])
  x2$pct[abs(x$aa) > 5] <- sum(x2$pct[abs(x$aa) > 5])
  if(any(abs(x$aa) > 5)){
    x <- x2[c(which(abs(x$aa) <= 5), which(abs(x$aa) > 5)[1]),] 
  }
  x
})


inframe_result_indelphi <- lapply(total_ids, function(id){
  tmp <- inframe_result_indelphi[[id]]
  if(nrow(tmp) == 0){
    tmp <- data.frame(indel=0,
                      fq=0,
                      pct=0,
                      aa=0, 
                      pct2 = 0)
  }
  tmp$pct2 <- tmp$pct / sum(tmp$pct) * 100
  tmp$id <- id
  tmp
})
names(inframe_result_indelphi) <- total_ids
inframe_result_indelphi <- data.frame(do.call(rbind, inframe_result_indelphi))



inframe_result_forecast <- lapply(total_ids, function(id){
  tmp <- forecast_118_res[[id]]
  info <- total_used_seq_tissue[total_used_seq_tissue$id == id,]
  total_reads <- sum(tmp$fq)
  
  if(info$mutType == "d"){
    tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 1, ]
    tmp$inframe <- as.numeric(tmp$indel_size) - 1
  } else {
    tmp <- tmp[as.numeric(tmp$indel_size) %% 3 == 2, ]
    tmp$inframe <- as.numeric(tmp$indel_size) + 1
  }
  
  tmp <- lapply(split(tmp$fq, tmp$inframe), sum)
  tmp <- data.frame(indel = as.numeric(names(tmp)), fq = unlist(tmp))
  tmp$pct <- tmp$fq / total_reads * 100
  res <- tmp
  res <- res[res$pct != 0,]
  res$aa <- as.numeric(res$indel) / 3
  res
})
names(inframe_result_forecast) <- total_ids
inframe_result_forecast <- lapply(inframe_result_forecast, function(x){
  if(nrow(x) == 0)
    return(x)
  x <- x[order(abs(x$aa)),]
  x2 <- x
  x2$aa[abs(x2$aa) > 5] <- "Others"
  x2$fq[abs(x$aa) > 5] <- sum(x2$fq[abs(x$aa) > 5])
  x2$pct[abs(x$aa) > 5] <- sum(x2$pct[abs(x$aa) > 5])
  if(any(abs(x$aa) > 5)){
    x <- x2[c(which(abs(x$aa) <= 5), which(abs(x$aa) > 5)[1]),] 
  }
  x
})


inframe_result_forecast <- lapply(total_ids, function(id){
  tmp <- inframe_result_forecast[[id]]
  if(nrow(tmp) == 0){
    tmp <- data.frame(indel=0,
                      fq=0,
                      pct=0,
                      aa=0, 
                      pct2 = 0)
  }
  tmp$pct2 <- tmp$pct / sum(tmp$pct) * 100
  tmp$id <- id
  tmp
})
names(inframe_result_forecast) <- total_ids
inframe_result_forecast <- data.frame(do.call(rbind, inframe_result_forecast))


print(load("~/Nutstore Files/Tobin/Previous/inframe_spec_result_cell_tissue.rda"))
inframe_result_processed$id2 <- str_remove(inframe_result_processed$id, "[-].*")
cell_inframe_result_processed$id2 <- str_remove(cell_inframe_result_processed$id, "[-].*")



total_used_seq_tissue$tissueInframe <- 0
total_used_seq_tissue$cellInframe <- 0
total_used_seq_tissue$forecastInframe <- 0
total_used_seq_tissue$indelphiInframe <- 0
total_used_seq_tissue$tissue0AA <- 0
total_used_seq_tissue$cell0AA<- 0
total_used_seq_tissue$forecast0AA <- 0
total_used_seq_tissue$indelphi0AA <- 0

for(i in 1 : nrow(total_used_seq_tissue)){
  id <- total_used_seq_tissue$id[i]
  total_used_seq_tissue$tissueInframe[i] <- sum(inframe_result_processed$pct[inframe_result_processed$id2 == id])
  total_used_seq_tissue$cellInframe[i] <- sum(cell_inframe_result_processed$pct[cell_inframe_result_processed$id2 == id])
  total_used_seq_tissue$forecastInframe[i] <- sum(inframe_result_forecast$pct[inframe_result_forecast$id == id])
  total_used_seq_tissue$indelphiInframe[i] <- sum(inframe_result_indelphi$pct[inframe_result_indelphi$id == id])
  total_used_seq_tissue$tissue0AA[i] <- sum(inframe_result_processed$pct2[inframe_result_processed$id2 == id & 
                                                                            inframe_result_processed$aa == "0"])
  total_used_seq_tissue$cell0AA[i] <- sum(cell_inframe_result_processed$pct2[cell_inframe_result_processed$id2 == id & 
                                                                               cell_inframe_result_processed$aa == "0"])
  total_used_seq_tissue$forecast0AA[i] <- sum(inframe_result_forecast$pct2[inframe_result_forecast$id == id & 
                                                                             inframe_result_forecast$aa == "0"])
  total_used_seq_tissue$indelphi0AA[i] <- sum(inframe_result_indelphi$pct2[inframe_result_indelphi$id == id & 
                                                                             inframe_result_indelphi$aa == "0"])
  
}



del1_info <- total_used_seq_tissue[total_used_seq_tissue$mutType == "d",]
plot_data <- hist(del1_info$tissueInframe, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)
plot_data$color <- "#B46F2877"
plot_data$color[plot_data$inframePct > 30] <-"#B46F28FF" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p1 <- ggplot(plot_data) + geom_bar(aes(y = x, x = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color, 
                                   orientation = "y") +
  geom_text(x = -15, y = 60, label = anno_info, data = data.frame(a=1)) +
  ylab("") + xlab("") +
  scale_y_continuous(position = "right") +
  scale_x_reverse(breaks = c(0, 5, 10, 15, 20, 25), limits = c(25, 0), expand = c(0,0)) +
  # scale_y_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), limits = c(0, 100)) +
  theme_bw() +
  theme(plot.margin =  unit(c(0,-0.4,0,1), "cm"), 
        axis.text.y = element_blank(),
        panel.grid = element_blank())
# dev.off()

total_plot_data <- list()
for(name in c("cellInframe", "forecastInframe", "indelphiInframe")){
  plot_data <- hist(del1_info[,name], breaks = seq(from=0, to=105, by=5))
  total_plot_data[[name]] <- data.frame(inframePct = plot_data$breaks[-1], 
                                        counts = plot_data$counts, type = name)
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data$type <- str_remove(total_plot_data$type, "Inframe")

max(plot_data$counts)
total_plot_data$alpha <- 0.45
total_plot_data$alpha[total_plot_data$inframePct > 30] <- 1 #B46F28
total_plot_data$x <- total_plot_data$inframePct - 2.5
total_info <- list()
for(name in unique(total_plot_data$type)){
  tmp <- total_plot_data[total_plot_data$type == name,]
  anno_info <- lapply(c(30, 50, 60, 90), function(pct){
    sum(tmp$counts[tmp$x > pct])
  })
  anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
  anno_info <- paste0(name, " >30:", anno_info$counts[1])
  
  
  
  total_info[[name]] <- anno_info
}
total_info <- data.frame(do.call(rbind, total_info))
total_info <- paste0(total_info[,1], collapse = "\n")


p2 <- ggplot(total_plot_data) + geom_bar(aes(y = x, x = counts, fill = type),
                                         color = "black",
                                         position = position_dodge(),
                                         orientation = "y", 
                                         stat="identity", alpha = total_plot_data$alpha) +
  geom_text(x = 15, y = 50, label = total_info, data = data.frame(a=1)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25), limits = c(0, 25), expand = c(0,0)) +
  scale_fill_manual(values = c("#A7B14E", "#4EA7B1", "#B14EA7")) + 
  ylab("Counts") + xlab("") + 
  theme_bw() + 
  theme(plot.margin =  unit(c(0,1,0,0), "cm"), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 0.5),
        panel.grid = element_blank())
pdf("~/Nutstore Files/Tobin/Previous/inframe_histgram_delMut_addAll_v2.pdf", width = 5, height = 5)
print(ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T))
dev.off()










####ins
ins1_info <- total_used_seq_tissue[total_used_seq_tissue$mutType == "i",]
plot_data <- hist(ins1_info$tissueInframe, breaks = seq(from=0, to=105, by=5))
plot_data <- data.frame(inframePct = plot_data$breaks[-1], 
                        counts = plot_data$counts)
max(plot_data$counts)
plot_data$color <- "#0A502877"
plot_data$color[plot_data$inframePct > 30] <-"#0A5028" 
plot_data$x <- plot_data$inframePct - 2.5
anno_info <- lapply(c(30, 50, 60, 90), function(pct){
  sum(plot_data$counts[plot_data$x > pct])
})
anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
anno_info <- paste(c(paste0("Total Counts:", sum(plot_data$counts)), paste("Counts of >", anno_info$pct, ":", anno_info$counts, sep="")), collapse = "\n")

p1 <- ggplot(plot_data) + geom_bar(aes(y = x, x = counts),
                                   color = "black",
                                   stat="identity", fill = plot_data$color, 
                                   orientation = "y") +
  geom_text(x = -15, y = 60, label = anno_info, data = data.frame(a=1)) +
  ylab("") + xlab("") +
  scale_y_continuous(position = "right") +
  scale_x_reverse(breaks = c(0, 5, 10, 15, 20, 25), limits = c(25, 0), expand = c(0,0)) +
  # scale_y_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90, 100), limits = c(0, 100)) +
  theme_bw() +
  theme(plot.margin =  unit(c(0,-0.4,0,1), "cm"), 
        axis.text.y = element_blank(),
        panel.grid = element_blank())
# dev.off()

total_plot_data <- list()
for(name in c("cellInframe", "forecastInframe", "indelphiInframe")){
  plot_data <- hist(ins1_info[,name], breaks = seq(from=0, to=105, by=5))
  total_plot_data[[name]] <- data.frame(inframePct = plot_data$breaks[-1], 
                                        counts = plot_data$counts, type = name)
}
total_plot_data <- data.frame(do.call(rbind, total_plot_data))
total_plot_data$type <- str_remove(total_plot_data$type, "Inframe")

max(plot_data$counts)
total_plot_data$alpha <- 0.45
total_plot_data$alpha[total_plot_data$inframePct > 30] <- 1 #B46F28
total_plot_data$x <- total_plot_data$inframePct - 2.5
total_info <- list()
for(name in unique(total_plot_data$type)){
  tmp <- total_plot_data[total_plot_data$type == name,]
  anno_info <- lapply(c(30, 50, 60, 90), function(pct){
    sum(tmp$counts[tmp$x > pct])
  })
  anno_info <- data.frame(pct = c(30, 50, 60, 90), counts = unlist(anno_info))
  anno_info <- paste0(name, " >30:", anno_info$counts[1])
  
  
  
  total_info[[name]] <- anno_info
}
total_info <- data.frame(do.call(rbind, total_info))
total_info <- paste0(total_info[,1], collapse = "\n")


p2 <- ggplot(total_plot_data) + geom_bar(aes(y = x, x = counts, fill = type),
                                         color = "black",
                                         position = position_dodge(),
                                         orientation = "y", 
                                         stat="identity", alpha = total_plot_data$alpha) +
  geom_text(x = 15, y = 50, label = total_info, data = data.frame(a=1)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25), limits = c(0, 25), expand = c(0,0)) +
  scale_fill_manual(values = c("#A7B14E", "#4EA7B1", "#B14EA7")) + 
  ylab("Counts") + xlab("") + 
  theme_bw() + 
  theme(plot.margin =  unit(c(0,1,0,0), "cm"), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = 0.5),
        panel.grid = element_blank())
pdf("~/Nutstore Files/Tobin/Previous/inframe_histgram_insMut_addAll_v2.pdf", width = 5, height = 5)
print(ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T))
dev.off()
