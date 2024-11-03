#Mapping Stat
setwd("~/data/project/ear_project/gene_therapy_ll/Second/Batch1/240318-A00599B/")
samples <- list.dirs(recursive = F)
sgRNA_tmp <- read.table("../sg_info.txt")
sgRNA_tmp <- sgRNA_tmp[-c(11, 12, 37 : 40),]
total_sample_stat <- list()
for(sample in samples){
  setwd(sample)
  total_sample_stat[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    setwd(paste0("CRISPResso_on_",sg))
    mapping_rate <- read.table("CRISPResso_mapping_statistics.txt", 
                               sep="\t", header = T)
    total_sample_stat[[sample]][[sg]] <- c(unlist(mapping_rate[1,c(1,3)]))
    setwd("..")
  }
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples, sg)
total_sample_stat <- lapply(total_sample_stat, function(x){
  x <- data.frame(do.call(rbind, x))
  x$unmapped <- x[,1] - x[,2]
  x$sg <- rownames(x)
  colnames(x) <- c("Total", "Mapped", "Unmapped", "ID")
  x
})

total_sample_stat_pct <- lapply(total_sample_stat, function(x){
  for(i in 3 : 1){
    x[,i] <- x[, i] / x[,1]
  }
  x
})

output_res<- data.frame(do.call(rbind, total_sample_stat))
output_res$Sample <- rownames(output_res)
output_res$Sample <- str_remove(str_remove(output_res$Sample, "./"),"[.].*")
output_res <- output_res[gtools::mixedorder(output_res$ID),]
output_res$Sample <- str_remove(output_res$Sample, "CRISPRessoPooled_on_")
output_res <- output_res[,c(4,5,1:3)]

output_res2<- data.frame(do.call(rbind, total_sample_stat_pct))
output_res2$Sample <- rownames(output_res2)
output_res2$Sample <- str_remove(str_remove(output_res2$Sample, "./"),"[.].*")
output_res2 <- output_res2[gtools::mixedorder(output_res2$ID),]
output_res2$Sample <- str_remove(output_res2$Sample, "CRISPRessoPooled_on_")
output_res2 <- output_res2[,c(4,5,1:3)]
colnames(output_res2) <- paste0(colnames(output_res2), "(Pct)")
output_res <- data.frame(cbind(output_res, output_res2[,c(-1,-2)]))
openxlsx::write.xlsx(list(Res = output_res), "../Result/batch1_mapping_stat1.xlsx",
                     colNames=T, rowNames=F)
rm(output_res, output_res2)
pdf("../Result/batch1_mapping_stat.pdf", width = 25, height = 8)
for(name in names(total_sample_stat)){
  plot_data <- total_sample_stat[[name]]
  plot_data <- plot_data[gtools::mixedorder(plot_data$ID),]
  plot_data$ID <- factor(plot_data$ID, levels = plot_data$ID)
  plot_data_tmp <- data.frame(counts= c(plot_data$Mapped, plot_data$Unmapped),
                              ID = c(plot_data$ID, plot_data$ID),
                              type = c(rep(c("Mapped","Unmapped"), c(nrow(plot_data),nrow(plot_data)))))
  p <- ggplot(plot_data) + geom_bar(aes(x = ID, y = Total), stat = "identity") + 
    geom_bar(data = plot_data_tmp, aes(x = ID, y = counts, fill = type), 
             stat = "identity", position = "stack") +
    geom_text(aes(x = ID, y = Total, label = Total), size = 4) + 
    theme_bw() + xlab("") + ylab("Counts") +  
    theme(panel.grid = element_blank(),text = element_text(size=15), 
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    ggsci::scale_fill_nejm() + 
    ggtitle(str_remove(name, "./"))
  print(p)
  rm(p)
  rm(plot_data, plot_data_tmp)
  rm(name)
}
dev.off()



###Modify stat

setwd("~/data/project/ear_project/gene_therapy_ll/Second/Batch1/240318-A00599B/")
samples <- list.dirs(recursive = F)
total_sample_stat <- list()
for(sample in samples){
  setwd(sample)
  total_sample_stat[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    setwd(paste0("CRISPResso_on_",sg))
    mapping_rate <- read.table("CRISPResso_quantification_of_editing_frequency.txt", 
                               sep="\t", header = T)
    total_sample_stat[[sample]][[sg]] <- c(unlist(mapping_rate[1,c(7:12)]))
    setwd("..")
  }
  setwd("..")
  rm(mapping_rate)
}
rm(sample, samples, sg)
total_sample_stat <- lapply(total_sample_stat, function(x){
  x <- data.frame(do.call(rbind, x))
  x$sg <- rownames(x)
  colnames(x) <- c("Unmodified", "Modified", "Discarded", "Insertions","Deletions", "Substitutions", "ID" )
  x
})
total_sample_stat <- lapply(total_sample_stat, function(x){
  x <- x[,c(-3)]
})

#output Excel
output_res<- data.frame(do.call(rbind, total_sample_stat))
output_res$Sample <- rownames(output_res)
output_res$Sample <- str_remove(str_remove(output_res$Sample, "./"),"[.].*")
output_res <- output_res[gtools::mixedorder(output_res$ID),]
output_res$Sample <- str_remove(output_res$Sample, "CRISPRessoPooled_on_")
output_res <- output_res[,c(6,7,1:5)]

total_sample_stat_pct <- lapply(total_sample_stat, function(x){
  total_counts <- unlist(rowSums(x[,-ncol(x)]))
  for(i in 1 : 5){
    x[,i] <- x[, i] / total_counts
  }
  x
})
output_res2<- data.frame(do.call(rbind, total_sample_stat_pct))
output_res2$Sample <- rownames(output_res2)
output_res2$Sample <- str_remove(str_remove(output_res2$Sample, "./"),"[.].*")
output_res2 <- output_res2[gtools::mixedorder(output_res2$ID),]
output_res2$Sample <- str_remove(output_res2$Sample, "CRISPRessoPooled_on_")
output_res2 <- output_res2[,c(6,7,1:5)]
colnames(output_res2) <- paste0(colnames(output_res2), "(Pct)")
output_res <- data.frame(cbind(output_res, output_res2[,c(-1,-2)]))
openxlsx::write.xlsx(list(Res = output_res), "../Result/batch1_edit_stat.xlsx",
                     colNames=T, rowNames=F)



pdf("../Result/batch1_edit_stat.pdf", width = 25, height = 8)
for(name in names(total_sample_stat)){
  plot_data <- total_sample_stat[[name]]
  plot_data <- plot_data[,c(-2)]
  plot_data <- reshape2::melt(plot_data, id.vars= 'ID')
  colnames(plot_data) <- c("ID", "type", "counts")
  plot_data <- data.frame(do.call(rbind, lapply(split(plot_data, plot_data$ID), function(x){
    x$pct <- x$counts/ sum(x$counts)
    x
  })))
  plot_data$ID <- factor(plot_data$ID, levels = gtools::mixedsort(unique(plot_data$ID)))
  
  p <- ggplot(plot_data, aes(x = ID, y = pct, fill = type, , label = round(pct, 2))) + 
    geom_bar(stat = "identity", position = "fill") +
    geom_text(size = 4, position = position_stack(vjust = 1)) + 
    theme_bw() + xlab("") + ylab("Percent") +  
    theme(panel.grid = element_blank(),text = element_text(size=15), 
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    ggsci::scale_fill_nejm() + 
    ggtitle(str_remove(name, "./"))
  print(p)
  rm(p)
  rm(plot_data)
  rm(name)
}
dev.off()

