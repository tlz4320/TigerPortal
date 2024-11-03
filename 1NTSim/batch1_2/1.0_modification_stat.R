setwd("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599/")
setwd("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/")
samples <- list.dirs(recursive = F)
sgRNA_tmp <- read.table("sg_info.txt")
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
openxlsx::write.xlsx(list(Res = output_res), "../Result/batch2_edit_stat1.xlsx",
                     colNames=T, rowNames=F)



pdf("../Result/batch2_edit_stat2.pdf", width = 25, height = 8)
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

