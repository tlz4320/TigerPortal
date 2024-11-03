setwd("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599/")
samples <- list.dirs(recursive = F)
sgRNA_tmp <- read.table("sg_info.txt")
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
openxlsx::write.xlsx(list(Res = output_res), "../Result/batch2_mapping_stat1.xlsx",
                     colNames=T, rowNames=F)
rm(output_res, output_res2)
pdf("../Result/batch2_mapping_stat2.pdf", width = 25, height = 8)
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



###unused
total_reads <- read.table("clean_data/reads_count.txt")
total_reads$V1 <- str_remove(total_reads$V1, "_merged.fq.gz")
total_reads$V2 <- total_reads$V2 / 4
total_bc_reads <- lapply(total_sample_stat, function(x){
  c(sum(x[,1]), sum(x[,2]))
})
total_bc_reads <- data.frame(do.call(rbind, total_bc_reads))
total_bc_reads$sample <- str_remove(rownames(total_bc_reads), "./")
total_bc_reads <- merge(total_reads, total_bc_reads, by.x="V1", by.y = "sample")
colnames(total_bc_reads) <- c("Sample", "Total", "haveBC", "Mapped")
total_bc_reads_pct <- total_bc_reads
for(i in ncol(total_bc_reads_pct) : 2){
  total_bc_reads_pct[,i] <- total_bc_reads_pct[,i] / total_bc_reads_pct[,2]
}
openxlsx::write.xlsx(list(Counts = total_bc_reads, Percent = total_bc_reads_pct), 
                     file="Barcode_stat.xlsx", rowNames=F, colNames=T)

plot_data <- reshape2::melt(total_bc_reads, id.vars=c("Sample"))
colnames(plot_data) <- c("Sample", "Type", "Counts")

plot_data_pct <- plot_data
plot_data_pct <- data.frame(do.call(rbind, lapply(split(plot_data, plot_data$Sample), function(x){
  x[,3] <- x[,3] / max(x[,3])
  x
})))

pdf("Result/sample_barcode_stat_max2.pdf", width = 12, height = 8)
ggplot(plot_data) + geom_bar(aes(x = Sample, y = Counts, fill = Type), 
                             stat = "identity", position = "dodge") + 
  geom_text(aes(x = Sample, y = Counts, group = Type, label = scales::scientific_format()(Counts)), 
            size = 4, position = position_dodge(width = 0.8)) + 
  theme_bw() + xlab("") + ylab("Counts") +  
  theme(panel.grid = element_blank(),text = element_text(size=20), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggsci::scale_fill_nejm()

ggplot(plot_data_pct) + geom_bar(aes(x = Sample, y = Counts, fill = Type), 
                                 stat = "identity", position = "dodge") + 
  geom_text(aes(x = Sample, y = Counts, group = Type, label = round(digits = 2, Counts)), 
            size = 5, position = position_dodge(width = 0.8)) + 
  theme_bw() + xlab("") + ylab("Percent") +  
  theme(panel.grid = element_blank(),text = element_text(size=20), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggsci::scale_fill_nejm()
dev.off()
