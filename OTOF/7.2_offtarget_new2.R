setwd("~/data/project/ear_project/gene_therapy_qsw/batch4/fastq/")
samples <- list.files(pattern = "CRISPRessoPooled_on")[seq(1,6,2)]
sgs <- paste0("OFF", 1 : 6)
offtarget_edit_table <- list()
for(sample in samples){
  setwd(sample)
  offtarget_edit_table[[sample]] <- list()
  for(sg in sgs){
    setwd(paste0("CRISPResso_on_",sg))
    filename <- list.files(pattern = "Alleles_frequency_table_around.*txt$")
    filename <- filename[which.max(str_length(filename))]
    edit_table <- read.table(filename, 
                             sep="\t", header = T, comment.char = "")
    offtarget_edit_table[[sample]][[sg]] <- edit_table
    setwd("..")
  }
  setwd("..")
}


offtarget_indel_pct <- lapply(offtarget_edit_table, function(x){
  lapply(names(x), function(n){
    y <- x[[n]]
    res <- data.frame(sample = n, total_reads = sum(y[,7]), 
                      indel_reads = sum(y[y$n_deleted + y$n_inserted != 0,7]))
    res$pct <- res$indel_reads / res$total_reads
    res
  })
})
offtarget_indel_pct <- data.frame(do.call(rbind, lapply(names(offtarget_indel_pct), function(n){
  x <- offtarget_indel_pct[[n]]
  x <- data.frame(do.call(rbind, x))
  x$batch <- n
  x
})))
offtarget_indel_pct$batch <- str_remove(offtarget_indel_pct$batch, "CRISPRessoPooled_on_")

setwd("~/data/project/ear_project/gene_therapy_qsw/batch5/fastq/")
samples <- list.files(pattern = "CRISPRessoPooled_on")[seq(1,6,2)]
sgs <- paste0("OFF", 1 : 6)
offtarget_edit_table <- list()
for(sample in samples){
  setwd(sample)
  offtarget_edit_table[[sample]] <- list()
  for(sg in sgs){
    setwd(paste0("CRISPResso_on_",sg))
    filename <- list.files(pattern = "Alleles_frequency_table_around.*txt$")
    filename <- filename[which.max(str_length(filename))]
    edit_table <- read.table(filename, 
                             sep="\t", header = T, comment.char = "")
    offtarget_edit_table[[sample]][[sg]] <- edit_table
    setwd("..")
  }
  setwd("..")
}


offtarget2_indel_pct <- lapply(offtarget_edit_table, function(x){
  lapply(names(x), function(n){
    y <- x[[n]]
    res <- data.frame(sample = n, total_reads = sum(y[,7]), 
                      indel_reads = sum(y[y$n_deleted + y$n_inserted != 0,7]))
    res$pct <- res$indel_reads / res$total_reads
    res
  })
})
offtarget2_indel_pct <- data.frame(do.call(rbind, lapply(names(offtarget2_indel_pct), function(n){
  x <- offtarget2_indel_pct[[n]]
  x <- data.frame(do.call(rbind, x))
  x$batch <- n
  x
})))
offtarget2_indel_pct$batch <- str_remove(offtarget2_indel_pct$batch, "CRISPRessoPooled_on_")


setwd("~/data/project/ear_project/gene_therapy_qsw/batch3/fastq/")
samples <- list.files(pattern = "CRISPResso_on_7L")[seq(1, 6, 2)]
otof_result <- list()
for(sample in samples){
  setwd(sample)
  filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
  filename <- filename[which.max(str_length(filename))]
  edit_table <- read.table(filename, 
                           sep="\t", header = T, comment.char = "")
  otof_result[[str_remove(sample, "CRISPResso_on_")]] <- edit_table
  setwd("..")
}
otof_indel_pct <-lapply(names(otof_result), function(n){
  y <- otof_result[[n]]
  res <- data.frame(sample = n, total_reads = sum(y[,7]), 
                    indel_reads = sum(y[y$n_deleted + y$n_inserted != 0,7]))
  res$pct <- res$indel_reads / res$total_reads
  res
})

otof_indel_pct <- data.frame(do.call(rbind, otof_indel_pct))
otof_indel_pct$batch <- paste0("Rep", 1 : 3)
otof_indel_pct$sample <- "7L"
plot_data <- data.frame(rbind(otof_indel_pct, offtarget_indel_pct))
offtarget2_indel_pct$sample <- paste0(offtarget_indel_pct$sample, "-NoViru")
plot_data <- data.frame(rbind(plot_data, offtarget2_indel_pct))

library(ggbreak)
pdf("~/Nutstore Files/Tobin/Merged1NT/OTOF_offtarger_result.pdf", width = 16, height = 5)
ggplot(plot_data, aes(x = sample, y = pct * 100)) + geom_boxplot() + 
  geom_point(position = position_jitter()) + ylab("Indel Percent in Total Reads")+
  scale_y_break(breaks = c(2, 65), space = 0)+
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  theme_classic2()
dev.off()
openxlsx::write.xlsx(plot_data, file="~/Nutstore Files/Tobin/Merged1NT/OTOF_offtarget_result.xlsx", colNames=T, rowNames=F)
