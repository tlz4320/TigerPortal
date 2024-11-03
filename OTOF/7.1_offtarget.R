setwd("~/data/project/ear_project/gene_therapy_qsw/batch2/fastq/")
sample_info <- data.frame(id = paste0("Q", 1  : 21), 
                          type = rep(c(paste0("Sg", 1 : 6), "Control"), rep(3, 7)))
offtarget_edit_table <- list()
for(sample in sample_info$id){
  setwd(paste0(sample, "_res"))
  setwd(paste0("CRISPResso_on_", sample))
  filename <- list.files(pattern = "Alleles_frequency_table_around.*txt$")
  filename <- filename[which.max(str_length(filename))]
  offtarget_edit_table[[sample]] <- read.table(filename, 
                           sep="\t", header = T, comment.char = "")
  setwd("../..")
}
offtarget_effic <- lapply(offtarget_edit_table, function(x){
  total_reads <- sum(x[,7])
  indel_reads <- sum(x[x$n_deleted != 0 | x$n_inserted != 0, 7])
  indel_reads / total_reads * 100
})
offtarget_effic <- data.frame(id = names(offtarget_effic), indel = unlist(offtarget_effic))
offtarget_effic <- merge(offtarget_effic, sample_info, by="id")

setwd("~/data/project/ear_project/gene_therapy_qsw/mOTOF/Result/")
ontarget_edit_table <- list()
samples <- list.files(pattern = "Rescue")
samples <- samples[-grep("html",samples)]
samples <- samples[-grep("zip",samples)]
for(sample in samples){
  setwd(sample)
  filename <- list.files(pattern = "Alleles_frequency_table_around.*txt$")
  filename <- filename[which.max(str_length(filename))]
  ontarget_edit_table[[sample]] <- read.table(filename, 
                                               sep="\t", header = T, comment.char = "")
  setwd("..")
}
ontarget_effic <- lapply(ontarget_edit_table, function(x){
  total_reads <- sum(x[,7])
  indel_reads <- sum(x[x$n_deleted != 0 | x$n_inserted != 0, 7])
  indel_reads / total_reads * 100
})

ontarget_effic <- data.frame(id = names(ontarget_effic), indel = unlist(ontarget_effic))
ontarget_effic$type <- "Rescue"

plot_data <- data.frame(rbind(ontarget_effic, offtarget_effic))
plot_data$type <- factor(plot_data$type, levels = c("Rescue", paste0("Sg", 1:6), "Control"))
ggplot(plot_data, aes(x = type, y = indel)) + geom_boxplot() + 
  geom_point(position = position_jitter()) + theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
