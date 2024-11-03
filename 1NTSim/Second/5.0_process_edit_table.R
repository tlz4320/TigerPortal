#在2.0_total_First_pop_peak_indel1.R
####Read batch1_1 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/Second/Batch1/240318-A00599B/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
sgRNA_tmp <- read.table("../sg_info.txt")
sgRNA_tmp <- sgRNA_tmp[c(19:30,32,33,35,36),]
total_edit_table <- list()
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
      setwd(paste0("CRISPResso_on_",sg))
      filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
      filename <- filename[which.max(str_length(filename))]
      edit_table <- read.table(filename, 
                               sep="\t", header = T, comment.char = "")
      total_edit_table[[sample]][[sg]] <- edit_table
      setwd("..")

  }
  setwd("..")
}
###Read batch1_2 edit table
setwd("~/data/project/ear_project/gene_therapy_ll/Second/Batch2/240426-A00599B/")
samples <- list.dirs(recursive = F)
samples <- str_remove(samples, "./")
sgRNA_tmp <- read.table("../sg_info.txt")
sgRNA_tmp <- sgRNA_tmp[c(1 : 18, 37 : 160), ]
rm_list <- c("169","170","177","178","207","208","209","210","221","224","225","240","250","256","258","268","279","281","284","285","299","311")
sgRNA_tmp <- sgRNA_tmp[!sgRNA_tmp$V1 %in% paste0("Sg_24_", rm_list),]
for(sample in samples){
  setwd(sample)
  total_edit_table[[sample]] <- list()
  for(sg in sgRNA_tmp$V1){
    setwd(paste0("CRISPResso_on_",sg))
    filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
    filename <- filename[which.max(str_length(filename))]
    edit_table <- read.table(filename, 
                             sep="\t", header = T, comment.char = "")
    total_edit_table[[sample]][[sg]] <- edit_table
    setwd("..")
  }
  setwd("..")
}
rm(sample, samples, sg, sgRNA_tmp, edit_table)


sgRNA <- load("~/data/project/ear_project/gene_therapy_ll/Result/second_sgCmp.rda")
sgRNA <- sgCmp2



tmp <- unique(unlist(lapply(total_edit_table, names)))
total_edit_table_rev <- lapply(tmp, function(x){
  res <- list()
  for(name in names(total_edit_table)){
    if(x %in% names(total_edit_table[[name]])){
      res[[name]] <- total_edit_table[[name]][[x]]
    }
  }
  res
})

names(total_edit_table_rev) <- tmp
total_edit_table_second <- total_edit_table
total_edit_table_second_rev <- total_edit_table_rev
save(total_edit_table_second, total_edit_table_second_rev, 
     file="~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_second.rda")
