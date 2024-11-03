setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep1-replace-by0829data/")
samples <- list.files(pattern = "^CRI")
total_insert_stat <- list()
for(sample in samples){
  setwd(sample)
  total_insert_stat[[sample]] <- list()
  filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
  edit_table <- read.table(filename, 
                             sep="\t", header = T, comment.char = "")
  total_insert_stat[[sample]] <- edit_table
  setwd("..")
  rm(edit_table)
}
rm(sample, samples)
old_rep1_edit <- total_insert_stat

setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep2-replace-by0829data/")
samples <- list.files(pattern = "^CRI")
total_insert_stat <- list()
for(sample in samples){
  setwd(sample)
  total_insert_stat[[sample]] <- list()
  filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
  edit_table <- read.table(filename, 
                           sep="\t", header = T, comment.char = "")
  total_insert_stat[[sample]] <- edit_table
  setwd("..")
  rm(edit_table)
}
rm(sample, samples)
old_rep2_edit <- total_insert_stat

setwd("~/data/project/ear_project/gene_therapy_ll/Previews/Rep3-replace-by0829data/")
samples <- list.files(pattern = "^CRI")
total_insert_stat <- list()
for(sample in samples){
  setwd(sample)
  total_insert_stat[[sample]] <- list()
  filename <- list.files(pattern = "Alleles_frequency_table_around_sgRNA.*txt$")
  edit_table <- read.table(filename, 
                           sep="\t", header = T, comment.char = "")
  total_insert_stat[[sample]] <- edit_table
  setwd("..")
  rm(edit_table)
}
rm(sample, samples)
old_rep3_edit <- total_insert_stat

save(old_rep1_edit, old_rep2_edit, old_rep3_edit, file="../Result/old_rep123_edit.rda")

