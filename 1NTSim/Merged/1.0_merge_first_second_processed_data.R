library(Biostrings)
library(stringr)
#这个文件就是把两次设计的结果都合并在一起，方便后续分析不需要再修改过多的代码
# source("~/code/earProject/gene_therapy/1NTSim/batch1_merge/0.2_AlignSeq.R")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_sgCmp.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/second_sgCmp.rda")
colnames(sgCmp2) <- colnames(sgCmp)
sgCmp <- data.frame(rbind(sgCmp, sgCmp2))
rm(sgCmp2)
save(sgCmp, file="~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_sgRNA_cmp_reDiff.txt", sep="\t")
bulge_pos2 <- read.table("~/data/project/ear_project/gene_therapy_ll/Second_sgRNA_cmp_reDiff.txt", sep="\t")

bulge_pos <- data.frame(rbind(bulge_pos, bulge_pos2))
rm(bulge_pos2)
ToNX::write_tb(bulge_pos, file="~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt")
load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_second.rda")

for(name in names(total_edit_table_second_rev)){
  total_edit_table_rev[[name]] <- total_edit_table_second_rev[[name]]
}
for(name in names(total_edit_table_second)){
  total_edit_table[[name]] <- total_edit_table_second[[name]]
}
save(total_edit_table, total_edit_table_rev, file="~/data/project/ear_project/gene_therapy_ll/Result/total_edit_table_first_second.Rda")
