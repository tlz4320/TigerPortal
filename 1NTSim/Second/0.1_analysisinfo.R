library(Biostrings)
library(ShortRead)
library(parallel)
library(stringr)
library(openxlsx)

###输出Second的信息，同时发现这个部分有些序列很接近的会被错误的pcr出来，所以改用genome序列
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/Second1NT/second_batch_design_total160_replace_bad_fix.xlsx")

sgRNA <- sgRNA[c(1 : 40),]
sgRNA$wt_seqs <- unlist(lapply(1 : nrow(sgRNA), function(i){
  seq <- sgRNA$genomeSeq[i]
  left <- sgRNA$L.Primer[i]
  right <- as.character(reverseComplement(DNAString(sgRNA$R.Primer[i])))
  leftpos <- str_locate(seq, left)[1]
  rightpos <- str_locate(seq, right)[2]
  right <- as.character(reverseComplement(DNAString(sgRNA$barcodeR[i])))
  str_sub(seq, leftpos, rightpos)
}))
sgRNA_tmp <- data.frame(seq = sgRNA$genomeSeq, sg = paste0(sgRNA$sgRNA, sgRNA$sgRNA_NGG))
sgRNA_tmp$id <- paste("Sg", sgRNA$ID_new, sep = "-")
sgRNA_tmp$id <- str_replace_all(sgRNA_tmp$id, "-", "_")
sgRNA_tmp <- sgRNA_tmp[,c(3,1,2)]

ToNX::write_tb(sgRNA_tmp[,c(1,2)], "~/data/project/ear_project/gene_therapy_ll/Second/Batch1/sgRNA_seq40.txt")


library(Biostrings)
library(ShortRead)
library(parallel)
library(stringr)
library(openxlsx)

###输出Second的信息，同时发现这个部分有些序列很接近的会被错误的pcr出来，所以改用genome序列
library(openxlsx)
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/Second1NT/second_batch_design_total160_replace_bad_fix.xlsx")

sgRNA$wt_seqs <- unlist(lapply(1 : nrow(sgRNA), function(i){
  seq <- sgRNA$genomeSeq[i]
  left <- sgRNA$L.Primer[i]
  right <- as.character(reverseComplement(DNAString(sgRNA$R.Primer[i])))
  leftpos <- str_locate(seq, left)[1]
  rightpos <- str_locate(seq, right)[2]
  right <- as.character(reverseComplement(DNAString(sgRNA$barcodeR[i])))
  str_sub(seq, leftpos, rightpos)
}))
sgRNA_tmp <- data.frame(seq = sgRNA$wt_seqs, sg = paste0(sgRNA$sgRNA, sgRNA$sgRNA_NGG))
sgRNA_tmp$id <- paste("Sg", sgRNA$ID_new, sep = "-")
sgRNA_tmp$id <- str_replace_all(sgRNA_tmp$id, "-", "_")
sgRNA_tmp <- sgRNA_tmp[,c(3,1,2)]

ToNX::write_tb(sgRNA_tmp[,c(1,2)], "~/data/project/ear_project/gene_therapy_ll/Second/Batch2/sgRNA_seq.txt")


sgRNA <- read.xlsx("~/Nutstore Files/Tobin/Second1NT/second_batch_design_total160_replace_bad_fix.xlsx")
sgRNA$wt_seqs <- unlist(lapply(1 : nrow(sgRNA), function(i){
  seq <- sgRNA$genomeSeq[i]
  left <- sgRNA$L.Primer[i]
  right <- as.character(reverseComplement(DNAString(sgRNA$R.Primer[i])))
  leftpos <- str_locate(seq, left)[1]
  rightpos <- str_locate(seq, right)[2]
  right <- as.character(reverseComplement(DNAString(sgRNA$barcodeR[i])))
  paste0(sgRNA$barcodeL[i], str_sub(seq, leftpos, rightpos), right)
}))
sgRNA_tmp <- data.frame(seq = sgRNA$wt_seqs, sg = paste0(sgRNA$sgRNA, sgRNA$sgRNA_NGG))
sgRNA_tmp$id <- paste("Sg", sgRNA$ID_new, sep = "-")
sgRNA_tmp$id <- str_replace_all(sgRNA_tmp$id, "-", "_")
sgRNA_tmp <- sgRNA_tmp[,c(3,1,2)]
ToNX::write_tb(sgRNA_tmp, "~/data/project/ear_project/gene_therapy_ll/Second/Batch2/sg_info.txt")
