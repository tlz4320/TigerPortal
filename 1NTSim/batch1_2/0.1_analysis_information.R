library(Biostrings)
library(ShortRead)
library(parallel)
library(stringr)
library(openxlsx)
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/2024_1_12_integrated_design_result.xlsx")
sgRNA <- sgRNA[c(49, 99 : 160),]
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
sgRNA_tmp$id <- paste("Sg", sgRNA$new_id, sep = "-")
sgRNA_tmp$id <- str_replace_all(sgRNA_tmp$id, "-", "_")
sgRNA_tmp <- sgRNA_tmp[,c(3,1,2)]
ToNX::write_tb(sgRNA_tmp, "~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/sg_info.txt")

ToNX::write_tb(sgRNA_tmp[,c(1,2)], "~/data/project/ear_project/gene_therapy_ll/batch2/sgRNA_seq2.txt")
