tmp1 <- read.table("~/result2.txt", header = F, sep = "\t")
tmp2 <- read.table("/home/bioinfo/data/share/extract_result2.txt", header = F, sep = "\t")
tmp1$id <- paste0(tmp1$V1, tmp1$V2)
tmp2$id <- paste0(tmp2$V1, tmp2$V2)
tmp <- merge(tmp1, tmp2, by="id")
tmp <- tmp[,-1]
tmp <- tmp[,c(-6,-7)]
tmp <- tmp[order(tmp[,4], decreasing = T),]
colnames(tmp) <- c("Chr", "Start", "End", "Similar", "Pattern", "Seq")
openxlsx::write.xlsx(tmp, file="/home/bioinfo/data/share/off_target_blast_result.xlsx", rowNames=F, colNames=T)
