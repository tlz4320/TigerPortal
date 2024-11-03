allsrr <- read.csv("~/data2/indelphi_data/all_info.csv")
table(unlist(lapply(allsrr$SampleName, function(x){unlist(strsplit(x,"[_]"))[1]})))
allsrr2 <- allsrr[grep("libb", str_to_lower(allsrr$SampleName)),]
ToNX::write_tb(allsrr2[,1], file="~/data2/indelphi_data/U2OS_libb_list.txt")

allsrr <- allsrr[grep("U2OS", allsrr$SampleName),]
ToNX::write_tb(allsrr[,1], file="~/data2/indelphi_data/U2OS_list.txt")

fqfile <- ShortRead::FastqStreamer("/home/bioinfo/data2/indelphi_data/SRR7536346/SRR7536346.fastq")
tmp_fq <- list()
while(length(fq <- yield(fqfile)@sread)){
  tmp_fq[[length(tmp_fq) + 1]] <- as.character(fq)
}
fq <- unlist(tmp_fq)
fq <- data.frame(table(fq))
fq <- fq[order(fq$Freq,decreasing = T),]
rm(tmp_fq)

allsrr <- allsrr[grep("LibA", allsrr$SampleName),]

liba <- read.table("~/data2/indelphi_data/sginfo.txt", sep="\t")
dup <- liba[duplicated(liba$V2),]
liba <- liba[!liba$V2 %in% dup$V2,]
ToNX::write_tb(liba, "~/data2/indelphi_data/sginfo.txt")


