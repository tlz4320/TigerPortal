#mOTOF CDS
#TCCCTCTGCTCGGTAAATTTTCACATAGAACCGTGCCCACTGCCGTTC{G}GGGGGCACGCCCTCGGGGAGCAGCAAGTTCCCTTCAATGTCGTCCTCATC
checked <- checkRS("rs397515581")
if(class(checked) == "logical"){
  if(!checked){
    shinyalert("Not a valid RS ID", "Not a valid RS ID", type="error")
    return()
  }
}
inputseq <- "GATGAGGACGACATTGAAGGGAACTTGCTGCTCCCCGAGGGCGTGCCCCC{C}GAACGGCAGTGGGCACGGTTCTATGTGAAAATTTACCGAGCAGAGGGA"
cds_region <- c(0, 0)
checked <- formatCheck(inputseq)
if(class(checked) == "logical"){
  if(!checked){
    shinyalert("Not a valid input", "Not a valid input", type="error")
    return()
  }
}
indel_pos <- checked$pos
seq <- checked$seq
isIndel <- checked$isIndel
indelNT <- checked$indelNT
possiblePAM <- findPAM(seq, pam = "NGG")
if(T){
  reverseSeq <- as.character(reverseComplement(DNAString(seq)))
  possiblePAM2 <- findPAM(reverseSeq, pam = "NGG")
} else {
  #只是为了方便下面判断而已
  possiblePAM2 <- possiblePAM
}
if(nrow(possiblePAM) != 0){
  possiblePAM$dis <- abs(possiblePAM$start - 4 - indel_pos)
  coloredSeq_f <- colorSeq(seq, possiblePAM, indel_pos, isIndel,indelNT)
  coloredSeq_f$strand <- "+"
} else {
  coloredSeq_f <- data.frame()
}
if(T & nrow(possiblePAM2) != 0){
  possiblePAM2$dis <- abs(possiblePAM2$start - 4 - (str_length(seq) - indel_pos + 1))
  coloredSeq_r <- colorSeq(reverseSeq, possiblePAM2, str_length(reverseSeq) - indel_pos + 1, isIndel,indelNT)
  coloredSeq_r$strand <- "-"
} else {
  coloredSeq_r <- data.frame()
}

if(nrow(coloredSeq_r) != 0 & nrow(coloredSeq_f) != 0){
  coloredSeq <- data.frame(rbind(coloredSeq_f, coloredSeq_r))
} else if(nrow(coloredSeq_r) == 0){
  coloredSeq <- coloredSeq_f
} else {
  coloredSeq <- coloredSeq_r
}

#f(isIndel != 0){
  coloredSeq <- coloredSeq[order(coloredSeq$dis),]
  coloredSeq <- coloredSeq[coloredSeq$dis < 20,]
#}
checkBoxs <- unlist(lapply(1 : nrow(coloredSeq), function(i){
  as.character(checkboxInput(paste0("id_", i),NULL, value = i <= 3, width = "10%"))
}))
df2 <- data.frame(Seq = coloredSeq$colored, 
                  Strand = coloredSeq$strand, 
                  choose=checkBoxs)
infos <- coloredSeq
seq <- seq
infos2 <- checked

sgInfos <- infos
sel <- unlist(lapply(1 : nrow(sgInfos), function(i){
  T
}))
sgInfos$id <- paste0("id_", 1 : nrow(sgInfos))
sgInfos_sel <- sgInfos[sel,]
for_indelphi <- lapply(1 : nrow(sgInfos_sel), function(i){
  
  strand <- sgInfos_sel$strand[i]
  start <- sgInfos_sel$start[i]
  if(strand == "+"){
    seq <- seq
  } else {
    seq <- reverseSeq
  }
  data.frame(id = sgInfos_sel$id[i],
             seq = seq, pos = start - 4)
  
})
for_indelphi <- data.frame(do.call(rbind, for_indelphi))
output_file <- tempfile(fileext = ".csv")
td <- tempfile(pattern = "Result")
dir.create(td, recursive = T)
write.table(for_indelphi, file = output_file, col.names=T, row.names=F, quote = T, sep=",")
system(paste("/home/bioinfo/micromamba/envs/bio/bin/python /home/bioinfo/code/inDelphi-model/script4R.py",
             "U2OS", output_file, td))

readPredict(td, sgInfos_sel, infos2, cds_region)
nowwd <- getwd()
setwd(tempdir())
system(paste0("tar -czf ", basename(td), ".tar.gz ", basename(td)))
setwd(nowwd)
