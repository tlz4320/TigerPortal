human_sg_res <- read.table("~/data/project/ear_project/gene_therapy/find1NTSim/threshold30000/human_sg_final_anno_threshold30000.txt")
human_sg_res[,16] <- str_sub(human_sg_res[,2], 1, 19)
human_sg_res[,17] <- str_sub(human_sg_res[,9], 1, 19)
table(human_sg_res[,16] == human_sg_res[,17])
human_sg_res <- human_sg_res[human_sg_res[,7] == human_sg_res[,14],]
human_sg_res <- human_sg_res[human_sg_res[,16] != human_sg_res[,17],]
human_sg_res$res <- ""
for(i in 1 :nrow(human_sg_res)){
  res <- paste(human_sg_res[i, 8], "\n", human_sg_res[i, 15], sep="")
  human_sg_res$res[i] <- res
}

human_sg_res <- human_sg_res[,-c(8,15,16,17)]
colnames(human_sg_res) <- c("Gene", "sgRNA1", "sgRNA1_Chr", "sgRNA1_Pos", "sgRNA1_Strand", 
                            "sgRNA1_Type","sgRNA1_NGG", 
                            "sgRNA2", "sgRNA2_Chr", "sgRNA2_Pos", "sgRNA2_Strand", 
                            "sgRNA2_Type","sgRNA2_NGG", "Compare" )

second_sg <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Second1NT/second_batch_design_total160_replace_bad_fix.xlsx")
second_sg <- second_sg[1 : 140,]

table(second_sg$sgRNA %in% c(human_sg_res$sgRNA1, human_sg_res$sgRNA2))


human_sg_res$id <- apply(human_sg_res[,c("sgRNA1", "sgRNA2")], 1, function(x){
  paste0(sort(x), collapse = "-")
})
human_sg_res <- human_sg_res[!duplicated(human_sg_res$id),]
tmp <- lapply(seq(1, nrow(second_sg), 2), function(i){
  paste0(sort(second_sg$sgRNA[c(i, i + 1)]), collapse = "-")
})
human_sg_res <- human_sg_res[human_sg_res$id %in% tmp,]
tmp <- c(unlist(lapply(human_sg_res$Compare, function(x){
  unlist(strsplit(x, "[\n]"))[1]
})), unlist(lapply(human_sg_res$Compare, function(x){
  unlist(strsplit(x, "[\n]"))[2]
})))

tmp <- data.frame(sg = c(human_sg_res$sgRNA1, human_sg_res$sgRNA2), 
                  cmp = tmp)
rownames(tmp) <- tmp$sg
second_sg$Cmp <- tmp[second_sg$sgRNA, 2]
second_sg$ID2 <- str_replace_all(paste0("Sg_", second_sg$ID_new), "-", "_")
save(second_sg, file="second_sg.rda")
