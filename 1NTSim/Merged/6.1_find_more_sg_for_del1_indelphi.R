###这个脚本是用来筛选一些新的数据用于indelphi预测
left_right_nt_table <- read.table("~/data/project/ear_project/gene_therapy_ll/left_right_nt_table_after_filtered.txt")
left_right_nt_table_withbulge <- left_right_nt_table[seq(1, nrow(left_right_nt_table), 2),]
left_right_nt_table_withoutbulge <- left_right_nt_table[seq(2, nrow(left_right_nt_table), 2),]
table(left_right_nt_table_withbulge$V2)
table(left_right_nt_table_withoutbulge$V2)

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
human_sg_res <- human_sg_res[abs(human_sg_res$sgRNA1_Pos - human_sg_res$sgRNA2_Pos) > 200,]
human_sg_res$append <- unlist(lapply(human_sg_res$Compare, function(x){
  x <- unlist(strsplit(x, "[\n]"))
  if(str_sub(x[1], -1, -1) == "-" | str_sub(x[2], -1, -1) == "-"){
    return("Last")
  }
  return("First")
}))
human_sg_res <- human_sg_res[human_sg_res$append == 'First',]

batch2_design_last_sgRNA <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/human_sgRNA_lastNTSim_filter_pos.xlsx")
batch2_design_last_sgRNA <- batch2_design_last_sgRNA[,-ncol(batch2_design_last_sgRNA)]
batch2_design_last_sgRNA$append <- "First"
human_sg_res <- data.frame(rbind(human_sg_res, batch2_design_last_sgRNA))

total_design <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/First1NT/2024_1_12_integrated_design_result.xlsx")
total_design2 <- openxlsx::read.xlsx("~/Nutstore Files/Tobin/Second1NT/second_batch_design_total160_replace_bad_fix.xlsx")
total_design <- c(total_design$sgRNA, total_design2$sgRNA)
human_sg_res <- human_sg_res[!(human_sg_res$sgRNA1 %in% total_design | 
                               human_sg_res$sgRNA2 %in% total_design),]
human_sg_res$included <- unlist(lapply(human_sg_res$Compare, function(x){
  x <- unlist(strsplit(x, "[\n]"))
  bulgeone <- x[1]
  theotherone <- x[2]
  if(str_sub(x[2], 1, 1) == "-"){
    bulgeone <- x[2]
    theotherone <- x[1]
  }
  bulgeone <- str_remove_all(bulgeone, "[-]")
  theotherone <- str_remove_all(theotherone, "[-]")
  cutnt_withbulge <- str_sub(bulgeone, 16, 19)
  cutnt_withoutbulge <- str_sub(theotherone, 16, 19)
  cutnt_withbulge <- cutnt_withbulge %in% left_right_nt_table_withbulge$V2
  cutnt_withoutbulge <- cutnt_withoutbulge %in% left_right_nt_table_withoutbulge$V2
  return(2 - sum(c(cutnt_withbulge, cutnt_withoutbulge)))
}))
human_sg_res$Compare <- str_replace(human_sg_res$Compare, "\n", " ")
human_sg_res$Compare <- str_replace(human_sg_res$Compare, "\r", "")
human_sg_res <- human_sg_res[human_sg_res$included !=0,]
human_sg_res <- human_sg_res[human_sg_res$sgRNA1_Chr != "chrY",]
ToNX::write_tb(human_sg_res, file="~/data/project/ear_project/gene_therapy/for_find_del1_with2NT.txt")

sg_genome <- read.table("~/data/project/ear_project/gene_therapy_ll/for_find_del1_with2NT_genome.txt", sep="\t")
sg_genome <- sg_genome[-grep("(AAAAAAAAAA)|(TTTTTTTTTT)|(CCCCCCCCCC)|(GGGGGGGGGG)|(NNNNNNNNNN)]",sg_genome$V2),]

human_sg_res <- human_sg_res[human_sg_res$sgRNA1 %in% sg_genome$V1 & human_sg_res$sgRNA2 %in% sg_genome$V1,]
human_sg_res <- human_sg_res[order(human_sg_res$included,decreasing = T),]


human_sg_res$withBulge4NT <- unlist(lapply(human_sg_res$Compare, function(x){
  x <- unlist(strsplit(x, "[ ]"))
  bulgeone <- x[1]
  theotherone <- x[2]
  if(str_sub(x[2], 1, 1) == "-"){
    bulgeone <- x[2]
    theotherone <- x[1]
  }
  bulgeone <- str_remove_all(bulgeone, "[-]")
  theotherone <- str_remove_all(theotherone, "[-]")
  str_sub(bulgeone, 16, 19)
}))
human_sg_res$withoutBulge4NT <- unlist(lapply(human_sg_res$Compare, function(x){
  x <- unlist(strsplit(x, "[ ]"))
  bulgeone <- x[1]
  theotherone <- x[2]
  if(str_sub(x[2], 1, 1) == "-"){
    bulgeone <- x[2]
    theotherone <- x[1]
  }
  bulgeone <- str_remove_all(bulgeone, "[-]")
  theotherone <- str_remove_all(theotherone, "[-]")
  str_sub(theotherone, 16, 19)
}))
human_sg_res$seq <- unlist(lapply(1 : nrow(human_sg_res), function(x){
  paste(sort(c(human_sg_res$sgRNA1[x], human_sg_res$sgRNA2[x])),collapse = "-")
}))
human_sg_res <- human_sg_res[!duplicated(human_sg_res$seq),]

human_sg_res$id <- paste0("Pred_", 1 : nrow(human_sg_res))


dup_sg <- c(human_sg_res$sgRNA1, human_sg_res$sgRNA2)
dup_sg <- dup_sg[duplicated(dup_sg)]
human_sg_res_uniq <- human_sg_res[!(human_sg_res$sgRNA1 %in% dup_sg | human_sg_res$sgRNA2 %in% dup_sg ),]



for_indelphi <- data.frame(id = c(human_sg_res_uniq$id, human_sg_res_uniq$id),
                           sg = c(human_sg_res_uniq$sgRNA1, human_sg_res_uniq$sgRNA2))
sg_genome <- sg_genome[sg_genome$V1 %in% for_indelphi$sg,]
sg_genome <- sg_genome[!duplicated(sg_genome$V1),]
rownames(sg_genome) <- sg_genome$V1
for_indelphi$genomeSeq <- sg_genome[for_indelphi$sg, 2]


for_indelphi$cutpos <- unlist(lapply(1 : nrow(for_indelphi), function(i){
  sg <- for_indelphi$sg[i]
  genome <- for_indelphi$genomeSeq[i]
  unlist(str_locate(genome, sg)[1]) + 16
}))
for_indelphi <- for_indelphi[c(for_indelphi$cutpos %in% c(245 : 255)),]
for_indelphi$id2 <- paste0("Pred_", 1 : nrow(for_indelphi))
for_indelphi <- for_indelphi[for_indelphi$id %in% unlist(
  lapply(split(for_indelphi$id, for_indelphi$id), function(x){
    if(length(unlist(x)) == 2){
      return(unlist(x))
    }
    return("")
  })
),]
human_sg_res_uniq <- human_sg_res_uniq[human_sg_res_uniq$id %in% for_indelphi$id,]
for_indelphi$bulge <- unlist(lapply(for_indelphi$sg, function(x){
  cmp <- human_sg_res_uniq$Compare[human_sg_res_uniq$sgRNA1 == x | human_sg_res_uniq$sgRNA2 == x]
  cmp <- unlist(strsplit(cmp, "[ ]"))
  bulgeone <- cmp[1]
  theotherone <- cmp[2]
  if(str_sub(cmp[2], 1, 1) == "-"){
    bulgeone <- cmp[2]
    theotherone <- cmp[1]
  }
  bulgeone <- str_remove_all(bulgeone, "[-]")
  return(x == bulgeone)
}))

openxlsx::write.xlsx(for_indelphi, file="~/Nutstore Files/Tobin/del1_indelphi_predict.xlsx")
for_indelphi <- for_indelphi[,c("id2", "genomeSeq", "cutpos")]

for_indelphi$genomeSeq <- unlist(lapply(1 : nrow(for_indelphi), function(i){
  str_sub(for_indelphi$genomeSeq[i], for_indelphi$cutpos[i] - 49, for_indelphi$cutpos[i] + 50)
}))
for_indelphi$cutpos <- 50
write.table(for_indelphi, "~/for_del1_predict_indelphi.csv", col.names = T, row.names = F,quote = T,sep = ",")

for_indelphi <- read.csv("~/for_del1_predict_indelphi.csv")
forcast_file <- data.frame(ID = for_indelphi$id2,
                           Target = for_indelphi$genomeSeq,
                           "PAM Index" = for_indelphi$cutpos+ 3)
ToNX::write_tb(forcast_file, file="~/for_del1_predict_forecast.txt", col.names = F)
