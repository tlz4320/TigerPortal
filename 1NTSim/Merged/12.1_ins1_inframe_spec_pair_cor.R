load("~/Nutstore Files/Tobin/Merged1NT/ins1_inframe_spec_result.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/first_second_sgCmp.rda")
bulge_pos <- read.table("~/data/project/ear_project/gene_therapy_ll/First_Second_sgRNA_cmp_reDiff.txt", sep="\t")
tmp <- bulge_pos[,c(4,5)]
tmp <- merge(tmp, sgCmp, by.x="V4", by.y="id")
tmp <- tmp[,c("V4","id2", "V5")]
inframe_result_add_pos <- merge(inframe_result_processed, tmp, by.x="id", by.y="id2")
paired_id <-inframe_result_add_pos[!duplicated(inframe_result_add_pos$id),]
paired_id <- split(paired_id$id, paired_id$V4)
inframe_cor <- lapply(paired_id, function(id){
  id <- unlist(id)
  one_inframe <- inframe_result_add_pos[inframe_result_add_pos$id == id[1],]
  two_inframe <- inframe_result_add_pos[inframe_result_add_pos$id == id[2],]
  
  one_inframe_array <- rep(0, 12)
  names(one_inframe_array) <- as.character(c(-5: 5, "Other"))
  two_inframe_array <- rep(0, 12)
  names(two_inframe_array) <- as.character(c(-5: 5, "Other"))
  for(i in 1 : nrow(one_inframe)){
    one_inframe_array[one_inframe$aa[i]]<-one_inframe$pct[i]
  }
  for(i in 1 : nrow(two_inframe)){
    two_inframe_array[two_inframe$aa[i]]<-two_inframe$pct[i]
  }
  cor(one_inframe_array, two_inframe_array)
})
inframe_cor <- data.frame(pair = names(inframe_cor), cor = unlist(inframe_cor))
inframe_cor <- merge(inframe_cor, tmp, by.x="pair", by.y="V4")
inframe_cor <- inframe_cor[!duplicated(inframe_cor$pair),]
pdf("~/Nutstore Files/Tobin/Merged1NT/ins1_pair_inframe_spec_cor.pdf", width = 7, height = 4)
ggplot(inframe_cor) + geom_boxplot(aes(x = as.character(V5), y = cor)) + 
  xlab("Bulge Pos") + ylab("Cor of ins1 inframe spec") + 
  theme_bw()
dev.off()

