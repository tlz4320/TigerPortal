setwd("~/data/project/ear_project/gene_therapy_ll/batch2/")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/indel_stat2.rda")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599//indel_stat.rda")
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/First1NT/2024_1_12_integrated_design_result.xlsx")
sgRNA <- sgRNA[c(29 : nrow(sgRNA)),]
sgRNA$id2 <- paste("Sg", sgRNA$new_id, sep = "-")
sgRNA$id2 <- str_replace_all(sgRNA$id2, "-", "_")
max_region <- c(-200, 30)
batch2 <- batch2_1
for(name in names(batch2_2)){
  batch2[[name]] <- batch2_2[[name]]
}
total_sg <- unique(unlist(lapply(batch2, names)))
batch2_rev <- lapply(total_sg, function(sg){
  res <- list()
  for(name in names(batch2)){
    n <- str_remove(name, "CRISPRessoPooled_on_")
    n <- str_remove(n, "[0-9]+")
    if(sg %in% names(batch2[[name]])){
      res[[n]] <- batch2[[name]][[sg]]
    }
  }
  res
})
names(batch2_rev) <- total_sg

###计算各个sgRNA相关性在ABC样本两两之间的结果
max_region <- c(-200, 30)
total_sg_cor <- lapply(total_sg, function(y){
  sel_sg <- batch2_rev[[y]]
  cmb <- combn(length(sel_sg), 2)
  cor_name <- lapply(1 : ncol(cmb), function(x){
    paste(names(sel_sg)[unlist(cmb[,x])], collapse = "-")
  })
  cor_res <- lapply(1 : ncol(cmb), function(x){
    ids <- unlist(cmb[,x])
    sel_sg1 <- sel_sg[[ids[1]]]
    sel_sg2 <- sel_sg[[ids[2]]]
    max_del <- min(c(sel_sg1[,1], sel_sg2[,1]))
    max_del <- max(c(max_del, max_region[1]))
    max_in <- max(c(sel_sg1[,1], sel_sg2[,1]))
    max_in <- min(c(max_in, max_region[2]))
    region <- max_del : max_in
    tmp_region <- region[!region %in% sel_sg1[,1]]
    if(length(tmp_region) != 0){
      sel_sg1 <- data.frame(rbind(sel_sg1, data.frame(indel_size = tmp_region, fq = 0)))
    }
    tmp_region <- region[!region %in% sel_sg2[,1]]
    if(length(tmp_region) != 0){
      sel_sg2 <- data.frame(rbind(sel_sg2, data.frame(indel_size = tmp_region, fq = 0)))
    }
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    sel_sg2[sel_sg2[,1] == 0, 2] <- 0
    common_indel <- intersect(sel_sg1[,1], sel_sg2[,1])
    sel_sg1 <- sel_sg1[sel_sg1[,1] %in% common_indel,]
    sel_sg2 <- sel_sg2[sel_sg2[,1] %in% common_indel,]
    sel_sg1 <- sel_sg1[order(sel_sg1[,1]),]
    sel_sg2 <- sel_sg2[order(sel_sg2[,1]),]
    # rm0 <- sel_sg1[,2] == 0 & sel_sg2[,2] == 0
    #cor.test(sel_sg1[!rm0,2], sel_sg2[!rm0,2])
    cor(sel_sg1[,2], sel_sg2[,2])
  })
  data.frame(name = unlist(cor_name), cor = unlist(cor_res), sg = y)
})
total_sg_cor <- data.frame(do.call(rbind, total_sg_cor))
total_sg_cor$sg <- factor(total_sg_cor$sg, 
                          levels = gtools::mixedsort(unique(total_sg_cor$sg)))

pdf("Result/cor_of_each_sg_in_each_two_sample.pdf", width = 30, height = 6)
ggplot(total_sg_cor, aes(x = sg, y = cor)) + 
  geom_boxplot() +
  geom_point(aes(color = name), position = position_dodge2(width = 0.2)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        text = element_text(size=15))
dev.off()


