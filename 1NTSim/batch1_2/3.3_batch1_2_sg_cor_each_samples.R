setwd("~/data/project/ear_project/gene_therapy_ll/batch2/")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240226-A00133B/indel_stat2.rda")
load("~/data/project/ear_project/gene_therapy_ll/batch2/240131-A00599//indel_stat.rda")
sgRNA <- read.xlsx("~/Nutstore Files/Tobin/First1NT/2024_1_12_integrated_design_result.xlsx")
sgRNA <- sgRNA[c(29 : nrow(sgRNA)),]
sgRNA$id2 <- paste("Sg", sgRNA$new_id, sep = "-")
sgRNA$id2 <- str_replace_all(sgRNA$id2, "-", "_")
max_region <- c(-200, 30)
batch2 <- list()
for(name in names(batch2_1)){
  n <- str_remove(name, "CRISPRessoPooled_on_")
  n <- str_remove(n, "[0-9]+")
  if(!n %in% names(batch2)){
    batch2[[n]] <- list()
  }
  for(sg in names(batch2_1[[name]])){
    batch2[[n]][[sg]] <- batch2_1[[name]][[sg]]
  }
}
for(name in names(batch2_2)){
  n <- str_remove(name, "CRISPRessoPooled_on_")
  n <- str_remove(n, "[0-9]+")
  if(!n %in% names(batch2)){
    batch2[[n]] <- list()
  }
  for(sg in names(batch2_2[[name]])){
    batch2[[n]][[sg]] <- batch2_2[[name]][[sg]]
  }
}
sgRNA_pair <- split(sgRNA$id2, sgRNA$id)
###计算各个sgRNA相关性在ABC样本两两之间的结果
total_sg_cor_each_sample <- list()

for(name in names(batch2)){
  edit_data <- batch2[[name]]
  sgRNA_pair_remain <- sgRNA_pair[unlist(lapply(sgRNA_pair, function(x){
    x <- unlist(x)
    return(sum(x %in% names(edit_data)) == 2)
  }))]
  tmp_cor <- lapply(sgRNA_pair_remain, function(y){
    ids <- unlist(y)
    sel_sg1 <- edit_data[[ids[1]]]
    sel_sg2 <- edit_data[[ids[2]]]
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
    cor.test(sel_sg1[,2], sel_sg2[,2])
  })
  names(tmp_cor) <- unlist(lapply(sgRNA_pair_remain, function(x){
    paste(x, collapse = "-")
  }))
  total_sg_cor_each_sample[[name]] <- data.frame(pair = names(tmp_cor), 
                                  cor = unlist(lapply(tmp_cor, function(x){x$estimate})),
                                  pval = unlist(lapply(tmp_cor, function(x){x$p.value})), 
                                  sample = name)
}

plot_data <- data.frame(do.call(rbind, total_sg_cor_each_sample))
plot_data$pair <- factor(plot_data$pair, levels = gtools::mixedsort(unique(plot_data$pair)))
pdf("Result/sgRNA_pair_edit_pattern_each_sample_cor.pdf", width = 30, height = 6)
ggplot(plot_data, aes(x = pair, y = cor)) + 
  geom_boxplot() +
  geom_point(aes(color = sample), position = position_dodge2(width = 0.2)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        text = element_text(size=15))
dev.off()


